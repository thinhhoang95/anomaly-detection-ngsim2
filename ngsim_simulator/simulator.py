import numpy as np
import pandas as pd

# add radia library to the search directory
import sys
sys.path.append("C:\\Users\\nxf67027\\Documents\\anomaly-detection-ngsim")
# sys.path.append("C:/Users/nxf67027/Documents/anomaly-detection-ngsim/otqd")

from radia.GlobalMapObject import GlobalMapObject as gmo
from radia import raytrace as rt
from otqd.otqd import OTQD

import copy
from radia.visualize_traffic_situation import plot_ngsim_memory

class Simulator:
    def __init__(self, path_to_csv, gtime_start, gtime_duration, offset_x, offset_y, make_first_global_time_zero = True, otqd_on = True):
        print('NGSIM Simulator: NGSIM-CSV file load')
        self.df = pd.read_csv(path_to_csv)
        self.gtime_start = gtime_start
        self.gtime_duration = gtime_duration
        self.gtime = gtime_start
        self.otqd_on = otqd_on # whether to trigger OTQD score calculation or not
        if make_first_global_time_zero:
            # Subtract the global time with the first entry of Global_Time in the DataFrame
            self.df['Global_Time'] = self.df['Global_Time'] - self.df['Global_Time'][0]
        self.vehicle_array_memory = []

        # Offset trajectory at the beginning (for making the trajectories starting at the same point - visualize the FDGRX graph).
        # Usually the offset_y and offset_x will coincide with the patch coordinates
        self.trajectory_offset_y = offset_y
        self.trajectory_offset_x = offset_x

        # Flags for initialization of other modular functionalities
        self.otqd_initialized = False
        self.in_patch_filter_initialized = False # to ensure the patch parameters must be initialized

        # Patch parameters
        self.x_min = 0
        self.x_max = 0
        self.y_min = 0
        self.y_max = 0

    def initialize_otqd(self, info_a, mu_a, info_e2, pca_mean, pca_components, i_max = 1, i_spacing = 10):
        self.info_a = info_a.copy()
        self.mu_a = mu_a.copy()
        # Measurement noise covariance
        self.info_e2 = info_e2
        # FPCA components
        self.pca_mean = pca_mean
        self.pca_components = pca_components
        # Hypothesis on when the observation started
        self.i_max = i_max
        self.i_spacing = i_spacing
        print('OTQD module initialized successfully')

        self.otqd_initialized = True

    def initialize_in_patch_filter(self, x_min, x_max, y_min, y_max, override_offsets = False):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.in_patch_filter_initialized = True
        if override_offsets:
            self.trajectory_offset_y = y_min
            self.trajectory_offset_x = x_min
        else:
            print('Warning: the trajectory offset values are set manually and could be different from the patch coordinates. This might lead to unexpected behaviours.')

    def get_df_trajs_in_patch(self, df):
        if not self.in_patch_filter_initialized:
            raise Exception("In patch filter was called but was not initialized. Did you forget to call initialize_in_patch_filter?")
        return df[(df['Local_X'] > self.x_min) &
                       (df['Local_X'] < self.x_max) &
                       (df['Local_Y'] > self.y_min) &
                       (df['Local_Y'] < self.y_max)]

    def get_vehicle_list_at_gtime(self, gtime, filter_in_patch = True):
        df_at_gtime = self.df[self.df['Global_Time']==gtime]
        if filter_in_patch:
            df_at_gtime = self.get_df_trajs_in_patch(df_at_gtime)
        vehicle_ids = df_at_gtime['Vehicle_ID'].unique()
        vehicle_array = []
        for vehicle_id in vehicle_ids:
            vehicle_row = df_at_gtime[df_at_gtime['Vehicle_ID'] == vehicle_id].iloc[0]
            vehicle = gmo(vehicle_row.Vehicle_ID, (vehicle_row.Local_X, vehicle_row.Local_Y), vehicle_row.Vehicle_EKF_Velocity, vehicle_row.Vehicle_EKF_Accel, vehicle_row.Vehicle_EKF_Theta + np.pi/2)
            vehicle_array.append(vehicle)
        return vehicle_array

    def find_state_in_vehicle_array(self, vehicle_id, vehicle_array):
        for vehicle in vehicle_array:
            if vehicle_id == vehicle.id:
                return vehicle.posx

    def trace_at_gtime(self, vehicle_array_with_mem, vehicle_array):
        # vehicle_array_with_mem = vehicle_array.copy()
        # perform ray tracing
        for vehicle in vehicle_array:
            rt.frontRadarRayTrace(vehicle, vehicle_array)

        # compare the state before and after trace. Vehicles that did not disappear will have its state added to the vehicle's memory
        vehicle_array_with_mem_ids = [v.id for v in vehicle_array_with_mem]
        vehicle_array_ids = [v.id for v in vehicle_array]

        vehicles_disappeared = [x for x in vehicle_array_with_mem_ids if x not in vehicle_array_ids] # in memory, but not now so the vehicle has disappeared
        vehicles_appeared = [x for x in vehicle_array_ids if x not in vehicle_array_with_mem_ids] # now, but not in memory so these vehicles are new

        # firstly, remove the vehicles disappeared from the memory
        for idx, vehicle in enumerate(vehicle_array_with_mem.copy()):
            # remove the vehicle from the vehicle list
            if vehicle.id in vehicles_disappeared:
                vehicle_array_with_mem.pop(idx)
            # consider remove this vehicle from all vehicles' memory?

        # secondly, add the vehicles freshly appeared in this scan
        for vehicle in vehicle_array:
            if vehicle.id in vehicles_appeared:
                vehicle_array_with_mem.append(vehicle)

        # thirdly, update the state of the vehicles in the vehicle_array_with_mem
        for vehicle_mem in vehicle_array_with_mem:
            for vehicle in vehicle_array: 
                if vehicle.id == vehicle_mem.id:
                    vehicle_mem.posx = vehicle.posx 
                    vehicle_mem.orientation = vehicle.orientation
                    vehicle_mem.spd = vehicle.spd
                    vehicle_mem.acc = vehicle.acc 
                    vehicle_mem.neighbors = vehicle.neighbors
                    # Do NOT transfer the memory!

        # fourthly, append the state of neighbors acquired from the current scan
        for index, vehicle in enumerate(vehicle_array_with_mem):
            for neighbor_id in vehicle_array[index].neighbors:
                neighbor_state = self.find_state_in_vehicle_array(neighbor_id, vehicle_array)
                if neighbor_id in vehicle.memory:
                    # This is an already existed vehicle
                    vehicle.memory[neighbor_id].append(neighbor_state)
                else:
                    # This is a freshly added vehicle
                    vehicle.memory[neighbor_id] = [neighbor_state]
                # Trigger OTQD calculation if configured
                if self.otqd_on:
                    self.otqd_trigger(vehicle.id, neighbor_id)

        return vehicle_array_with_mem

    def restart_simulation(self):
        self.gtime = self.gtime_start

    def step(self):
        if self.gtime == self.gtime_start:
            # Initialization
            self.vehicle_array_memory = self.get_vehicle_list_at_gtime(self.gtime)
            print('Initialized vehicle list for the first time')
            self.gtime += 100
            return True

        if self.gtime < self.gtime_start + self.gtime_duration:
            vehicle_array = self.get_vehicle_list_at_gtime(self.gtime)
            self.trace_at_gtime(self.vehicle_array_memory, vehicle_array)
            self.gtime += 100
            return True
        else:
            print('The simulation has ended')
            return False

    def get_vehicle_by_id(self, vehicle_id):
        """
        Return the vehicle whose id is vehicle_id from vehicle_array_memory
        :param vehicle_id: ID of vehicle
        :return: vehicle object from vehicle_array_memory
        """
        for vehicle in self.vehicle_array_memory:
            if vehicle.id == vehicle_id:
                return vehicle
        raise Exception("Vehicle lookup unsuccessful: " + str(vehicle_id))

    def otqd_trigger(self, vehicle_id, neighbor_id):
        """
        Perform OTQD on the memory of all objects
        :return:
        """

        if not self.otqd_initialized:
            raise Exception("OTQD was not initialized for the simulator.")

        vehicle = self.get_vehicle_by_id(vehicle_id)
        trajectory = vehicle.memory[neighbor_id] # lookup the trajectory by the key, which is neighbor_id

        # Initialize an OTQD for each individual trajectory
        # For Y components
        otqd = OTQD(self.info_a, mu_a=self.mu_a, info_e2=10., pca_mean=self.pca_mean, pca_components=self.pca_components, i_max=self.i_max, i_spacing=self.i_spacing)
        likelihood_curve = np.zeros((len(trajectory),))
        for k in range(len(trajectory)):
            otqd.new_measurement(trajectory[k][1] - self.trajectory_offset_y) # for Y component
            likelihood_curve[k] = np.max(otqd.calculate_log_likelihood()) # best trajectory initial time
        # TODO: For X components

        # Save the estimated likelihood curve
        vehicle.otqd_entropy[neighbor_id] = likelihood_curve