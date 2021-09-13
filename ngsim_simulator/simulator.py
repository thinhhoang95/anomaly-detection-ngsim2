import numpy as np
import pandas as pd

# add radia library to the search directory
import sys
sys.path.append("~/Documents/anomaly-detection-ngsim/radia")

from radia.GlobalMapObject import GlobalMapObject as gmo
from radia import raytrace as rt

import copy

class Simulator:
    def __init__(self, path_to_csv, gtime_start, gtime_duration, make_first_global_time_zero = True):
        print('NGSIM Simulator: NGSIM-CSV file load')
        self.df = pd.read_csv(path_to_csv)
        self.gtime_start = gtime_start
        self.gtime_duration = gtime_duration
        self.gtime = gtime_start
        if make_first_global_time_zero:
            self.df['Global_Time'] = self.df['Global_Time'] - self.df['Global_Time'][0]
        self.vehicle_array_memory = []

    def get_vehicle_list_at_gtime(self, gtime):
        df_at_gtime = self.df[self.df['Global_Time']==gtime]
        vehicle_ids = df_at_gtime['Vehicle_ID'].unique()
        vehicle_array = []
        for vehicle_id in vehicle_ids:
            vehicle_row = df_at_gtime[df_at_gtime['Vehicle_ID'] == vehicle_id].iloc[0]
            vehicle = gmo(vehicle_row.Vehicle_ID, (vehicle_row.Local_X, vehicle_row.Local_Y), vehicle_row.Vehicle_EKF_Velocity, vehicle_row.Vehicle_EKF_Accel, vehicle_row.Vehicle_EKF_Theta)
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

        # thirdly, append the state of neighbors acquired from the current scan
        for index, vehicle in enumerate(vehicle_array_with_mem):
            for neighbor_id in vehicle_array[index].neighbors:
                neighbor_state = self.find_state_in_vehicle_array(neighbor_id, vehicle_array)
                if neighbor_id in vehicle.memory:
                    # This is an already existed vehicle
                    vehicle.memory[neighbor_id].append(neighbor_state)
                else:
                    # This is a freshly added vehicle
                    vehicle.memory[neighbor_id] = [neighbor_state]

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