import random

import pandas as pd


class NGSIM_Processor:
    def __init__(self, df):
        print('NGSIM Processor is initialized')
        self.x_min = 0
        self.x_max = 0
        self.y_min = 0
        self.y_max = 0
        self.df = df
        self.max_trace_steps = 200

    def set_patch(self, x_min, x_max, y_min, y_max):
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max

    def get_df_trajs_in_patch(self):
        return self.df[(self.df['Local_X'] > self.x_min) &
                       (self.df['Local_X'] < self.x_max) &
                       (self.df['Local_Y'] > self.y_min) &
                       (self.df['Local_Y'] < self.y_max)]

    def is_veh_in_patch(self, lx, ly):
        if lx < self.x_min or lx > self.x_max or ly < self.y_min or ly > self.y_max:
            return False
        else:
            return True

    def get_trajectory_random_in_patch(self):
        """
        To return part of the trajectory where a random vehicle stays in the predefined road patch
        :return: traj_t, traj_x, traj_y
        """
        df_patch_trajectories = self.get_df_trajs_in_patch()
        num_of_rows = df_patch_trajectories.shape[0]
        random_row_index = random.randrange(0, num_of_rows - 1)
        vehicle_id_of_the_random_row = df_patch_trajectories.iloc[random_row_index].Vehicle_ID
        gtime_at_random_row = df_patch_trajectories.iloc[random_row_index].Global_Time
        # State variables for saving the trajectory
        traj_x = []
        traj_y = []
        traj_t = []
        traj_v = []

        gtime = gtime_at_random_row
        # Trace backward in time until the vehicle is out of the patch
        while True:
            df_query = self.df[(self.df['Vehicle_ID'] == vehicle_id_of_the_random_row) & (self.df['Global_Time'] == gtime)]
            if df_query.shape[0] == 0:
                break
            row_of_interest = df_query.iloc[0]
            lx = row_of_interest.Local_X
            ly = row_of_interest.Local_Y
            v = row_of_interest.Vehicle_EKF_Velocity
            if not self.is_veh_in_patch(lx, ly):
                break
            else:
                traj_x.append(lx)
                traj_y.append(ly)
                traj_t.append(gtime)
                traj_v.append(v)
                gtime -= 100

        gtime = gtime_at_random_row
        # Trace forward in time until the vehicle is out of the patch
        while True:
            gtime += 100
            df_query = self.df[
                (self.df['Vehicle_ID'] == vehicle_id_of_the_random_row) & (self.df['Global_Time'] == gtime)]
            if df_query.shape[0] == 0:
                break
            row_of_interest = df_query.iloc[0]
            lx = row_of_interest.Local_X
            ly = row_of_interest.Local_Y
            v = row_of_interest.Vehicle_EKF_Velocity
            if not self.is_veh_in_patch(lx, ly):
                break
            else:
                traj_x.append(lx)
                traj_y.append(ly)
                traj_t.append(gtime)
                traj_v.append(v)

        # Perform a final sort
        traj_x = [x for _, x in sorted(zip(traj_t, traj_x))]
        traj_y = [x for _, x in sorted(zip(traj_t, traj_y))]
        traj_v = [x for _, x in sorted(zip(traj_t, traj_v))]
        traj_t.sort()

        # print("Trajectory length " + str(len(traj_t)) + " generated")

        return traj_t, traj_x, traj_y

    def get_trajectory_random_in_patch_for_a_duration(self, duration):
        """
        To return part of the trajectory where a random vehicle stays in the predefined road patch
        :return: traj_t, traj_x, traj_y
        """
        df_patch_trajectories = self.get_df_trajs_in_patch()
        num_of_rows = df_patch_trajectories.shape[0]
        random_row_index = random.randrange(0, num_of_rows - 1)
        vehicle_id_of_the_random_row = df_patch_trajectories.iloc[random_row_index].Vehicle_ID
        gtime_at_random_row = df_patch_trajectories.iloc[random_row_index].Global_Time
        # State variables for saving the trajectory
        traj_x = []
        traj_y = []
        traj_t = []
        traj_v = []

        # print('HEYEHE')

        gtime = gtime_at_random_row
        # Trace backward in time until the vehicle is out of the patch
        while True:
            df_query = self.df[(self.df['Vehicle_ID'] == vehicle_id_of_the_random_row) & (self.df['Global_Time'] == gtime)]
            if df_query.shape[0] == 0:
                break
            row_of_interest = df_query.iloc[0]
            lx = row_of_interest.Local_X
            ly = row_of_interest.Local_Y
            v = row_of_interest.Vehicle_EKF_Velocity
            if not self.is_veh_in_patch(lx, ly):
                break
            else:
                traj_x.append(lx)
                traj_y.append(ly)
                traj_t.append(gtime)
                traj_v.append(v)
                gtime -= 100

        relative_gtime = gtime

        gtime = gtime_at_random_row
        # Trace forward in time for a duration specified by the parameter "duration"

        while gtime <= relative_gtime + duration:
            gtime += 100
            # Query the vehicle in question at selected gtime
            df_query = self.df[
                (self.df['Vehicle_ID'] == vehicle_id_of_the_random_row) & (self.df['Global_Time'] == gtime)]
            if df_query.shape[0] == 0:
                print('The vehicle {} disappeared before completing the trj of duration {}'.format(vehicle_id_of_the_random_row, duration))
                break
            row_of_interest = df_query.iloc[0]
            lx = row_of_interest.Local_X
            ly = row_of_interest.Local_Y
            v = row_of_interest.Vehicle_EKF_Velocity
            traj_x.append(lx)
            traj_y.append(ly)
            traj_t.append(gtime)
            traj_v.append(v)

        #if len(traj_t) < int(duration/100):
        #    raise Exception('Trajectory length is {}'.format(len(traj_t)))

        # print(len(traj_t))

        # Perform a final sort
        traj_x = [x for _, x in sorted(zip(traj_t, traj_x))]
        traj_y = [x for _, x in sorted(zip(traj_t, traj_y))]
        traj_v = [x for _, x in sorted(zip(traj_t, traj_v))]
        traj_t.sort()

        # print("Trajectory length " + str(len(traj_t)) + " generated")

        return traj_t, traj_x, traj_y

    def get_trajectory_random_in_patch_with_velocity(self):
            """
            To return part of the trajectory where a random vehicle stays in the predefined road patch
            :return: traj_t, traj_x, traj_y
            """
            df_patch_trajectories = self.get_df_trajs_in_patch()
            num_of_rows = df_patch_trajectories.shape[0]
            random_row_index = random.randrange(0, num_of_rows - 1)
            vehicle_id_of_the_random_row = df_patch_trajectories.iloc[random_row_index].Vehicle_ID
            gtime_at_random_row = df_patch_trajectories.iloc[random_row_index].Global_Time
            # State variables for saving the trajectory
            traj_x = []
            traj_y = []
            traj_t = []
            traj_v = []

            gtime = gtime_at_random_row
            # Trace backward in time until the vehicle is out of the patch
            while True:
                df_query = self.df[(self.df['Vehicle_ID'] == vehicle_id_of_the_random_row) & (self.df['Global_Time'] == gtime)]
                if df_query.shape[0] == 0:
                    break
                row_of_interest = df_query.iloc[0]
                lx = row_of_interest.Local_X
                ly = row_of_interest.Local_Y
                v = row_of_interest.Vehicle_EKF_Velocity
                if not self.is_veh_in_patch(lx, ly):
                    break
                else:
                    traj_x.append(lx)
                    traj_y.append(ly)
                    traj_t.append(gtime)
                    traj_v.append(v)
                    gtime -= 100

            gtime = gtime_at_random_row
            # Trace forward in time until the vehicle is out of the patch
            while True:
                gtime += 100
                df_query = self.df[
                    (self.df['Vehicle_ID'] == vehicle_id_of_the_random_row) & (self.df['Global_Time'] == gtime)]
                if df_query.shape[0] == 0:
                    break
                row_of_interest = df_query.iloc[0]
                lx = row_of_interest.Local_X
                ly = row_of_interest.Local_Y
                v = row_of_interest.Vehicle_EKF_Velocity
                if not self.is_veh_in_patch(lx, ly):
                    break
                else:
                    traj_x.append(lx)
                    traj_y.append(ly)
                    traj_t.append(gtime)
                    traj_v.append(v)

            # Perform a final sort
            traj_x = [x for _, x in sorted(zip(traj_t, traj_x))]
            traj_y = [x for _, x in sorted(zip(traj_t, traj_y))]
            traj_v = [x for _, x in sorted(zip(traj_t, traj_v))]
            traj_t.sort()

            # print("Trajectory length " + str(len(traj_t)) + " generated")

            return traj_t, traj_x, traj_y, traj_v

    def get_trajectory_random_in_patch_with_velocity_and_acceleration(self):
            """
            To return part of the trajectory where a random vehicle stays in the predefined road patch
            :return: traj_t, traj_x, traj_y
            """
            df_patch_trajectories = self.get_df_trajs_in_patch()
            num_of_rows = df_patch_trajectories.shape[0]
            random_row_index = random.randrange(0, num_of_rows - 1)
            vehicle_id_of_the_random_row = df_patch_trajectories.iloc[random_row_index].Vehicle_ID
            gtime_at_random_row = df_patch_trajectories.iloc[random_row_index].Global_Time
            # State variables for saving the trajectory
            traj_x = []
            traj_y = []
            traj_t = []
            traj_v = []
            traj_a = []

            gtime = gtime_at_random_row
            # Trace backward in time until the vehicle is out of the patch
            while True:
                df_query = self.df[(self.df['Vehicle_ID'] == vehicle_id_of_the_random_row) & (self.df['Global_Time'] == gtime)]
                if df_query.shape[0] == 0:
                    break
                row_of_interest = df_query.iloc[0]
                lx = row_of_interest.Local_X
                ly = row_of_interest.Local_Y
                v = row_of_interest.Vehicle_EKF_Velocity
                a = row_of_interest.Vehicle_EKF_Accel
                if not self.is_veh_in_patch(lx, ly):
                    break
                else:
                    traj_x.append(lx)
                    traj_y.append(ly)
                    traj_t.append(gtime)
                    traj_v.append(v)
                    traj_a.append(a)
                    gtime -= 100

            gtime = gtime_at_random_row
            # Trace forward in time until the vehicle is out of the patch
            while True:
                gtime += 100
                df_query = self.df[
                    (self.df['Vehicle_ID'] == vehicle_id_of_the_random_row) & (self.df['Global_Time'] == gtime)]
                if df_query.shape[0] == 0:
                    break
                row_of_interest = df_query.iloc[0]
                lx = row_of_interest.Local_X
                ly = row_of_interest.Local_Y
                v = row_of_interest.Vehicle_EKF_Velocity
                a = row_of_interest.Vehicle_EKF_Accel
                if not self.is_veh_in_patch(lx, ly):
                    break
                else:
                    traj_x.append(lx)
                    traj_y.append(ly)
                    traj_t.append(gtime)
                    traj_v.append(v)
                    traj_a.append(a)

            # Perform a final sort
            traj_x = [x for _, x in sorted(zip(traj_t, traj_x))]
            traj_y = [x for _, x in sorted(zip(traj_t, traj_y))]
            traj_v = [x for _, x in sorted(zip(traj_t, traj_v))]
            traj_a = [x for _, x in sorted(zip(traj_t, traj_a))]
            traj_t.sort()

            # print("Trajectory length " + str(len(traj_t)) + " generated")

            return traj_t, traj_x, traj_y, traj_v, traj_a