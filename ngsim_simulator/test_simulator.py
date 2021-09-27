from simulator import Simulator as NGSIM
import numpy as np
import pickle as pkl
import sys

from radia.visualize_traffic_situation import plot_ngsim_memory

# sys.path.append("~/Documents/anomaly-detection-ngsim/otqd")

ngsim = NGSIM('kalmaned2.csv', 290000, 3000, offset_x=0.0, offset_y=500.0, make_first_global_time_zero=True) # the offset x, y is the patch start coordinate
# Initialize the parameters for OTQD module of NGSIM
gmy = pkl.load(open("y.pkl", "rb"))

ngsim.initialize_otqd(info_a = np.linalg.inv(gmy['cov']), mu_a = gmy['mu'].transpose(), info_e2 = 10., pca_mean = gmy['mean_curve'], pca_components = gmy['friendly_basis'], i_max = 10, i_spacing = 5)
ngsim.initialize_in_patch_filter(0.0, 75.0, 500.0, 625.0)
while ngsim.step():
    print('Sim time ' + str(ngsim.gtime))
    if (ngsim.gtime == ngsim.gtime_start + ngsim.gtime_duration - 100):
        plot_ngsim_memory(ngsim.vehicle_array_memory)

# Some statistical information for debugging
print('Vehicles ranking that have neighbors: ')
print([id for id, x in enumerate(ngsim.vehicle_array_memory) if len(x.neighbors) > 0])


print('Simulation completed')