import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# add radia library to the search directory
import sys
sys.path.append("~/Documents/anomaly-detection-ngsim/radia")
sys.path.append("~/Documents/anomaly-detection-ngsim/otqd")

def plot_ngsim_memory(memory):
    """
    Plot the current traffic situation to verify the correctness of the radar ray tracing
    :param memory: from ngsim.vehicle_array_memory
    :return: None
    """

    # all_x = [x.posx[0] for x in memory]
    # all_y = [x.posx[1] for x in memory]
    #
    # max_x = max(all_x)
    # min_x = min(all_x)
    # max_y = max(all_y)
    # min_y = min(all_y)

    for vehicle in memory:
        plt.scatter(vehicle.posx[0], vehicle.posx[1])
        plt.text(vehicle.posx[0] - 2, vehicle.posx[1] - 2, vehicle.id)

    plt.show()