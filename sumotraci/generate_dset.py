 # File: generate_dset.py 
 # Author: Hoang Dinh Thinh 
 # Date: 20/04/2022 
 # Description: Generate a trajectory set from a SUMO simulation scenario 
 # Params: the sumocfg file in the sumoCmd variable

import os
import sys
if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools') # for Macbook Pro machine
    sys.path.append(tools)
else:
    sys.exit("please declare environment variable 'SUMO_HOME'")

import traci
sumoBinary = "sumo-gui" # for Macbook Pro machine
sumoCmd = [sumoBinary, "-c", "highway/highway.sumocfg"] # for Macbook Pro machine
traci.start(sumoCmd)

# Delete the file if exists
try:
    os.remove('trajectory.csv')
except OSError:
    print('Cannot remove the file trajectory.csv')
    print(OSError.strerror)
    pass
with open('trajectory.csv', 'a') as csvFile:
    step = 0
    while step < 300:
        traci.simulationStep()
        listAllVehiclesId = traci.vehicle.getIDList()
        for vehId in listAllVehiclesId:
            rowToWrite = str(step) + ","
            rowToWrite += str(vehId) + "," # vehicle ID
            rowToWrite += str(traci.vehicle.getRoadID(vehId)) + ","
            rowToWrite += str(traci.vehicle.getLaneID(vehId)) + ","
            pos = traci.vehicle.getPosition(vehId)
            rowToWrite += str(pos[0]) + ","
            rowToWrite += str(pos[1]) + ","
            csvFile.write(rowToWrite + "\n")
        step += 1
traci.close()