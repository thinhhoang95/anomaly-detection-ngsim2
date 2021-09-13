from simulator import Simulator as NGSIM

ngsim = NGSIM('../kalmaned2.csv', 400000, 400, make_first_global_time_zero=True)

while ngsim.step():
    print('Sim time ' + str(ngsim.gtime))

print('Simulation completed')