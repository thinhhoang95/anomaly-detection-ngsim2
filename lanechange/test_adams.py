import sys
sys.path.append('/Users/thinhhoang/Documents/anomaly-detection-ngsim/lanechange')
import numpy as np
import onlinebayesian as olbayes
import pickle


amean = np.zeros((2,1))
acov = np.array([[4027.36937562,  -10.59873913],[ -10.59873913,   79.46485122]])
info_e2 = 1e-4


fpca_discretized = pickle.load(open('/Users/thinhhoang/Documents/anomaly-detection-ngsim/lanechange/fpca.pyo', 'rb'))
basis = fpca_discretized.components_.data_matrix
basis = np.reshape(basis[:2], (-1,basis.shape[1])).transpose()
mean = fpca_discretized.mean_.data_matrix.reshape((-1,1))
bayesFactorizerParams = {
    'info_a' : np.linalg.inv(acov) * 1.,
    'mu_a': amean, 
    'info_e2': info_e2, 
    'pca_mean': mean,
    'pca_components': basis, 
    'i_max': 1
}
detector = olbayes.OnlineBayesian(bayesFactorizerParams)
sample_traj = pickle.load(open('/Users/thinhhoang/Documents/anomaly-detection-ngsim/lanechange/sample.pyo', 'rb'))

# Start our shiny detector
rt_posterior = np.zeros((sample_traj.shape[0],sample_traj.shape[0]))
for t in range(sample_traj.shape[0]):
    detector.new_datum(sample_traj[t])
    for rt in range(detector.t):
        rt_posterior[detector.t, rt] = detector.calculate_rl_posterior(rt)
    detector.shift_time()


print('Done')