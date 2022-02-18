# %% [markdown]
# # ASPTest Script
# 
# We will test our Online Bayesian Detector customized for detection of the presence of a stochastic process.

# %% [markdown]
# # Load data

# %%
import pickle
import numpy as np
from matplotlib import pyplot as plt 
import sys 

sys.path.append('/Users/thinhhoang/Documents/anomaly-detection-ngsim/lanechange')

pfile = pickle.load(open( "lanechange/lane.pyo", "rb" ))
tvec = pfile['t_vec'][0]
xvec = pfile['x_vec']

# for i in range(xvec.shape[0]):
#     plt.plot(tvec, xvec[i,:])

# plt.title('Lateral trajectories')

# %% [markdown]
# # Load FPCA

# %%
import pickle
fpca_discretized = pickle.load(open('lanechange/fpca.pyo', 'rb'))
basis = fpca_discretized.components_.data_matrix.reshape((-1,fpca_discretized.components_.data_matrix.shape[1])).transpose()
basis = basis[:,:2]
mean = fpca_discretized.mean_.data_matrix.reshape((-1,1))

# %%
# plt.subplot(2,1,1)
# plt.plot(basis)
# plt.title('Basis functions')
# plt.subplot(2,1,2)
# plt.plot(mean)
# plt.title('Mean')

# %% [markdown]
# # Get a sample trajectory

# %%
lcdts = pickle.load(open('lanechange/lcsp.pyo','rb'))

# %%
# sample_traj = xvec[16,:]
sample_traj = np.zeros_like(mean)
sample_traj[:24] = (lcdts[0].reshape((-1,1)) - mean)[0]
sample_traj[24:] = (lcdts[0].reshape((-1,1)) - mean)[:-24]
# sample_traj = xvec[16,:]
plt.plot(sample_traj)
plt.show()

# %% [markdown]
# # Load up the detector

# %%
from ASPDetector import ASPDetector

# %%
asp = ASPDetector(
    np.array([[0],[0]]),
    np.array([[30,0],[0,15]]),
    1e-2,
    mean,
    basis,
    75
)

# %%
# UNIT TEST 1: MARGINALIZATION MODULE ===

from scipy.stats import norm

if False:
    sample_traj = asp.kl_basis @ np.array([[40],[50]])
    llh_log = np.zeros_like(mean)
    x = []
    for t in range(mean.shape[0]):
        datum = sample_traj[t]
        x.append(datum)
        llh, S, m, _ = asp.marginalize_repr_a(np.array(x), asp.mu_a, asp.sigma_a, asp.sigma_e)
        llh_log[t] = llh
        # === debugging code ===
        if t == 25:
            print('m at {} is: {}'.format(t,m))
            # reconstructing the trajectory from m
            plt.plot(asp.kl_basis @ m, label='Bayesian')
            plt.plot(np.array(x), label='True data')
            plt.show()

    sys.exit(0)

# %%
# UNIT TEST 2: PREDICTION MODULE ===

if False:
    sample_traj = asp.kl_basis @ np.array([[40],[50]])
    llh_log = np.zeros_like(mean)
    x = []
    for t in range(mean.shape[0]):
        datum = sample_traj[t]
        x.append(datum)
        llh, S, m, _ = asp.marginalize_repr_a(np.array(x), asp.mu_a, asp.sigma_a, asp.sigma_e)
        llh_log[t] = llh
        # === debugging code ===
        if t == 5:
            print('m at {} is: {}'.format(t,m))
            print('Last datum: {}'.format(datum))
            # reconstructing the trajectory from m
            # now let's predict the value at t+1
            x_samples = norm.rvs(size=1000, loc=datum.item(), scale=1)
            llh_pred = np.zeros_like(x_samples)
            
            for i, sample in enumerate(x_samples):
                llhh, SS, mm, _ = asp.marginalize_repr_a(np.array([[sample]]), m, S, asp.sigma_e)
                llh_pred[i] = llhh.item()
            
            plt.hist(llh_pred, bins=100)
            plt.show()

    sys.exit(0)

# %%
# UNIT TEST 3: OBSERVATION LIKELIHOOD TEST ===
if False:
    sample_traj = asp.kl_basis @ np.array([[40],[50]])
    llh_log = np.zeros_like(mean)
    x = []
    for t in range(mean.shape[0]):
        datum = sample_traj[t]
        x.append(datum)
        llh, S, m, _ = asp.marginalize_repr_a(np.array(x), asp.mu_a, asp.sigma_a, asp.sigma_e)
        llh_log[t] = llh
        # === debugging code ===
        if t == 30:
            print('m at {} is: {}'.format(t,m))
            print('Previous datum: {}'.format(x[-2]))
            print('Latest datum: {}'.format(x[-1]))
            # reconstructing the trajectory from m
            # now let's predict the value at t+1
            asp.t = t
            asp.x = np.array(x)
            print('Likelihood of the latest datum: ', asp.p_prediction(t-1)) # use all the past data!

    sys.exit(0)

# %%
for t in range(50):
    asp.add_datum(sample_traj[t])
    print('Datum ', t, ' added')

plt.imshow(asp.prt_given_xt)# %%
plt.show()