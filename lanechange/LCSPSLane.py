from tokenize import Number
import numpy as np 
from scipy.stats import norm
from scipy.stats import expon
from BasisExtender import extend_basis
from LCSPS import LCSPS

class LCSPSLane(LCSPS): # inherits from LCSPS
    # Linear Combination of Stochastic Process Segmentation in Time Series (LC-SPS)
    def __init__(self, mu_a, sigma_a, sigma_e, kl_mean, kl_basis, lanes) -> None:
        super().__init__(mu_a, sigma_a, sigma_e, kl_mean, kl_basis)
        self.lanes = lanes
        
    def marginalize_repr_all_subseqs(self, plot_x = False):
        # calculate the marginal distribution p(x_{ij}) of every single substrings i, j
        # setup the data vector
        x_vec = np.array(self.x)
        x_vec_len = np.shape(self.x)[0]
        # p will store all marginal proba of subsequence x_{i:j}
        p = np.zeros((len(self.lanes), x_vec_len, x_vec_len))
        # this pattern will avoid confusions about the way numpy array indices work
        i_from = 1
        i_to = x_vec_len
        for lane, lane_start in enumerate(self.lanes):
            for i in range(i_from - 1, i_to): # the first index can be anything from the first datum to the last
                j_from = i 
                j_to = x_vec_len
                for j in range(j_from, j_to):
                    x = x_vec[i:j+1].reshape((-1,1)) - lane_start
                    # next we are going to marginalize this x vector against the kl_basis
                    p[lane, i,j], _, __, ___, ____ = self.marginalize_repr_a(x, self.mu_a, self.sigma_a, self.sigma_e, log_result=True, plot_x=plot_x)
                    # print(i,j,p[i,j])
        return np.amax(p, axis=0)