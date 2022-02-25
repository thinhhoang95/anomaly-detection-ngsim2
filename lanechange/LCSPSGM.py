from tokenize import Number
import numpy as np 
from scipy.stats import norm
from scipy.stats import expon
from BasisExtender import extend_basis
from LCSPS import LCSPS

class LCSPSGM(LCSPS): # inherits from LCSPS
    # Linear Combination of Stochastic Process Segmentation in Time Series (LC-SPS)
    def __init__(self, mu_a, sigma_a, sigma_e, kl_mean, kl_basis, kl_weight) -> None:
        super().__init__(mu_a, sigma_a, sigma_e, kl_mean, kl_basis)
        self.kl_weight = kl_weight
        self.mixture_components = kl_weight.shape[0]
        
    def marginalize_repr_all_subseqs(self, plot_x = False):
        # calculate the marginal distribution p(x_{ij}) of every single substrings i, j
        # setup the data vector
        x_vec = np.array(self.x)
        x_vec_len = np.shape(self.x)[0]
        # p will store all marginal proba of subsequence x_{i:j}
        p = np.zeros((x_vec_len, x_vec_len))
        # this pattern will avoid confusions about the way numpy array indices work
        i_from = 1
        i_to = x_vec_len
        for i in range(i_from - 1, i_to): # the first index can be anything from the first datum to the last
            j_from = i 
            j_to = x_vec_len
            for j in range(j_from, j_to):
                x = x_vec[i:j+1].reshape((-1,1)) - x_vec[i]
                for mixture_id in range(self.mixture_components):
                    mixture_weight = self.kl_weight[mixture_id] 
                    mixture_mean = self.mu_a[mixture_id,:] # associated mean of the p(ai)
                    mixture_cov = self.sigma_a[mixture_id,:] # associated covariance of the p(ai)
                    # print('---')
                    # print(mixture_mean)
                    # print(mixture_cov)
                # next we are going to marginalize this x vector against the kl_basis
                pm, _, __, ___, ____ = self.marginalize_repr_a(x, mixture_mean, mixture_cov, self.sigma_e, log_result=True, plot_x=plot_x)
                p[i,j] += mixture_weight * pm
                # p[i,j], _, __, ___, ____ = self.marginalize_repr_a(x, self.mu_a, self.sigma_a, self.sigma_e, log_result=True, plot_x=plot_x)
                # print(i,j,p[i,j])
        return p