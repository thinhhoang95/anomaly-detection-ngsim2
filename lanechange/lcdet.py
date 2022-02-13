from tokenize import Number
import numpy as np 
from scipy.stats import norm

class LCDetector:
    def __init__(self, mu_a, sigma_a, sigma_e, kl_mean, kl_basis) -> None:
        self.mu_a = mu_a 
        self.sigma_a = sigma_a 
        self.sigma_e = sigma_e 
        self.info_e = 1/sigma_e
        self.info_a = np.linalg.inv(sigma_a)
        self.kl_mean = kl_mean
        self.kl_basis = kl_basis
        self.kl_length = np.shape(kl_basis)[0]
        self.x = []
        self.t = -1
        print('Initialization of LCDetector')

    def add_datum(self, x) -> None:
        self.t += 1
        self.x.append(x)

    def p_representation_marginalized(self, t_from, t_to, include_t_to = False) -> Number:
        if include_t_to:
            x_tilde = np.array(self.x[t_from : t_to])
        else:
            x_tilde = np.array(self.x[t_from : t_to + 1])
        x_length = np.shape(x_tilde)[0]
        if x_length > self.kl_length:
            # We will have to construct a sliding window in this case
            max_window_offset = x_length - self.kl_length
            for offset in range(1, max_window_offset):
                if include_t_to:
                    x_tilde = np.array(self.x[t_from + offset : t_to + offset])
                else:
                    x_tilde = np.array(self.x[t_from + offset : t_to + offset + 1])
                # We solve the marginalization then picks the highest likelihood
        elif x_length == self.kl_length:
            F_tilde = self.kl_basis[:,:]
        else: # x_length < self.kl_length
            F_tilde = self.kl_basis[:x_length,:]
        S_mat = np.linalg.inv(self.info_a + F_tilde.T @ F_tilde * 1/(self.sigma_e**2))
        m_vec = S_mat @ (x_tilde.T @ F_tilde * 1/(self.sigma_e**2) + self.mu_a.T @ self.info_a).T
        