import numpy as np 

class OTQD:
    def __init__(self, info_a, mu_a, info_e2, pca_mean, pca_components, i_max = 1, i_spacing = 10):
        # Hypothesis of nominal traffic: a Gaussian mixture
        self.info_a = info_a.copy()
        self.mu_a = mu_a.copy()
        # Measurement noise covariance
        self.info_e2 = info_e2
        # FPCA components
        self.pca_mean = pca_mean
        self.pca_components = pca_components
        # Last observation time
        self.t = 0

        # Hypothesis on when the observation started
        self.i_max = i_max
        self.i_spacing = i_spacing

        # State variables
        self.mu_pre_prime = np.zeros((i_max,2,1)) 
        self.info_prime = np.zeros((i_max,2,2))
        self.preresidual = np.zeros((i_max))
        
        # Initialize state variables
        for i in range(i_max):
            self.mu_pre_prime[i,:,:] = (info_a @ mu_a).reshape((-1,1))
            self.info_prime[i,:,:] = info_a
            self.preresidual[i] = mu_a.transpose() @ info_a @ mu_a 

    def get_pca_at_t(self, t):
        return (self.pca_mean[t], self.pca_components[t,:])
    
    def new_measurement(self, x):
        for i in range(self.i_max): # i: delay for the PCA values
            mean_t, pca_components_t = self.get_pca_at_t(self.t + i * self.i_spacing)
            pca_components_t = pca_components_t.reshape((-1,1))
            CX = self.info_e2 * np.matmul(pca_components_t, pca_components_t.transpose()) # should yield a n x n matrix, NOT a scalar
            CD = self.info_e2 * (-mean_t + x) * pca_components_t
            iota_t = self.info_e2 * np.square(x - mean_t)

            # Update state vector
            self.info_prime[i,:,:] += CX 
            self.mu_pre_prime[i,:,:] += CD
            self.preresidual[i] += iota_t 
        self.t += 1

    def calculate_likelihood(self):
        covar_bare = np.linalg.inv(self.info_a)
        # print('DBG: Covar bare: ', covar_bare)
        likelihood = np.zeros((self.i_max,))
        likelihood_w = np.zeros((self.i_max,))
        for i in range(self.i_max):
            print('Delay i = ', i)
            covar_prime = np.linalg.inv(self.info_prime[i,:,:])
            print('DBG: Covar prime ', covar_prime)
            mu_prime = covar_prime @ self.mu_pre_prime[i,:,:]
            print('mu pprime', self.mu_pre_prime[i,:,:])
            print('mu prime ', mu_prime)
            residual = self.preresidual[i] - mu_prime.transpose() @ self.info_prime[i,:,:] @ mu_prime
            print('DBG: Preresidual ', self.preresidual[i])
            print('DBG: Preresidual minus ', mu_prime.transpose() @ self.info_prime[i,:,:] @ mu_prime)
            # print('DBG: Residual ', residual)
            likelihood[i] = np.power(self.info_e2 * np.sqrt(2*np.pi), self.t / 2.) * np.sqrt(np.linalg.det(covar_prime)/np.linalg.det(covar_bare)) * np.exp(-.5 * residual) 
            likelihood_w[i] = self.preresidual[i]
            print('---')
        return likelihood, likelihood_w

    def calculate_log_likelihood(self):
        covar_bare = np.linalg.inv(self.info_a)
        # print('DBG: Covar bare: ', covar_bare)
        likelihood = np.zeros((self.i_max,))
        likelihood_w = np.zeros((self.i_max,))
        for i in range(self.i_max):
            #print('Delay i = ', i)
            covar_prime = np.linalg.inv(self.info_prime[i,:,:])
            #print('DBG: Covar prime ', covar_prime)
            mu_prime = covar_prime @ self.mu_pre_prime[i,:,:]
            #print('mu pprime', self.mu_pre_prime[i,:,:])
            #print('mu prime ', mu_prime)
            residual = self.preresidual[i] - mu_prime.transpose() @ self.info_prime[i,:,:] @ mu_prime
            #print('DBG: Preresidual ', self.preresidual[i])
            #print('DBG: Preresidual minus ', mu_prime.transpose() @ self.info_prime[i,:,:] @ mu_prime)
            # print('DBG: Residual ', residual)
            likelihood[i] = -self.t * np.log(self.info_e2) + .5 * np.log(np.linalg.det(covar_prime)/np.linalg.det(covar_bare)) -.5 * residual
            likelihood_w[i] = self.preresidual[i]
            # print('---')
        return likelihood