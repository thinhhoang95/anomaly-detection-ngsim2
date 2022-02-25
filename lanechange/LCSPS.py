from tokenize import Number
import numpy as np 
from scipy.stats import norm
from scipy.stats import expon
from BasisExtender import extend_basis
from matplotlib import pyplot as plt

class LCSPS:
    # Linear Combination of Stochastic Process Segmentation in Time Series (LC-SPS)
    def __init__(self, mu_a, sigma_a, sigma_e, kl_mean, kl_basis) -> None:
        self.mu_a = mu_a 
        self.sigma_a = sigma_a 
        self.sigma_e = sigma_e 
        self.info_e = 1/sigma_e
        self.info_a = np.linalg.inv(sigma_a)
        self.kl_mean = kl_mean
        self.kl_basis = kl_basis # should be like [20,140] where 20 is the number of components, 140 are timesteps
        self.kl_length = np.shape(kl_basis)[0]
        self.x = []
        self.t = -1

        self.prtxt = np.zeros((55,55)) # should be large enough to contain the whole time series!
        self.prtxt[0,0] = 1 # run = 0, at no data, proba is 1
        self.prt_given_xt = np.zeros((55,55)) 
        print('Initialization of LCSPS completed')

    def add_datum(self, x):
        self.x.append(x)
        self.t += 1

    def marginalize_repr_a(self, x, mu_a, sigma_a, sigma_e, kl_t = 0, kl_length = 0, log_result = False, plot_x = False) -> Number:
        # note that x is subtracted by mean_curve already
        info_a = np.linalg.inv(sigma_a)
        info_e_squared = np.square(1/sigma_e)
        x_length = np.shape(x)[0]
        # the x vector is longer than the basis temporal size, we must extend the basis to match the x vector
        if x_length > self.kl_length:
            F = extend_basis(self.kl_basis, x_length)
            M = extend_basis(self.kl_mean, x_length)
        else: # x_length is <= basis_length, we just trim the basis!
            if kl_t > 0:
                # we do not start from the "beginning" of the KL basis, but rather from kl_t
                F = self.kl_basis[kl_t:kl_t+kl_length,:]
                M = self.kl_mean[kl_t:kl_t+kl_length,:]
            else:
                F = self.kl_basis[0:x_length,:] # so now it will be like [140, 20]
                M = self.kl_mean[0:x_length]
        x = x - M # removing the mean
        if plot_x:
            plt.figure()
            plt.plot(x)
        # new covariance matrix
        info_S = info_a + F.T @ F * info_e_squared
        S = np.linalg.inv(info_S)
        # new mean value
        m = S @ (x.T @ F * info_e_squared + mu_a.T @ info_a).T
        # and the residual is
        r = info_e_squared * x.T @ x + mu_a.T @ info_a @ mu_a - m.T @ info_S @ m
        kl_next = x_length
        if not log_result: # log_result = False
            llh = (1/(np.sqrt(2 * np.pi) * sigma_e) ** x_length) * np.sqrt(np.linalg.det(S) / np.linalg.det(sigma_a)) * np.exp(-r/2)
            return llh, S, m, r, kl_next
        else: # log_result = True
        # loglikelihood for better numerical stability?
            logllh = -x_length * np.log(np.sqrt(2 * np.pi) * sigma_e) + np.log(np.sqrt(np.linalg.det(S) / np.linalg.det(sigma_a))) - r/2
            return logllh, S, m, r, kl_next

    def hazard(self, t) -> Number: # use Planck (i.e., discrete exponential distribution with parameter lambda as interval prior)
        return expon.pdf(t, scale=1/self.h_lambda) / (1 - expon.pdf(t, scale=1/self.h_lambda))

    def marginalize_repr_all_subseqs(self, plot_x=False):
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
                x = x_vec[i:j+1].reshape((-1,1))
                # next we are going to marginalize this x vector against the kl_basis
                p[i,j], _, __, ___, ____ = self.marginalize_repr_a(x, self.mu_a, self.sigma_a, self.sigma_e, log_result=True, plot_x=plot_x)
                # print(i,j,p[i,j])
        return p

    def marginalize_repr_subseq(self, i, j, plot_x=False):
        x_vec = np.array(self.x)
        x_vec_len = np.shape(self.x)[0]
        x = x_vec[i:j+1].reshape((-1,1)) - x_vec[i]
        # next we are going to marginalize this x vector against the kl_basis
        p, _, __, ___, ____ = self.marginalize_repr_a(x, self.mu_a, self.sigma_a, self.sigma_e, log_result=True, plot_x=plot_x)
        # print(i,j,p[i,j])
        return p

    def get_change_point(self, p, K = 1):
        x_vec_len = np.shape(self.x)[0]
        # we sequentially find all changepoints, limited by the upper variable K 
        ck = np.zeros((K+2,)) # all changepoints are 1-based, not zero-based
        ck_llh = np.zeros((K+2,))
        ck[0] = -1
        ck[-1] = self.t - 1
        for k in range(1,K+1): # k from 1 to K
            # although written here as k we find the index for changepoint k+1 (since k starts from 0)
            hypotheses = []
            hypotheses_indices = []
            hypotheses_t1 = []
            hypotheses_t2 = []
            for i in range(int(ck[k-1]+1), x_vec_len - 1): # the changepoint cannot coincide with the last datum
                first_seq_id_from = int(ck[k-1]+1)
                first_seq_id_to = i # inclusive
                second_seq_id_from = i+1 
                second_seq_id_to = x_vec_len-1 # inclusive
                hypotheses_t1.append([first_seq_id_from, first_seq_id_to])
                hypotheses_t2.append([second_seq_id_from, second_seq_id_to])
                p_first_seq = p[first_seq_id_from, first_seq_id_to]
                p_second_seq = p[second_seq_id_from, second_seq_id_to]
                hypotheses.append(p_first_seq + p_second_seq)
                hypotheses_indices.append(i)

                # first_seq = self.x[first_seq_id_from:first_seq_id_to+1] #inclusive
                # second_seq = self.x[second_seq_id_from:second_seq_id_to+1]
            arghypot = np.argmax(np.array(hypotheses))

            ck[k] = hypotheses_indices[arghypot]
            ck_llh[k] = hypotheses[arghypot]
            first_t = hypotheses_t1[arghypot]
            second_t = hypotheses_t2[arghypot]
            print('Changepoint placement k={}: {}'.format(k, first_t))
            print('Changepoint placement k={}: {}'.format(k, second_t))
        return ck, ck_llh, first_t, second_t