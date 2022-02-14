from tokenize import Number
import numpy as np 
from scipy.stats import norm
from scipy.stats import planck

class ASPDetector:
    def __init__(self, mu_a, sigma_a, sigma_e, kl_mean, kl_basis, h_lambda) -> None:
        self.mu_a = mu_a 
        self.sigma_a = sigma_a 
        self.sigma_e = sigma_e 
        self.info_e = 1/sigma_e
        self.info_a = np.linalg.inv(sigma_a)
        self.kl_mean = kl_mean
        self.kl_basis = kl_basis # should be like [20,140] where 20 is the number of components, 140 are timesteps
        self.kl_length = np.shape(kl_basis)[1]
        self.x = []
        self.t = -1
        self.h_lambda = h_lambda

        self.prtxt = np.zeros((1000,1000)) # should be large enough to contain the whole time series!
        print('Initialization of ASPDetector completed')

    def marginalize_repr_a(self, x, mu_a, sigma_a, sigma_e) -> Number:
        # note that x is subtracted by mean_curve already
        info_a = np.linalg.inv(sigma_a)
        info_e_squared = np.square(1/sigma_e)
        x_length = np.shape(x)[0]
        # the x vector is longer than the basis temporal size, we just take the last kl_length entries of the x_vector:
        if x_length > self.kl_length:
            x = x[-self.kl_length:]
        F = self.kl_basis[:,0:x_length].T # so now it will be like [140, 20]
        # new covariance matrix
        S = np.linalg.inv(info_a + F.T @ F * info_e_squared)
        info_S = np.linalg.inv(S)
        # new meann value
        m = S @ (x.T @ F * info_e_squared + mu_a.T @ info_a).T 
        # and the residual is
        r = info_e_squared * x.T @ x + mu_a.T @ info_a @ mu_a - m.T @ info_S @ m
        llh = (1/(np.sqrt(2 * np.pi) * sigma_e) ** x_length) * np.sqrt(np.linalg.det(S) / np.linalg.det(sigma_a)) * np.exp(-r/2)
        return llh, S, m, r

    def hazard(self, t) -> Number: # use Planck (i.e., discrete exponential distribution with parameter lambda as interval prior)
        return planck.pmf(t, self.h_lambda) / (1 - planck.cdf(t, self.h_lambda))

    def p_rt_transit(self, rtm1, rt):
        if rt == 0:
            return self.hazard(rtm1 + 1)
        elif rt == rtm1 + 1:
            return 1 - self.hazard(rtm1 + 1)
        else:
            return 0

    def p_prediction(self, rtm1):
        x_infl_lbi = min(rtm1, self.kl_length) # x_influential_lower_bound_index
        if x_infl_lbi == 0:
            # use p(a) instead of p(a|x_infl), in other words, S and m are the same as mu_a and sigma_a
            x_now = self.x[-1:] # latest value in the array 
            llh, SS, mm, __ = self.marginalize_repr_a(x_now, self.mu_a, self.sigma_a, self.sigma_e)
            return llh
        else:
            llh_candidates = np.zeros((x_infl_lbi,))
            for l in range(x_infl_lbi):
                x_infl = self.x[-l-1:-1] # the "previous" value is x[-2], not x[-1]!
                _, S, m, __ = self.marginalize_repr_a(x_infl, self.mu_a, self.sigma_a, self.sigma_e) # we just need S and m to plug in to the formula on page 6
                x_now = self.x[-1:] # latest value in the array 
                llh, SS, mm, __ = self.marginalize_repr_a(x_now, m, S, self.sigma_e)
                llh_candidates[l] = llh
            return np.max(llh_candidates)

    def add_datum(self, x) -> None:
        self.t += 1
        self.x.append(x)

        for rtm1 in range(self.t): # note that self.t is excluded in this loop
            p_xt = self.p_prediction(rtm1)
            # case: rt = rtm1 + 1
            p_rt_equals_rtm1_plus_1 = self.p_rt_transit(rtm1, rtm1+1)
            # case: rt = 0
            p_rt_equals_zero = self.p_rt_transit(rtm1, 0)
        
        for rt in range(self.t+1):
            if rt == 0:
                for rtm1 in range(self.t): # if rt=0 then rtm1 can be anything from 0 to t-1
                    pass
            else: # rt>0
                rtm1 = rt-1
                p_rt = self.p_rt_transit(rtm1, rt)
                p_xt = self.p_prediction(rtm1)
                if self.t==0:
                    p_rtxt_prev = 1 # basically p(r0) = 1 when r0=0 and 0 everywhere else
                else:
                    p_rtxt_prev = self.prtxt[self.t-1, rtm1]
                
