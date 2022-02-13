import sys
from numbers import Number
from typing import Dict
sys.path.append('/Users/thinhhoang/Documents/anomaly-detection-ngsim/otqd')
from otqd import OTQD as BayesFactorizer
import numpy as np
from scipy.stats import geom

class OnlineBayesian:
    def __init__(self, bayesFactorizerParams: Dict) -> None:
        self.bayesFactorizerParams = bayesFactorizerParams
        self.info_e2 = bayesFactorizerParams['info_e2']
        self.prx = [] # to store the p(rt, x_{1:t}) variables
        self.x_vec = [] # to store the history of observations
        self.t = 0
    
    def new_datum(self, x) -> None:
        self.x_vec.append(x)
    
    def calculate_upm(self, r, t) -> Number: # will calculate UPM up to the entry self.x_vec[t]
        self.bayesFactorizer = BayesFactorizer(**self.bayesFactorizerParams)
        if r>t:
            raise Exception("Why should run r = {} > t = {}?".format(r,t))
        # This method will calculate the UPM distribution with run length r given
        x_vec_of_this_run = self.x_vec[t-r+1:t+1] # get r values, ending with x_vec[t] (not x_vec[t+1])
        if r == 0:
            return 1 # initial conditions for UPM, this is the simple model
            # denominator_llh = 1.
            # nominator_llh = 1
        else:
            if r==1:
                denominator_llh = 1
            for i in range(r-1):
                # we gradually "feed" the BayesFactorizer with values from x_vec_of_this_run, up
                # to the next-to-last datum (to retrieve the denominator first) 
                self.bayesFactorizer.new_measurement(x_vec_of_this_run[i])
            log_llh, cov_llh = self.bayesFactorizer.calculate_log_likelihood_with_covar()
            denominator_llh = np.exp(log_llh)
        self.bayesFactorizer.new_measurement(x_vec_of_this_run[-1])
        log_llh, cov_llh = self.bayesFactorizer.calculate_log_likelihood_with_covar()
        nominator_llh = np.exp(log_llh)
        return 1/(np.sqrt(2 * np.pi / self.info_e2)) * nominator_llh / denominator_llh

    def hazard(self, n, distrib = 'geometric') -> Number:
        if distrib == 'geometric':
            return geom.pmf(n, p=0.5)/(1-geom.cdf(n, p=0.5))
        else:
            raise Exception('Distribution not supported')

    def calculate_prx(self, t, l):
        if l>0:
            prxm1 = self.calculate_prx(t=t-1, l=l-1)
            return prxm1 * self.calculate_upm(r=l-1,t=t-1) * (1-self.hazard(l-1))
        elif l==0:
            result = 0
            for rtm1 in range(t): # rtm1 from 0 to t-1
                prxm1 = self.calculate_prx(t-1, rtm1)
                result += prxm1 * self.calculate_upm(r=rtm1, t=t-1) * self.hazard(rtm1)
            return result

    def calculate_normalization_term(self) -> Number:
        result = 0
        for rt in range(self.t):
            result += self.calculate_prx(self.t, rt)
        return result

    def calculate_rl_posterior(self, rt) -> Number:
        return self.calculate_prx(self.t, rt) / self.calculate_normalization_term()

    def shift_time(self):
        self.t += 1 # should be called last
        print('Time shift: {:d}'.format(self.t))

    # I particularly find UPM formula extremely confusing, it's unclear what l refers to

    # def calculate_prx(self, t, rt, rtm1 = 0) -> Number:
    #     # params >> rtm1: r_{t-1}, prxm1: p(r_{t-1}, x_{1:t-1})
    #     # return >> prx: p(r_t, x{1:t})
    #     if (t==1):
    #         # This is the beginning of the time series (NOT the beginning of the run)
    #         # but of course as this is our first measurement, the run length is also zero 
    #         # note that the vice versa is not neccessarily true!
    #         return 1. * self.calculate_upm(0) * (1.-self.hazard(1))
    #     else:
    #         # t>1, this is NOT our first measurement
    #         if rt == rtm1 + 1:
    #             prxm1 = self.prx[-1]
    #             result = prxm1 * self.calculate_upm(rtm1) * (1.-self.hazard(rtm1))
    #             self.prx.append(result)
    #             return result 
    #         if rt == 0:
    #             result = 0
    #             for i in range(t):
    #                 result += self.prx[i] * self.calculate_upm(i) * self.hazard(i)
    #             return result
    #         else:
    #             return 0 # because from r_{t-1} there are only 2 cases: r_t = r_{t-1} + 1 or r_t = 0
    #         else:


