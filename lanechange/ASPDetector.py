from tokenize import Number
import numpy as np 
from scipy.stats import norm
from scipy.stats import expon

class ASPDetector:
    def __init__(self, mu_a, sigma_a, sigma_e, kl_mean, kl_basis, h_lambda) -> None:
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
        self.h_lambda = h_lambda

        self.prtxt = np.zeros((55,55)) # should be large enough to contain the whole time series!
        self.prtxt[0,0] = 1 # run = 0, at no data, proba is 1
        self.prt_given_xt = np.zeros((55,55)) 
        print('Initialization of ASPDetector completed')

    def marginalize_repr_a(self, x, mu_a, sigma_a, sigma_e, kl_t = 0, kl_length = 0) -> Number:
        # note that x is subtracted by mean_curve already
        info_a = np.linalg.inv(sigma_a)
        info_e_squared = np.square(1/sigma_e)
        x_length = np.shape(x)[0]
        # the x vector is longer than the basis temporal size, we just take the last kl_length entries of the x_vector:
        if x_length > self.kl_length:
            x = x[-self.kl_length:]
        if kl_t > 0:
            # we do not start from the "beginning" of the KL basis, but rather from kl_t
            F = self.kl_basis[kl_t:kl_t+kl_length,:]
        else:
            F = self.kl_basis[0:x_length,:] # so now it will be like [140, 20]
        # new covariance matrix
        info_S = info_a + F.T @ F * info_e_squared
        S = np.linalg.inv(info_S)
        # new mean value
        m = S @ (x.T @ F * info_e_squared + mu_a.T @ info_a).T
        # and the residual is
        r = info_e_squared * x.T @ x + mu_a.T @ info_a @ mu_a - m.T @ info_S @ m
        llh = (1/(np.sqrt(2 * np.pi) * sigma_e) ** x_length) * np.sqrt(np.linalg.det(S) / np.linalg.det(sigma_a)) * np.exp(-r/2)
        kl_next = x_length
        return llh, S, m, r, kl_next

    def hazard(self, t) -> Number: # use Planck (i.e., discrete exponential distribution with parameter lambda as interval prior)
        return expon.pdf(t, scale=1/self.h_lambda) / (1 - expon.pdf(t, scale=1/self.h_lambda))

    def p_rt_transit(self, rtm1, rt):
        if rt == 0:
            return self.hazard(rtm1 + 1)
        elif rt == rtm1 + 1:
            return 1 - self.hazard(rtm1 + 1)
        else:
            return 0

    def p_prediction(self, rtm1):
        x_infl_lbi = min(rtm1, self.kl_length) # x_influential_lower_bound_index
        # if x_infl_lbi == 0:
        #     # use p(a) instead of p(a|x_infl), in other words, S and m are the same as mu_a and sigma_a
        #     x_now = np.array(self.x[-1:]) # latest value in the array 
        #     llh, SS, mm, __ = self.marginalize_repr_a(x_now, self.mu_a, self.sigma_a, self.sigma_e)
        #     print('No influential history1!')
        #     return llh
        # else:
        llh_candidates = np.zeros((x_infl_lbi + 1,))
        # print("LBI: ", x_infl_lbi)
        for l in range(x_infl_lbi + 1):
            x_infl = np.array(self.x[-l-2:-1]) # the "previous" value is x[-2], not x[-1]!
            # print("SP's beginning hypothesis: last {} entries > Influential history: {}".format(l, x_infl))
            if x_infl.shape[0] == 0: # no influential past history
                S = self.sigma_a
                m = self.mu_a
                # print('No influential history!')
            else:
                _, S, m, __, kl_next = self.marginalize_repr_a(x_infl, self.mu_a, self.sigma_a, self.sigma_e) # we just need S and m to plug in to the formula on page 6
            # print("Mean after influential history fusion: {}".format(m))
            x_now = np.array(self.x[-1]) # latest value in the array 
            llh, SS, mm, __, ___ = self.marginalize_repr_a(x_now, m, S, self.sigma_e, kl_t = kl_next, kl_length = 1)
            # print("Datum likelihood: {}".format(llh))
            # print("-")
            llh_candidates[l] = llh
        result = np.max(llh_candidates)
        if result > 0:
            # print("Best Datum Likelihood: {}".format(result))
            pass
        return result

    def add_datum(self, x) -> None:
        self.t += 1
        self.x.append(x)
        
        for rt in range(self.t+1): # the current run (run at time t) could either continue or restart from zero
        # hence it could take any value from 0 to t
            if rt == 0:
                if self.t == 0:
                    p_rtxt = 1 # basically this only executes when rt=0, t=0
                else: # self.t >= 1
                    p_rtxt = 0
                    for rtm1 in range(self.t): # if rt=0 then rtm1 can be anything from 0 to t-1
                        p_rt = self.p_rt_transit(rtm1, rt)
                        p_xt = self.p_prediction(rtm1)
                        p_rtxt_prev = self.prtxt[rtm1, self.t-1] # correspond of x: self.t-1 because rtm1 can be from 0 to t-1
                        p_rtxt += p_rt * p_xt * p_rtxt_prev
                # print('rt=0')
            else: # rt>=1 so rtm1 must be rt-1, no other values are possible
                rtm1 = rt-1
                if self.t==30 and rtm1==29:
                    print('Breakpoint')
                p_rt = self.p_rt_transit(rtm1, rt)
                p_xt = self.p_prediction(rtm1)
                p_rtxt_prev = self.prtxt[rtm1, self.t-1]
                p_rtxt = p_rt * p_xt * p_rtxt_prev
                # print('rt={}'.format(rt))
            
            if np.isnan(p_rtxt):
                raise Exception("NaN in calculation")

            if p_rtxt < 1e-60:
                p_rtxt = 1e-60 # for numerical stability

            self.prtxt[rt, self.t] = p_rtxt

        for rt in range(self.t+1):
            self.prt_given_xt[rt, self.t] = self.p_rt_given_xt(rt)
            if np.isnan(self.p_rt_given_xt(rt)):
                raise Exception("NaN detected")


    def p_xt(self):
        return np.sum(self.prtxt[:, self.t])

    def p_rt_given_xt(self, rt):
        if np.isnan(self.p_xt()):
            raise Exception("NaN in calculation")
        return self.prtxt[rt, self.t] / self.p_xt()
