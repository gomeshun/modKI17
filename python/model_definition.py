# Experimental code to demonstrate 'biased' KI17 analysis
# 
# 

import numpy as np
import pandas as pd
from scipy.stats import norm as gauss
from scipy.integrate import dblquad 
import MCgenerator

# definitions of true parameters
# r_e: typical scale of the member distribution
# s: contamination ratio
# v_mem, v_fg: mean velocity of the member/fg stars
# dv_fg: standard deviation of foreground stars
# a,b : parameters of dispersiom curve: \sigma_{los}(R) = a + b(R/r_e) 

initial_parameters_mem = {
    'r_e' : 100, 'v_mem' : 50, 'a' : 10,'b' : 10
    }
initial_parameters_fg = {
    'v_fg' : 200, 'dv_fg' : 200
    }
initial_parameters = {
    **initial_parameters_mem,**initial_parameters_fg,
    's' : 0.5
    }

RoI_R = 200
RoI_v_lo = 0
RoI_v_hi = 100

def dist_func_mem(v,R,r_e,v_mem,a,b):
    '''
    distribution function of the member stars.
        f(v,R) = prob_R(R) * prob_v_R(v,R)
            = 2\pi R * (Plummer model for R) * (Gaussian for v) /norm
        where
        norm = \int\int_{RoI}{dv}{dR} (2\pi R) * (Plummer model) * (Gaussian for v)
    '''
    prob_R = lambda R: 2*np.pi*R * np.power(1+(R/r_e)*(R/r_e),-2) # Plummer model, NOTE: no r_e^-2 but it is not required because the normalization factor also contain r_e^-2
    prob_v_R = lambda v,R: gauss.pdf(v,loc=v_mem,scale=(a+b*(R/r_e)))
    norm = quad(
        func = lambda R: (gauss.cdf(RoI_v_hi,loc=v_mem,scale=(a+b*(R/r_e)))-gauss.cdf(RoI_v_lo,loc=v_mem,scale=(a+b*(R/r_e))))*prob_R(R),
        a = 0, b = RoI_R
    )
    #norm,err_norm = dblquad(
    #    func = lambda v,R: prob_v_R(v,R)*prob_R(R), 
    #    a = 0, b = RoI_R, # for R
    #    gfun = lambda R: RoI_v_lo, hfun = lambda R: RoI_v_hi # for v
    #    )
    # return prob_v_R*prob_R for (v,R) in RoI and 0 (False) for not in RoI
    isin_RoI = np.prod((RoI_v_lo<v,v<RoI_v_hi,0<R,R<RoI_R),axis=0).astype(np.bool) # astype(np.bool) is required because fancy_index is available only for bool but (0,1)
#    return prob_v_R(v,R)*prob_R(R)/norm * isin_RoI # This may be return NaN because NaN * 0 = NaN
    if np.array(v).ndim == 0: # return scaler
        return prob_v_R(v,R)*prob_R(R)/norm if isin_RoI else 0
    else: # return array
        ret = np.zeros(v.shape)
        ret[isin_RoI] = (prob_v_R(v,R)*prob_R(R)/norm)[isin_RoI] 
        return ret

def dist_func_fg(v,R,v_fg,dv_fg):
    '''
    distribution function of the foreground stars.
        f(v,R) = prob_R(R) * prob_v_R(v,R)
            = 2\pi R * (constant) * (Gaussian for v) /norm
        where
        norm = \int\int_{RoI}{dv}{dR} (2\pi R) * (constant) * (Gaussian for v)
        We note that the constant factor is not required because it is divided
    '''
    prob_R = lambda R: 2*np.pi*R # constant model
    prob_v_R = lambda v,R: gauss.pdf(v,loc=v_fg,scale=dv_fg)
    norm = quad(
        func = lambda R: (gauss.cdf(RoI_v_hi,loc=v_fg,scale=dv_fg)-gauss.cdf(RoI_v_lo,loc=v_fg,scale=dv_fg))*prob_R(R),
        a = 0, b = RoI_R
    )
    #norm,err_norm = dblquad(
    #    func = lambda v,R: prob_v_R(v,R)*prob_R(R), 
    #    a = 0, b = RoI_R, # for R
    #    gfun = lambda R: RoI_v_lo, hfun = lambda R: RoI_v_hi # for v
    #    )
    # return prob_v_R*prob_R for (v,R) in RoI and 0 (False) for not in RoI
    isin_RoI = np.prod((RoI_v_lo<v,v<RoI_v_hi,0<R,R<RoI_R),axis=0).astype(np.bool)
#    return prob_v_R(v,R)*prob_R(R)/norm * isin_RoI # This may be return NaN because NaN * 0 = NaN
    if np.array(v).ndim == 0: # return scaler
        return prob_v_R(v,R)*prob_R(R)/norm if isin_RoI else 0
    else: # return array
        ret = np.zeros(v.shape)
        ret[isin_RoI] = (prob_v_R(v,R)*prob_R(R)/norm)[isin_RoI] 
        return ret

def dist_func_tot(v,R,r_e,v_mem,a,b,v_fg,dv_fg,s):
    return s*dist_func_mem(v,R,r_e,v_mem,a,b)+(1-s)*dist_func_fg(v,R,v_fg,dv_fg)

def generate_mem(r_e,v_mem,a,b,N_burnin=1000, N_sampling=1000, N_stride=5,**options):
    logfunc = lambda v,R: np.log(dist_func_mem(v,R,r_e,v_mem,a,b))
    gen = MCgenerator.MCgenerator(logfunc,args_logpdf_init={'v':25,'R':30})
    gen.update_MCparameter(dargs_logpdf={'v':15,'R':30},**options)
    # sampling
    print("burnin start"), gen.generate_tuned(N_burnin,repeat=5) # for burn in
    print("sampling start"), gen.generate(N_sampling*N_stride) # for sampling
    args_logpdf_chain = gen.args_logpdf_chain[-N_sampling*N_stride::N_stride] 
    return args_logpdf_chain

def generate_fg(v_fg,dv_fg,N_burnin=1000, N_sampling=1000, N_stride=5,**options):
    logfunc = lambda v,R: np.log(dist_func_fg(v,R,v_fg,dv_fg))
    gen = MCgenerator.MCgenerator(logfunc,args_logpdf_init={'v':25,'R':30})
    gen.update_MCparameter(dargs_logpdf={'v':25,'R':30},**options)
    # sampling
    print("bunnin start"),gen.generate_tuned(N_burnin,repeat=5) # for burn in
    print("sampling start"),gen.generate(N_sampling*N_stride) # for sampling
    args_logpdf_chain = gen.args_logpdf_chain[-N_sampling*N_stride::N_stride] 
    return args_logpdf_chain

if __name__ == "__main__":
    N_mem,N_fg = 20000,20000
    push_time = 5
    
    yn_create_mock = 'y' #input("create mocks? (y/n)")
    iscreate_mock = (yn_create_mock=='y')
    
    if iscreate_mock:
        print("create mock start")
        import MCgenerator
        import matplotlib.pyplot as plt
        print("mock_mem creation start.")
        mem = generate_mem(**initial_parameters_mem,N_sampling=N_mem,push_time=push_time)
        print("modk_fg createion start.")
        fg = generate_fg(**initial_parameters_fg,N_sampling=N_fg,push_time=push_time)
        print("mem,fg mocks has been created.")
        plt.hist2d(mem.v,mem.R)
        plt.hist2d(fg.v,fg.R)
        plt.xlabel("v")
        plt.ylabel("R")
        plt.show()
        yn_savemock = 'y' #input("save the created mock? (y/n)")
        issavemock = (yn_savemock=='y')
        if issavemock:
            mem.to_csv("mock_mem_v002.csv",index=None)
            fg.to_csv("mock_fg_v002.csv",index=None)
            print("mocks saved.")
    else:
        print("nooutput.")





