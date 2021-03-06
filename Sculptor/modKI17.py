# version 1.0.1
# 
# update of version 1.0.1
# - Use SkyCoord.radial_velocity_err as the err of the input array.
import MCgenerator, dSph_model, coord
from numpy import power,sqrt
from scipy.special import logsumexp
from scipy.stats import norm
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

#update likelihood
DEBUG = False

dSph_property = pd.read_csv("dSph_property.csv",index_col=0)
prior_lim = []
draco_prop = dSph_property.loc["Draco"]
sculptor_prop = dSph_property.loc["Sculptor"]
#RA0 = draco_prop.RAdeg
#DE0 = draco_prop.DEdeg
#DIST = draco_prop.DIST
#err_DIST = draco_prop.err_DIST
#dSph_property

def is_positive(*args):
    return np.array(args)>0

class modKI17:
    def __init__(self,sc_obsdata,dsph_name,paramlims_fname,prior_norm_fname):
        """
        sc_obsdata: SkyCoord of observed data
        sc_center: SkyCoord of ad hoc center of dSph
        paramlims_fname: filename of prior configuration
        """
        # limit of parameter
        self.prior_lim = pd.read_csv(paramlims_fname,index_col=0)
        self.prior_lim_param_names = self.prior_lim.index
        
        # prior_norm of parameter
        _df_prior_norm = pd.read_csv(prior_norm_fname,index_col=0)
        self.prior_norm = norm(loc=_df_prior_norm["loc"],scale=_df_prior_norm["scale"])
        self.prior_norm_param_names = _df_prior_norm.index
        
        # index for the likelihood, prior and posterior
        #self.param_names = ["re_pc","odds","dra0","dde0",
        #    "log10_rs_pc","log10_rhos_Msunpc3","a","b","g",
        #    "mlog10_1manib",
        #    "vmem","vfg0","vfg1","dvfg0","dvfg1", "sfg0","dist"]
        self.param_names = self.prior_lim.index
        
        #if not np.any(self.param_names == self.prior_lim.index):
        #    raise TypeError("parameter order is not matched!")
        
        #for val in vs,dRAs,dDEs:
        #    if isinstance(val,pd.DataFrame):
        #        val = val.values
        self.sc_obsdata = sc_obsdata
        if hasattr(self.sc_obsdata,"radial_velocity_err"):
            print("sc_obsdata has radial_velocity_err. We use the likelihood function with velocity error:\n{}".format(self.sc_obsdata.radial_velocity_err))
        #self.sc_center0 = sc_center0
        
        # read property
        RA0,DE0 = dSph_property.loc[dsph_name][["RAdeg","DEdeg"]]
        DIST,err_DIST = dSph_property.loc[dsph_name][["DIST","err_DIST"]]
        self.sc_center0 = SkyCoord(ra=RA0*u.deg,dec=DE0*u.deg,distance=DIST*u.pc)
        self.sc_center0.distance_err = err_DIST
        self.RoI_R = np.max(self.Rs(0,0,DIST)) # use Rs.max as the RoI
        
        self.beta = 1/np.log(len(self.sc_obsdata))
        
        #print(Rs.describe())
        #print("beta: {}".format(self.beta)) if beta != 1 else None
        mem = dSph_model.plummer_model(re_pc=200)
        dm = dSph_model.NFW_model(
            a=2.78,b=7.78,g=0.675,
            rhos_Msunpc3=np.power(10,-2.05),rs_pc=np.power(10,3.96),
            R_trunc_pc=2000
        )
        self.dsph = dSph_model.dSph_model(anib=-0.5,submodels_dict={"stellar_model":mem,"DM_model":dm},show_init=True)
        
        self.fg = dSph_model.uniform2d_model(Rmax_pc=self.RoI_R,show_init=True)
    
    @staticmethod
    def params_to_series(index,**kwargs):
        """
        return the series whoes elements are "kwargs" ordered by "index".
        it is useful to avoid to mistake the order of the parameters
        """
        sr = pd.Series(index=index)
        for key in kwargs:
            sr[key] = kwargs[key]
            #display(sr)
        display("params_to_ser",sr) if DEBUG else None
        return sr
    
    @staticmethod
    def params_to_array(index,**kwargs):
        """
        return the series whoes elements are "kwargs" ordered by "index".
        it is useful to avoid to mistake the order of the parameters
        """
        ret = modKI17.params_to_series(index=index,**kwargs).values
        display("params_to_arr",ret) if DEBUG else None
        return ret
    
    @staticmethod
    def array_to_series(index,params):
        ret = pd.Series(params,index=index)
        display("array_to_ser",ret) if DEBUG else None
        return ret
    
    def sc_center(self,dra0,dde0,dist):
        return SkyCoord(
            ra=self.sc_center0.ra+dra0*u.deg,
            dec=self.sc_center0.dec+dde0*u.deg,
            distance=dist*u.pc)
        
    def Rs(self,dra0,dde0,dist):
        c = self.sc_center(dra0,dde0,dist)
        return c.distance.value*np.sin(self.sc_obsdata.separation(c).rad)

    def lnprior(self,params):
        if not self.is_parameters_in_domain(params):
            display("out ouf dom:",self.is_parameters_in_domain(params)) if DEBUG else None
            return -np.inf
        else:
            #note that the order of the 
            # _prior_param_names = ["re_pc","odds","dra0","dde0","dist"]
            sr_params = self.array_to_series(index=self.param_names,params=params)
            #args_prior = self.params_to_array(index=self.prior_norm_param_names,
            #                re_pc=sr_params.re_pc, odds=sr_params.odds,
            #                dra0=sr_params.dra0, dde0=sr_params.dde0,
            #                dist=sr_params.dist
            #               )
            args_prior = np.array([sr_params[p_name] for p_name in self.prior_norm_param_names])
            logGs = self.prior_norm.logpdf(args_prior)
            #logGs = []
            #logGs.append(norm.logpdf(re,  loc=191,      scale=5.7)     )
            #logGs.append(norm.logpdf(odds,loc=8.794,    scale=0.5107)  )
            #logGs.append(norm.logpdf(dra0,loc=4.212e-3, scale=7.052e-3))
            #logGs.append(norm.logpdf(dde0,loc=-1.991e-3,scale=3.302e-3))
            #logGs.append(norm.logpdf(dist,loc=self.sc_center0.distance,scale=self.sc_center0.distance_err))
            return np.sum(logGs)

    def __call__(self,**args):
        params = self.params_to_array(index=self.param_names,**args)
        ret = self.lnprior(params)
        if ret == -np.inf:
            return ret
        else:
            return ret + np.sum(self.lnlikeli(params)) 
        
    def lnprob(self,params):
        sr_params = self.array_to_series(index=self.param_names,params=params)
        return self.__call__(**sr_params)
    
    def lnposterior_general(self,params):
        lnp = self.lnprior(params)
        if lnp > -np.inf:
            lnl = self.lnlikeli(params) 
            return (self.beta*lnl+lnp, lnl)
        else:
            return (-np.inf, np.nan)
    
    def is_parameters_in_domain(self,params):
        """
        check the parameters are in valid domain.
        note that the params is ordered.
        """
        sr_params = self.array_to_series(index=self.param_names,params=params)
        prms_df = pd.DataFrame({**self.prior_lim,"inprms":sr_params})
        display(prms_df) if DEBUG else None
        is_in_minmax = np.all((prms_df.prms_min < prms_df.inprms) & (prms_df.inprms < prms_df.prms_max))
        #is_ordered = prms_df.inprms.vfg0 < prms_df.inprms.vfg1
        if hasattr(prms_df,"sfg1"):
            sfg0,sfg1 = prms_df.inprms.sfg0,prms_df.inprms.sfg1
            sfg2 = 1-sfg0-sfg1
            is_ordered = (sfg0 > sfg1) & (sfg1 > sfg2) & (sfg2 > 0)
        elif hasattr(prms_df,"sfg0"):
            sfg0 = prms_df.inprms.sfg0
            sfg1 = 1-sfg0
            is_ordered = (sfg0 > sfg1) & (sfg1 > 0)
        else:
            is_ordered = True
        display("is_in_minmax",(prms_df.prms_min < prms_df.inprms) & (prms_df.inprms < prms_df.prms_max),"is_ordered",is_ordered) if DEBUG else None
        return is_in_minmax & is_ordered
    """
    def is_parameters_in_domain(self,re,odds,dra0,dde0,
            log10_rs_pc,log10_rhos_Msunpc3,a,b,g,
            mlog10_1manib,
            vmem,vfg0,vfg1,dvfg0,dvfg1, sfg0):
        is_positive_ = np.all(is_positive(re,odds,dvfg0,dvfg1))
        is_vfg_ordered_ = vfg0 < vfg1
        is_ffg_normalized_ = 0<sfg0<1
        is_in_domain_ = -1<mlog10_1manib<1 and -4<log10_rhos_Msunpc3<4 and 0<log10_rs_pc<5 and 0.5<a<3 and 3<b<10 and 0<g<1.2
        return (is_positive_ and is_vfg_ordered_ and is_ffg_normalized_ and is_in_domain_)
    """
    
    def lnfmems(self,re_pc,dra0,dde0,
            log10_rs_pc,log10_rhos_Msunpc3,a,b,g,
            mlog10_1manib,
            vmem,dist):

        # for debug
        prms = pd.Series(locals()).drop("self")
        display("args of loglikeli:",prms) if DEBUG else None
        
        vs=self.sc_obsdata.radial_velocity.value
        mem,dm= self.dsph.submodels["stellar_model"],self.dsph.submodels["DM_model"]
        
        #update parameters
        mem.update({"re_pc":re_pc*(dist/self.sc_center0.distance.value)}) # Note that re_pc given by the stelar fit is just the angle (re_rad), not re_pc !!!
        dm.update({"rs_pc":power(10,log10_rs_pc),"rhos_Msunpc3":power(10,log10_rhos_Msunpc3),"a":a,"b":b,"g":g})
        self.dsph.update({"anib":1-power(10,-mlog10_1manib)})
        ref_R = mem.half_light_radius() # 1.67834699001666*re
        
        Rs = self.Rs(dra0,dde0,dist) # here 
        
        sigmalos = self.dsph.sigmalos_dequad_interp1d_downsampled(Rs)
        #sigmalos = self.dsph.sigmalos_dequad(Rs)
        
        vobs_err = (self.sc_obsdata.radial_velocity_err.value if hasattr(self.sc_obsdata,"radial_velocity_err") else 0)
        ret = norm.logpdf(vs,loc=vmem,scale=sqrt(sigmalos**2+vobs_err**2))
        return ret 
    
    def lnlikeli(self,params):
        return self.loglikelihood(*params)
    
    def loglikelihood(
        self,re_pc,odds,dra0,dde0,
            log10_rs_pc,log10_rhos_Msunpc3,a,b,g,
            mlog10_1manib,
            vmem,vfg0,vfg1,dvfg0,dvfg1, sfg0, dist
    ): # R_trunc_pc is fixed 2000 pc but its not affect to the result (truncation radius is not used in the calculation of sigmalos in "dSph_Model")

        prms = pd.Series(locals()).drop("self")
        display("args of loglikeli:",prms) if DEBUG else None
        
        vs=self.sc_obsdata.radial_velocity.value
        mem,dm,fg = self.dsph.submodels["stellar_model"],self.dsph.submodels["DM_model"],self.fg
        
        #update parameters
        mem.update({"re_pc":re_pc*(dist/self.sc_center0.distance.value)}) # Note that re_pc given by the stelar fit is just the angle (re_rad), not re_pc !!!
        dm.update({"rs_pc":power(10,log10_rs_pc),"rhos_Msunpc3":power(10,log10_rhos_Msunpc3),"a":a,"b":b,"g":g})
        self.dsph.update({"anib":1-power(10,-mlog10_1manib)})
        ref_R = mem.half_light_radius() # 1.67834699001666*re
        
        Rs = self.Rs(dra0,dde0,dist) # here 
        
        s = 1/(1+ 1/(odds * mem.density_2d_normalized_re(Rs)))
        sigmalos = self.dsph.sigmalos_dequad_interp1d_downsampled(Rs)
        #sigmalos = self.dsph.sigmalos_dequad(Rs)
        sfg1 = 1-sfg0
        
        vobs_err = (self.sc_obsdata.radial_velocity_err.value if hasattr(self.sc_obsdata,"radial_velocity_err") else 0)
        logfmem = norm.logpdf(vs,loc=vmem,scale=sqrt(sigmalos**2+vobs_err**2))
        
        logffg0 = norm.logpdf(vs,loc=vfg0,scale=sqrt(dvfg0**2+vobs_err**2))
        logffg1 = norm.logpdf(vs,loc=vfg1,scale=sqrt(dvfg1**2+vobs_err**2))
        
        logfs = [logfmem,logffg0,logffg1]
        ss = [s,(1-s)*sfg0,(1-s)*sfg1]
                         
        display("sigmalos:{}".format(sigmalos)) if DEBUG else None
        print("fmem:{}".format(fmem)) if DEBUG else None
        print("s*fmem+(1-s)*ffg:{}".format(s*fmem+(1-s)*ffg)) if DEBUG else None
        
        ret = np.sum(logsumexp(a=logfs,b=ss,axis=0)) # note that logsumexp must be used to aviod over/underflow of the likelihood
        
        return ret

    
    
class modKI17_1gauss(modKI17):
    def loglikelihood(self,re_pc,odds,dra0,dde0,
            log10_rs_pc,log10_rhos_Msunpc3,a,b,g,
            mlog10_1manib,
            vmem,vfg0,dvfg0,dist
    ):
        mem = self.dsph.submodels["stellar_model"]
        mem.update({"re_pc":re_pc*(dist/self.sc_center0.distance.value)}) # TODO: "update" is required. This is very confusing and may cause bugs.
        
        vs = self.sc_obsdata.radial_velocity.value
        vobs_err = (self.sc_obsdata.radial_velocity_err.value if hasattr(self.sc_obsdata,"radial_velocity_err") else 0)
        s = 1/(1+ 1/(odds * mem.density_2d_normalized_re(self.Rs(dra0,dde0,dist)))) # In this line, we use "mem.foo", so we should update "mem" in advance.
        
        lnfmems = self.lnfmems(re_pc,dra0,dde0,log10_rs_pc,log10_rhos_Msunpc3,a,b,g,mlog10_1manib,vmem,dist)
        lnffg0s = norm.logpdf(vs,loc=vfg0,scale=sqrt(dvfg0**2+vobs_err**2))
    
        logfs = [lnfmems,lnffg0s]
        ss = [s,(1-s)]
        
        ret = np.sum(logsumexp(a=logfs,b=ss,axis=0)) # note that logsumexp must be used to aviod over/underflow of the likelihood
        
        return ret
    
class modKI17_memonly:
    def __init__(self,vs,dRAs,dDEs,dsph_name,paramlims_fname,beta=1):
        """
        vs: velosity data
        dRAs, dDEs: position data
        RA0,DE0: adhoc central position of dSph 
        paramlims_fname: filename of prior configuration
        """
        self.params_name = ["re","dra0","dde0","log10_rs_pc","log10_rhos_Msunpc3","a","b","g","mlog10_1manib","vmem","dist"]
        self.prior_lim = pd.read_csv(paramlims_fname,index_col=0).loc[self.params_name]
        
        #for val in vs,dRAs,dDEs:
        #    if isinstance(val,pd.DataFrame):
        #        val = val.values
        self.vs = vs
        self.dRAs,self.dDEs = dRAs,dDEs
        
        # read property
        self.RA0,self.DE0 = dSph_property.loc[dsph_name][["RAdeg","DEdeg"]]
        self.DIST,self.err_DIST = dSph_property.loc[dsph_name][["DIST","err_DIST"]]

        self.RoI_R = np.max(self.Rs(0,0,self.DIST)) # use Rs.max as the RoI
        self.beta = beta
        #print(Rs.describe())
        print("beta: {}".format(self.beta)) if beta != 1 else None
        mem = dSph_model.plummer_model(re_pc=200)
        dm = dSph_model.NFW_model(
            a=2.78,b=7.78,g=0.675,
            rhos_Msunpc3=np.power(10,-2.05),rs_pc=np.power(10,3.96),
            R_trunc_pc=2000
        )
        self.dsph = dSph_model.dSph_model(anib=-0.5,submodels_dict={"stellar_model":mem,"DM_model":dm},show_init=True)
        
        #self.fg = dSph_model.uniform2d_model(Rmax_pc=self.RoI_R,show_init=True)

    def Rs(self,dra0,dde0,dist):
        return coord.projected_distance(
            dist=dist,
            ra_center = self.RA0+dra0,
            de_center = self.DE0+dde0,
            ra = self.RA0+self.dRAs,
            de = self.DE0+self.dDEs,
            dtype="deg")

    def lnprior(self,**prms):
        if not self.is_parameters_in_domain(**prms):
            display("out ouf dom:",self.is_parameters_in_domain(**prms)) if DEBUG else None
            return -np.inf
        else:
            #G_re = norm.pdf(re,loc=191,scale=5.7)
            #G_odds = norm.pdf(odds,loc=8.794,scale=0.5107)
            #G_dra0 = norm.pdf(dra0,loc=4.212e-3,scale=7.052e-3)
            #G_dde0 = norm.pdf(dde0,loc=-1.991e-3,scale=3.302e-3)
            G_dist = norm.pdf(prms["dist"],loc=self.DIST,scale=self.err_DIST)
        
            return np.sum(np.log([G_dist,1]))

    def __call__(self,**args):
        ret = self.lnprior(**args)
        if ret == -np.inf:
            return ret
        else:
            return ret + np.sum(self.lnlikelis(**args)) 
        
    def lnprob(self,params):
        prms = {key:val for key,val in zip(self.params_name,params)}
        return self.__call__(**prms)
    
    def is_parameters_in_domain(self,**prms):
        """
        check the parameters are in valid domain.
        """
        prms_df = pd.DataFrame({**self.prior_lim,"inprms":prms})
        display(prms_df) if DEBUG else None
        is_in_minmax = np.all((prms_df.prms_min < prms_df.inprms) & (prms_df.inprms < prms_df.prms_max))
        #is_ordered = prms_df.inprms.vfg0 < prms_df.inprms.vfg1
        display("is_in_minmax",(prms_df.prms_min < prms_df.inprms) & (prms_df.inprms < prms_df.prms_max)) if DEBUG else None
        return is_in_minmax
    """
    def is_parameters_in_domain(self,re,odds,dra0,dde0,
            log10_rs_pc,log10_rhos_Msunpc3,a,b,g,
            mlog10_1manib,
            vmem,vfg0,vfg1,dvfg0,dvfg1, sfg0):
        is_positive_ = np.all(is_positive(re,odds,dvfg0,dvfg1))
        is_vfg_ordered_ = vfg0 < vfg1
        is_ffg_normalized_ = 0<sfg0<1
        is_in_domain_ = -1<mlog10_1manib<1 and -4<log10_rhos_Msunpc3<4 and 0<log10_rs_pc<5 and 0.5<a<3 and 3<b<10 and 0<g<1.2
        return (is_positive_ and is_vfg_ordered_ and is_ffg_normalized_ and is_in_domain_)
    """
    
    def lnlikelis(
        self,re,dra0,dde0,
            log10_rs_pc,log10_rhos_Msunpc3,a,b,g,
            mlog10_1manib,
            vmem, dist
    ): # R_trunc_pc is fixed 2000 pc
        prms = pd.Series(locals()).drop("self")
        display("args of loglikeli:",prms) if DEBUG else None
        
        vs=self.vs
        mem,dm = self.dsph.submodels["stellar_model"],self.dsph.submodels["DM_model"]
        
        #update parameters
        mem.update({"re_pc":re*(dist/self.DIST)}) # Note that re_pc given by the stelar fit is just the angle (re_rad), not re_pc !!!
        dm.update({"rs_pc":power(10,log10_rs_pc),"rhos_Msunpc3":power(10,log10_rhos_Msunpc3),"a":a,"b":b,"g":g})
        self.dsph.update({"anib":1-power(10,-mlog10_1manib)})
        ref_R = mem.half_light_radius() # 1.67834699001666*re
        
        Rs = self.Rs(dra0,dde0,dist) # here 
        
        #s = 1/(1+ 1/(odds * mem.density_2d_normalized_re(Rs)))
        sigmalos = self.dsph.sigmalos_dequad(Rs)
        if np.any(np.isnan(sigmalos)):
            raise TypeError("sigmalos is nan! {}".format(sigmalos))
        #sigmalos = self.dsph.sigmalos_dequad(Rs)
        #sfg1 = 1-sfg0
        
        fmem = norm.pdf(vs,loc=vmem,scale=sigmalos)
        
        #ffg0 = sfg0 * norm.pdf(vs,loc=vfg0,scale=dvfg0)
        #ffg1 = sfg1 * norm.pdf(vs,loc=vfg1,scale=dvfg1)
        #ffg = ffg0+ffg1
        
        display("sigmalos:{}".format(sigmalos)) if DEBUG else None
        print("fmem:{}".format(fmem)) if DEBUG else None
        #print("s*fmem+(1-s)*ffg:{}".format(s*fmem+(1-s)*ffg)) if DEBUG else None
        #ret = np.log(s*fmem+(1-s)*ffg)
        ret = np.log(fmem)
        
        return self.beta * ret
 