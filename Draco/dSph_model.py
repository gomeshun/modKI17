import pandas as pd
import numpy as np
import multiprocessing as multi

from numpy import array,pi,sqrt,exp,power,log,log10,cos,tan,sin
from scipy.stats import norm
from scipy.special import k0, betainc, beta, hyp2f1, erf
from scipy import integrate
from scipy.constants import parsec, degree # parsec in meter, degree in radian
from scipy.integrate import quad
from scipy.interpolate import interp1d

GMsun_m3s2 = 1.32712440018e20
R_trunc_pc = 1866.

class model:
    #params, required_params_name = pd.Series(), ['',]
    '''
    params, required_prams_name is undefined; must be defined in child class
    '''
    def __init__(self,submodels_dict=None,**params):
        self.params = pd.Series(params)
        self.required_params_name = list(self.required_params_name)
        self.submodels = submodels_dict
        if 'submodels' in self.params.index:
            raise TypeError('submodels -> submodels_dict?')
        if set(self.params.index) != set(self.required_params_name):
            raise TypeError(self.name+' has the paramsters: '+str(self.required_params_name)+" but in put is "+str(self.params.index))
        if self.submodels != None:
            self.name = ': ' + ' and '.join((model.name for model in self.submodels.values()))

    def __repr__(self):
        ret = self.name + ":\n" + self.show_params().__repr__()
        #if len(self.submodels) > 0:
        #    ret += '\n'
        #    for model in self.submodels.values():
        #        ret += model.__repr__() + '\n'
        return ret

    def show_params(self,models='all'):
        ret = self.params
        if self.submodels != None and models=='all':
            ret = pd.concat([model.show_params('all') for model in self.submodels.values()])
        return ret

    def show_required_params_name(self,models='all'):
        ret = self.required_params_name[:] # need copy because we must keep self.required_params_name
        if self.submodels != None and models=='all':
            [ ret.extend(model.show_required_params_name('all')) for model in self.submodels.values() ]
        return ret

    def is_required_params_name(self,params_name_candidates):
        return [ (p in self.required_params_name) for p in params_name_candidates ]

    def update(self,new_params_dict,target='all'):
        new_params = pd.Series(new_params_dict) if type(new_params_dict)==dict else None
        #print(np.isin(new_params.index, self.show_required_params_name('all')+['this','all']))
        if not np.any(np.isin(new_params.index, self.show_required_params_name('all')+['this','all'])):
            raise TypeError("new params has no required parameters")
        if target in ('this','all'):
            self.params.update(new_params)
            [model.params.update(new_params) for model in self.submodels.values()] if target in ('all',) else None
        else:
            [(model.params.update(new_params) if target==model.name else None) for model in self.submodels.values()]


class stellar_model(model):
    name = "stellar model"
    def density(self,distance_from_center,dimension):
        if dimension == "2d":
            return self.density_2d(distance_from_center)
        elif dimension == "3d":
            return self.density_3d(distance_from_center)


class plummer_model(stellar_model):
    name = "Plummer model"
    required_params_name = ['re_pc',]
    def density_2d(self,R_pc):
        re_pc= self.params.re_pc
        return 1/(1+(R_pc/re_pc)**2)**2 /np.pi/re_pc**2
    def density_3d(self,r_pc):
        re_pc= self.params.re_pc
        return (3/4/np.pi/re_pc**3)/np.sqrt(1+(r_pc/re_pc)**2)**5
    def cdf_R(self,R_pc):
        re_pc= self.params.re_pc
        return 1/(1+(re_pc/R_pc)**2)


class exp2d_model(stellar_model):
    name = "exp2d model"
    required_params_name = ['re_pc',]
    def density_2d(self,R_pc):
        re_pc = self.params.re_pc
        return (1./2/pi/re_pc**2)*exp(-R_pc/re_pc) 
    def density_3d(self,r_pc):
        re_pc = self.params.re_pc
        return (1./2/pi**2/re_pc**3)*k0(r_pc/re_pc)

class uniform2d_model(stellar_model):
    name = "uniform model"
    required_params_name = ['Rmax_pc',]
    def density_2d(self,R_pc):
        return 1./(pi*self.params.Rmax_pc**2)
    def cdf_R(self,R_pc):
        return (R_pc/self.params.Rmax_pc)**2    

class DM_model(model):
    name = "DM model"

class NFW_model(DM_model):
    name = "NFW model"
    required_params_name = ['rs_pc','rhos_Msunpc3','a','b','g','R_trunc_pc']
    
    #def __init__(self,params):
    #    super().__init__(params)
    #    if set(self.params.index) != set(NFW_params_name):
    #        raise TypeError('NFW_model has the paramsters: '+str(NFW_params_name))

    def density_3d(self,r_pc):
        pass
    def enclosure_mass(self,r_pc):
        #ret = (array(r_pc.shape) if len(r_pc)>1 else 0)
        rs_pc, rhos_Msunpc3,a,b,g = self.params.rs_pc, self.params.rhos_Msunpc3, self.params.a, self.params.b,self.params.g
        
        is_in_Rtrunc = r_pc<self.params.R_trunc_pc
        is_outof_Rtrunc = np.logical_not(is_in_Rtrunc)
        
        x = power(r_pc*is_in_Rtrunc/rs_pc,a)
        x_truncd = power(R_trunc_pc/rs_pc,a)
        argbeta0 = (3-g)/a
        argbeta1 = (b-3)/a
        
        #ret[is_in_Rtrunc] = (4.*pi*rs_pc**3*rhos_Msunpc3/a)*beta(argbeta0,argbeta1)*betainc(argbeta0,argbeta1,x/(1+x))
        #ret[is_outof_Rtrunc] = (4.*pi*rs_pc**3*rhos_Msunpc3/a)*beta(argbeta0,argbeta1)*betainc(argbeta0,argbeta1,x_truncd/(1+x_truncd))
        return is_in_Rtrunc * (4.*pi*rs_pc**3*rhos_Msunpc3/a)*beta(argbeta0,argbeta1)*betainc(argbeta0,argbeta1,x/(1+x)) + is_outof_Rtrunc * (4.*pi*rs_pc**3*rhos_Msunpc3/a)*beta(argbeta0,argbeta1)*betainc(argbeta0,argbeta1,x_truncd/(1+x_truncd))
        
class dSph_model(model):
    name = 'dSph_model'
    required_params_name = ['anib','ra_center_deg','de_center_deg','dist_center_pc']
#    def __init__(self,stellar_model,DM_model,**params_dSph_model):
#        """
#        params_dSph_model: pandas.Series, index = (params_stellar_model,params_DM_model,center_of_dSph)
#        """
#        # NOTE IT IS NOT COM{PATIBLE TO THE CODE BELOW!!!
#        super().__init__(**params_dSph_model)
#        self.submodels = (stellar_model,DM_model)
#        self.name = ' and '.join((model.name for model in self.submodels))
    def sigmar2(self,r_pc):
        RELERROR_INTEG = 1e-6
        anib = self.params.anib
        integrand = lambda r,r1: self.submodels["stellar_model"].density_3d(r)*np.power(r/r1,-2*anib)*GMsun_m3s2*self.submodels["DM_model"].enclosure_mass(r)/r**2/self.submodels["stellar_model"].density_3d(r_pc)*1e-6/parsec
        integ, abserr = integrate.quad(integrand,r_pc,np.inf,args=(r_pc,))
        return integ
    
    def naive_sigmalos2(self,R_pc):
        RELERROR_INTEG = 1e-6
        anib = self.params.anib
        integrand = lambda r: (1-anib*(R_pc/r)**2)*self.submodels["stellar_model"].density_3d(r)*self.sigmar2(r)/np.sqrt(1-(R_pc/r)**2)
        rs_interp = np.logspace(-2,6,51)
        integrand_interp = interp1d(rs_interp,[integrand(r) for r in rs_interp],kind="quadratic") 
        integ, abserr = integrate.quad(integrand_interp,R_pc,np.inf)
        return 2*integ/self.submodels["stellar_model"].density_2d(R_pc)

    def integrand_sigmalos2(self,t,arg_R_pc):
        params_dSph,params_stellar,params_DM = self.params,self.submodels["stellar_model"].params,self.submodels["DM_model"].params
        anib = params_dSph.anib
        #print(params_dSph,params_stellar,params_DM)
        rhos_Msunpc3,rs_pc,a,b,g = [params_DM[key] for key in ('rhos_Msunpc3','rs_pc','a','b','g')]
        re_pc = params_stellar.re_pc
        z2 = 1./t/t
        #####weight = sqrt(z2-1)*((1.5-anib)*power(z2,anib)*hyp2f1(0.5,anib,1.5,1-z2)-0.5)
        #weight = sqrt(z2-1)*((1.5-anib)*z2*hyp2f1(1.0,1.5-anib,1.5,1-z2)-0.5) # looks like faster
        ####weight = sqrt(z2-1)*((1.5-anib)*hyp2f1(1,anib,1.5,1-1./z2)-0.5) #NOTE: It looks simple but slow because of the calculation of 2f1
        return sqrt(z2-1)*((1.5-anib)*z2*hyp2f1(1.0,1.5-anib,1.5,1-z2)-0.5) * self.submodels["stellar_model"].density_3d(arg_R_pc/t)*GMsun_m3s2*self.submodels["DM_model"].enclosure_mass(arg_R_pc/t)/parsec # NOTE: without z^-2 !! 1/parsec because we divide it by Sigmaast [pc^-2] so convert one of [m] -> [pc]     

    def sigmalos2_scaler(self,R_pc): # sigmalos2[km^2 s^-2] 
        params_dSph,params_stellar,params_DM = self.params,self.submodels["stellar_model"].params,self.submodels["DM_model"].params
        anib = params_dSph.anib
        #print(params_dSph,params_stellar,params_DM)
        rhos_Msunpc3,rs_pc,a,b,g = [params_DM[key] for key in ('rhos_Msunpc3','rs_pc','a','b','g')]
        re_pc = params_stellar.re_pc

        RELERR_INTEG = 1e-6
        #print("R_pc:",R_pc)
        integ, abserr =  integrate.quad(self.integrand_sigmalos2, 0,1, args=(R_pc,),points=(0.5,0.9,0.99,0.999,0.9999,0.99999))
        return 2*integ/self.submodels["stellar_model"].density_2d(R_pc)*1e-6
    
    def sigmalos2_vector(self,Rs_pc):
        sigmalos2_scaler = self.sigmalos2_scaler
        return [sigmalos2_scaler(R_pc) for R_pc in Rs_pc]

    def sigmalos2(self,R_pc):
        return (self.sigmalos2_scaler(R_pc) if len(R_pc)==1 else self.sigmalos2_vector(R_pc))
    
    def distribution_func(self,vs,Rs):
        pass

class KI17_model:
    def __init__(self,params_KI17_model):
        """
        params_KI17_model: pandas.Series, index = (params_dSph_model,params_FG_model,s)
        """
        pass


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    dm_model = NFW_model(a=2.78,b=7.78,g=0.675,rhos_Msunpc3=np.power(10,-2.05),rs_pc=np.power(10,3.96),R_trunc_pc=2000)
    mystellar_model = plummer_model(re_pc=221)
    draco_model = dSph_model(submodels_dict={"DM_model":dm_model,"stellar_model":mystellar_model},anib=1-np.power(10,0.13),dist_center_pc=76e3,ra_center_deg=0,de_center_deg=0)
    
    Rs = np.logspace(-1,9,100)
    ss = draco_model.sigmalos2(Rs)
    ss2 = [draco_model.sigmar2(R) for R in Rs]
    ss3 = [draco_model.naive_sigmalos2(R) for R in Rs]
    print(draco_model)
    print(draco_model.integrand_sigmalos2(1,1))
    plt.plot(Rs,np.sqrt(ss))
    plt.plot(Rs,np.sqrt(ss2))
    plt.plot(Rs,np.sqrt(ss3))
    plt.xscale("log")
    plt.yscale("log")
    #plt.ylim(0,40)
    plt.show()
    input("press any key")







