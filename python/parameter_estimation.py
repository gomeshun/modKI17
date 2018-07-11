import pandas as pd
import numpy as np
import sys
import model_definition as model
from model_definition import initial_parameters_mem, initial_parameters_mem, initial_parameters_fg
import MCgenerator

argvs = sys.argv
argc = len(argvs)
####    ####
mess = ['member stars filename: ','foreground stars filename: ','putput filename: ','log filename: ']
mem_fname = argvs[1] if argc>1 else input(mess[0])
fg_fname = argvs[2] if argc>2 else input(mess[1])
out_fname = argvs[3] if argc>3 else input(mess[2])
log_fname = argvs[4] if argc>4 else input(mess[3])
fnames = (mem_fname,fg_fname,out_fname,log_fname)
####    ####
[print(mes+fname) for mes,fname in zip(mess,fnames)]
if "--no-confirmation" in argvs:
    print("--no-confirmation option called. do without confirmation.")
else:
    yn = input("OK? (y/n)")
    if yn not in ("yes","y"):
        print("you typed: "+yn+" so exited.")
        quit()

print("parameter_estimation start.")

# define the MC parameter here
n_burnin = 10000 if ("--burnin" not in argvs) else int(argvs[argvs.index("--burnin")+1])
burnin_repeat = 2 if ("--burnin-repeat" not in argvs) else int(argvs[argvs.index("--burnin-repeat")+1])
n_sampling = 20000 if ("--sampling" not in argvs) else int(argvs[argvs.index("--sampling")+1])

mcsteps = [n_burnin,]*burnin_repeat
mcsteps.append(n_sampling)
print("mcsteps:\n",mcsteps)

# import unbiased observed data
mem = pd.read_csv(mem_fname)
fg = pd.read_csv(fg_fname)
if "kind" in mem.columns:
    mem = mem[mem.kind=="mem"]
if "kind" in fg.columns:
    fg = fg[fg.kind=="fg"]
print("read mocks...completed.")

# impose some bias here
# mem =
# fg = 
tot = mem.append(fg)
print("impose biases...completed.")

class KI17_loglikelihood:
    def __init__(self,df):
        self.vs = df.v
        self.Rs = df.R

    def __call__(self,v_mem,a,b,s):
        r_e,v_fg,dv_fg = initial_parameters_mem['r_e'],initial_parameters_fg['v_fg'],initial_parameters_fg['dv_fg']
        if a<0 or a+b*(model.RoI_R/r_e) < 0 or s<0 or 1<s:
            return -np.inf
        else:
            fmem,ffg = model.dist_func_mem,model.dist_func_fg
            loglikelis = np.log(s*fmem(self.vs,self.Rs,r_e,v_mem,a,b)+(1-s)*ffg(self.vs,self.Rs,v_fg,dv_fg))
            return np.sum(loglikelis)

loglikeli = KI17_loglikelihood(tot)
args_logpdf_init = {key:model.initial_parameters[key] for key in ("v_mem","a","b","s")}
dargs_logpdf_init = {'v_mem':0.2,'a':0.2,'b':0.2,'s':0.01}
print("likelihood has been defined.")

gen = MCgenerator.MCgenerator(
    logpdf_func = loglikeli,
    args_logpdf_init = args_logpdf_init,
    dargs_logpdf = dargs_logpdf_init,
    push_time=5)
print("MC generator setup.")


gen.generate_tuned(mcsteps)
gen.to_csv(out_fname)
gen.to_csv(log_fname,output_log=True)





