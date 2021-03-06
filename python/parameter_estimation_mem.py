import pandas as pd
import numpy as np
import sys
import glob
import model_definition as model
from model_definition import initial_parameters_mem, initial_parameters_mem, initial_parameters_fg
import MCgenerator

RoI_R = 200
RoI_vs = (0,100)

argvs = sys.argv
argc = len(argvs)
####    ####
mess = ['member stars filename: ','foreground stars filename: ','putput filename: ','log filename: ']
mem_fname = argvs[1] if argc>1 else input(mess[0])
out_fname = argvs[2] if argc>3 else input(mess[1])
log_fname = argvs[3] if argc>4 else input(mess[2])
fnames = (mem_fname,out_fname,log_fname)
####    ####
[print(mes+fname) for mes,fname in zip(mess,fnames)]
if "--no-confirmation" in argvs:
    print("--no-confirmation option called. do without confirmation.")
else:
    yn = input("OK? (y/n)")
    if yn not in ("yes","y"):
        print("you typed: "+yn+" so exited.")
        quit()

for fname in (out_fname,log_fname):
    is_already_exist = (fname in glob.glob("*"))
    if is_already_exist:
        print(fname+' is already exist!')
        
        
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
#mem = mem[(mem.R<RoI_R).values * (RoI_vs[0]<mem.R).values * (mem.R<RoI_vs[1]).values]
if "kind" in mem.columns:
    mem = mem[mem.kind=="mem"]
print("read mocks...completed.")

# impose some bias here
# mem =
# fg = 
tot = mem #mem.append(fg)
print("impose biases...completed.")

class KI17_loglikelihood_memonly:
    def __init__(self,df):
        self.vs = df.v
        self.Rs = df.R

    def __call__(self,v_mem,a,b):
        r_e = initial_parameters_mem['r_e']
        if a<0 or a+b*(model.RoI_R/r_e) < 0:
            return -np.inf
        else:
            fmem = model.dist_func_mem
            loglikelis = np.log(fmem(self.vs,self.Rs,r_e,v_mem,a,b))
            return np.sum(loglikelis)

loglikeli = KI17_loglikelihood_memonly(tot)
args_logpdf_init = {key:model.initial_parameters[key] for key in ("v_mem","a","b")}
dargs_logpdf_init = {'v_mem':0.2,'a':0.2,'b':0.2}
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





