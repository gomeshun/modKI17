from model_definition import generate_mem,generate_fg,initial_parameters_mem,initial_parameters_fg
from MCgenerator import MCgenerator

initial_parameters_mem["a"] = 30
initial_parameters_mem["b"] = -10

fname_mem = "mock_mem_narrowed_v002.csv"
fname_fg = "mock_fg_narrowed_v002.csv"
#fname_mem = "mock_mem_v005.csv"
#fname_fg = "mock_fg_v005.csv"

N_mem,N_fg = 25000,25000
N_stride = 30
push_time = 5

print("mock_mem creation start.")
mem = generate_mem(**initial_parameters_mem,N_sampling=N_mem,N_stride=N_stride,push_time=push_time)
print("modk_fg createion start.")
fg = generate_fg(**initial_parameters_fg,N_sampling=N_fg,N_stride=N_stride,push_time=push_time)
print("mem,fg mocks has been created.")
#plt.hist2d(mem.v,mem.R)
#plt.hist2d(fg.v,fg.R)
#plt.xlabel("v")
#plt.ylabel("R")
#plt.show()
yn_savemock = 'y' #input("save the created mock? (y/n)")
issavemock = (yn_savemock=='y')
if issavemock:
    mem.to_csv(fname_mem,index=None)
    fg.to_csv(fname_fg,index=None)
    print("mocks saved.")
