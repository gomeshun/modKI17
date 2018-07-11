from model_definition import *
from plot_setup import *
from copy import deepcopy

plot_steps = 16
alpha = 0.2

sig_only = deepcopy(initial_parameters)
sig_only['s'] = 1
fg_only = deepcopy(initial_parameters)
fg_only['s'] = 0

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

Rs = np.linspace(0-20,RoI_R+20,plot_steps)
vs = np.linspace(RoI_v_lo-20,RoI_v_hi+20,plot_steps)

vs_meshed,Rs_meshed = np.meshgrid(vs,Rs)

fs_meshed = dist_func_tot(v=vs_meshed,R=Rs_meshed,**initial_parameters)
fs_meshed_sig = initial_parameters['s']*dist_func_tot(v=vs_meshed,R=Rs_meshed,**sig_only)
fs_meshed_fg = (1-initial_parameters['s'])*dist_func_tot(v=vs_meshed,R=Rs_meshed,**fg_only)

ax.plot_surface(vs_meshed,Rs_meshed,fs_meshed,alpha=alpha,rstride=1,cstride=1)
ax.plot_surface(vs_meshed,Rs_meshed,fs_meshed_sig,alpha=alpha,rstride=1,cstride=1)
ax.plot_surface(vs_meshed,Rs_meshed,fs_meshed_fg,alpha=alpha,rstride=1,cstride=1)
ax.plot_wireframe(vs_meshed,Rs_meshed,fs_meshed,alpha=alpha,rstride=1,cstride=1,color='C0')
ax.plot_wireframe(vs_meshed,Rs_meshed,fs_meshed_sig,alpha=alpha,rstride=1,cstride=1,color='C1')
ax.plot_wireframe(vs_meshed,Rs_meshed,fs_meshed_fg,alpha=alpha,rstride=1,cstride=1,color='C2')
ax.scatter3D((),(),(),label='tot',color='C0')
ax.scatter3D((),(),(),label='mem',color='C1')
ax.scatter3D((),(),(),label='fg',color='C2')
ax.set_xlabel("vs")
ax.set_ylabel("Rs")
ax.set_zlabel("fs")
ax.legend()

fig.show()
input("press any key")
