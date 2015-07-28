# Python code to generate LBM input file
import math
import argparse
import numpy as np
import scipy.io

"""
Input file generator for the LBM code for 3D flow past an obstacle whose geometry is described
by 'turbine_blade.mat'

Num_ts = total number of LBM time steps to execute
ts_rep_frequency = time-step reporting frequency. (LBM code gives periodic updates on progress)
Warmup_ts = number of time steps to take before recording data
plot_freq = number of time-steps between data dumps on pressure and velocity
Re = Flow Reynolds number
Cs = Turbulence model parameter

Some suggested problem inputs:

To compare with non-turbulent version
Num_ts = 2000
ts_rep_freq = 1000
Warmup_ts = 0
plot_freq = 500
Re = 67
dt = .0025
Cs = 0

For some initial problem testing:
Num_ts = 30000
ts_rep_freq = 1000
Warmup_ts = 0
plot_freq = 10000
Re = 5000
dt = 0.0005
Cs = 2

For longer testing with turbulent flow simulations that you expect to run stably,
conduct an appropriate warmup period (you should determine this with preliminary runs)
followed by a comparatively short number of time steps during which you collect data frequently.


"""


# build input parser
parser = argparse.ArgumentParser(prog='genInput.py',description='Process LBM input arguments')
parser.add_argument('Num_ts',type=int)
parser.add_argument('ts_rep_freq',type=int)
parser.add_argument('Warmup_ts',type=int)
parser.add_argument('plot_freq',type=int)
parser.add_argument('Re',type=float)
parser.add_argument('dt',type=float)
parser.add_argument('Cs',type=float)
# parse input arguments
args = parser.parse_args()

# assign to required variables
Num_ts = args.Num_ts
ts_rep_freq = args.ts_rep_freq
Warmup_ts = args.Warmup_ts
plot_freq = args.plot_freq
Re = args.Re
dt = args.dt
Cs = args.Cs




#----You should not have to edit anything below this line -------------------

turb_input = scipy.io.loadmat('turbine_blade.mat')
# overall domain dimensions
Lx_p = float(turb_input['Lx_p'])
Ly_p = float(turb_input['Ly_p'])
Lz_p = float(turb_input['Lz_p'])
Lo = float(turb_input['Lo'])
Ny_divs = int(turb_input['Ny_divs'])
obstList = list((turb_input['gnn']).flatten())

Ny = math.ceil((Ny_divs-1)*(Ly_p/Lo))+1
Nx = math.ceil((Ny_divs-1)*(Lx_p/Lo))+1
Nz = math.ceil((Ny_divs-1)*(Lz_p/Lo))+1
nnodes = Nx*Ny*Nz

# compute geometric data only once
x = np.linspace(0.,Lx_p,Nx).astype(np.float32);
y = np.linspace(0.,Ly_p,Ny).astype(np.float32);
z = np.linspace(0.,Lz_p,Nz).astype(np.float32);
numEl = Nx*Ny*Nz
Y,Z,X = np.meshgrid(y,z,x);

XX = np.reshape(X,numEl)
YY = np.reshape(Y,numEl)
ZZ = np.reshape(Z,numEl)

print 'There are %d nodes in the obstacle'%len(obstList)
print 'Writing those nodes to file'
# now write this obstList to file.
obstFilename = 'obst_file.lbm'
obstFile = open(obstFilename,'w')
obstFile.write('%i \n'%len(obstList))
for i in range(len(obstList)):
    obstFile.write('%i \n'%obstList[i])
obstFile.close()


rho_p = 1000. # physical density
nu_p = 0.001/rho_p #physical kinematic viscosity

# non-dimensionalization
Uo = nu_p*Re/Lo
To = Lo/Uo
Uavg = Uo

Ld = 1.; Td = 1.; Ud = (To/Lo)*Uavg;
nu_d = 1./Re
dx = 1./(Ny_divs - 1.)
u_lbm = (dt/dx)*Ud
nu_lbm = (dt/(dx**2))*nu_d
omega = 1./(3.*nu_lbm+0.5)

u_conv_fact = (dx/dt)*(Lo/To)
t_conv_fact = (dt*To)
l_conv_fact = dx*Lo
p_conv_fact = (((l_conv_fact/t_conv_fact)**2)*(1./3.))/(l_conv_fact**3)

rho_lbm = rho_p*(l_conv_fact**3)

print 'l_conv_fact = %g.\n'%l_conv_fact
print 'p_conv_fact = %g.\n'%p_conv_fact


print 'Number of lattice points = %d.' % nnodes
print 'Number of time-steps = %d.' % Num_ts
print 'LBM viscosity = %g.' % nu_lbm
print 'LBM relaxation parameter (omega) = %g.' % omega
print 'LBM flow Mach number = %g. ' % u_lbm
print 'Nx = %d' % Nx
print 'Ny = %d' % Ny
print 'Nz = %d' % Nz

#run_dec = raw_input('Would you like to continue? [Y/n]: ')
run_dec = 'y' # just let it run

if run_dec!='n' and run_dec!='N':
    print 'Ok! Cross your fingers!!'
    # write the input file
    params = open('params.lbm','w')
    params.write('%d \n'%1) # lattice selection (keep.  We might actually use this)
    
    params.write('%d \n'%Num_ts)
    params.write('%d \n'%ts_rep_freq)
    params.write('%d \n'%Warmup_ts)
    params.write('%d \n'%plot_freq)
    params.write('%g \n'%Cs)
    params.write('%g \n'%rho_lbm) # density
    params.write('%g \n'%u_lbm) # scaled maximum velocity
    params.write('%g \n'%omega) # relaxation parameter
    params.write('%d \n'%Nx) # number of nodes in the x, y and z direction
    params.write('%d \n'%Ny)
    params.write('%d \n'%Nz)
    
    # the following will not be used by the MPI code, but will be available
    # for the post-processing script
    
    params.write('%f \n'%Lx_p) # physical dimensions in the x,y and z dimensions
    params.write('%f \n'%Ly_p)
    params.write('%f \n'%Lz_p)
    
    params.write('%f \n'%t_conv_fact)  # time, length and pressure conversion factors
    params.write('%f \n'%l_conv_fact)
    params.write('%g \n'%p_conv_fact)
    
    params.close()
    
else:
    print 'Run aborted.  Better luck next time!'
    
    
    
