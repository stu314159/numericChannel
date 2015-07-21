# Python code to generate LBM input file
import math
import argparse
import numpy as np

"""
Input file generator for the LBM code for 3D flow past a wall mounted brick.

Num_ts = total number of LBM time steps to execute
ts_rep_frequency = time-step reporting frequency. (LBM code gives periodic updates on progress)
plot_freq = number of time-steps between data dumps on pressure and velocity
Ny_divs = number of lattice points across the entrance of the channel

The code is set below to provide inputs that will help ensure stable LBM execution,
but reliable data will be produced only if the reported LBM flow mach number
is less than 0.1.  I recommend you do not procede with actually performing
the simulation unless this is the case.

Otherwise, more accurate results will be obtained as you increase Ny_divs.
You can get more detailed data if you increase the plot_freq.

Some suggested problem inputs:

# debug
Num_ts = 2
ts_rep_freq = 1
plot_freq = 1
Ny_divs = 5

# Small -- to make sure everything is running right:
Num_ts = 2000
ts_rep_freq = 100
plot_freq = 500
Ny_divs = 15


# Medium - to get a feel for how long codes will take to run as the problem size is increased:
Num_ts = 10000
ts_rep_freq = 500
plot_freq = 5000
Ny_divs = 25

# Large - for benchmarking performance scaling
Num_ts = 30000
ts_rep_freq = 1000
plot_freq = 29999
Ny_divs = 40




"""


# build input parser
parser = argparse.ArgumentParser(prog='genInput.py',description='Process LBM input arguments')
parser.add_argument('Num_ts',type=int)
parser.add_argument('ts_rep_freq',type=int)
parser.add_argument('plot_freq',type=int)
parser.add_argument('Ny_divs',type=int)
# parse input arguments
args = parser.parse_args()

# assign to required variables
Num_ts = args.Num_ts
ts_rep_freq = args.ts_rep_freq
plot_freq = args.plot_freq
Ny_divs = args.Ny_divs




#----You should not have to edit anything below this line -------------------

Re = 67.
# overall domain dimensions
Lx_p = 4. # "thickness"
Ly_p = 3. # "height"
Lz_p = 14. # "length"

# describe brick dimensions and location
h_brick = 1.
z_brick = 6.
x_brick = 2.


# for now, just doing Poiseuille flow.
Lo = h_brick;#characteristic length is block height
R=0; x_c = 0; z_c = 0;# not used.

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


print 'Computing global node numbers within the wall mounted brick'
# to do:  identify the elements that lie within the obstacle
# bone-headed way:
obstList = []
for gid in range(int(numEl)):
    xG = XX[gid]; yG = YY[gid]; zG = ZZ[gid];
    if ((xG>= (x_brick - h_brick/2.)) and (xG <= (x_brick + h_brick/2.)) \
        and (yG <= h_brick)  and (zG >= (z_brick-h_brick/2.)) and \
        (zG <= (z_brick + h_brick/2.))):
        obstList.append(gid)

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
# set omega as a constant at 1.6
omega = 1.74935
nu_lbm = ((1./omega) - 0.5)/3.
dt = nu_lbm*(dx**2)/nu_d
u_lbm = (dt/dx)*Ud


u_conv_fact = (dx/dt)*(Lo/To)
t_conv_fact = (dt*To)
l_conv_fact = dx*Lo
p_conv_fact = ((l_conv_fact/t_conv_fact)**2)*(1./3.)

rho_lbm = rho_p





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
    
    
    #params.write('%d \n'%1) # dynamics (not needed - add later if used)
    #params.write('%d \n'%0) # entropic (boolean) ( not needed )
    #params.write('%d \n'%0) # initialization (not needed)
    
    
    params.write('%d \n'%Num_ts)
    params.write('%d \n'%ts_rep_freq)
    params.write('%d \n'%plot_freq)
   
    params.write('%f \n'%rho_lbm) # density
    params.write('%f \n'%u_lbm) # scaled maximum velocity
    params.write('%f \n'%omega) # relaxation parameter
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
    
    
    
