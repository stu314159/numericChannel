WMBrick3D - 3 Diminsional Lattice Boltzmann Simulation of a Wall Mounted Brick.

As the title suggests, this code solves the isothermal incompressible Navier Stokes Equations
via the lattice Boltzmann method. 

The problem domain is a 3D rectangular channel with prescribed constant-velocity inlet and constant
pressure outlet boundary conditions.  The channel is aligned with the Z axis with no-slip boundaries on the 
Z-X plane at Y = 0 and Y = Ly_p.  Periodic boundary conditions are used on the Z-Y plane at X = 0 and X = Lx_p
so really, this code simulates a periodically repeating channel.

The name suggests that inside the channel there is a wall mounted brick, but indeed, the code can be 
modified to include any obstruction in the channel.  That said, the current input file generator -- genInput.py --
produces a cube of size and specified location in the channel with a user-adjustable Reynolds number and fluid properties.

The LBM uses a 15-speed 3D lattice (D3Q15) with regularized boundary conditions and LBGK bulk dynamics.  

Parallelism is provided via MPI and OpenACC.  By adjusting the MPI_FLAGS in the Makefile, the code can, alternatively,
be compiled to run with MPI and OpenMP hybrid parallelism.

To generate inputs, use genInput.py
>>python genInput.py [Num_ts] [ts_reporting frequency] [plot frequency] [Ny_divs]

For initial testing it is recommended to use:
Num_ts = 2000
ts_reporting frequency = 500
plot frequency = 500
Ny_divs = 25


