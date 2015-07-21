#vtkHelper.py

import numpy as np
from array import array

def saveStructuredPointsVTK_ascii(scalar,scalar_name,filename,dims,origin,spacing):
    """
    scalar = iterable data structure containing scalar data (assume floating point)
    scalar_name = string holding name of data
    dims = iterable containing [Nx, Ny, Nz]
    origin = iterable containing the origin
    spacing = iterable containing the space size between points
    
    """
    numEl = dims[0]*dims[1]*dims[2]
    file = open(filename,'w')
    file.write('# vtk DataFile Version 3.0\n')
    file.write('VTK file for data post-processed with Python\n')
    file.write('ASCII\n\n')
    file.write('DATASET STRUCTURED_POINTS\n')
    file.write('DIMENSIONS %d %d %d \n'%(dims[0],dims[1],dims[2]))
    file.write('ORIGIN %g %g %g \n '%(origin[0],origin[1],origin[2]))
    file.write('SPACING %g %g %g \n'%(spacing[0],spacing[1],spacing[2]))
    file.write('POINT_DATA %d \n'%numEl)
    file.write('SCALARS %s float 1 \n'%scalar_name)
    file.write('LOOKUP_TABLE default \n')
    for i in range(numEl):
        file.write('%g \n'%scalar[i])
    file.close()
    


def saveScalarStructuredGridVTK_binary(scalar,scalar_name,x,y,z,filename,dims):
    """
    save a VTK 3.0 legacy STRUCTURED_GRID data set in binary format.
    scalar = scalar data to be saved
    scalar_name = string with the name of the variable (e.g. pressure)
    x,y,z = numpy array with the geometry data
    filename = string with desired file name
    dims = 3-tuple with number of rank of each dimension
    
    
    
    """
    
    numEl_size = x.size; numEl = np.prod(numEl_size);
    # open the file and write the ASCII header:
    file = open(filename,'w')
    file.write('# vtk DataFile Version 3.0\n')
    file.write('VTK file for data post-processed with Python\n')
    file.write('Binary\n\n')
    file.write('DATASET STRUCTURED_GRID\n')
    file.write('DIMENSIONS %d %d %d \n'%(dims[0],dims[1],dims[2]))
    file.write('POINTS %d float\n'%(numEl))
    file.close()
    
    # append binary x,y,z data
    file = open(filename,'ab')
    for i in range(len(x)): # there really needs to be a better way.
        pt = [x[i],y[i],z[i]]
        pt_buf = array('f',pt)
        pt_buf.byteswap()
        file.write(pt_buf)
         
    
    file.close()
    
     # append another ASCII sub header for the scalar data
    file = open(filename,'a')
    file.write('\nSCALARS %s int\n'%scalar_name)
    file.write('LOOKUP_TABLE default\n')
    file.close()
    
    # append binary scalar data
    file = open(filename,'ab')
    p_buf = array('f',scalar); p_buf.byteswap()
    file.write(p_buf)
    file.close()
    
    

def saveVelocityAndPressureVTK_binary(pressure,u,v,w,x,y,z,filename,dims):
    """
    save a VTK 3.0 legacy STRUCTURED_GRID data set in binary format.
    pressure = numpy array with pressure values (data type??)
    u,v,w = numpy array with velocity data
    x,y,z = numpy array with geometry data
    filename = string with desired file name
    dims = 3-tuple with the number of rank of each dimension
    
    """
    numEl_size = u.size; numEl = np.prod(numEl_size);
    # open the file and write the ASCII header:
    file = open(filename,'w')
    file.write('# vtk DataFile Version 3.0\n')
    file.write('VTK file for data post-processed with Python\n')
    file.write('Binary\n\n')
    file.write('DATASET STRUCTURED_GRID\n')
    file.write('DIMENSIONS %d %d %d \n'%(dims[0],dims[1],dims[2]))
    file.write('POINTS %d float\n'%(numEl))
    file.close()
    
    # append binary x,y,z data
    file = open(filename,'ab')
    for i in range(len(x)): # there really needs to be a better way.
        pt = [x[i],y[i],z[i]]
        pt_buf = array('f',pt)
        pt_buf.byteswap()
        file.write(pt_buf)
         
    
    file.close()
    
    # append an ASCII sub header
    file = open(filename,'a')
    file.write('POINT_DATA %d \n'%numEl)
    file.write('VECTORS velocity_vectors float\n')
    file.close()
    
    # append binary u,v,w data
    file = open(filename,'ab')
    for i in range(len(u)):
        pt = [u[i],v[i],w[i]]
        pt_buf = array('f',pt)
        pt_buf.byteswap()
        file.write(pt_buf)
   
    file.close()
    
    # append ASCII sub header for scalar velocity magnitude data
    file = open(filename,'a')
    file.write('SCALARS VelocityMagnitude float\n')
    file.write('LOOKUP_TABLE default\n')
    
    file.close()
    
    file = open(filename,'ab')
    v_mag = np.sqrt(u**2+v**2+w**2)
    file = open(filename,'ab')
    p_buf = array('f',v_mag); p_buf.byteswap()
    file.write(p_buf)
    file.close()
    
    
    # append another ASCII sub header for the scalar pressure data
    file = open(filename,'a')
    file.write('SCALARS Pressure float\n')
    file.write('LOOKUP_TABLE default\n')
    file.close()
    
    # append binary pressure data
    file = open(filename,'ab')
    p_buf = array('f',pressure); p_buf.byteswap()
    file.write(p_buf)
    file.close()