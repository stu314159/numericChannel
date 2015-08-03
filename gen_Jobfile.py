#gen_Jobfile.py
"""
more automated generation of PBS jobfile 

"""
jobfileName = 'wmb_jobfile.pbs'
executableName = 'WMBrick3D'

jobName = 'N4P8_intel' # must be 15 characters or less
nnodes = 4
ppn = 8 # for LBM jobs on the Cray XC30 machines, 8 ppn seems to saturate memory bandwidth.
mpi_procs_per_node = 1 # only change this if you need more memory per process.
walltime='04:00:00'
platform='SHEPARD'
queue = 'standard'

proc_script='process_lbm_data.py'

filesToCopy=[ executableName, 'obst_file.lbm','params.lbm',proc_script,'vtkHelper.py','validate_r3.py','gold_standard.npy']



#--------- more-or-less fixed code below -----------------

proj_id = 'USNAM37752431' # <-- don't post this anywhere public
mpi_procs = mpi_procs_per_node*nnodes


# open the file
jf = open(jobfileName,'w')

# essential PBS directives
jf.write('#!/bin/bash \n') # the shell
jf.write('#PBS -A %s \n'%proj_id) # project identifier
jf.write('#PBS -q %s \n'%queue) # specify queue
jf.write('#PBS -l select=%d:ncpus=%d:mpiprocs=%d \n'% \
         (nnodes,ppn,mpi_procs_per_node))
jf.write('#PBS -l walltime=%s \n'%walltime)
jf.write('#PBS -l ccm=1 \n') # specify cluster compatibility mode.  Why wouldn't you?

#optional PBS directives
jf.write('#PBS -N %s \n'%jobName)
jf.write('#PBS -j oe \n')
jf.write('#PBS -V \n')
jf.write('#PBS -S /bin/bash \n')


# Execution block
jf.write('cd $WORKDIR\n')
jf.write("JOBID=`echo $PBS_JOBID | cut -d '.' -f 1` \n")
jf.write('if [ ! -d $JOBID ]; then \n')
jf.write('  mkdir -p $JOBID \n')
jf.write('fi \n')
jf.write('cd $JOBID \n')
# copy files 
for s in filesToCopy:
    jf.write('cp $PBS_O_WORKDIR/%s . \n'% s)

## move to the $JOBDIR
#jf.write('cd $JOBDIR \n')  #<--- this was an error

# invoke execution
jf.write('export OMP_NUM_THREADS=%d \n'%ppn)
jf.write('aprun -B ./%s \n'%(executableName))
#jf.write('python ./%s \n'%proc_script)

# create job to cleanup and archive data
#jf.write('cd $WORKDIR \n')
#jf.write('rm -f archive_job \n')
#jf.write('cat > archive_job << END \n')
#jf.write('#!/bin/bash \n')
#jf.write('#PBS -l walltime=06:00:00 \n')
#jf.write('#PBS -q transfer \n')
#jf.write('#PBS -A %s \n'% proj_id)
#jf.write('#PBS -l select=1:ncpus=1 \n')
#jf.write('#PBS -j oe \n')
#jf.write('#PBS -S /bin/bash \n')
#jf.write('cd $WORKDIR \n')
#jf.write('rsh $ARCHIVE_HOST mkdir $ARCHIVE_HOME/$JOBID \n')
#jf.write('rcp -r $JOBID $ARCHIVE_HOST:$ARCHIVE_HOME/ \n')
#jf.write('rsh $ARCHIVE_HOST ls -l $ARCHIVE_HOME/$JOBID \n')
#jf.write('rm -rf $JOBID \n')
#jf.write('END \n')
#jf.write('qsub archive_job \n')

# end of archive script portion.

# close the file
jf.close()
