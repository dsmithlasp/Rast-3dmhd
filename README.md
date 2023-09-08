# Rast-3dmhd
 THREE-D MAGNETOHYDRODYNAMICS SIMULATION
 written by Dr. Rast Laboratory for Atmospheric & Space Physics @ University of Colorado Boulder

git clone https://github.com/dsmithlasp/Rast-3dmhd

Feel free to reach out to me for questions Doug.Smith@lasp.colorado.edu

Source code consists of 5 files
The fortran code files
  3dmhd.f
  3dmhdset.f
  3dmhdsub.f
  3dmhdparam.f

Plus a makefile

BUILD HINTS
-----------
you will need to edit the makefile to suit your environment and compilers
The makefile is fairly straighforward you will need to know where your fortran compileres are 
as well as MPI resources.

All the build, run and data dumps all happen in this directory.

Running 'make' in this directory should build the 3dmhd.exe executable in accordance with your makefile.

This code has been run extensively with Intel & GNU Compilers. 
This code has been built and run with OpenMPI, MPICH & Intel MPI libraries.

'make clean' will clean up past object files and you should run this between builds.

RUN HINTS
---------
A typical run command may look like the following:
  'mpirun -np 24 ./3dmhd.exe'
Where -np is the number of processes to use for the job.

This simulation setup as the default in the code should take less than 10 minutes using 24 processes.

OUTPUT
------
THe simulation will create a number of output files.
There should be one output data file for each process containing timestep dumps plus one .par file.
There should be a .lis file which contains simulation summary information. This file is important since 
it will let you know your simulation ran to completion with the following output at the bottom of the file:

"------------------------------------------------------------------------------

     Iteration       1 completed
     -----------------------------------
     Total simulation time:   5.5363D-04
     Present simulation time: 5.5363D-04
     Maximum Mach number:     0.0000D+00
     Present time step:       5.5363D-04
------------------------------------------------------------------------------
               SIMULATION COMPLETE
------------------------------------------------------------------------------
------------------------------------------------------------------------------"


NOTE:
Take some time to peruse the source files there are many usefule tidbits of information in the 
code comments.
