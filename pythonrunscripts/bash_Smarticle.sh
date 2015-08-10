#!/bin/sh

#PBS -N SmarticleSimulation
#PBS -l nodes=1:ppn=1,walltime=4:00:00

cd $PBS_O_WORKDIR

/home/wsavoie/ChronoSrc/SmarticlesBuild/SmarticlesSystem > bash_Smarticle.out.$PBS_JOBID.`hostname`
