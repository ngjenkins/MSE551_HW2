#!/bin/bash

# set working directory
WD=$PWD
# set base file name
name=1-HW2_1_
# remove previous files
rm -r ${name}_out
rm -r ${name}_in
rm -r ${name}_job
rm -r ${name}_data
rm *.e*
rm *.o*
#make destination directories
mkdir ${name}_out
mkdir ${name}_in
mkdir ${name}_job
mkdir ${name}_data

export N=25

for i in `seq 0 $N`
do

#################################################input_file###########################################################
cat >${name}_in/$name$i.in <<!

units lj
atom_style atomic
boundary p p p
lattice sc 1.0

# create simulation cell
region r1 block -5.0 5.0 -5.0 5.0 -5.0 5.0
create_box 1 r1

# required must come after box is created
mass 1 1.0

# create two atoms,
variable N equal "$N" 
variable ri equal "0.9"
variable rf equal "4.0"
variable rs equal "(v_rf-v_ri)/v_N"
variable r equal "(v_rs*$i)+v_ri"

create_atoms 1 single  0.0 0.0 0.0
create_atoms 1 single  \$r 0.0 0.0

# set non-bonded potential
pair_style lj/cut 5.0
pair_coeff 1 1 1.0 1.0

run_style verlet
timestep  0.005   

minimize 1.0e-10 1.0e-10 0 0

write_data ${name}_data/LJoptimized$name$i.data
!
#################################################end_input_file###########################################################



###################################################Run_job################################################################
cat >${name}_job/lammps$i.job <<!

#!/bin/bash

#PBS -N $name$i
#PBS -W group_list=mse551-2017
#PBS -q windfall
### "pcmem=6gb" is the memory attribute for all of the standard nodes
#PBS -l select=1:ncpus=1:mem=6gb:pcmem=6gb
#PBS -l place=free:shared
#PBS -l walltime=1:00:00
#PBS -l cput=1:00:00

### cd: set directory for job execution, ~netid = home directory path
cd $WD

### Load required modules/libraries
module load lammps/gcc/17Nov16 

export MPI_DSM_DISTRIBUTE
export OMP_NUM_THREADS=1

mpirun -np 1 lmp_mpi-gcc -sf opt < ${name}_in/$name$i.in > ${name}_out/$name$i.out

rm *.e*
rm *.o*

!
#################################################End_Run_job###############################################################

qsub ${name}_job/lammps$i.job

done

 



