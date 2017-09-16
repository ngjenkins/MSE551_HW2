#!/bin/bash

# set working directory
WD=$PWD
# set base file name
name=1-HW2_2_
# remove previous files
if [ -d ${name}_out ]
then
    rm -r ${name}_*
    rm *.e*
    rm *.o*

    #make destination directories
    mkdir ${name}_out
    mkdir ${name}_in
    mkdir ${name}_job
    mkdir ${name}_data
else
    #make destination directories
    mkdir ${name}_out
    mkdir ${name}_in
    mkdir ${name}_job
    mkdir ${name}_data
fi

export N=10

lat_types=('sc' 'fcc' 'bcc')


for i in ${lat_types[@]}
do

#################################################input_file###########################################################
cat >${name}_in/$name$i.in <<!

# Find minimum energy fcc configuration
# Mark Tschopp, 2010
# Modified by Abduljabar Alsayoud

# ---------- Initialize Simulation --------------------- 
clear 
units metal 
dimension 3 
boundary p p p 
atom_style atomic 

# ---------- Variables --------------------- 
variable lat equal 4


# ---------- Create Atoms --------------------- 
#lattice has to be specified first -> all geometry commands are based on it
lattice 	$i \${lat}
#region ID style args keyword (0 1 means 0 lat) (specifies the simulation cell)
region	box block 0 1 0 1 0 1 units lattice
#create_box N region-ID (N=# of atom types)
create_box	1 box

lattice	$i 4 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1  
#create_atoms type style
create_atoms 1 box
replicate 1 1 1

# ---------- Define Interatomic Potential --------------------- 
pair_style eam/alloy 
pair_coeff * * LAMMPS_POTENTIALS/AlCu.eam.alloy Al
neighbor 2.0 bin 
neigh_modify delay 10 check yes 
 
# ---------- Define Settings --------------------- 
#compute ID group-ID style 
#potentail energy per atom
compute poteng all pe/atom
#the sum of all poteng 
compute eatoms all reduce sum c_poteng 


# ---------- Run Minimization --------------------- 
#So timestep start at 0
reset_timestep 0 
fix 1 all box/relax iso 0.0 vmax 0.001
thermo 10 
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms 
min_style cg 
minimize 1e-25 1e-25 5000 10000
write_data ${name}_data/Aloptimized$name$i.data

variable natoms equal "count(all)" 
variable teng equal "c_eatoms"
variable length equal "lx"
variable ecoh equal "v_teng/v_natoms"

print "Total energy (eV) = \${teng};"
print "Number of atoms = \${natoms};"
print "Lattice constant (Angstoms) = \${length};"
print "Cohesive energy (eV) = \${ecoh};"

print "All done!" 
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
export LAMMPS_POTENTIALS=/extra/njenkins/LAMMPS/potentials
export LAMMPS_DATA=/extra/njenkins/LAMMPS/data

mpirun -np 1 lmp_mpi-gcc -sf opt < ${name}_in/$name$i.in > ${name}_out/$name$i.out

rm *.e*
rm *.o*

!
#################################################End_Run_job###############################################################

qsub ${name}_job/lammps$i.job

done

 



