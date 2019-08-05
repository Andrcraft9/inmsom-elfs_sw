MAKE = make

#Translator for Cluster
FC = mpif90

# Sequental options
FCFLAGS = -cpp -dM -ffree-line-length-0 -O3 -fopenmp -I./Inc

#OMP parrallel options
#FCFLAGS = -O3 -openmp -auto -stack_temps -I./Inc

RM = rm -rf
