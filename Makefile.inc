MAKE = make

#Translator for Cluster
FC = mpif90

# Sequental options
FCFLAGS = -w -ffree-line-length-0 -O3 -openmp -I./Inc

#OMP parrallel options
#FCFLAGS = -O3 -openmp -auto -stack_temps -I./Inc

RM = rm -rf

