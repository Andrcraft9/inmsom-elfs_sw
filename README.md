### Shallow water model, 2018. Version with explicit leap-frog scheme.
### High-Communication No-Extra Computation parallel algorithm.
### Fiji problem

### Run guide:

1. Configure basinpar.fi (base area, grid config, nonlinear terms)

2. Configure atmforcing.fi and atmpars.fi if you need (atmosphere data config)

3. Configure reclen.fi (i/o config) and locout.fi (output config)

4. Compile: 
```

make inmsom
```

5. Config phys_proc.par and ocean_run.par

6. Run, for example: 
```

export OMP_NUM_THREADS=1
mpirun -n 4 inmsom
```

7. See results with using GrADS

Note:
By default, model starts from some initial sea level (slf.dat in CP0 folder) with zero velocities. 
