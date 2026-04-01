Code to calibrate model parameters using 1D EAT framework and eatpy library (Bruggeman et. al. 2024, https://gmd.copernicus.org/articles/17/5619/2024/). 
The code assumes that observation files have variables structured in gridded time x depth format, where depth index corresponding to the real depth in m.
The basic functionality is in calibrate_models.py and the different tasks can be run from run.py, run_seasons.py and run_obs_reductions.py.

This is a fork of the original repo at https://github.com/JOZSKA/Parameter_calibration_1D_models.
The new features include:
- separate configuration file
- parallelization with MPI

### Quick start
- Prepare your configuration as derived class from BaseConfiguration. See L4 class in configurations.py as example.
- Prepare a python file including the desired calibrations. See main.py as an example.
- install and activate a conda env with EAT:
```
conda activate eat
```
- Run with MPI with the available number of processes:
```
AVAIL_PROCS=100
mpiexec -n $AVAIL_PROCS python main.py
```
