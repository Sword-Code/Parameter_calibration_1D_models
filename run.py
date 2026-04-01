# identifies the best matching parameters to the observations across an ensemble of runs

import csv

import numpy as np

from mpi import mpi
import eatpy
import calibrate_models as cm
from configurations import chosen_conf


# #below you need to set-up as appropriate
# 
# path_observations = "./L4_observations.nc"  #supply path to observational file. The observational file is assumed to be containing 1D, or 2D arrays with (time) or (time, depth) dimensions, where time is always labeled as number of days from 01/01/1998 (it is daily resolution) and depth is always spaced by 1m (so N vertical layers means going N meters deep). For simplification no fluctuations in the sea level height are considered...
# 
# model_directory = "/work/jos/eat/1D_configuration/L4/large_ensemble"    # path to folder with the model ensemble simulations - the model outputs are picked across the ensemble from there
# 
# perturbed_parameters_listed = ["B1_rR2", "B1_sR1", "P1_alpha", "B1_rR3", "B1_frR3", "B1_srs", "P2_alpha", "P1_sum", "B1_q10", "B1_pu", "P2_srs", "P1_srs", "P1_xqcn", "P2_xqcn", "P2_xqn", "P1_xqn"]   # parameters calibrated..
# 
# observed_types = ["ChlTot", "O2", "chl_f", "nit", "phos", "sil", "amm"]   # observed types of variables
# 
# model_types = ["total_chlorophyll_calculator_result_1", "O2_o", "total_chlorophyll_calculator_result_2", "N3_n", "N1_p", "N5_s", "N4_n"]  # model outputs corresponding to the observed variables
# 
# start_obs_index = 9*365+2   # since the observational files start on the 01/01/1998, to match up the model simulation period, the observations need to be extracted in a certain range. This is the start of the simulation period after spin up in number of days since 01/01/1998
# 
# length_of_data = 14*365-100  # length of simulation period in number of days
# 
# model_start = 365-31+4*365+1    # this is the start of simulation period after spin up 
# 
# n_ens_members = 5000 # number of model ensemble members


def run(conf):
    mpi.print(f"routine: run, conf: {conf.name}. Start!")
    
    path_observations = conf.path_observations
    model_directory = conf.model_directory
    perturbed_parameters_listed = conf.perturbed_parameters_listed
    observed_types = conf.observed_types
    model_types = conf.model_types
    start_obs_index = conf.start_obs_index
    length_of_data = conf.length_of_data
    model_start = conf.model_start
    n_ens_members = conf.n_ens_members
    
    mpi.print(f"routine: run, conf: {conf.name}. n_ens_members={n_ens_members}.")
    
    RMSE = np.zeros((n_ens_members//mpi.size + (n_ens_members%mpi.size >0)))  # the skill score is stored here

    calibration_init = cm.calibrate_model(length_period = length_of_data, obs_types = observed_types, path_obs = path_observations,  start_period_obs = start_obs_index, mask = 0.)     # initialize to get the observations out

    (observations, observations_depths, observations_times, observations_n_depths) = calibration_init.provide_observations() 

    # a loop to match ensemble members to the observations and calculate RMSE

    for rank_count, member in enumerate(range(1+mpi.rank,n_ens_members+1, mpi.size)): 
        
        print(member, flush=True)

        calibration_init = cm.calibrate_model(
            length_period = length_of_data,
            obs_types = observed_types,
            path_mod = model_directory+"/result_"+str(member).zfill(4)+".nc",
            mod_types = model_types,
            start_period_mod = model_start,
            observations = observations,
            observations_depths = observations_depths,
            observations_times = observations_times,
            n_depths_obs = observations_n_depths,
            )
        RMSE[rank_count] = calibration_init.RMSE_metric()
        
    # RMSE reconstruction at rank 0
    
    RMSEs=mpi.gather(RMSE)
    if mpi.rank!=0:
        return
    RMSE=np.array(RMSEs).transpose().reshape((-1))
    assert len(RMSE)>=n_ens_members
    assert RMSE[n_ens_members:].sum()==0
    RMSE=RMSE[:n_ens_members]
    
    # minimize RMSE to find the index of the best performing ensemble member

    best_ensemble_member = np.argwhere(RMSE==np.amin(RMSE))[0][0]+1  # identify the ensemble member number corresponding to minimum RMSE

    optimal_parameter_values = {}

    # identify and save the parameters belonging to the best performing member

    for parameter in perturbed_parameters_listed:
        with eatpy.models.gotm.YAMLEnsemble(model_directory+"/fabm_"+str(best_ensemble_member).zfill(4)+".yaml", 1) as fabm:
            parameter_value = fabm["instances/"+parameter[:2]+"/parameters/"+parameter[3:]]
            optimal_parameter_values.update({parameter:parameter_value})
            
            
    np.savetxt(conf.save_dir/f"{conf.name}_RMSE_all.txt", RMSE)
    
    with open(conf.save_dir/f"{conf.name}_best_parameters.csv", "w") as csv_file:
        w = csv.writer(csv_file)
        
        # loop over dictionary keys and values
        for key, val in optimal_parameter_values.items():
            # write every key and value to file
            w.writerow([key, val])
        
    mpi.print(f"routine: run, conf: {conf.name}. Done!")
        

def main():
    run(chosen_conf)

if __name__=="__main__":
    main()

     





