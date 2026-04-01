# identifies how parameter optimization changes with season - much of the code has similar structure than run.py (so look for comments there), with some modifications..

import csv

import numpy as np

from mpi import mpi
import eatpy
import calibrate_models as cm
from configurations import chosen_conf


# #below you need to set-up as appropriate
# 
# path_observations = "./L4_observations.nc"  
# 
# model_directory = "/work/jos/eat/1D_configuration/L4/large_ensemble"
# 
# perturbed_parameters_listed = ["B1_rR2", "B1_sR1", "P1_alpha", "B1_rR3", "B1_frR3", "B1_srs", "P2_alpha", "P1_sum", "B1_q10", "B1_pu", "P2_srs", "P1_srs", "P1_xqcn", "P2_xqcn", "P2_xqn", "P1_xqn"] 
# 
# observed_types = ["ChlTot", "O2", "chl_f", "nit", "phos", "sil", "amm"]   
# 
# model_types = ["total_chlorophyll_calculator_result_1", "O2_o", "total_chlorophyll_calculator_result_2", "N3_n", "N1_p", "N5_s", "N4_n"]  # model outputs corresponding to the observed variables
# 
# start_obs_index = 9*365+2   
# 
# length_of_data = 14*365-100  
# 
# model_start = 365-31+4*365+1     
# 
# n_ens_members = 5000

def run_seasons(conf):
    mpi.print(f"routine: run_seasons, conf: {conf.name}. Start!")
    
    path_observations = conf.path_observations
    model_directory = conf.model_directory
    perturbed_parameters_listed = conf.perturbed_parameters_listed
    observed_types = conf.observed_types
    model_types = conf.model_types
    start_obs_index = conf.start_obs_index
    length_of_data = conf.length_of_data
    model_start = conf.model_start
    n_ens_members = conf.n_ens_members
    
    mpi.print(f"routine: run_seasons, conf: {conf.name}. n_ens_members={n_ens_members}.")
    
    # RMSE is defined as a dictionary for each season
    
    seasons = ["winter", "spring", "summer", "autumn"]
    n_RMSE=n_ens_members//mpi.size + (n_ens_members%mpi.size >0)
    RMSE = {season:np.zeros((n_RMSE)) for season in seasons}

    # this is to get the full observations across all seasons

    init_obs = cm.calibrate_model(length_period = length_of_data, obs_types = observed_types, path_obs = path_observations,  start_period_obs = start_obs_index, mask = 0.)

    (observations_full, observations_depths_full, observations_times_full, observations_n_depths_full) = init_obs.provide_observations() 

    # loop through model ensemble members
    
    for rank_count, member in enumerate(range(1+mpi.rank,n_ens_members+1, mpi.size)): 
        
        print(member, flush=True)     
        
        # this is to get the full model data matching the full observations (across all seasons)
        
        init_mod = cm.calibrate_model(
            length_period = length_of_data,
            obs_types = observed_types,
            path_mod = model_directory+"/result_"+str(member).zfill(4)+".nc",
            mod_types = model_types,
            start_period_mod = model_start,
            observations = observations_full,
            observations_depths = observations_depths_full,
            observations_times = observations_times_full,
            n_depths_obs = observations_n_depths_full,
            ) 
        model_full = init_mod.provide_matching_model()           

        # loop through seasons

        for season in seasons:
        
            # this is to select a specific season from obseervations and model data
        
            season_init = cm.calibrate_model(
                length_period = length_of_data,
                obs_types = observed_types,
                model=model_full,
                mod_types = model_types,
                observations = observations_full,
                observations_times = observations_times_full,
                season = season,
                )       
            (observations, model) = season_init.select_seasons()
            
            # use the season-specific data to calibrate the model and store in RMSE dictionary
            
            calibration_init = cm.calibrate_model(
                length_period = length_of_data,
                obs_types = observed_types,
                mod_types = model_types,
                model=model,
                observations = observations,
                )
            RMSE[season][rank_count] = calibration_init.RMSE_metric()
    
    # RMSE reconstruction at rank 0
    
    RMSEs=mpi.gather(RMSE)
    if mpi.rank!=0:
        return
    RMSE={season:np.array([RMSE_part[season] for RMSE_part in RMSEs]).transpose().reshape((-1)) for season in seasons}
    for season in seasons:
        assert len(RMSE[season])>=n_ens_members
        assert RMSE[season][n_ens_members:].sum()==0
        RMSE[season]=RMSE[season][:n_ens_members]

    best_ensemble_member = {}
    optimal_parameter_values = {}

    # run through seasons, identify the index of the best performing ensemble member and store the parameter values
            
    for season in seasons:
        best_ensemble_member.update({season:np.argwhere(RMSE[season]==np.amin(RMSE[season]))[0][0]+1})  # identify the ensemble member number corresponding to minimum RMSE

        for parameter in perturbed_parameters_listed:
            with eatpy.models.gotm.YAMLEnsemble(model_directory+"/fabm_"+str(best_ensemble_member[season]).zfill(4)+".yaml", 1) as fabm:
                parameter_value = fabm["instances/"+parameter[:2]+"/parameters/"+parameter[3:]]
                optimal_parameter_values.update({season+"_"+parameter:parameter_value})
            
            
        np.savetxt(conf.save_dir/f"{conf.name}_RMSE_{season}.txt", RMSE[season])

    with open(conf.save_dir/f"{conf.name}_best_parameters_seasons.csv", "w") as csv_file:
        w = csv.writer(csv_file)
            
        # loop over dictionary keys and values
        for key, val in optimal_parameter_values.items():
            # write every key and value to file
            w.writerow([key, val])
            
    mpi.print(f"routine: run_seasons, conf: {conf.name}. Done!")

def main():
    run_seasons(chosen_conf)

if __name__=="__main__":
    main()

     

   
