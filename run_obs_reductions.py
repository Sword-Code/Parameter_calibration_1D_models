# This is to optimize parameters across different observational types and a number of model ensemble members (with perturbed parameters). The parameters are optimized across a randomly selected subset of the available observations of size N (this is typically a fixed ratio of the original number of observations for each type, determined by the variable "obs_ratio"). The random selection of the sub-set of observations is done across 50 random realizations (set by the variable "n_red_ens_members"). Across those realizations a mean value of optimized parameter, as well as the standard deviation of the calibrated parameter values are calculated. This tells us about sensitivity of parameter calibration to spatio-temporal intermittency of observations (e.g. their number is vaguely related to obs sampling frequency). A lot of the code is analogous to run.py, so look for comments there.

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
# observed_types = ["ChlTot", "O2", "chl_f", "nit", "phos", "sil", "amm"]   # observed types of variables
# 
# model_types = ["total_chlorophyll_calculator_result_1", "O2_o", "total_chlorophyll_calculator_result_2", "N3_n", "N1_p", "N5_s", "N4_n"]  # model outputs corresponding to the observed variables
# 
# start_obs_index = 9*365+2  
# 
# length_of_data = 14*365-100  
# 
# model_start = 365-31+4*365+1    
# 
# n_mod_ens_members = 5000  # this is the number of model ensemble members throughought which the parameters are estimated
# 
# n_red_ens_members = 50  # this is the number of random realizations of sub-sampling of the observation data
# 
# obs_ratio = 0.25  # this is a fixed ratio between the reduced number of observations and the total number of observations for each type - it is predefined here


def run_obs_reductions(conf, obs_ratio=0.25):
    mpi.print(f"routine: run_obs_reductions, obs_ratio: {obs_ratio}, conf: {conf.name}. Start!")
    
    path_observations = conf.path_observations
    model_directory = conf.model_directory
    perturbed_parameters_listed = conf.perturbed_parameters_listed
    observed_types = conf.observed_types
    model_types = conf.model_types
    start_obs_index = conf.start_obs_index
    length_of_data = conf.length_of_data
    model_start = conf.model_start
    n_mod_ens_members = conf.n_ens_members
    n_red_ens_members = conf.n_red_ens_members
    
    mpi.print(f"routine: run_obs_reductions, obs_ratio: {obs_ratio}, conf: {conf.name}. n_mod_ens_members={n_mod_ens_members}.")
    
    
    n_RMSE=n_mod_ens_members//mpi.size + (n_mod_ens_members%mpi.size >0)
    RMSE = np.zeros((n_RMSE, n_red_ens_members))   # here the final metric is stored as the number of ensemble members x number of sub-samplings

    # read initial set of observations...

    init_obs = cm.calibrate_model(length_period = length_of_data, obs_types = observed_types, path_obs = path_observations,  start_period_obs = start_obs_index, mask = 0.)

    (observations_full, observations_depths_full, observations_times_full, observations_n_depths_full) = init_obs.provide_observations()   


    num_obs = {}  # these are total numbers of observations in the sub-sampled set for each type - to be determined
    reduced_obs={}  # this is the set of sub-sampled observations (of size N) randomly drawn in a n_red_ensemble_mem=50 member ensemble - the method is to randomly shuffle observations (indexes) and then take first N members 
    total_indexes={}  # these are corresponding shuffled indexes - they are stored, so corresponding model data can be derived for the sub-sampled observations

    # create a dictionary storing number of observations for each type in the reduced data-set 

    for observation_type in observed_types:
        num_obs.update({observation_type:int(obs_ratio*len(observations_times_full[observation_type]))})

    # here there is an initialization to generate the sub-sets of observations corresponding to num_obs

    init_subsamp = cm.calibrate_model(obs_types = observed_types, observations=observations_full, num_obs=num_obs)


    # in the loop below the observations are drawn across the n_red_ensemble_mem=50 realizations
    
    if mpi.rank==0:
        for reduction in range(0, n_red_ens_members):

            observations, indexes = init_subsamp.select_random_obs()
            reduced_obs.update({str(reduction):observations})
            total_indexes.update({str(reduction):indexes})

    # Broadcast only the shuffled indexes; each rank will reconstruct reduced_obs locally
    total_indexes = mpi.bcast(total_indexes)

    # Reconstruct reduced_obs on all ranks from observations_full and total_indexes
    reduced_obs = {}
    for reduction_key, indexes in total_indexes.items():
        obs_subset = {}
        for observation_type in observed_types:
            # Use only the first num_obs indices to match the reduced observation subset size
            obs_indices = indexes[observation_type][:num_obs[observation_type]]
            obs_subset[observation_type] = observations_full[observation_type][obs_indices]
        reduced_obs[reduction_key] = obs_subset
    
    # this is a loop through the model ensemble members, first extracting the full model data and then compiling an equivalent model set to the sub-sampled set of observations. It is done in two steps, so the netCDF files with model data are opened only once, which should save time..

    for rank_count, member in enumerate(range(1+mpi.rank,n_mod_ens_members+1, mpi.size)): 
        
        print(member, flush=True)                

        # read the full model data

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

        # run through the sub-samplings

        for reduction  in range(0,n_red_ens_members):
        
            # for each sub-sampling all the model-equivalent values to the sub-sampled observations are stored in "model"
        
            model = {}
        
            for observation_type in observed_types:
            
                model.update({observation_type:model_full[observation_type][total_indexes[str(reduction)][observation_type]][:num_obs[observation_type]]}) 
        
            # here the RMSE metric is calculated from the sub-sampled observations and model equivalents. It is then stored in the RMSE 2D array..
        
            calibration_init = cm.calibrate_model(obs_types = observed_types, model=model, mod_types = model_types, observations = reduced_obs[str(reduction)])       

            RMSE[rank_count, reduction] = calibration_init.RMSE_metric()

        
    # RMSE reconstruction at rank 0
    
    RMSEs=mpi.gather(RMSE)
    if mpi.rank!=0:
        return
    RMSE=np.array(RMSEs).transpose((1,0,2)).reshape((-1, n_red_ens_members))
    assert len(RMSE)>=n_mod_ens_members
    assert RMSE[n_mod_ens_members:].sum()==0
    RMSE=RMSE[:n_mod_ens_members]

    # derive the best performing model ensemble index and corresponding parameters for each sub-sampling 

    best_ensemble_member = {}
    parameters = np.zeros((len(perturbed_parameters_listed), n_red_ens_members))
            
    # this disctionary stores the mean and std across the sub-samplings and outputs in .csv file        
    optimal_parameter_values = {}        
            
    for reduction in range(0,n_red_ens_members):
        best_ensemble_member.update({str(reduction):np.argwhere(RMSE[:,reduction]==np.amin(RMSE[:,reduction]))[0][0]+1})  # identify the ensemble member number corresponding to minimum RMSE

        for par_index, parameter in enumerate(perturbed_parameters_listed):
            with eatpy.models.gotm.YAMLEnsemble(model_directory+"/fabm_"+str(best_ensemble_member[str(reduction)]).zfill(4)+".yaml", 1) as fabm:
                parameters[par_index, reduction] = fabm["instances/"+parameter[:2]+"/parameters/"+parameter[3:]]
                

    for par_index, parameter in enumerate(perturbed_parameters_listed):   
        optimal_parameter_values.update({parameter+"_mean":np.mean(parameters[par_index,:])})
        optimal_parameter_values.update({parameter+"_std":np.std(parameters[par_index,:])})
            

    with open(conf.save_dir/f"{conf.name}_best_parameters_reductions_{obs_ratio}.csv", "w") as csv_file:
        w = csv.writer(csv_file)
            
        # loop over dictionary keys and values
        for key, val in optimal_parameter_values.items():
            # write every key and value to file
            w.writerow([key, val])

 
    mpi.print(f"routine: run_obs_reductions, obs_ratio: {obs_ratio}, conf: {conf.name}. Done!")


def main():
    run_obs_reductions(chosen_conf)

if __name__=="__main__":
    main()
     

   
