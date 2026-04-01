# This is to optimize parameters across different observational types and a number of model ensemble members (with perturbed parameters). The impact of introducing additional white multiplicative Gaussian noise to the observations is tested. The noise is defined with a "noise_scale" parameter,which amounts to the standard deviation of the Gauss distribution. Then optimization across 50 different realizations of the noise is performed (set by the variable "n_red_ens_members") and the spread of the optimizations is evaluated with standard deviation, which is saved. A lot of the code is analogous to run.py and run_obs_reductions.py, so look for comments there.

import csv

import numpy as np

from mpi import mpi
import eatpy
import calibrate_models as cm
from configurations import chosen_conf


# #below you need to set-up as appropriate
# 
# path_observations = "./BATS_observations_w_chl_all_2sigma.nc"  
# 
# model_directory = "/work/jos/eat/1D_configuration/L4/large_ensemble"  
# 
# perturbed_parameters_listed = ["B1_rR2", "B1_sR1", "P1_alpha", "B1_rR3", "B1_frR3", "B1_srs", "P2_alpha", "P1_sum", "B1_q10", "B1_pu", "P2_srs", "P1_srs", "P1_xqcn", "P2_xqcn", "P2_xqn", "P1_xqn"] 
# 
# observed_types = ['Satellite_chlorophyll', 'Oxygen', 'Insitu_chlorophyll', 'Nitrate', 'Phosphate', 'Silicate']# ["ChlTot", "O2", "chl_f", "nit", "phos", "sil", "amm"]   # observed types of variables
# 
# 
# 
# model_types = ["total_chlorophyll_calculator_result_1", "O2_o", "total_chlorophyll_calculator_result_2", "N3_n", "N1_p", "N5_s"]#, "N4_n"]  # model outputs corresponding to the observed variables
# 
# start_obs_index = 9*365+2  
# 
# length_of_data = 14*365-100  
# 
# model_start = 365-31+4*365+1    
# 
# n_mod_ens_members = 5000  # this is the number of model ensemble members throughought which the parameters are estimated
# 
# n_red_ens_members = 50  # this is the number of random realizations of the noise
# 
# noise_scale = 0.275

def run_obs_perturbations(conf, noise_scale=0.275):
    mpi.print(f"routine: run_obs_perturbations, noise_scale: {noise_scale}, conf: {conf.name}. Start!")
    
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
    
    mpi.print(f"routine: run_obs_perturbations, noise_scale: {noise_scale}, conf: {conf.name}. n_mod_ens_members={n_mod_ens_members}.")

    n_RMSE=n_mod_ens_members//mpi.size + (n_mod_ens_members%mpi.size >0)
    RMSE = np.zeros((n_RMSE, n_red_ens_members))   # here the final metric is stored as the number of ensemble members x number of sub-samplings


    # read initial set of observations...

    init_obs = cm.calibrate_model(length_period = length_of_data, obs_types = observed_types, path_obs = path_observations,  start_period_obs = start_obs_index, mask = 1.e+20)

    (observations_full, observations_depths_full, observations_times_full, observations_n_depths_full) = init_obs.provide_observations()   

    obs_with_noise={}  # this is the set of observations with noise added 



    # create a dictionary storing number of observations for each type in the reduced data-set 
    randomization = None
    if mpi.rank==0:
        randomization=np.random.randint(10**6)
    randomization=mpi.bcast(randomization)
    for noise in range(0,n_red_ens_members):
        
        # Use a noise-index–seeded RNG so all MPI ranks get identical noise realizations
        rng = np.random.RandomState(noise+randomization)
        obs={}
        for observation_type in observed_types:
            noise_multipliers = rng.normal(
                loc=1.0,
                scale=noise_scale,
                size=len(observations_full[observation_type]),
            )
            obs[observation_type] = np.clip(
                observations_full[observation_type] * noise_multipliers,
                0,
                None,
            )
        obs_with_noise[str(noise)] = obs
    
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

        # run through the noise realizations

        for noise  in range(0,n_red_ens_members):    
        
            # here the RMSE metric is calculated from the sub-sampled observations and model equivalents. It is then stored in the RMSE 2D array..
        
            calibration_init = cm.calibrate_model(obs_types = observed_types, model=model_full, mod_types = model_types, observations = obs_with_noise[str(noise)])       

            RMSE[rank_count, noise] = calibration_init.RMSE_metric()
            
    # RMSE reconstruction at rank 0
    
    RMSEs=mpi.gather(RMSE)
    if mpi.rank!=0:
        return
    RMSE=np.array(RMSEs).transpose((1,0,2)).reshape((-1, n_red_ens_members))
    assert len(RMSE)>=n_mod_ens_members
    assert RMSE[n_mod_ens_members:].sum()==0
    RMSE=RMSE[:n_mod_ens_members]

    # derive the best performing model ensemble index and corresponding parameters for each noise realization 

    best_ensemble_member = {}
    parameters = np.zeros((len(perturbed_parameters_listed), n_red_ens_members))
            
    # this disctionary stores the mean and std across the noise realizations and outputs in .csv file        
    optimal_parameter_values = {}        
            
    for noise in range(0,n_red_ens_members):
        best_ensemble_member.update({str(noise):np.argwhere(RMSE[:,noise]==np.amin(RMSE[:,noise]))[0][0]+1})  # identify the ensemble member number corresponding to minimum RMSE

        for par_index, parameter in enumerate(perturbed_parameters_listed):
            with eatpy.models.gotm.YAMLEnsemble(model_directory+"/fabm_"+str(best_ensemble_member[str(noise)]).zfill(4)+".yaml", 1) as fabm:
                parameters[par_index, noise] = fabm["instances/"+parameter[:2]+"/parameters/"+parameter[3:]]
                

    for par_index, parameter in enumerate(perturbed_parameters_listed):   
        optimal_parameter_values.update({parameter+"_mean":np.mean(parameters[par_index,:])})
        optimal_parameter_values.update({parameter+"_std":np.std(parameters[par_index,:])})
            

    with open(conf.save_dir/f"{conf.name}_best_parameters_obs_noise_{noise_scale}.csv", "w") as csv_file:
        w = csv.writer(csv_file)
            
        # loop over dictionary keys and values
        for key, val in optimal_parameter_values.items():
            # write every key and value to file
            w.writerow([key, val])

    mpi.print(f"routine: run_obs_perturbations, noise_scale: {noise_scale}, conf: {conf.name}. Done!")

def main():
    run_obs_perturbations(chosen_conf)

if __name__=="__main__":
    main()
     

   
