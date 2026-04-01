# Code to derive metrics for model calibration in different sensitivity experiments, i.e. how calibrated parameters depend on observation data abundance, spatial location, season and the observed types. Jozef Skakala April 2025
 

from netCDF4 import Dataset
import numpy as np
import csv
# from matplotlib import pyplot as plt

def distance(a,b):
    return np.abs(a-b)


# function to read observations and store them in dictionaries as obs type vs numpy array. It separately stores the information on obs values, times and depths - this is then used to derive the corresponding model outputs for parameter calibration.

def read_observations(path_obs, observed_types, start_obs_index, end_obs_index, mask):

    observations = {}  # here the observations are stored as 1D arrays
    observations_depths = {}  # here the depths of observations are stored as 1D arrays
    observations_times = {}  # here the days of observations are stored as 1D arrays
    observations_n_data = {}  # number of observations for each type

    # obs = Dataset(path_obs)  # read in the observational file
    with Dataset(path_obs) as obs:

        for observed_type in observed_types:   # read in all the observations into "observations" dictionary
        
            if obs.variables[observed_type][:].ndim==2:  # if the variable is 2D (time x depth)
                observations_spec_inputs = obs.variables[observed_type][:]
                observations_spec_inputs = observations_spec_inputs[start_obs_index:end_obs_index,:] 
                relevant = (observations_spec_inputs != mask) & ~np.isnan(observations_spec_inputs) & observations_spec_inputs.mask==False
                observations_spec = observations_spec_inputs[relevant]  # flattened array with observations
                observations_spec_depth = np.argwhere(relevant)[:,1]  # flattened array with observational depths (adding 0.5 puts it in the centre of the 1m thick layer)  
                observations_spec_time = np.argwhere(relevant)[:,0]  # flattened array with times of observations
                n_depths_obs = observations_spec_inputs.shape[1]
    
            else:  # this means the variable is 1D (time)
                observations_spec_inputs = np.array(obs.variables[observed_type])[start_obs_index:end_obs_index]
                relevant = (observations_spec_inputs != mask) & ~np.isnan(observations_spec_inputs) 
                observations_spec = observations_spec_inputs[relevant]  # flattened array with observations 
                observations_spec_time = np.argwhere(relevant)[:,0]  # flattened array with times of observations
                observations_spec_depth = np.zeros((len(observations_spec_time))) #formally set "surface" data depths to 0.5m 
        
        # build the output dictionaries
            
            observations.update({observed_type:observations_spec})
            observations_depths.update({observed_type:observations_spec_depth}) 
            observations_times.update({observed_type:observations_spec_time}) 
            observations_n_data.update({observed_type:len(observations_spec)}) 

    # obs.close()
    
    
    return observations, observations_depths, observations_times, n_depths_obs #probably wrong here!
    # return observations, observations_depths, observations_times, observations_n_data # correction?


# match depth indexes between models and observations   
    
def match_depth_indexes(path_mod, n_depths_obs, model_start, model_end):

    # modinp = Dataset(path_mod)
    with Dataset(path_mod) as modinp:
        depths = np.abs(np.array(modinp.variables["z"])[model_start:model_end,:,0,0].mean(axis=0)) 
    n_depths_mod = len(depths)  # number of model vertical layers  
    indexes = np.zeros((n_depths_obs))  # records the model vertical index corresponding to each observational vertical layer
            
    for layer in range(0,n_depths_obs):
        indexes[layer] = np.argwhere(distance(layer+0.5,depths)==np.amin(distance(layer+0.5,depths)))[0][0]
        
    return indexes


# derive model variables corresponding to observations when the times and depths of observations are provided   
    
def match_model_with_observations(path_mod, model_types, observed_types, observations_depths, observations_times, model_start_index, model_end_index, indexes):

    # modinp = Dataset(path_mod)
    with Dataset(path_mod) as modinp:

        model_obs_equiv_output = {}

        for observed_type, model_type in zip(observed_types, model_types):   # run through all observed variables
            
            if ("total_chlorophyll_calculator_result" in model_type) and ("total_chlorophyll_calculator_result" not in modinp.variables.keys()):  # error?
            # if ("total_chlorophyll_calculator_result" == model_type) and ("total_chlorophyll_calculator_result" not in modinp.variables.keys()):  # is this the correct one?
                model = modinp.variables["P1_Chl"][:] + modinp.variables["P2_Chl"][:] + modinp.variables["P3_Chl"][:] + modinp.variables["P4_Chl"][:] # model data for specific type
            else:                       
                model = modinp.variables[model_type][:]  # model data for specific type
                
            model = model[model_start_index:model_end_index,:,0,0]

            observations_n_data = len(observations_times[observed_type])
                    
            model_obs_equiv = np.zeros((observations_n_data))  # this is the model output corresponding to the observational array  

            for datapoint in range(0, observations_n_data): # here you construct the model output that corresponds to the observations                
                                            
                model_obs_equiv[datapoint] = model[observations_times[observed_type][datapoint], int(indexes[int(observations_depths[observed_type][datapoint])])] 
            
            model_obs_equiv_output.update({observed_type : model_obs_equiv})
       
    return model_obs_equiv_output
   
   
   
# this is a core class for model calibration. It enables a number of functionalities and data manipulation, e.g. to provide observations, matching model outputs, calculate score metrics, derive seasonal calibration, or randomly reduce the observational data for sensitivity experiments.

class calibrate_model:

    def __init__(self, **kwargs):
    
        self.description = "Calibrate models using 2D formatted observations"
        self.author = "Jozef Skakala"
        
        if "length_period" in kwargs:
            self.length_period = kwargs["length_period"]
        if "obs_types" in kwargs:
            self.obs_types = kwargs["obs_types"]
        
        if "observations" in kwargs:
            self.observations = kwargs["observations"]
            if "observations_depths" in kwargs:
                self.observations_depths = kwargs["observations_depths"]
            if "observations_times" in kwargs:    
                self.observations_times = kwargs["observations_times"]
            if "n_depths_obs" in kwargs:
                self.n_depths_obs = kwargs["n_depths_obs"]
        else:        
            self.path_obs = kwargs["path_obs"]
            self.start_period_obs = kwargs["start_period_obs"]
            self.mask = kwargs["mask"]
            (self.observations, self.observations_depths, self.observations_times, self.n_depths_obs) = read_observations(self.path_obs, self.obs_types, self.start_period_obs, self.start_period_obs + self.length_period, self.mask)
            
            
        if "mod_types" in kwargs:                                        
            self.mod_types = kwargs["mod_types"]
            if "model" in kwargs:
                self.model_matching_obs = kwargs["model"]
            else:
                self.path_mod = kwargs["path_mod"]
                self.start_period_mod = kwargs["start_period_mod"]
                self.matching_depth_indexes = match_depth_indexes(self.path_mod, self.n_depths_obs, self.start_period_mod, self.start_period_mod + self.length_period)
                self.model_matching_obs = match_model_with_observations(self.path_mod, self.mod_types, self.obs_types, self.observations_depths, self.observations_times, self.start_period_mod, self.start_period_mod + self.length_period, self.matching_depth_indexes)
        if "num_obs" in kwargs:
            self.num_obs = kwargs["num_obs"]
        if "season" in kwargs:
            self.season = kwargs["season"]
        
    def provide_observations(self):   # outputs observations if you provide path to the file, obs types and the relevant period
        return self.observations, self.observations_depths, self.observations_times, self.n_depths_obs
        
    def provide_matching_model(self):  # matches model to the observations 
        return self.model_matching_obs
        
    def select_random_obs(self):  # selects a random subset of observations for each obs type - the size of the subset is predefined
        observations_out = {}
        indexes_out = {}
        for observed_type in self.obs_types:
            indexes_array = np.arange(0,len(self.observations[observed_type]))
            np.random.shuffle(indexes_array)
            observations_out.update({observed_type:self.observations[observed_type][indexes_array][:self.num_obs[observed_type]]})
            indexes_out.update({observed_type:indexes_array})
        return observations_out, indexes_out
        
    def select_seasons(self):  # selects seasonal data from the model and observations across 22 year period
        observations_out = {}
        model_out = {}
        time=[]
        if self.season == "winter":
            for year in range(0,22):
                time.append(np.arange(int(year*365.25), int(year*365.25) + 31 + 28))
                time.append(np.arange(int((year+1)*365.25)-31, int((year+1)*365.25)))
        elif self.season == "spring":
            for year in range(0,22):
                time.append(np.arange(int(year*365.25)+31+28, int(year*365.25)+31+28+31+30+31))
        elif self.season == "summer":
            for year in range(0,22):
                time.append(np.arange(int(year*365.25)+31+28+31+30+31, int(year*365.25)+31+28+31+30+31+30+31+31))
        elif self.season == "autumn":
            for year in range(0,22):
                time.append(np.arange(int(year*365.25)+31+28+31+30+31+30+31+31, int((year+1)*365.25)-31))
        time = np.concatenate(time)
        for observed_type in self.obs_types:
            observations_out.update({observed_type:self.observations[observed_type][np.isin(self.observations_times[observed_type], time)]})
            model_out.update({observed_type:self.model_matching_obs[observed_type][np.isin(self.observations_times[observed_type], time)]})
            
        return observations_out, model_out            
        
    def R_metric(self):  # calculates R (Pearson correlation) metric
        R_out = []
        for observed_type in self.obs_types:
             R_out.append(np.corrcoef(self.observations[observed_type], self.model_matching_obs[observed_type])[0][1])
        return np.mean(R_out)
        
    def RMSE_metric(self):  # calculates RMSE metric
        RMSE_out = []
        for observed_type in self.obs_types:
            RMSE_out.append(np.sqrt(np.mean((self.observations[observed_type] - self.model_matching_obs[observed_type])**2))/np.std(self.observations[observed_type]))
        return np.mean(RMSE_out)
        
        
    
    
    
    
    
