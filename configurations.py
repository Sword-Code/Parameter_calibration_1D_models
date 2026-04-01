from pathlib import Path


class BaseConfiguration:
    name='BaseConfiguration'
    
    ###########################################
    # below you need to set-up as appropriate #
    ###########################################
    
    path_observations = "BATS_observations_w_chl_all_2sigma.nc"  #supply path to observational file. The observational file is assumed to be containing 1D, or 2D arrays with (time) or (time, depth) dimensions, where time is always labeled as number of days from 01/01/1998 (it is daily resolution) and depth is always spaced by 1m (so N vertical layers means going N meters deep). For simplification no fluctuations in the sea level height are considered...

    model_directory = "BATS_results"    # path to folder with the model ensemble simulations - the model outputs are picked across the ensemble from there
    
    output_dir = "outputs"  # path to outputs (csv and txt) folder 

    perturbed_parameters_listed = [
        'P3_p_q10',
        'B1_p_pu_ra',
        'P1_p_sum',
        'P2_p_q10',
        'P1_p_alpha_chl',
        'P3_p_sum',
        'Z5_p_minfood',
        'Z5_p_chuc',
        'P1_p_q10',
        'P1_p_qup',
        'P1_p_qlcPPY',
        'P3_p_alpha_chl',
        'P2_p_alpha_chl',
        'P1_p_xqp',
        'P3_p_qlcPPY',
        ]   # parameters calibrated..
    
    
    observed_types = [
        'Oxygen',
        'Nitrate',
        'Phosphate',
        'Silicate',
        'Satellite_chlorophyll',
        'Insitu_chlorophyll',
        ]   # observed types of variables
    
    
    _dict_model_variables = {
        'Oxygen': 'O2_o',
        'Nitrate': 'N3_n',
        'Phosphate': 'N1_p',
        'Silicate': 'N5_s',
        'Satellite_chlorophyll': 'total_chlorophyll',
        'Insitu_chlorophyll': 'total_chlorophyll',
        }  # dictionary observation_variables -> model_variables
    

    model_types = None  # model outputs corresponding to the observed variables. If None, it is filled automatically using _dict_model_variables.

    start_obs_index = 9*365+2   # since the observational files start on the 01/01/1998, to match up the model simulation period, the observations need to be extracted in a certain range. This is the start of the simulation period after spin up in number of days since 01/01/1998

    length_of_data = 14*365+4  # length of simulation period in number of days

    model_start = 365    # this is the start of simulation period after spin up 

    n_ens_members = None # number of model ensemble members used. If None, it automatically takes the maximum available number.
    
    
    # only for obs_reduction and obs_perturbation: 
    
    n_red_ens_members = 50  # this is the number of random realizations of sub-sampling of the observation data / of the perturbation noise
    
    ########################
    # end of user settings #
    ########################
    
    _save_dir=None
    
    @property
    def save_dir(self):
        if self._save_dir is None:
            self._save_dir=Path(self.output_dir)/self.name
            self._save_dir.mkdir(parents=True, exist_ok=True)
        return self._save_dir
    
    def __init__(self):
        if self.model_types is None:
            self.model_types=[self._dict_model_variables[observed_type] for observed_type in self.observed_types]
        
        if self.n_ens_members is None:
            path_model=Path(self.model_directory)
            if not path_model.exists():
                raise FileNotFoundError(
                    f"Model directory '{self.model_directory}' does not exist; "
                    "expected to find ensemble files matching 'result_????.nc'."
                )
            self.n_ens_members=len(list(path_model.glob('result_????.nc')))
            if self.n_ens_members <= 0:
                raise ValueError(
                    f"No ensemble members found in directory '{self.model_directory}'. "
                    "Expected at least one file matching 'result_????.nc'."
                )
            # print(f'Found {self.n_ens_members} in model directory: {path_model}')


class Boussole(BaseConfiguration):
    name='Boussole'
    
    path_observations = "/g100_work/OGS_test2528/sspada00/SEAMLESS/BOUSSOLE_observations_w_sat.nc"  #supply path to observational file. The observational file is assumed to be containing 1D, or 2D arrays with (time) or (time, depth) dimensions, where time is always labeled as number of days from 01/01/1998 (it is daily resolution) and depth is always spaced by 1m (so N vertical layers means going N meters deep). For simplification no fluctuations in the sea level height are considered...

    model_directory = "/g100_scratch/userexternal/ateruzzi/WP6_ms/BOUSSOLE/LARGE_ENSEMBLE_SIMULATIONS/BOUSSOLE_allparameters"    # path to folder with the model ensemble simulations - the model outputs are picked across the ensemble from there
    
    observed_types = [
        'Oxygen',
        'Nitrate',
        'Phosphate',
        'Silicate',
        'Satellite_chlorophyll',
        # 'Insitu_chlorophyll',
        ]   # observed types of variables


class BoussoleSatOnly(Boussole):
    name='BoussoleSatOnly'
    
    observed_types = [
        # 'Oxygen',
        # 'Nitrate',
        # 'Phosphate',
        # 'Silicate',
        'Satellite_chlorophyll',
        # 'Insitu_chlorophyll',
        ]   # observed types of variables


class BATS(BaseConfiguration):
    name='BATS'
    
    path_observations = "/g100_work/OGS_test2528/sspada00/SEAMLESS/BATS_observations_w_chl_all_2sigma.nc"  #supply path to observational file. The observational file is assumed to be containing 1D, or 2D arrays with (time) or (time, depth) dimensions, where time is always labeled as number of days from 01/01/1998 (it is daily resolution) and depth is always spaced by 1m (so N vertical layers means going N meters deep). For simplification no fluctuations in the sea level height are considered...

    model_directory = "/g100_scratch/userexternal/ateruzzi/WP6_ms/BATS/LARGE_ENSEMBLE_SIMULATIONS/BATS_allparameters"    # path to folder with the model ensemble simulations - the model outputs are picked across the ensemble from there
    
    observed_types = [
        'Oxygen',
        'Nitrate',
        'Phosphate',
        'Silicate',
        'Satellite_chlorophyll',
        'Insitu_chlorophyll',
        ]   # observed types of variables
    

class BATSSatOnly(BATS):
    name='BATSSatOnly'
    
    observed_types = [
        # 'Oxygen',
        # 'Nitrate',
        # 'Phosphate',
        # 'Silicate',
        'Satellite_chlorophyll',
        # 'Insitu_chlorophyll',
        ]   # observed types of variables
    
    
class L4(BaseConfiguration):
    name='L4'
    
    #below you need to set-up as appropriate

    path_observations = "./L4_observations.nc"  #supply path to observational file. The observational file is assumed to be containing 1D, or 2D arrays with (time) or (time, depth) dimensions, where time is always labeled as number of days from 01/01/1998 (it is daily resolution) and depth is always spaced by 1m (so N vertical layers means going N meters deep). For simplification no fluctuations in the sea level height are considered...

    model_directory = "/work/jos/eat/1D_configuration/L4/large_ensemble"    # path to folder with the model ensemble simulations - the model outputs are picked across the ensemble from there

    perturbed_parameters_listed = ["B1_rR2", "B1_sR1", "P1_alpha", "B1_rR3", "B1_frR3", "B1_srs", "P2_alpha", "P1_sum", "B1_q10", "B1_pu", "P2_srs", "P1_srs", "P1_xqcn", "P2_xqcn", "P2_xqn", "P1_xqn"]   # parameters calibrated..

    observed_types = ["ChlTot", "O2", "chl_f", "nit", "phos", "sil", "amm"]   # observed types of variables

    model_types = ["total_chlorophyll_calculator_result_1", "O2_o", "total_chlorophyll_calculator_result_2", "N3_n", "N1_p", "N5_s", "N4_n"]  # model outputs corresponding to the observed variables

    start_obs_index = 9*365+2   # since the observational files start on the 01/01/1998, to match up the model simulation period, the observations need to be extracted in a certain range. This is the start of the simulation period after spin up in number of days since 01/01/1998

    length_of_data = 14*365-100  # length of simulation period in number of days

    model_start = 365-31+4*365+1    # this is the start of simulation period after spin up 

    n_ens_members = 5000 # number of model ensemble members
    
    
    # only for obs_reduction and obs_perturbation: 
    
    n_red_ens_members = 50  # this is the number of random realizations of sub-sampling of the observation data / of the perturbation noise
    
    
class TestConf(BATS):
    name='TestConf'
    n_ens_members=2


chosen_conf=TestConf()
