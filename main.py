from mpi import mpi
from run import run
from run_seasons import run_seasons
from run_obs_reductions import run_obs_reductions
from run_obs_perturbations import run_obs_perturbations
from run_pars_separately import run_pars_separately
from configurations import TestConf, L4

def main():
    # mpi.root = 0  # Uncomment to specify a preferred root rank.
    # Valid values are in the range [0, mpi.size); invalid values may raise ValueError.
    
    confs=[
        # TestConf(),
        L4(),
        ]
    
    for conf in confs:
        run(conf)
        run_seasons(conf)
        for i in range(1,5):
            run_obs_reductions(conf, obs_ratio=0.5**i)
        for noise_scale in [0.15, 0.275, 0.5]:
            run_obs_perturbations(conf, noise_scale=noise_scale)
        run_pars_separately(conf)

if __name__=="__main__":
    main()
     
