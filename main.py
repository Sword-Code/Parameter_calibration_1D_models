from run import run
from run_seasons import run_seasons
from run_obs_reductions import run_obs_reductions
from run_obs_perturbations import run_obs_perturbations
from configurations import TestConf, L4


def main():
    confs=[
        # TestConf(),
        L4(),
        ]
    
    for conf in confs:
        run(conf)
        run_seasons(conf)
        for i in range(1,5):
            run_obs_reductions(conf, obs_ratio=0.5**i)
        run_obs_perturbations(conf, noise_scale=0.275)

if __name__=="__main__":
    main()
     
