import argparse
import numpy as np

def sample_Bi_incubation_period(num_samples):
    # Samples incubation period lognormal params according to Bi et al, appendix Table S2
    # https://www.medrxiv.org/content/10.1101/2020.03.03.20028423v3.full.pdf
    # meanlog parameter estimate 1.57 (1.44, 1.69) 95%CI and sdlog estimate 0.65 (0.56, 0.73) 95%CI
    meanlogs = np.random.normal(1.57, (1.69 - 1.44)/2/1.96, num_samples)
    sdlogs = np.random.normal(0.65, (0.73-0.56)/2/1.96, num_samples)

    meanlogs[meanlogs < 0] = 1.44
    sdlogs[sdlogs < 0] = 0.56

    return np.round(meanlogs, 3), np.round(sdlogs, 3)

def sample_Bi_serial_interval(num_samples):
    # Samples serial interval gamma parameters according to Bi et al, appendix Table S2
    # https://www.medrxiv.org/content/10.1101/2020.03.03.20028423v3.full.pdf

    # Shape param A has mean 2.29 and 95% CI (1.77, 3.34)
    # To symmetrise the CI we use sample log(A - 1.26) as Gaussian before inverting
    # This gives symmetric CI with mean 0.03 and 95% CI (-0.67,0.73)
    logA = np.random.normal(0.03, (0.73+0.67)/2/1.96, num_samples)
    A = np.exp(logA) + 1.26
    
    # Rate param B has mean 0.36 and 95% CI (0.26, 0.57), 
    # To symmetrise the CI we use sample log(B - 0.16) as Gaussian before inverting
    # This gives symmetric CI with mean -1.6 and 95% CI (-2.3, -0.9), base e
    logB = np.random.normal(-1.6, (-0.9+2.3)/2/1.96, num_samples)
    B = np.exp(logB) + 0.16

    return np.round(A, 3), np.round(B, 3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--save_dir",
        type=str,
        help="dir where bootstrap samples are saved",
    )
    parser.add_argument(
        "--num_samples",
        type=int,
        default=10,
    )
    args = parser.parse_args()
    n_samples = args.num_samples
    Aip, Bip = sample_Bi_serial_interval(n_samples)

    Adp, Bdp = sample_Bi_incubation_period(n_samples)

    fname = open(f"{args.save_dir}/bootstrap_params.txt", "w")
    for i in range(n_samples):
        fname.write(f"--Aip {Aip[i]} --Bip {Bip[i]} --Adp {Adp[i]} --Bdp {Bdp[i]}")
        fname.write("\n")    
