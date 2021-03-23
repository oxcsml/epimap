from joblib import Parallel, delayed
import queue
import os
import subprocess
from sklearn.model_selection import ParameterGrid
from pathlib import Path
import numpy as np
import time
import shutil
import sys
import itertools
import fire

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.
    https://stackoverflow.com/a/3041990/3160671
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).
    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

def dict_as_tuple(x):
    if isinstance(x, dict):
        return tuple(sorted(x.items()))
    elif isinstance(x, list):
        return [dict_as_tuple(d) for d in x]


def tuple_as_dict(x):
    if isinstance(x, tuple):
        return dict(x)
    elif isinstance(x, list):
        return [tuple_as_dict(d) for d in x]


def dict2args(dct):
    "Turns a dictionary of parameter values into a string of args"
    args_str = ""
    for param, val in dct.items():
        args_str += f"--{param} {val} "
    return args_str


# %%

# ---------------------------------------------------------------------------- #
#                               group comparison                               #
# ---------------------------------------------------------------------------- #

# Rscript covidmap/stage2_run.r $4 \
#             --region_id \$SLURM_ARRAY_TASK_ID \
#             && echo run: DONE
    # "#SBATCH --job-name=Rmap_regional_start_", 
    # "#SBATCH --output="$3"/regional/output/run_%A_%a.out", 
# DATA_NAME = "spring_dynamics"

RUN_NAME = "Rmap_backtest_run_regions" 
max_parallel_jobs = 100
script_cmd = f"Rscript covidmap/stage2_run.r "

sbatch_file_preamble = (
    "#SBATCH --mail-user=$USER@stats.ox.ac.uk",
    "#SBATCH --mail-type=ALL", 
    "#SBATCH --partition=ziz-large", 
    "#SBATCH --ntasks=1", 
    "#SBATCH --cpus-per-task=1", 
    "#SBATCH --mem-per-cpu=10G", 
)

# space_scale_list = ["1", "2"]
# first_day_modelled_list = ["0", "4"]

def get_grid(space_scale_list, first_day_modelled_list, backtest_directory="fits/backtests"):
    first_day_modelled_list = first_day_modelled_list.split(",") # list of strings
    space_scale_list = list(space_scale_list) # list of floats/integers

    param_ranges = {
        "approximation": ["regional"],
        "globalkernel": ["none"],
        "spatialkernel": ["matern12"],
        # "results_directory": results_directory_list,
        "fixed_gp_time_length_scale": ["100.0"],
        "fixed_gp_space_length_scale": space_scale_list,
        "first_day_modelled": first_day_modelled_list,
        "weeks_modelled": ["15"],
        "region_id": list(range(1, 10)),
    }

    grid = list(ParameterGrid(param_ranges))

    for g in grid:
        g["results_directory"] = f"{backtest_directory}/space_{g['fixed_gp_space_length_scale']}/start_{g['first_day_modelled']}_weeks_{g['weeks_modelled']}"
        if g["fixed_gp_space_length_scale"] == 0:
            del g["fixed_gp_space_length_scale"]
            g['spatialkernel'] = "none"

    return grid

grid = fire.Fire(get_grid)
num_jobs = len(grid)

# raise ValueError
# %%
# ---------------------------------------------------------------------------- #
#                   Scheduler (shouldn't need to touch this)                   #
# ---------------------------------------------------------------------------- #

# Make directory for logging
base_dir = f'./slurm/backtest_evaluation/outputs/{RUN_NAME}/{time.strftime("%d-%m-%Y-%H%M%S")}'
Path(base_dir).mkdir(parents=True, exist_ok=True)
Path(os.path.join(base_dir, "scripts")).mkdir(parents=True, exist_ok=True)

parallel_launcher_file = Path(__file__)

shutil.copyfile(
    parallel_launcher_file,
    os.path.join(base_dir, "scripts", os.path.basename(parallel_launcher_file)),
)

# Compose the sbatch file
sbatch_file_lines = [
    "#!/bin/bash",
    f"#SBATCH --job-name={RUN_NAME}",
    f"#SBATCH --array=0-{len(grid) - 1}%{max_parallel_jobs}",
    f"#SBATCH --wait", # don't exit until job finishes
]
sbatch_file_lines += list(sbatch_file_preamble)
sbatch_file_lines += [
    f"#SBATCH --output={Path(base_dir, 'job_%a', 'output.o')}",
    f"#SBATCH --error={Path(base_dir, 'job_%a', 'output.e')}",
]

for i, params in enumerate(grid):

    args_str = dict2args(params)
    full_cmd = script_cmd + args_str

    Path(base_dir, f"job_{i}").mkdir(parents=False, exist_ok=False)

    sbatch_file_lines += [
        f'commands[{i}]="{full_cmd}" '
        f"\n"
        f'logdirs[{i}]="{Path(base_dir, f"job_{i}").absolute()}"'
    ]


sbatch_file_lines += [
    "set -e",
    "cd /data/ziz/not-backed-up/scratch/szaidi/Rmap/",
    'echo "Running on $(hostname)"',
    'echo "$(hostname) ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} ${SLURM_JOB_ID} ${commands[$SLURM_ARRAY_TASK_ID]}" > '
    + "${logdirs[$SLURM_ARRAY_TASK_ID]}/slurm_jobid.txt",
    "set +e",
    "srun ${commands[$SLURM_ARRAY_TASK_ID]}",
]

# Save the sbatch file
sbatch_filename = Path(base_dir, "scripts", f"slurm_{RUN_NAME}.sh")
with open(sbatch_filename, "w") as file:
    file.write("\n".join(sbatch_file_lines))

print()
print(
    f"Running a total of {num_jobs} jobs using at most {max_parallel_jobs} GPUs in parallel:"
)
for parameters_dict in grid:
    print(parameters_dict)
    print()

if query_yes_no("Type 'yes' to execute the above runs, or type 'no' to abort."):
    # Execute the sbatch file
    os.system(f"sbatch {sbatch_filename}")
