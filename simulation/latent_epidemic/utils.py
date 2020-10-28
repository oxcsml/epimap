import os
import shutil


def save_simulation_with_data(X, C, out_dir, base_data_dir="data"):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    shutil.copytree(base_data_dir, out_dir)

    X.to_csv(os.path.join(out_dir, "latent.csv"))
    C.to_csv(os.path.join(out_dir, "cases.csv"))
