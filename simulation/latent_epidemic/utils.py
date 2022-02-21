import os
import csv
import json
import shutil


def save_simulation_with_data(X, C, R, parameters, out_dir, base_data_dir="data"):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    shutil.copytree(base_data_dir, out_dir)

    X.to_csv(os.path.join(out_dir, "latent.csv"))
    C.to_csv(os.path.join(out_dir, "cases.csv"))
    R.to_csv(os.path.join(out_dir, "R.csv"))
    Rt = (
        R.reset_index()
        .melt(id_vars="area", var_name="Date", value_name="Rt")
        .sort_values(by=["area", "Date"])
    )
    Rt.update('"' + Rt[["area"]].astype(str) + '"')
    print(Rt)
    Rt.to_csv(os.path.join(out_dir, "Rt.csv"), index=False, quoting=csv.QUOTE_NONE)
    with open(os.path.join(out_dir, "parameters.json"), "w") as f:
        json.dump(parameters, f)
