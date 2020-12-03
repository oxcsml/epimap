import pandas as pd
import numpy as np

traffic_raw = pd.read_csv("data/uk_traffic_mergedflows.csv", index_col=0)
name_code = pd.read_csv("data/metadata.csv")
rad_flux = pd.read_csv("data/radiation_flux_ls=0.1.csv", index_col=0)
nhs_scot = pd.read_csv("data/nhs_scotland_health_boards.csv")

traffic = pd.DataFrame(0.0, index=rad_flux.index, columns=rad_flux.columns)

BUCKINGHAM_SUBREGIONS = {
    "_B4": "E07000004",
    "_B5": "E07000005",
    "_B6": "E07000006",
    "_B7": "E07000007",
}


def get_flow(code_from, code_to):
    flow = traffic_raw[
        (traffic_raw["From"] == code_from) & (traffic_raw["To"] == code_to)
    ]["Flow"]

    if len(flow) == 0:
        if (code_from not in list(traffic_raw["From"])) or (
            code_to not in list(traffic_raw["To"])
        ):
            print(f"Didn't find one of the codes at all in the table.") # shouldn't happen
        return 0
    else:
        return flow.sum()


def region2code(region):
    if region in BUCKINGHAM_SUBREGIONS.keys():
        return BUCKINGHAM_SUBREGIONS[region]
    elif region == "Cornwall and Isles of Scilly":
        return "E06000052,E06000053"
    elif region == "Hackney and City of London":
        return "E09000012"  # this is just Hackney
    elif region == "Westminster":
        return "E09000001,E09000033"  # westminster and city of london

    codes = name_code[name_code["AREA"] == region]["CODE"]

    if len(codes) == 1:
        return codes.item()
    else:
        assert region == "Highland" # should only happen for this?
        return list(codes)[0]


def get_subregions(region):
    if region in list(nhs_scot["NHS Scotland Health Board"]):
        subregions = list(
            nhs_scot[nhs_scot["NHS Scotland Health Board"] == region]["area"]
        )
    elif region == "Buckinghamshire":  # special case
        subregions = list(BUCKINGHAM_SUBREGIONS.keys())
    else:
        subregions = [region]

    return subregions


if __name__ == "__main__":
    NN = len(traffic.index) * len(traffic.columns)
    k = 0
    for i, region_from in enumerate(traffic.index):
        for j, region_to in enumerate(traffic.columns):

            subregions_from = get_subregions(region_from)
            subregions_to = get_subregions(region_to)

            flow = 0
            for subregion_from in subregions_from:
                for subregion_to in subregions_to:
                    code_from = region2code(subregion_from)
                    code_to = region2code(subregion_to)

                    flow += get_flow(code_from, code_to)

            traffic.loc[region_from, region_to] = flow

            k += 1

            if k % 1000 == 0:
                print(f"Done {k}/{NN} entries.")

    traffic.to_csv("data/uk_traffic.csv")
