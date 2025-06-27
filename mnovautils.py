import json
import pathlib
import pandas as pd
from shutil import copyfile
import datetime
import jinja2
import webbrowser

# from qtutils import warning_dialog


def warning_dialog(warning_message):
    print(warning_message)
    # read in template html file
    with open(r"html/warning_template.html", "r") as f:
        html_template = f.read()

    # create the html file
    environment = jinja2.Environment()

    template = environment.from_string(html_template)

    warning_html = template.render(warning_message=warning_message)

    warning_html_path = pathlib.Path("html", "warning_dialog.html")

    with open(warning_html_path, "w") as f:
        f.write(warning_html)

    ret = webbrowser.open("file://" + warning_html_path.absolute().__str__())


pulseSequences_to_experimentType = {
    "reset_psychetse.ptg": "H1_pureshift",
    "zg": "H1_1D",
    "cosygpppqf_ptype": "COSY",
    "zgpg30": "C13_1D",
    "hsqcetgpsisp2.2": "HSQC",
    "hmbcetgpl3nd.ptg": "HMBC",
    "hsqc_clip_cosy_mc_notation.eeh": "HSQC_CLIPCOSY",
    "hcdeptedetgpzf": "DDEPT_CH3_only",
    "noesygpphppzs": "NOESY",
}


def make_excel_backup(fn_excel: pathlib.Path) -> bool:

    ret = False
    if fn_excel.exists():
        fn_excel_bckup = (
            fn_excel.parent
            / "excel_backup"
            / (
                fn_excel.stem
                + "_"
                + datetime.datetime.now().strftime("%d%b%Y_%H%M%S")
                + ".xlsx"
            )
        )
        fn_excel_bckup.parents[0].mkdir(parents=True, exist_ok=True)
        ret = copyfile(fn_excel, fn_excel_bckup)
        if len(str(ret)) > 0:
            ret = True
    return ret


def return_nonempty_mnova_datasets(data: dict) -> dict:
    # remove datasets where multiplet counts is zero and peaks counts is zero and integrals counts is zero

    dicts_to_keep = {}
    for k, v in data.items():
        print(k)
        if k == "nmrAssignments":
            dicts_to_keep[k] = v

        elif isinstance(v, dict):
            if (
                (v["multiplets"]["count"] > 0)
                or (v["peaks"]["count"] > 0)
                or (v["integrals"]["count"] > 0)
            ):
                dicts_to_keep[k] = v
        else:
            dicts_to_keep[
                k
            ] = v  # must be smiles string, maybe need to check for this at some point

    return dicts_to_keep


def add_technique(key: str, technique_keys: dict, technique_counts: dict, expt: str):

    if expt in technique_keys.values():
        technique_counts[expt] += 1
        technique_keys[key] = f"{expt}_{technique_counts[expt]}"
    else:
        technique_keys[key] = expt


def read_in_mesrenova_json(fn: pathlib.Path) -> dict:
    """Read in the JSON file exported from MestReNova."""

    # check the type of fn is a pathlib.Path
    if isinstance(fn, pathlib.Path):

        with open(fn, "r") as file:
            data_orig = json.load(file)
        print(data_orig.keys())
    elif isinstance(fn, dict):
        data_orig = fn

    data = return_nonempty_mnova_datasets(data_orig)
    print(data.keys())

    # Identify the technique keys present in the JSON data
    technique_keys = {}
    technique_counts = {
        "HSQC": 0,
        "HMBC": 0,
        "COSY": 0,
        "C13_1D": 0,
        "H1_1D": 0,
        "NOESY": 0,
        "H1_pureshift": 0,
        "HSQC_CLIPCOSY": 0,
        "DDEPT_CH3_only": 0,
    }
    for key in data:
        if isinstance(data[key], dict):

            subtype = data[key].get("subtype", "")
            pulse_sequence = data[key].get("pulsesequence", "")

            if pulse_sequence in pulseSequences_to_experimentType:
                print("pulse_sequence found: ", pulse_sequence)
                add_technique(
                    key,
                    technique_keys,
                    technique_counts,
                    pulseSequences_to_experimentType[pulse_sequence],
                )
            else:

                if (
                    subtype.lower().find("hsqc") != -1
                    and data[key].get("type", "").lower() == "2d"
                ):
                    add_technique(key, technique_keys, technique_counts, "HSQC")

                elif (
                    subtype.lower().find("hmbc") != -1
                    and data[key].get("type", "").lower() == "2d"
                ):
                    add_technique(key, technique_keys, technique_counts, "HMBC")
                elif (
                    subtype.lower().find("cosy") != -1
                    and data[key].get("type", "").lower() == "2d"
                ):
                    add_technique(key, technique_keys, technique_counts, "COSY")
                elif (
                    subtype.lower().find("13c") != -1
                    and data[key].get("type", "").lower() == "1d"
                ):
                    add_technique(key, technique_keys, technique_counts, "C13_1D")
                elif (
                    subtype.lower().find("1h") != -1
                    and data[key].get("type", "").lower() == "1d"
                    and "psyche" in data[key].get("pulsesequence", "").lower()
                ):
                    add_technique(key, technique_keys, technique_counts, "H1_pureshift")
                elif (
                    subtype.lower().find("1h") != -1
                    and data[key].get("type", "").lower() == "1d"
                ):
                    add_technique(key, technique_keys, technique_counts, "H1_1D")
                elif (
                    subtype.lower().find("1h") != -1
                    and data[key].get("type", "").lower() == "2d"
                ):
                    add_technique(key, technique_keys, technique_counts, "NOESY")

    for k, v in technique_keys.items():
        data[v] = data[k]
        del data[k]

    return data


def get_2D_dataframe_from_json(json_data: dict, technique: str) -> pd.DataFrame:
    """
    Returns a pandas dataframe from the json_data dictionary for the specified technique.
    """
    print("technique: ", technique, json_data[technique]["peaks"]["count"])
    df_data = []
    for i in range(json_data[technique]["peaks"]["count"]):
        # check if the peak exists in the json data
        if str(i) not in json_data[technique]["peaks"]:
            print(f"{technique} :: Peak {i} not found in json data")
            continue
        df_data.append(
            [
                json_data[technique]["peaks"][str(i)]["delta2"],
                json_data[technique]["peaks"][str(i)]["delta1"],
                json_data[technique]["peaks"][str(i)]["intensity"],
                json_data[technique]["peaks"][str(i)]["type"],
            ]
        )

    df = pd.DataFrame(df_data, columns=["f2 (ppm)", "f1 (ppm)", "Intensity", "Type"])

    # sort the dataframe by f2 (ppm), descending order, reset the index and start the index at 1
    df = df.sort_values(by=["f2 (ppm)"], ascending=False).reset_index(drop=True)
    df.index += 1

    return df


def get_1d_dataframe_from_json(json_data: dict, technique: str) -> pd.DataFrame:
    df_data = []
    if json_data[technique]["multiplets"]["count"] == 0:
        # find peaks from  from peaks key
        for i in range(json_data[technique]["peaks"]["count"]):
            if str(i) in json_data[technique]["peaks"]:
                df_data.append(
                    [
                        json_data[technique]["peaks"][str(i)]["delta1"],
                        json_data[technique]["peaks"][str(i)]["intensity"],
                        json_data[technique]["peaks"][str(i)]["type"],
                    ]
                )

        df = pd.DataFrame(df_data, columns=["ppm", "Intensity", "Type"])

    else:
        # find peaks from  from multiplets key
        # Name	Shift	Range	H's	Integral	Class	J's	Method

        count = json_data[technique]["multiplets"]["count"]
        normValue = json_data[technique]["multiplets"]["normValue"]
        for i in [str(i) for i in range(count)]:
            if str(i) in json_data[technique]["multiplets"]:
                row = [
                    json_data[technique]["multiplets"][i]["delta1"],
                    json_data[technique]["multiplets"][i]["integralValue"],
                    json_data[technique]["multiplets"][i]["nH"],
                    json_data[technique]["multiplets"][i]["category"],
                ]

                # create a string from the list of J values and add it to df_data
                j_values = json_data[technique]["multiplets"][i]["jvals"]
                j_string = ", ".join([f"{j:1.3}" for j in j_values])
                j_string = f"{j_string}"
                row.append(j_string)
                df_data.append(row)

        df = pd.DataFrame(df_data, columns=["ppm", "Integral", "H's", "Class", "J's"])
        df["Integral"] = df["Integral"] / normValue

    # sort the dataframe by f2 (ppm), descending order, reset the index and start the index at 1
    df = df.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
    df.index += 1

    return df


def create_dataframes_from_mresnova_json(data: dict) -> dict:
    """
    Returns a dictionary of pandas dataframes for each technique in the data dictionary.
    """
    _dataframes = {}
    for k, v in data.items():
        print(k)
        if k in [
            "H1_1D",
            "C13_1D",
            "HSQC",
            "HMBC",
            "COSY",
            "NOESY",
            "H1_pureshift",
            "HSQC_CLIPCOSY",
            "DDEPT_CH3_only",
        ]:
            if v["type"].lower() == "2d":
                df = get_2D_dataframe_from_json(data, k)
                _dataframes[k] = df
            elif v["type"].lower() == "1d":
                df = get_1d_dataframe_from_json(data, k)
                _dataframes[k] = df
        elif k in ["smiles"]:
            _dataframes["molecule"] = pd.DataFrame([data["smiles"]], columns=["smiles"])

        elif k in ["nmrAssignments"]:
            _dataframes["nmrAssignments"] = pd.DataFrame.from_dict(
                data["nmrAssignments"], orient="index"
            )

    return _dataframes


def write_excel_file_from_mresnova_df_dict(
    df_frames: dict, excel_path: pathlib.Path, qtstarted: bool = False
) -> bool:

    # check if path is valid
    if not excel_path.parent.exists():
        return False
    else:
        try:
            with pd.ExcelWriter(excel_path) as writer:
                # check if path is valid
                for k, df in df_frames.items():
                    df.to_excel(writer, sheet_name=k)
        except OSError:
            # print exception
            warning_dialog(
                f"Exception occurred attempting to write excel file\n{str(excel_path)}"
            )
            return False
        return True

    # def check_for_multiple_HSQC_expts(data: dict) -> dict:
    #     """Check for multiple HSQC experiments and rename them."""

    #     hsqc_expts = [e for e in data.keys() if e.find("HSQC") != -1]
    #     print(hsqc_expts)

    #     if len(hsqc_expts) == 2:
    #         # compare the number of peaks in each experiment
    #         if (
    #             data[hsqc_expts[0]]["peaks"]["count"]
    #             >= data[hsqc_expts[1]]["peaks"]["count"]
    #         ):
    #             # rename the 2nd experiment
    #             data["HSQC_CH"] = data.pop(hsqc_expts[1])
    #         else:
    #             # rename the 1st experiment
    #             data["HSQC_CH"] = data.pop(hsqc_expts[0])
    #             # rename the 2nd experiment
    #             data["HSQC"] = data.pop(hsqc_expts[1])

    #     return data


if __name__ == "__main__":
    fn = pathlib.Path(
        "exampleProblems/EVB_330b/EVB_330b_annotated_assignments_mresnova.json"
    )
    fn = pathlib.Path(
        r"exampleProblems\2-ethyl-1-indanone_mnova_assign\2-ethyl-1-indanone_assignments_mresnova.json"
    )
    # fn = pathlib.Path( "exampleProblems/EVB_330b/EVB_330b_mresnova.json")
    print(fn, fn.exists())
    if fn.exists():
        data = read_in_mesrenova_json(fn)
        print(data.keys())
    else:
        print("File not found: ", fn)

    # print(data["nmrAssignments"])

    # df = pd.DataFrame.from_dict(data["nmrAssignments"], orient="index")
    # print(df)

    dataframes = create_dataframes_from_mresnova_json(data)
    print(dataframes.keys())
    print(dataframes["nmrAssignments"])

    # drop duplicates based on f1_ppm and f2_ppm1 column in the dataframe dataframes["nmrAssignments"]

    dataframes["nmrAssignments"] = dataframes["nmrAssignments"].drop_duplicates(
        subset=["f1_ppm", "f2_ppm1"], keep="first"
    )
    print(dataframes["nmrAssignments"])

    # sort the dataframe by f1_ppm, descending order, reset the index and start the index at 1
    # save the index first as mol_idx
    dataframes["nmrAssignments"]["mol_idx"] = dataframes["nmrAssignments"].index
    dataframes["nmrAssignments"] = (
        dataframes["nmrAssignments"]
        .sort_values(by=["f1_ppm"], ascending=False)
        .reset_index(drop=True)
    )
    dataframes["nmrAssignments"].index += 1
    print(dataframes["nmrAssignments"])

    print("\nC13_1D\n", dataframes["C13_1D"])
    print("\nH1_1D\n", dataframes["H1_1D"])
    print("\nHSQC\n", dataframes["HSQC"])
    print("\nHMBC\n", dataframes["HMBC"])

    # create pandas dataframe from "v"
    # print(data["H1_1D"]["peaks"]["count"])
    # print(data["C13_1D"]["peaks"]["count"])
    # print(data["HSQC"]["peaks"]["count"])
    # print(data["HMBC"]["peaks"]["count"])
    # print(data["COSY"]["peaks"]["count"])
    # print(data["NOESY"]["peaks"]["count"])
    # print(data["H1_pureshift"]["peaks"]["count"])
    # print(data["HSQC_CLIPCOSY"]["peaks"]["count"])
    # print(data["DDEPT_CH3_only"]["peaks"]["count"])
    # print(data["NOESY"]["peaks"]["count"])

    # df_frames = create_dataframes_from_mresnova_json(data)
    # print(df_frames.keys())
