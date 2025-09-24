import webbrowser
import platform
import itertools

from pathlib import Path
import pandas as pd
import numpy as np
import json

from scipy import stats
import networkx as nx
from networkx.readwrite import json_graph
from rdkit import Chem

# from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import jinja2

from excelheaders import EXCEL_ORIG_DF_COLUMNS as excel_orig_df_columns

from globals import SVG_DIMENSIONS as svgDimensions

from globals import CARBONSEPARATION
from globals import PROTONSEPARATION
from globals import NMREXPERIMENTS


def find_nearest(true_values: list, value: float):
    arraynp = np.asarray(true_values)
    idx = (np.abs(arraynp - value)).argmin()
    return arraynp[idx]


def tidyup_ppm_values(
    df: pd.DataFrame, true_values: list, column_name: str, ppm_tolerance=0.005
) -> pd.DataFrame:

    # make a copy of the column_name adding a suffix orig
    df[f"{column_name}_orig"] = df[column_name]

    # make a probability column to see how far replacement is from original
    df[f"{column_name}_prob"] = 0

    # create dataframe with ppm values and their nearest true value
    # dfnew = df.assign(nearest_true_value = df[column_name].apply(lambda x: find_nearest(true_values, x)))
    df[column_name] = df[column_name].apply(lambda x: find_nearest(true_values, x))

    # calculate probabilities
    for idx in df.index:
        df.loc[idx, f"{column_name}_prob"] = stats.norm.pdf(
            df.loc[idx, column_name],
            loc=df.loc[idx, f"{column_name}_orig"],
            scale=ppm_tolerance,
        )

    return df


def return_nonempty_mnova_datasets(data: dict) -> dict:
    # remove datasets where multiplet counts is zero and peaks counts is zero and integrals counts is zero

    print("\nreturn_nonempty_mnova_datasets(data: dict) -> dict:\n")
    print("data.keys()\n", data.keys(), "\n")
    dicts_to_keep = {}
    for k, v in data.items():
        print(f"Processing key: {k}")
        print(f"\tdata[{k}].keys()\n", v.keys(), "\n")
        print(f"v.items():\n", v.items(), "\n")
        if v["datatype"] == "nmrspectrum":
            if (
                (v["multiplets"]["count"] > 0)
                or (v["peaks"]["count"] > 0)
                or (v["integrals"]["count"] > 0)
            ):
                dicts_to_keep[k] = v
        else:
            dicts_to_keep[k] = v

    return dicts_to_keep


def add_technique(key: str, technique_keys: dict, technique_counts: dict, expt: str):

    if expt in technique_keys.values():
        technique_counts[expt] += 1
        technique_keys[key] = f"{expt}_{technique_counts[expt]}"
    else:
        technique_keys[key] = expt


def read_in_mesrenova_json(fn: Path) -> dict:
    """Read in the JSON file exported from MestReNova."""

    print("\nread_in_mesrenova_json(fn: Path) -> dict:\n")

    # check the type of fn is a pathlib.Path
    if isinstance(fn, Path):
        with open(fn, "r") as file:
            data_orig = json.load(file)
    elif isinstance(fn, dict):
        data_orig = fn

    print("data_orig keys: ", data_orig.keys())

    data = return_nonempty_mnova_datasets(data_orig)

    # # create a dictionary of filename to experiment type ie  HSQC: filename, HMBC: filename, etc
    # chosen_spectra = { v.split(" ")[-1]: v.split(" ")[-2] for _,v in  data["chosenSpectra"]["data"].items()}
    # # replace the nmrdata key names with the experiment type ie HSQC, HMBC, etc
    # for k, v in chosen_spectra.items():
    #     data[k] = data[v]
    #     del data[v]

    return data


def read_in_mesrenova_json_multiple(fn: Path) -> dict:
    """Read in the JSON file exported from MestReNova."""

    print("\nread_in_mesrenova_json(fn: Path) -> dict:\n")

    # check the type of fn is a pathlib.Path
    if isinstance(fn, Path):
        with open(fn, "r") as file:
            data_orig = json.load(file)
    elif isinstance(fn, dict):
        data_orig = fn

    data = return_nonempty_mnova_datasets(data_orig)

    # # create a dictionary of filename to experiment type ie  HSQC: filename, HMBC: filename, etc
    # chosen_spectra = { v.split(" ")[-1]: v.split(" ")[-2] for _,v in  data["chosenSpectra"]["data"].items()}
    # # replace the nmrdata key names with the experiment type ie HSQC, HMBC, etc
    # for k, v in chosen_spectra.items():
    #     data[k] = data[v]
    #     del data[v]

    return data


def get_2D_dataframe_from_json(json_data: dict, technique: str) -> pd.DataFrame:
    """
    Returns a pandas dataframe from the json_data dictionary for the specified technique.
    """
    print("technique: ", technique)
    df_data = []

    signaltype = "peaks"

    df_data = []
    for i in range(json_data[technique][signaltype]["count"]):
        # check if the peak exists in the json data
        if str(i) not in json_data[technique][signaltype]["data"]:
            continue
        df_data.append(
            [
                json_data[technique]["peaks"]["data"][str(i)]["delta2"],
                json_data[technique]["peaks"]["data"][str(i)]["delta1"],
                json_data[technique]["peaks"]["data"][str(i)]["intensity"],
                json_data[technique]["peaks"]["data"][str(i)]["type"],
                json_data[technique]["peaks"]["data"][str(i)]["annotation"],
            ]
        )

    df = pd.DataFrame(
        df_data, columns=["f2 (ppm)", "f1 (ppm)", "Intensity", "Type", "Annotation"]
    )

    # sort the dataframe by f2 (ppm), descending order, reset the index and start the index at 1
    df = df.sort_values(by=["f2 (ppm)"], ascending=False).reset_index(drop=True)
    df.index += 1

    return df


def get_1d_dataframe_from_json(json_data: dict, technique: str) -> pd.DataFrame:
    df_data = []
    if json_data[technique]["multiplets"]["count"] == 0:
        if json_data[technique]["peaks"]["count"] == 0:
            print(f"no peaks found in {technique}")
            return pd.DataFrame()
        else:
            pks_data = json_data[technique]["peaks"]["data"]
            for pk in pks_data:
                df_data.append(
                    [
                        pks_data[pk]["delta1"],
                        pks_data[pk]["intensity"],
                        pks_data[pk]["type"],
                    ]
                )

        df = pd.DataFrame(df_data, columns=["ppm", "Intensity", "Type"])

        print(f"\nget_1d_dataframe_from_json {technique}:\n\n {df}\n")

    else:
        # find peaks from  from multiplets key
        # Name	Shift	Range	H's	Integral	Class	J's	Method

        count = json_data[technique]["multiplets"]["count"]
        normValue = json_data[technique]["multiplets"]["normValue"]
        for i in [str(i) for i in range(count)]:
            if str(i) in json_data[technique]["multiplets"]["data"]:
                row = [
                    json_data[technique]["multiplets"]["data"][i]["delta1"],
                    json_data[technique]["multiplets"]["data"][i]["integralValue"],
                    json_data[technique]["multiplets"]["data"][i]["nH"],
                    json_data[technique]["multiplets"]["data"][i]["category"],
                ]

                # create a string from the list of J values and add it to df_data
                j_values = json_data[technique]["multiplets"]["data"][i]["jvals"]
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
    print("\ncreate_dataframes_from_mresnova_json(data: dict) -> dict:\n")
    print("***** EEH ****")
    _dataframes = {}
    for k, v in data.items():
        print(f"Processing key: {k}")
        # split the key by "_" and check if the first part is in NMREXPERIMENTS
        k1 = k.split("_")
        if len(k1) == 2:
            k1 = k1[0]
        elif len(k1) == 3:
            k1 = k1[0] + "_" + k1[1]
        elif len(k1) == 4:
            k1 = k1[0] + "_" + k1[1] + "_" + k1[2]
        else:
            k1 = k1[0]

        print(f"k1: {k1}")
        if k1 in NMREXPERIMENTS:
            if v["type"].lower() == "2d":
                df = get_2D_dataframe_from_json(data, k)
                _dataframes[k] = df
            elif v["type"].lower() == "1d":
                df = get_1d_dataframe_from_json(data, k)
                _dataframes[k] = df

        

        elif k1 in [
            "allAtomsInfo",
            "nmrAssignments",
            "carbonAtomsInfo",
            "c13predictions",
        ]:
            _dataframes[k1] = pd.DataFrame.from_dict(data[k1]["data"], orient="index")

        elif k1 in [
            "smiles",
            "molfile",
            "carbonCalcPositionsMethod",
            "MNOVAcalcMethod",
            "workingDirectory",
            "workingFilename",
            "simulatedAnnealing",
            "randomizeStart",
            "startingTemperature",
            "endingTemperature",
            "coolingRate",
            "numberOfSteps",
            "ppmGroupSeparation",
        ]:
            _dataframes[k1] = pd.DataFrame(v["data"], index=[0])
            _dataframes[k1].columns = [k1]

        elif k1 in ["spectraWithPeaks", "chosenSpectra"]:
            _dataframes[k1] = pd.DataFrame.from_dict(v["data"], orient="index")
            _dataframes[k1].columns = [k1]

            if k1 == "chosenSpectra":
                expts = [s.split()[-1] for s in v["data"].values()]
                filenames = [s.split()[-2] for s in v["data"].values()]
                _dataframes[k1]["filename"] = filenames
                _dataframes[k1]["expt"] = expts

                _dataframes[k1]["NMRdimensions"] = "unknown"
                for idx, row in _dataframes[k1].iterrows():
                    kwds = row[k1].split()
                    if len(kwds) == 5:
                        nmr_nuclei = kwds[0]
                    else:
                        nmr_nuclei = "".join(kwds[:2])

                    if len(nmr_nuclei) in [2, 3]:
                        expt_type = "1D"
                    elif len(nmr_nuclei) in [7, 8, 9]:
                        expt_type = "2D"
                    else:
                        expt_type = "unknown"

                    _dataframes[k1].loc[idx, "NMRdimensions"] = expt_type
            else:
                filenames = [s.split()[-1] for s in v["data"].values()]
                _dataframes[k1]["filename"] = filenames
        else:
            print(f"Skipping {k1}")

    print("_dataframes keys: ", _dataframes.keys())

    return _dataframes


def standardize_column_headings(dataframes: dict):
    for k in dataframes.keys():
        dataframes[k].rename(
            columns={
                "H's": "numProtons",
                "Integral": "integral",
                "J's": "jCouplingVals",
                "Class": "jCouplingClass",
                "Intensity": "intensity",
                "Shift": "ppm",
                "Range": "range",
                "f2 (ppm)": "f2_ppm",
                "f1 (ppm)": "f1_ppm",
                "Type": "signaltype",
            },
            inplace=True,
        )

    return dataframes


def add_missing_columns_to_nmrExpt_dataframes(dataframes: dict):

    nmrExpts = NMREXPERIMENTS

    for k in dataframes:
        if k in nmrExpts:
            if "Type" in dataframes[k]:
                continue
            if "signaltype" not in dataframes[k]:
                dataframes[k]["signaltype"] = 0

    return dataframes


def attach_symmetry_atom_index(nmrAssignments: pd.DataFrame) -> pd.DataFrame:
    nmrAssignments["sym_atom_idx"] = ""
    nmrAssignments["sym_atomNumber"] = ""

    for idx, row in nmrAssignments.iterrows():
        df = nmrAssignments[nmrAssignments.f1_ppm == row["f1_ppm"]]

        if len(df) > 1:
            # remove row from df that has the same index as the current row
            df = df.drop(idx)

            for idx2, row2 in df.iterrows():
                if nmrAssignments.at[idx2, "sym_atom_idx"] == "":
                    nmrAssignments.at[idx2, "sym_atom_idx"] = f"{row['atom_idx']}"
                    nmrAssignments.at[idx2, "sym_atomNumber"] = f"{row['atomNumber']}"
                else:
                    nmrAssignments.at[
                        idx2, "sym_atom_idx"
                    ] = f"{nmrAssignments.at[idx2, 'sym_atom_idx']}, {row['atom_idx']}"
                    nmrAssignments.at[
                        idx2, "sym_atomNumber"
                    ] = f"{nmrAssignments.at[idx2, 'sym_atomNumber']}, {row['atomNumber']}"

    return nmrAssignments


def add_jCouplingVals_jCouplingClass_to_h1(
    h1: pd.DataFrame, h1_1D: pd.DataFrame
) -> pd.DataFrame:
    """
    Add jCouplingVals and jCouplingClass to h1 from h1_1D
    """
    h1["jCouplingClass"] = ""
    h1["jCouplingVals"] = ""

    if len(h1_1D) == len(h1):
        # sort the dataframes by ppm in descending order
        h1 = h1.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
        h1_1D = h1_1D.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)

        h1["jCouplingVals"] = h1_1D["jCouplingVals"]
        h1["jCouplingClass"] = h1_1D["jCouplingClass"]
    else:
        # loop through h1 and find the row in h1_1D with the same ppm value
        for idx, row in h1.iterrows():
            # find the row in h1_1D with the same ppm value
            h1_1D_row = h1_1D[h1_1D["ppm"] == row["ppm"]]
            if len(h1_1D_row) > 0:
                h1.loc[idx, "jCouplingVals"] = h1_1D_row["jCouplingVals"].values[0]
                h1.loc[idx, "jCouplingClass"] = h1_1D_row["jCouplingClass"].values[0]

    return h1


def add_jCouplingVals_jCouplingClass_to_nmrAssignments(
    nmrAssignments: pd.DataFrame, h1: pd.DataFrame
) -> pd.DataFrame:
    """
    Add jCouplingVals and jCouplingClass to nmrAssignments from h1
    """
    nmrAssignments["jCouplingClass"] = ""
    nmrAssignments["jCouplingVals"] = ""

    for idx, row in nmrAssignments.iterrows():
        # search h1_1d based on mol_idx
        df = h1[h1.atom_idx == row["atom_idx"]]
        # remove rows where jCouplingVals is empty or ""
        df = df[df["jCouplingVals"] != ""]

        if len(df) == 0:
            continue
        if len(df) > 1:
            jcouplingVals = " & ".join(df["jCouplingVals"].values)
        elif len(df) == 1:
            jcouplingVals = df["jCouplingVals"].values[0]

        nmrAssignments.at[idx, "jCouplingVals"] = jcouplingVals

    # repeat for jCouplingClass
    for idx, row in nmrAssignments.iterrows():
        # search h1_1d based on mol_idx
        df = h1[h1.atom_idx == row["atom_idx"]]
        # remove rows where jCouplingClass is empty or ""
        df = df[df["jCouplingClass"] != ""]

        if len(df) == 0:
            continue
        if len(df) > 1:
            jcouplingClass = " & ".join(df["jCouplingClass"].values)
        elif len(df) == 1:
            jcouplingClass = df["jCouplingClass"].values[0]

        nmrAssignments.at[idx, "jCouplingClass"] = jcouplingClass

    return nmrAssignments


def create_svg_string(mol, molWidth=1000, molHeight=600, svgWidth=1200, svgHeight=700):

    translateWidth = int((svgWidth - molWidth) / 2)
    translateHeight = int((svgHeight - molHeight) / 2)

    # AllChem.Compute2DCoords(mol)

    d2d = Draw.rdMolDraw2D.MolDraw2DSVG(molWidth, molHeight)
    d2d.drawOptions().minFontSize = 18

    d2d.DrawMolecule(mol)
    d2d.TagAtoms(mol)
    d2d.FinishDrawing()
    sss = d2d.GetDrawingText().replace(
        f"width='{molWidth}px' height='{molHeight}px'",
        f"width={molWidth} height={molHeight}",
    )
    sss = sss.replace("fill:#FFFFFF", "fill:none").replace(
        "<svg", '<svg class="center"'
    )

    sss = sss.replace(
        f"<!-- END OF HEADER -->",
        f"<!-- END OF HEADER -->\n<g transform='translate({translateWidth}, {translateHeight})'>",
    )
    sss = sss.replace("</svg>", "</g>\n</svg>")
    sss = sss.replace(
        f"width={molWidth} height={molHeight} viewBox='0 0 {molWidth} {molHeight}'",
        f"width='78%' height='100%' viewBox='0 0 {svgWidth} {svgHeight}'",
    )

    idx_list = []
    xxx = []
    yyy = []

    new_xy3 = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C":
            idx = atom.GetIdx()
            point = d2d.GetDrawCoords(idx)
            idx_list.append(idx)
            xxx.append(point.x / molWidth)
            yyy.append(point.y / molHeight)
            new_xy3[idx] = (point.x / molWidth, point.y / molHeight)

    return sss, new_xy3


def h1_from_nmrAssignments(nmrAssignments: pd.DataFrame) -> pd.DataFrame:
    """_summary_

    Args:
        nmrAssignments (_type_): _description_
    """
    f2_ppm1 = nmrAssignments[nmrAssignments.f2_ppm1 != "undefined"]
    f2_ppm2 = nmrAssignments[nmrAssignments.f2_ppm2 != "undefined"]
    f2_ppm3 = nmrAssignments[nmrAssignments.f2_ppm3 != "undefined"]  # not used

    # loop through f2_ppm2 and if the value of in column f2_ppm

    rows_to_add = {}
    for idx, row in f2_ppm2.iterrows():
        # if row["f2_ppm2"] =="undefined" then continue
        if row["f2_ppm2"] == "undefined":
            continue

        if row["f2_ppm2"] != f2_ppm1.loc[idx, "f2_ppm1"]:
            # #swap columns f2_ppm1 with f2_ppm2 in row
            row["f2_ppm1"] = row["f2_ppm2"]
            rows_to_add[idx] = row

    f2_ppm1 = pd.concat([f2_ppm1, pd.DataFrame(rows_to_add).T], ignore_index=True)

    h1 = f2_ppm1.sort_values(by=["f2_ppm1"], ascending=False).reset_index(drop=True)

    # change column name f2_ppm1 to ppm
    h1 = h1.rename(columns={"f2_ppm1": "ppm"})
    # drop columns f2_ppm2 and f2_ppm3
    h1 = h1.drop(columns=["f2_ppm2", "f2_ppm3"])

    # set signaltype to 0
    h1["signaltype"] = 0

    return h1


def orientation(p, q, r):
    val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
    if val == 0:
        return 0  # Collinear
    elif val > 0:
        return 1  # Clockwise
    else:
        return 2  # Counterclockwise


def on_segment(p, q, r):
    if min(p[0], r[0]) <= q[0] <= max(p[0], r[0]) and min(p[1], r[1]) <= q[1] <= max(
        p[1], r[1]
    ):
        return True
    return False


def do_intersect(p1, q1, p2, q2):
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if o1 != o2 and o3 != o4:
        return True

    if o1 == 0 and on_segment(p1, p2, q1):
        return True

    if o2 == 0 and on_segment(p1, q2, q1):
        return True

    if o3 == 0 and on_segment(p2, p1, q2):
        return True

    if o4 == 0 and on_segment(p2, q1, q2):
        return True

    return False


def create_network_graph(c13, h1):

    # create nodes of graph
    G2 = nx.Graph()
    # for i in c13.mol_idx:
    #     print(i)
    #     G2.add_node(i)
    #     G2.nodes[i]["symbol"] = "C"

    for (
        i,
        ppm,
        ppm_calculated,
        numProtons,
        atomNumber,
        iupacLabel,
        jcouplingVals,
        jcouplingClass,
        x,
        y,
    ) in zip(
        c13.atom_idx,
        c13.f1_ppm,
        c13.ppm_calculated,
        c13.numProtons,
        c13.atomNumber,
        c13.iupacLabel,
        c13.jCouplingVals,
        c13.jCouplingClass,
        c13.x,
        c13.y,
    ):
        idx = i
        G2.add_node(idx)
        G2.nodes[idx]["symbol"] = "C"
        G2.nodes[idx]["atomNumber"] = atomNumber
        G2.nodes[idx]["ppm"] = ppm
        G2.nodes[idx]["ppm_calculated"] = ppm_calculated
        G2.nodes[idx]["numProtons"] = numProtons
        G2.nodes[idx]["iupacLabel"] = iupacLabel
        G2.nodes[idx]["jCouplingVals"] = jcouplingVals
        G2.nodes[idx]["jCouplingClass"] = jcouplingClass
        G2.nodes[idx]["x"] = x
        G2.nodes[idx]["y"] = y

    # Add proton ppm information to the graph nodes
    for c in h1.atom_idx.unique():
        idx = c
        G2.nodes[idx]["H1_ppm"] = list(h1[h1.atom_idx == c].ppm.values)

    return G2


def create_cosy_combinations(cosy: pd.DataFrame, h1: pd.DataFrame) -> dict:
    cosy_combinations = {}
    for i, row in cosy.iterrows():
        # find the rows in h1 that have the same ppm vaslue as the f1_ppm value in the row

        cosy_node = [row.f1_ppm, row.f2_ppm]
        cosy_node.sort()
        f1_ppm = cosy_node[0]
        f2_ppm = cosy_node[1]

        cosy_combinations[(f1_ppm, f2_ppm)] = []

        df1 = h1[h1.ppm == f1_ppm]
        df2 = h1[h1.ppm == f2_ppm]

        for j, row1 in df1.iterrows():
            for k, row2 in df2.iterrows():
                cosy_combinations[(f1_ppm, f2_ppm)].append(
                    (
                        row1.atom_idx,
                        row2.atom_idx,
                        row1.sym_atom_idx,
                        row2.sym_atom_idx,
                        ((row1.x, row1.y), (row2.x, row2.y)),
                    )
                )

    return cosy_combinations


def create_clipcosy_combinations(cosy: pd.DataFrame, h1: pd.DataFrame) -> dict:
    cosy_combinations = {}
    for i, row in cosy.iterrows():
        # find the rows in h1 that have the same ppm vaslue as the f1_ppm value in the row
        cosy_node = [row.f1_ppm, row.f2_ppm]
        cosy_node.sort()
        f1_ppm = cosy_node[0]
        f2_ppm = cosy_node[1]  # hopefully after sorting this is still the C13 ppm value

        cosy_combinations[(f1_ppm, f2_ppm)] = []

        df1 = h1[h1.ppm == f1_ppm]
        df2 = h1[h1.f1_ppm == f2_ppm]

        for j, row1 in df1.iterrows():
            for k, row2 in df2.iterrows():
                cosy_combinations[(f1_ppm, f2_ppm)].append(
                    (
                        row1.atom_idx,
                        row2.atom_idx,
                        row1.sym_atom_idx,
                        row2.sym_atom_idx,
                        ((row1.x, row1.y), (row2.x, row2.y)),
                    )
                )

    return cosy_combinations


def remove_duplicate_cosy_combinations(cosy_combinations: dict) -> dict:
    # remove any cosy edges that cross in groups of 4
    keep = {}

    for k, v in cosy_combinations.items():
        keep[k] = []
        if len(v) == 4:
            for ipos in itertools.combinations([0, 1, 2, 3], 2):
                if do_intersect(
                    v[ipos[0]][-1][0],
                    v[ipos[0]][-1][1],
                    v[ipos[1]][-1][0],
                    v[ipos[1]][-1][1],
                ):
                    continue
                else:
                    # keep the combination
                    keep[k].append(v[ipos[0]])
                    keep[k].append(v[ipos[1]])
        else:
            keep[k] = v

    return keep


def add_cosy_edges_to_graph(G2: nx.Graph, keep: dict) -> nx.Graph:

    # add cosy edges to graph
    for k, v in keep.items():
        for vv in v:
            node1 = vv[0]
            node2 = vv[1]
            if G2.has_edge(node1, node2):
                # set cosy attribute to True
                G2.get_edge_data(node1, node2)["cosy"] = True
            else:
                G2.add_edge(node1, node2)
                # set cosy attribute to True
                G2.get_edge_data(node1, node2)["cosy"] = True

    return G2


def add_clipcosy_edges_to_graph(G2: nx.Graph, keep: dict) -> nx.Graph:

    # add cosy edges to graph
    for k, v in keep.items():
        for vv in v:
            node1 = vv[0]
            node2 = vv[1]
            if G2.has_edge(node1, node2):
                # set cosy attribute to True
                G2.get_edge_data(node1, node2)["cosy"] = True
            else:
                G2.add_edge(node1, node2)
                # set cosy attribute to True
                G2.get_edge_data(node1, node2)["cosy"] = True

    return G2


def add_hmbc_edges_to_graph(
    G2: nx.Graph, hmbc: pd.DataFrame, h1: pd.DataFrame, c13: pd.DataFrame
) -> nx.Graph:
    # add HMBC edges to G2
    print("Adding HMBC edges 2")
    count = 0
    for i, row in hmbc.iterrows():
        # find the rows in h1 that have the same ppm vaslue as the f1_ppm value in the row
        df1 = h1[h1.ppm == row.f2_ppm]
        df2 = c13[c13.f1_ppm == row.f1_ppm]

        for j, row1 in df1.iterrows():
            for k, row2 in df2.iterrows():
                e1 = row1.atom_idx
                e2 = row2.atom_idx
                ppm1 = row1.f1_ppm
                ppm2 = row2.f1_ppm
                if (e1 == e2) or (ppm1 == ppm2):
                    continue
                if not G2.has_edge(e1, e2):
                    G2.add_edge(e1, e2, weight=1)
                    count += 1
                G2.get_edge_data(e1, e2)["hmbc"] = True
    return G2


def create_htmlpage_from_graph(
    G2: nx.Graph, fn_html: str, fn_html_template: str, jinjadata: dict
):

    # open html template file
    html_fn = Path(fn_html_template)
    with open(html_fn, "r") as f:
        html_template = f.read()

    environment = jinja2.Environment()
    template = environment.from_string(html_template)
    htmlpage = template.render(jinjadata)

    # replace True with true and False with false in the htmlpage
    htmlpage = htmlpage.replace("True", "true")
    htmlpage = htmlpage.replace("False", "false")

    # write the html page to a file
    html_fn = Path(fn_html)
    with open(html_fn, "w") as f:
        f.write(htmlpage)

    # open the html page in a browser
    webbrowser.open(html_fn.absolute().as_uri())





def calc_minimum_ppm_separation(ppm_pks, ppmSeparation):

    if len(ppm_pks) < 2:
        return ppmSeparation
    # Calculate pairwise differences
    differences = np.abs(ppm_pks[:, None] - ppm_pks)
    np.fill_diagonal(differences, np.inf)  # Ignore self-comparisons

    # Find the indices of the minimum difference
    i, j = np.unravel_index(np.argmin(differences), differences.shape)

    # Get the two closest values
    closest_values = (ppm_pks[i], ppm_pks[j])

    print("The two closest values are:", closest_values)

    value1, value2 = closest_values
    scale = 1.0  # Start with an initial scale value
    while True:
        pdf_result = stats.norm.pdf(value1, loc=value2, scale=scale)
        if pdf_result == 0.0:
            break
        scale /= 2  # Decrease the scale value

    if scale < ppmSeparation:
        return scale
    else:
        return ppmSeparation

def combine_multiple_nmrExpt_dataframes(dataframes):
    """
    Combine multiple NMR experiment dataframes into a single dataframe.
    This is useful when there are multiple 1D or 2D spectra for the same experiment.
    """

    print("\ncombine_multiple_nmrExpt_dataframes(dataframes: dict) -> dict:\n")

    print("dataframes keys: ", dataframes.keys())

    # group the expt dataframes by their experiment type
    nmrExptsDataframes = {}
    nmrExptsToDelete = []
    for k, v in dataframes.items():
        # split the key by "_" once from the right
        print(f"Processing key: {k}")
        nameSplit = k.rsplit("_", 1)
        k_number = ""
        k1 = nameSplit[0]
        if len(nameSplit) == 2:
            k_number = nameSplit[1]
        if k1 in NMREXPERIMENTS:
            print( f"{k} shape: {v.shape}"  )
            if k1 not in nmrExptsDataframes:
                nmrExptsDataframes[k1] = []
            nmrExptsDataframes[k1].append(v)
            if k_number.isdigit():
                nmrExptsToDelete.append(k)

    # # print the number of dataframes in each group
    # print("\nNumber of dataframes in each NMR experiment group:")
    # for k, v in nmrExptsDataframes.items():
    #     print(f"{k}: {len(v)} dataframes")
    
    # combine the dataframes in each group
    for k, v in nmrExptsDataframes.items(): 
        if len(v) > 1:
            # combine the dataframes in v
            combined_df = pd.concat(v, ignore_index=True)
            # remove duplicate rows based on ppm or f1_ppm and f2_ppm columns
            if "ppm" in combined_df.columns:
                combined_df = combined_df.drop_duplicates(subset=["ppm"])
            elif "f1_ppm" in combined_df.columns and "f2_ppm" in combined_df.columns:
                combined_df = combined_df.drop_duplicates(subset=["f1_ppm", "f2_ppm"])
            # sort the dataframe by ppm in descending order
            if "ppm" in combined_df.columns:
                combined_df = combined_df.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
            elif "f1_ppm" in combined_df.columns:
                combined_df = combined_df.sort_values(by=["f1_ppm"], ascending=False).reset_index(drop=True)
            else:
                print(f"No ppm column found in {k}, skipping sorting")
            # set the index to start at 1
            combined_df.index += 1
            dataframes[k] = combined_df
        elif len(v) == 1:
            dataframes[k] = v[0]
        else:
            print(f"No data for {k}")
    # print out the keys of the combined dataframes and the number of rows in each dataframe
    print("\nCombined NMR experiment dataframes:")
    for k, v in dataframes.items():
        print(f"{k}: {v.shape[0]} rows")

    # remove dataframes that were found in nmrExptsToDelete
    for k in nmrExptsToDelete:
        if k in dataframes:
            del dataframes[k]
            print(f"Removed {k} from dataframes")

    print("\nFinal dataframes keys:")
    for k in dataframes.keys():
        print(f"{k}: {dataframes[k].shape[0]} rows")

    return dataframes


class NMRProblem:
    def __init__(self, dataframes: dict):
        self.dataframes = dataframes
        print("\nNMRProblem.__init__(self, dataframes: dict):\n")
        print("*********************************************")
        print("self.dataframes.keys():", self.dataframes.keys() )
        # data is a dictionary of pandas dataframes
        self.add_missing_spectra()
        print("\nadd_missing_spectra()")
        print("915 self.dataframes.keys():\n", self.dataframes.keys() )
        self.dataframes = standardize_column_headings(self.dataframes)
        print("\nstandardize_column_headings(self.dataframes)")
        print("917 self.dataframes.keys():\n", self.dataframes.keys() )

        # print(self.dataframes["C13_1D_0"].head())

        self.dataframes = combine_multiple_nmrExpt_dataframes(self.dataframes)
        print("\ncombine_multiple_nmrExpt_dataframes(self.dataframes)")
        print("\n923 self.dataframes.keys():", self.dataframes.keys() )

        # print(self.dataframes["C13_1D"].head())

        self.dataframes = add_missing_columns_to_nmrExpt_dataframes(self.dataframes)
        print("\nadd_missing_columns_to_nmrExpt_dataframes(self.dataframes)")
        print("\n927 self.dataframes.keys():", self.dataframes.keys() )

        print(self.dataframes["C13_1D"].head())

        hsqc_H1_pks = self.dataframes["HSQC"].f2_ppm.values
        c13_pks = self.dataframes["C13_1D"].ppm.values
        # self.carbonSeparation = calc_minimum_ppm_separation(c13_pks, CARBONSEPARATION)
        # self.protonSeparation = calc_minimum_ppm_separation(hsqc_H1_pks, PROTONSEPARATION)

        self.carbonSeparation = CARBONSEPARATION
        self.protonSeparation = PROTONSEPARATION

        print("*********************************************")

        print("self.carbonSeparation", self.carbonSeparation)
        print("self.protonSeparation", self.protonSeparation)
        print("*********************************************")

        self.exact_ppm_values = self.exact_ppm_values_only()
        print("self.exact_ppm_values:", self.exact_ppm_values)

    @classmethod
    def from_mnova_json_file(cls, fn: Path):
        """
        Create an NMRProblem object from a json file.
        """
        with open(fn, "r") as f:
            cls.json_data = json.load(f)

        # convert the json data into a dictionary of pandas dataframes

        data = read_in_mesrenova_json(cls.json_data)
        dataframes = create_dataframes_from_mresnova_json(data)
        print("dataframes keys: ", dataframes.keys())
        return cls(dataframes)
    
    # @classmethod
    # def from_mnova_json_file_multiple(cls, fn: Path):
    #     """
    #     Create an NMRProblem object from a json file.
    #     """
    #     with open(fn, "r") as f:
    #         cls.json_data = json.load(f)

    #     # convert the json data into a dictionary of pandas dataframes

    #     data = read_in_mesrenova_json_multiple(cls.json_data)
    #     dataframes = create_dataframes_from_mresnova_json(data)
    #     return cls(dataframes)

    @classmethod
    def from_mnova_dict(cls, json_data: dict):
        """
        Create an NMRProblem object from a json file.
        """
        print("\nfrom_mnova_dict(cls,json_data: dict):\n")
        cls.json_data = json_data
        # convert the json data into a dictionary of pandas dataframes
        data = read_in_mesrenova_json(cls.json_data)
        dataframes = create_dataframes_from_mresnova_json(data)
        print("\nfrom_mnova_dict\n")
        print("dataframes keys: ", dataframes.keys())
        return cls(dataframes)
    
        # @classmethod
    # def from_mnova_dict_multiple(cls, json_data: dict):
    #     """
    #     Create an NMRProblem object from a json file.
    #     """
    #     print("\nfrom_mnova_dict(cls,json_data: dict):\n")
    #     cls.json_data = json_data
    #     # convert the json data into a dictionary of pandas dataframes
    #     data = read_in_mesrenova_json_multiple(cls.json_data)
    #     dataframes = create_dataframes_from_mresnova(data)
    #     return cls(dataframes)

    @classmethod
    def from_excel_file(cls, fn: Path):
        """
        Create an NMRProblem object from an excel file.
        excel files contains a number of sheets
        """
        dataframes = pd.read_excel(fn, sheet_name=None)
        return cls(dataframes)
    

    def exact_ppm_values_only(self):

        expts_with_peaks = self.dataframes["chosenSpectra"][~self.dataframes["chosenSpectra"]["expt"].str.contains("SKIP")].expt.values

        if 'C13_1D' in expts_with_peaks and 'HSQC' in expts_with_peaks:
            f1_ppm = set(self.dataframes["C13_1D"]["ppm"].values)
            f2_ppm = set(self.dataframes["HSQC"]["f2_ppm"].values)

            f2_ppm_expts = set()
            f1_ppm_expts = set()

            for expt in expts_with_peaks:
                if expt in ["HSQC", "HSQC_CLIPCOSY", "HMBC", "DDEPT_CH3_ONLY"]:
                    # add f2_ppm values to f2_ppm_expts set
                    print(f"1036 Processing expt: {expt}")
                    print(self.dataframes[expt].columns)
                    f2_ppm_set = set(self.dataframes[expt]["f2_ppm"].values)
                    f2_ppm_expts.update(f2_ppm_set)
                    
                    f1_ppm_set = set(self.dataframes[expt]["f1_ppm"].values)
                    f1_ppm_expts.update(f1_ppm_set)

                    print("\n****************************\n\n")
                    print(f"{expt} f1_ppm: {len(f1_ppm_set - f1_ppm)}, {f1_ppm_set - f1_ppm}, f2_ppm: {len(f2_ppm_set - f2_ppm)}, {f2_ppm_set - f2_ppm}")
                    print("\n****************************\n\n")


                elif expt in ["COSY"]:
                    f2_ppm_set = set(self.dataframes[expt]["f2_ppm"].values)
                    f2_ppm_expts.update(f2_ppm_set)
                    f1_ppm_set = set(self.dataframes[expt]["f1_ppm"].values)
                    f2_ppm_expts.update(f1_ppm_set)

                    
                    print("\n****************************\n\n")
                    print(f"{expt} f1_ppm: {len(f1_ppm_set - f1_ppm)}, {f1_ppm_set - f1_ppm}, f2_ppm: {len(f2_ppm_set - f2_ppm)}, {f2_ppm_set - f2_ppm}")
                    print("\n****************************\n\n")

                    


            print("f1_ppm_expts:", f1_ppm_expts)
            print("f2_ppm_expts:", f2_ppm_expts)

            print("extra f1_ppm:", f1_ppm - f1_ppm_expts)
            print("extra f2_ppm:", f2_ppm - f2_ppm_expts)

            #  check if f1_ppm - f1_ppm_expts is empty
            if len(f1_ppm - f1_ppm_expts) > 0 or len(f2_ppm - f2_ppm_expts) > 0:
                print("There are extra ppm values in f1_ppm and f2_ppm that are not in the expts_with_peaks.")
                return False
            else:
                print("All ppm values in f1_ppm and f2_ppm are accounted for in the expts_with_peaks.")
                return True
        else:
            return False
    

    def is_prediction(self):
        """
        Check if the data is a prediction
        """

        print("self.dataframes[\"MNOVAcalcMethod\"].loc[0, \"MNOVAcalcMethod\"]\n", self.dataframes["MNOVAcalcMethod"].loc[0, "MNOVAcalcMethod"])

        if self.dataframes["MNOVAcalcMethod"].loc[0, "MNOVAcalcMethod"] in [
            "MNOVA Predict",
            "NMRSHIFTDB2 Predict",
            "JEOL Predict"
        ]:
            return True
        else:
            return False


    def prediction_from_nmrshiftdb2(self):
        """
        Check if the data is a prediction from NMRSHIFTDB2
        """

        if (
            self.dataframes["MNOVAcalcMethod"].loc[0, "MNOVAcalcMethod"]
            == "NMRSHIFTDB2 Predict"
        ):
            print("Prediction from NMRSHIFTDB2")
            return True
        elif (
            self.dataframes["MNOVAcalcMethod"].loc[0, "MNOVAcalcMethod"]
            == "JEOL Predict"
        ):
            print("Prediction from JEOL")
            return True
        else:
            print("Prediction from MNOVA")
            return False  
              
    def JEOL_prediction_used(self):
        """
        Check if the data is a prediction from JEOL
        """

        if (
            self.dataframes["MNOVAcalcMethod"].loc[0, "MNOVAcalcMethod"]
            == "JEOL Predict"
        ):
            print("Prediction from JEOL")
            return True
        elif (
            self.dataframes["MNOVAcalcMethod"].loc[0, "MNOVAcalcMethod"]
            == "MNOVA Predict"
        ):
            print("Prediction from MNOVA")
            return True
        else:
            print("Prediction from NMRSHIFTDB2")
            return False

    def add_missing_spectra(self):
        """
        Add missing spectra to the dataframes dictionary
        """
        exceptions_names = [
            "molfile",
            "carbonAtomsInfo",
            "nmrAssignments",
            "c13predictions",
            "molecule",
        ]

        for k in excel_orig_df_columns.keys():
            if k in exceptions_names:
                if k not in self.dataframes:
                    self.dataframes[k] = pd.DataFrame(columns=excel_orig_df_columns[k])
            else:
                if k + "_0" not in self.dataframes:
                    self.dataframes[k+"_0"] = pd.DataFrame(columns=excel_orig_df_columns[k])

    def get_molfile_string(self):
        molfile_string = self.dataframes["molfile"].loc[0, "molfile"]
        return molfile_string

    def prepare_network_graph(self):
        # copy the dataframes for convienence
        nmrAssignments = self.dataframes["nmrAssignments"]
        self.carbonAtomsInfo = self.dataframes["carbonAtomsInfo"]
        cosy = self.dataframes["COSY"]
        hmbc = self.dataframes["HMBC"]
        hsqc = self.dataframes["HSQC"]
        h1_1d = self.dataframes["H1_1D"]
        clipcosy = self.dataframes["HSQC_CLIPCOSY"]
        ddept = self.dataframes["DDEPT_CH3_ONLY"]

        nmrAssignments = attach_symmetry_atom_index(nmrAssignments)

        # create svg image of background molecule
        molWidth = svgDimensions.MOLWIDTH
        molHeight = svgDimensions.MOLHEIGHT

        svgWidth = svgDimensions.SVGWIDTH
        svgHeight = svgDimensions.SVGHEIGHT

        # fix coordinates of molecule before creating png

        self.mol = Chem.MolFromMolBlock(self.dataframes["molfile"].loc[0, "molfile"])

        # mol = Chem.AddHs(mol)
        # AllChem.EmbedMolecule(mol, randomSeed=3)
        # mol = Chem.RemoveHs(mol)

        # mol.Compute2DCoords()

        svg_str, new_xy3 = create_svg_string(
            self.mol,
            molWidth=molWidth,
            molHeight=molHeight,
            svgWidth=svgWidth,
            svgHeight=svgHeight,
        )
        new_xy3_plus1 = {key + 1: value for key, value in new_xy3.items()}

        # add x and y coordinates to nmrAssignments
        for idx, row in nmrAssignments.iterrows():
            nmrAssignments.at[idx, "x"] = new_xy3_plus1[row["atom_idx"]][0]
            nmrAssignments.at[idx, "y"] = new_xy3_plus1[row["atom_idx"]][1]

        print("\n994 nmrAssignments\n", nmrAssignments, "\n")

        print("Adding x and y coordinates to carbonAtomsInfo")
        print("\n997 carbonAtomsInfo\n", self.carbonAtomsInfo, "\n")

        print("\n999 new_xy3\n")
        for key, value in new_xy3.items():
            print(key, value)

        # add x and y coordinates to catomAtomsInfo
        for idx, row in self.carbonAtomsInfo.iterrows():
            print("idx", idx, row["atom_idx"])
            self.carbonAtomsInfo.at[idx, "x"] = new_xy3[row["atom_idx"]][0]
            self.carbonAtomsInfo.at[idx, "y"] = new_xy3[row["atom_idx"]][1]

        # change atom_idx to id in carbonAtomsInfo
        self.carbonAtomsInfo = self.carbonAtomsInfo.rename(columns={"atom_idx": "id"})

        # add additional columns to self.carbonAtomsInfo
        self.carbonAtomsInfo["ppm"] = ""
        self.carbonAtomsInfo["ppm_calculated"] = ""
        self.carbonAtomsInfo["iupacLabel"] = ""
        self.carbonAtomsInfo["jCouplingVals"] = ""
        self.carbonAtomsInfo["jCouplingClass"] = ""
        self.carbonAtomsInfo["H1_ppm"] = ""
        self.carbonAtomsInfo["visible"] = ""

        # create h1 dataframe from nmrAssignments dataframe
        h1 = h1_from_nmrAssignments(nmrAssignments)

        # drop duplicates from h1 based on ppm
        h1_ppm_unique_list = h1.drop_duplicates(subset=["ppm"])["ppm"].to_list()
        h1_f1_ppm_unique_list = h1.drop_duplicates(subset=["f1_ppm"])[
            "f1_ppm"
        ].to_list()  # used for hsqc, ddept and hsqc_clipcosy  c13
        c13_ppm_unique_list = nmrAssignments.drop_duplicates(subset=["f1_ppm"])[
            "f1_ppm"
        ].to_list()

        if not self.exact_ppm_values:
            # tidy up the ppm values in the dataframes
            hsqc = tidyup_ppm_values(
                hsqc, h1_ppm_unique_list, "f2_ppm", ppm_tolerance=self.protonSeparation
            )
            hsqc = tidyup_ppm_values(
                hsqc, h1_f1_ppm_unique_list, "f1_ppm", ppm_tolerance=self.carbonSeparation
            )

            h1_1d = tidyup_ppm_values(
                h1_1d, h1_ppm_unique_list, "ppm", ppm_tolerance=self.protonSeparation
            )

            cosy = tidyup_ppm_values(
                cosy, h1_ppm_unique_list, "f2_ppm", ppm_tolerance=self.protonSeparation
            )
            cosy = tidyup_ppm_values(
                cosy, h1_ppm_unique_list, "f1_ppm", ppm_tolerance=self.protonSeparation
            )

            hmbc = tidyup_ppm_values(
                hmbc, h1_ppm_unique_list, "f2_ppm", ppm_tolerance=self.protonSeparation
            )
            hmbc = tidyup_ppm_values(
                hmbc, c13_ppm_unique_list, "f1_ppm", ppm_tolerance=self.carbonSeparation
            )

            hmbc.drop(hmbc[hmbc.f1_ppm_prob == 0].index, inplace=True)
            hmbc.drop(hmbc[hmbc.f2_ppm_prob == 0].index, inplace=True)

            clipcosy = tidyup_ppm_values(
                clipcosy, h1_ppm_unique_list, "f2_ppm", ppm_tolerance=self.protonSeparation
            )
            clipcosy = tidyup_ppm_values(
                clipcosy,
                h1_f1_ppm_unique_list,
                "f1_ppm",
                ppm_tolerance=self.carbonSeparation,
            )

            ddept = tidyup_ppm_values(
                ddept, h1_ppm_unique_list, "f2_ppm", ppm_tolerance=self.protonSeparation
            )
            ddept = tidyup_ppm_values(
                ddept, h1_f1_ppm_unique_list, "f1_ppm", ppm_tolerance=self.carbonSeparation
            )
        else:
            print("Exact ppm values only, no tidy up required")

        # add jCouplingVals and jCouplingClass to h1 from h1_1D
        h1 = add_jCouplingVals_jCouplingClass_to_h1(h1, h1_1d)

        # add jCouplingVals and jCouplingClass to nmrAssignments from h1
        nmrAssignments = add_jCouplingVals_jCouplingClass_to_nmrAssignments(
            nmrAssignments, h1
        )

        c13_sorted = nmrAssignments.sort_values(
            by=["atom_idx"], ascending=True
        ).reset_index(drop=True)

        # add dummy ppm_calculated to c13_sorted
        c13_sorted["ppm_calculated"] = ""

        G2 = create_network_graph(c13_sorted, h1)

        cosy_combinations = create_cosy_combinations(cosy, h1)

        clipcosy_combinations = create_clipcosy_combinations(clipcosy, h1)

        kept_cosy_combinations = remove_duplicate_cosy_combinations(cosy_combinations)
        kept_clipcosy_combinations = remove_duplicate_cosy_combinations(
            clipcosy_combinations
        )

        G2 = add_cosy_edges_to_graph(G2, kept_cosy_combinations)
        G2 = add_cosy_edges_to_graph(G2, kept_clipcosy_combinations)

        G2 = add_hmbc_edges_to_graph(G2, hmbc, h1, c13_sorted)

        self.G2 = G2

        # jsonGraphData = json_graph.node_link_data(G2, edges="links")  # added edges='links' for compatibility
        jsonGraphData = json_graph.node_link_data(G2)  # reverted back to original as pythonanywhere lower version

        self.jinjadata = {
            "svg_container": svg_str,
            "graph_edges": jsonGraphData["links"],
            "graph_nodes": jsonGraphData["nodes"],
            "translateX": int((svgWidth - molWidth) / 2),
            "translateY": int((svgHeight - molHeight) / 2),
            "title": "dummy_title",
            "catoms": json.dumps(
                self.carbonAtomsInfo.to_dict(orient="records"), indent=4
            ),
        }


if __name__ == "__main__":

    if platform.system() == "Windows":
        fn_json = Path(
            r"C:\Users\vsmw51\Dropbox\PC\Downloads\exampleProblems\dimethylbenzene\dimethylbenzene_assignments_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/predictC13test/predictC13test_predictions_mresnova.json"
        )
        fn_json = Path(
            r"C:\Users\vsmw51\Dropbox\PC\Downloads\exampleProblems\EVB_330b_annotated\EVB_330b_annotated_assignments_mresnova.json"
        )
        fn_json = Path(
            r"C:\Users\vsmw51\OneDrive - Durham University\projects\programming\2022\python\python\rdkit37\projects\simpleNMR\exampleProblems\santonin\santonin_assignments_mresnova.json"
        )

    elif platform.system() == "Darwin":
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/dimethylbenzene/dimethylbenzene_assignments_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/predictC13test/predictC13test_predictions_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/EVB_330b_annotated/EVB_330b_annotated_assignments_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2022/python/python/rdkit37/projects/simpleNMR/exampleProblems/santonin/santonin_assignments_mresnova.json"
        )

    if fn_json.exists():
        problemdata_json = NMRProblem.from_mnova_json_file(fn_json)
    else:
        print(f"{fn_json.name} not found")

    problemdata_json.prepare_network_graph()

    create_htmlpage_from_graph(
        problemdata_json.G2,
        "test.html",
        r"templates/d3molplotmnova_template.html",
        problemdata_json.jinjadata,
    )
