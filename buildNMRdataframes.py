import sys
import numpy as np
import pandas as pd
from scipy import stats

from networkx.readwrite import json_graph

from excelheaders import EXCEL_DF_COLUMNS as excel_df_columns
from excelheaders import EXCEL_ORIG_DF_COLUMNS as excel_orig_df_columns

from html_from_assignments import create_cosy_combinations
from html_from_assignments import create_clipcosy_combinations
from html_from_assignments import remove_duplicate_cosy_combinations

from html_from_assignments import add_cosy_edges_to_graph
from html_from_assignments import add_hmbc_edges_to_graph

from html_from_assignments import create_svg_string
from html_from_assignments import NMRProblem
from html_from_assignments import create_network_graph

import expectedmolecule

from globals import SVG_DIMENSIONS as svgDimensions

from globals import CARBONSEPARATION
from globals import PROTONSEPARATION


def find_and_group_CH2s(df1):

    # get a list of the CH2 resonances
    ch2_vals = df1[df1.CH2].f1_ppm.tolist()
    ch2_idx_vals = df1[df1.CH2].index.tolist()

    unique_ch2s = []
    unique_idxs = []

    # start with first hmbc value in the hmbc list
    # create a probability distribution around it and obtain the probability of all the other values to
    # to see if they are close to the first value.
    # all values that have a +ve probability are close to the first value
    # all values with a zero probability are not.
    # add the +ve to a saved list of lists "similar_hmbcs"
    # then remove them from the original list and repeat until original list length is zero

    while len(ch2_vals):
        # choose first from the list
        p0 = ch2_vals[0]
        # find list of hmbc values that are similar to the first in the list
        similar_ch2s = [
            p for p in ch2_vals if stats.norm.pdf(p, loc=p0, scale=CARBONSEPARATION) > 0
        ]
        similar_idxs = [
            i
            for i, p in zip(ch2_idx_vals, ch2_vals)
            if stats.norm.pdf(p, loc=p0, scale=CARBONSEPARATION) > 0
        ]
        # save the list
        unique_ch2s.append(similar_ch2s)
        unique_idxs.append(similar_idxs)

        # keep only hmbc values that were not similar
        ch2_idx_vals = [
            i
            for i, p in zip(ch2_idx_vals, ch2_vals)
            if stats.norm.pdf(p, loc=p0, scale=CARBONSEPARATION) == 0
        ]
        ch2_vals = [
            p
            for p in ch2_vals
            if stats.norm.pdf(p, loc=p0, scale=CARBONSEPARATION) == 0
        ]

    return unique_idxs, unique_ch2s


def find_nearest(true_values: list, value: float):
    arraynp = np.asarray(true_values)
    idx = (np.abs(arraynp - value)).argmin()
    return arraynp[idx]


def tidyup_ppm_values(
    df: pd.DataFrame, true_values: list, column_name: str, ppm_tolerance=0.005
) -> pd.DataFrame:
    """_summary_

    Args:
        df (pd.DataFrame): nmr experiment dataframe
        true_values (list): list of true values to replace column_name values with
        column_name (str): column name to replace values
        ppm_tolerance (float, optional): grouping tolerance. Defaults to 0.005.

    Returns:
        pd.DataFrame: orignal dataframe with column_name values replaced with nearest true_values
    """

    # make a copy of the column_name adding a suffix orig
    df[f"{column_name}_orig"] = df[column_name]

    # make a probability column to see how far replacement is from original
    df[f"{column_name}_prob"] = 0.0

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


def init_c13_from_hsqc_and_hmbc(
    hsqc_df: pd.DataFrame, hmbc_df: pd.DataFrame
) -> pd.DataFrame:
    """works out the unique C13 values from the HSQC and HMBC dataframes

    Args:
        hsqc_df (pd.DataFrame): hsqc dataframe
        hmbc_df (pd.DataFrame): hmbc dataframe

    Returns:
        pd.DataFrame: c13 dataframe
    """

    # if self.C13_data_present:
    #     return self.c13_df, False

    # find CH2 groups in hsqc_df
    # CH2 groups where intensity is < 0.0

    hsqc_df["CH2"] = False
    hsqc_df.loc[hsqc_df.intensity < 0, "CH2"] = True

    unique_idxs, unique_ch2s = find_and_group_CH2s(hsqc_df)

    # replace hsqc CH2 values in f1_ppm of HSQC
    for idx, ch2 in zip(unique_idxs, unique_ch2s):
        hsqc_df.loc[idx, "f1_ppm"] = np.mean(ch2)

    # replace values for f2_ppm in hmbc starting from f2_ppm hsqc
    hmbc_df = tidyup_ppm_values(
        hmbc_df,
        sorted(hsqc_df.f2_ppm.unique().tolist(), reverse=True),
        "f2_ppm",
        ppm_tolerance=PROTONSEPARATION,
    )

    # if any hmbc f2_ppm_probs == 0 then drop the row
    hmbc_df.drop(hmbc_df[hmbc_df.f2_ppm_prob == 0].index, inplace=True)

    # find all f1_ppm HMBC idx resonances that are not showing up in the HSQC f2_ppm
    iii = []
    for i in hmbc_df.index:
        prob_vals = []
        for c in hsqc_df.f1_ppm.unique():
            prob_vals.append(
                stats.norm.pdf(hmbc_df.loc[i, "f1_ppm"], loc=c, scale=CARBONSEPARATION)
            )
        if np.array(prob_vals).sum() == 0:
            iii.append(i)
        else:
            pass

    # keep only the unique hmbc resonances not in f1_ppm HSQC
    # get a list of the hmbc resonances
    hmbcs = hmbc_df.loc[iii, "f1_ppm"].tolist()
    unique_hmbc = []

    # start with first hmbc value in the hmbc list
    # create a probability distribution around it and obtain the probability of all the other values to
    # to see if they are close to the first value.
    # all values that have a +ve probability are close to the first value
    # all values with a zero probability are not.
    # add the +ve to a saved list of lists "similar_hmbcs"
    # then remove them from the original list and repeat until original list length is zero

    while len(hmbcs):
        # choose first from the list
        p0 = hmbcs[0]
        # find list of hmbc values that are similar to the first in the list
        similar_hmbcs = [
            p for p in hmbcs if stats.norm.pdf(p, loc=p0, scale=CARBONSEPARATION) > 0
        ]
        # save the list
        unique_hmbc.append(similar_hmbcs)

        # keep only hmbc values that were not similar
        hmbcs = [
            p for p in hmbcs if stats.norm.pdf(p, loc=p0, scale=CARBONSEPARATION) == 0
        ]

    # create an array of mean values for hmbc f1_ppm not found in hsqc f1_ppm
    mean_unique_hmbc_vals = [np.mean(h) for h in unique_hmbc]

    # tidyup f1_ppm values in hmbc that are not in f1_ppm HSQC
    hmbc_1 = tidyup_ppm_values(
        hmbc_df.loc[iii],
        mean_unique_hmbc_vals,
        "f1_ppm",
        ppm_tolerance=CARBONSEPARATION,
    )

    # tidyup f1_ppm values in hmbc that are in f1_ppm HSQC
    hmbc_2 = tidyup_ppm_values(
        hmbc_df.drop(iii),
        hsqc_df.f1_ppm.unique(),
        "f1_ppm",
        ppm_tolerance=CARBONSEPARATION,
    )

    # rejoin two parts of HMBC data
    hmbc_df = pd.concat([hmbc_1, hmbc_2])
    hmbc_df.sort_index(inplace=True)

    # add f2p_ppm column to HSQC and HMBC tables
    # f2p_ppm is C13 one to one relation between f1_ppm and f2_ppm in HSQC
    hmbc_df["f2p_ppm"] = 0.0
    for idx, f2ppmHSQC, f1ppmHSQC in zip(hsqc_df.index, hsqc_df.f2_ppm, hsqc_df.f1_ppm):
        hmbc_df.loc[hmbc_df.f2_ppm == f2ppmHSQC, "f2p_ppm"] = f1ppmHSQC

    # return list of C13 values
    c13_list = sorted(
        set(hmbc_df.f1_ppm).union(set(hmbc_df.f2p_ppm), set(hsqc_df.f1_ppm)),
        reverse=True,
    )

    # create a dataframe of C13 values
    c13_df = pd.DataFrame(c13_list, columns=["ppm"])
    c13_df["signaltype"] = "Compound"
    c13_df.index = c13_df.index + 1

    # label c13 values as CH2 or not using the hsqc_df
    c13_df["CH2"] = False
    for idx, ppm in zip(c13_df.index, c13_df.ppm):
        if ppm in hsqc_df.f1_ppm.values:
            c13_df.loc[idx, "CH2"] = hsqc_df[hsqc_df.f1_ppm == ppm]["CH2"].values[0]

    #  label c13 values as quaternary if not present in hsqc_df
    c13_df["quaternary"] = False
    c13_df["CH0"] = False
    for idx, ppm in zip(c13_df.index, c13_df.ppm):
        if ppm not in hsqc_df.f1_ppm.values:
            c13_df.loc[idx, "quaternary"] = True
            c13_df.loc[idx, "CH0"] = True

    # label c13 values as CH3CH if quaternary and CH2 both False
    c13_df["CH3CH1"] = False
    for idx, ppm in zip(c13_df.index, c13_df.ppm):
        if not c13_df.loc[idx, "quaternary"] and not c13_df.loc[idx, "CH2"]:
            c13_df.loc[idx, "CH3CH1"] = True

    return c13_df, True


def init_h1_from_hsqc(hsqc_df):

    # if self.pureshift_data_present and self.H1_data_present:
    #     return self.h1_df, self.pureshift_df

    # if (not self.h1_df.empty) and (not self.pureshift_df.empty):
    #     return self.h1_df, self.pureshift_df

    # make a list of the columns to copy over from hsqc_df

    h1_cols = [
        "f2_ppm",
        "integral",
        "numProtons",
        "jCouplingClass",
        "jCouplingVals",
        "atom_idx",
        "sym_atom_idx",
        "atomNumber",
        "sym_atomNumber",
        "CH3",
        "CH2",
        "CH1",
        "CH0",
        "CH3CH1",
        "x",
        "y",
        "Compound",
        "intensity",
    ]
    h1_cols_keep = []
    for colname in hsqc_df.columns:
        if colname in h1_cols:
            h1_cols_keep.append(colname)

    h1_df = hsqc_df[h1_cols_keep].copy()
    h1_df["ppm"] = h1_df["f2_ppm"].copy()
    # h1_df["signaltype"] = self.hsqc_df["signaltype"].copy()

    # keep only rows where type is 'Compound' or  "" (empty)
    # self.h1_df = self.h1_df[self.h1_df["signaltype"] == "Compound"]

    if "numProtons" not in h1_df.columns:
        h1_df["numProtons"] = -1
    if "integral" not in h1_df.columns:
        h1_df["integral"] = -1.0
        if "numprotons" in h1_df.columns:
            h1_df["integral"] = h1_df["numprotons"].copy()

    # self.h1_df["jCouplingVals"] = -1
    # self.h1_df["jCouplingClass"] = "u"
    # self.h1_df["intensity"] = self.hsqc_df["intensity"].copy()
    # self.h1_df["range"] = 0.0

    # order h1_df by ppm in place highest to lowest, reset index to 1,2,3,4,5...
    h1_df = h1_df.sort_values("ppm", ascending=False, ignore_index=True)
    h1_df.index = h1_df.index + 1

    # create pureshift_df from h1_df
    # self.pureshift_df = self.h1_df.copy()

    return h1_df, True


def assign_rows(df1, df2, df1_ppm, df2_ppm, df2_CHn, *columns):
    """
    Assign values from df2 to df1 based on the ppm values in df1_ppm
    """

    df2_rows = (
        df2[df2[df2_CHn]]
        .sort_values(by=[df2_ppm], ascending=False)
        .reset_index(drop=True)
    )
    df1_rows = (
        df1[df1[df2_CHn]]
        .sort_values(by=[df1_ppm], ascending=False)
        .reset_index(drop=True)
    )
    df1_ppm_unique = df1_rows[df1_ppm].unique()

    for ppm, idx in zip(df1_ppm_unique, df2_rows.index):
        for col in columns:
            df1.loc[(df1[df1_ppm] == ppm), col] = df2_rows.loc[idx, col]


def complete_dataframes_c13_hsqc(
    mol,
    hsqc,
    c13_1d,
    ddept=pd.DataFrame(columns=excel_df_columns["ddept_ch3_only"]),
    hmbc=pd.DataFrame(columns=excel_df_columns["hmbc"]),
    cosy=pd.DataFrame(columns=excel_df_columns["cosy"]),
    hsqc_clipcosy=pd.DataFrame(columns=excel_df_columns["hsqc_clip_cosy"]),
):

    # create a dataframe of C13 values from c13_1d
    c13 = c13_1d.copy()

    # check if signaltype in dataframe and if it is keep only rows where signaltype is Compound or 0
    if "signaltype" in c13.columns:
        c13 = c13[(c13["signaltype"] == "Compound") | (c13["signaltype"] == 0)]
        if len(c13) == 0:
            c13 = c13_1d.copy()
    else:
        c13["signaltype"] = 0

    c13["numProtons"] = -1

    h1, ok = init_h1_from_hsqc(hsqc)

    # tidy up the f1_ppm values in the hsqc dataframe using the c13 dataframe ppm values
    hsqc = tidyup_ppm_values(hsqc, c13.ppm.values, "f1_ppm")

    # prepare hsqc dataframe
    hsqc["CH3CH1"] = False
    hsqc["CH3"] = False
    hsqc["CH2"] = False
    hsqc["CH1"] = False
    hsqc["CH0"] = False
    hsqc["numProtons"] = -1
    hsqc.loc[hsqc.intensity < 0, "CH2"] = True

    unique_idxs, unique_ch2s = find_and_group_CH2s(hsqc)

    # replace hsqc CH2 values in f1_ppm of HSQC
    for idx, ch2 in zip(unique_idxs, unique_ch2s):
        hsqc.loc[idx, "f1_ppm"] = np.mean(ch2)

    # set numProtons to 2 for CH2 values in hsqc dataframe
    hsqc.loc[hsqc.CH2, "numProtons"] = 2

    hsqc.loc[hsqc.intensity > 0, "CH3CH1"] = True

    # label c13 values as CH2 or not using the hsqc_df
    c13["CH2"] = False
    for idx, ppm in zip(c13.index, c13.ppm):
        if ppm in hsqc.f1_ppm.values:
            c13.loc[idx, "CH2"] = hsqc[hsqc.f1_ppm == ppm]["CH2"].values[0]

    #  label c13 values as quaternary if not present in hsqc_df
    c13["quaternary"] = False
    c13["CH0"] = False
    for idx, ppm in zip(c13.index, c13.ppm):
        if ppm not in hsqc.f1_ppm.values:
            c13.loc[idx, "quaternary"] = True
            c13.loc[idx, "CH0"] = True
            c13.loc[idx, "numProtons"] = 0

    # label c13 values as CH3CH if quaternary and CH2 both False
    c13["CH3CH1"] = False
    for idx, ppm in zip(c13.index, c13.ppm):
        if not c13.loc[idx, "CH0"] and not c13.loc[idx, "CH2"]:
            c13.loc[idx, "CH3CH1"] = True

    # # set the CH3Ch1 values in the hsqc dataframe intensity greater than 0
    #

    h1, c13, hsqc = produce_standardized_h1_c13_hsqc_dataframes(
        mol, h1, c13, hsqc, ddept, hmbc, cosy, hsqc_clipcosy
    )

    return h1, c13, hsqc


def complete_dataframes_hsqc_only(
    mol,
    hsqc,
    ddept=pd.DataFrame(columns=excel_df_columns["ddept_ch3_only"]),
    hmbc=pd.DataFrame(columns=excel_df_columns["hmbc"]),
    cosy=pd.DataFrame(columns=excel_df_columns["cosy"]),
    hsqc_clipcosy=pd.DataFrame(columns=excel_df_columns["hsqc_clip_cosy"]),
):

    c13, ok = init_c13_from_hsqc_and_hmbc(hsqc, hmbc)
    h1, ok = init_h1_from_hsqc(hsqc)

    h1, c13, hsqc = produce_standardized_h1_c13_hsqc_dataframes(
        mol, h1, c13, hsqc, ddept, hmbc, cosy, hsqc_clipcosy
    )

    return h1, c13, hsqc


def produce_standardized_h1_c13_hsqc_dataframes(
    mol, h1, c13, hsqc, ddept, hmbc, cosy, hsqc_clipcosy
):

    for df in [h1, c13, hsqc, hmbc, cosy, ddept, hsqc_clipcosy]:
        if "atom_idx" not in df.columns:
            df["atom_idx"] = -1
        if "sym_atom_idx" not in df.columns:
            df["sym_atom_idx"] = ""
        if "atomNumber" not in df.columns:
            df["atomNumber"] = -1
        if "sym_atomNumber" not in df.columns:
            df["sym_atomNumber"] = ""
        if "CH3" not in df.columns:
            df["CH3"] = False
        if "CH2" not in df.columns:
            df["CH2"] = False
        if "CH1" not in df.columns:
            df["CH1"] = False
        if "CH0" not in df.columns:
            df["CH0"] = False
        if "CH3CH1" not in df.columns:
            df["CH3CH1"] = False
        if "x" not in df.columns:
            df["x"] = -1
        if "y" not in df.columns:
            df["y"] = -1
        if "numProtons" not in df.columns:
            df["numProtons"] = -1
        if "jCouplingVals" not in df.columns:
            df["jCouplingVals"] = ""
        if "jCouplingClass" not in df.columns:
            df["jCouplingClass"] = ""

    true_c13_values = c13["ppm"].values
    true_h1_values = h1["ppm"].values

    # tidy up the C13 ppm values in the dataframes
    hmbc = tidyup_ppm_values(hmbc, true_c13_values, "f1_ppm")
    ddept = tidyup_ppm_values(ddept, true_c13_values, "f1_ppm")
    hsqc_clipcosy = tidyup_ppm_values(hsqc_clipcosy, true_c13_values, "f1_ppm")

    # tidy up the H1 ppm values in the dataframes
    hmbc = tidyup_ppm_values(hmbc, true_h1_values, "f2_ppm")
    ddept = tidyup_ppm_values(ddept, true_h1_values, "f2_ppm")
    hsqc_clipcosy = tidyup_ppm_values(hsqc_clipcosy, true_h1_values, "f2_ppm")
    cosy = tidyup_ppm_values(cosy, true_h1_values, "f2_ppm")
    cosy = tidyup_ppm_values(cosy, true_h1_values, "f1_ppm")

    print("\nddept\n", ddept)
    # add atom_idx and atomNumber to ddept_ch3_only using the values in mol.molprops_df
    # extract the CH3 rows from mol.molprops_df and sort based on descending ppm values
    if len(ddept) > 0:
        ch3_rows = (
            mol.molprops_df[mol.molprops_df["CH3"] == True]
            .sort_values(by=["ppm"], ascending=False)
            .reset_index(drop=True)
        )

        # sort ddept_ch3_only by ppm in descending order
        ddept = ddept.sort_values(by=["f1_ppm"], ascending=False).reset_index(drop=True)

        # copy over the atom_idx and atomNumber values from ch3_rows to ddept_ch3_only
        ddept["atom_idx"] = ch3_rows["atom_idx"]
        ddept["atomNumber"] = ch3_rows["atomNumber"]
        ddept["x"] = ch3_rows["x"]
        ddept["y"] = ch3_rows["y"]
        ddept["sym_atom_idx"] = ch3_rows["sym_atom_idx"]
        ddept["sym_atomNumber"] = ch3_rows["sym_atomNumber"]

        ddept["CH3"] = True
        ddept["CH2"] = False
        ddept["CH1"] = False
        ddept["CH0"] = False
        ddept["CH3CH1"] = True

        ddept["numProtons"] = 3

    # assign CH2, CH1, CH0 numprotons values in the hsqc dataframe

    # CH2s should have already been set in the hsqc dataframe
    hsqc.loc[hsqc.CH2, "numProtons"] = 2

    # set the CH3CH1 values in the hsqc dataframe
    hsqc.loc[~hsqc.CH2, "CH3CH1"] = True

    # if no CH3 values in expected molecule then set CH values to True based on CH1CH3 values being true_c13_values
    if mol.num_CH3_carbon_atoms == 0:
        hsqc.loc[hsqc.CH3CH1, "CH1"] = True
        hsqc.loc[hsqc.CH3CH1, "numProtons"] = 1

        # assign atom_idx and atomNumber values to the hsqc dataframe
    elif mol.num_CH1_carbon_atoms == 0:
        hsqc.loc[hsqc.CH3CH1, "CH3"] = True
        hsqc.loc[hsqc.CH3CH1, "numProtons"] = 3

    elif len(ddept) > 0:
        # use the ddpt_ch3_only dataframe to update the CH3 column in the hsqc dataframe by matching the f1_ppm values
        for idx, row in ddept.iterrows():
            hsqc.loc[hsqc.f1_ppm == row["f1_ppm"], "CH3"] = True
            hsqc.loc[hsqc.f1_ppm == row["f1_ppm"], "numProtons"] = 3

        # set the rest of the CH1CH3 values that are True and CH3 is False to CH
        hsqc.loc[(hsqc.CH3CH1) & (~hsqc.CH3), "CH1"] = True
        hsqc.loc[(hsqc.CH3CH1) & (~hsqc.CH3), "numProtons"] = 1

        print("hsqc")
        print(
            hsqc[["f1_ppm", "f2_ppm", "atom_idx", "CH3", "CH1", "CH3CH1", "numProtons"]]
        )

    else:  # else attempt to assign the CH3 values via the expected molecule and matching the f1_ppm values

        hsqc_CH3CH1_df = hsqc[hsqc.CH3CH1]

        molprops_CH3CH1_df = mol.molprops_df[mol.molprops_df.CH3CH1]
        molprops_sym_CH3CH1_df = mol.sym_molprops_df[mol.sym_molprops_df.CH3CH1]

        if len(hsqc_CH3CH1_df) == len(molprops_sym_CH3CH1_df):
            molprops_df = molprops_sym_CH3CH1_df
        elif len(hsqc_CH3CH1_df) == len(molprops_CH3CH1_df):
            molprops_df = molprops_CH3CH1_df
        else:
            molprops_combo = []
            len_molprops_CH3 = len(molprops_CH3CH1_df[molprops_CH3CH1_df.CH3])
            len_molprops_sym_CH3 = len(
                molprops_sym_CH3CH1_df[molprops_sym_CH3CH1_df.CH3]
            )
            len_molprops_CH1 = len(molprops_CH3CH1_df[molprops_CH3CH1_df.CH1])
            len_molprops_sym_CH1 = len(
                molprops_sym_CH3CH1_df[molprops_sym_CH3CH1_df.CH1]
            )

            molprops_CH3 = molprops_CH3CH1_df[molprops_CH3CH1_df.CH3]
            molprops_sym_CH3 = molprops_sym_CH3CH1_df[molprops_sym_CH3CH1_df.CH3]

            molprops_CH1 = molprops_CH3CH1_df[molprops_CH3CH1_df.CH1]
            molprops_sym_CH1 = molprops_sym_CH3CH1_df[molprops_sym_CH3CH1_df.CH1]

            if len(hsqc_CH3CH1_df) == len_molprops_CH3 + len_molprops_sym_CH1:
                molprops_combo.append(molprops_CH3)
                molprops_combo.append(molprops_sym_CH1)

            elif len(hsqc_CH3CH1_df) == len_molprops_sym_CH3 + len_molprops_CH1:
                molprops_combo.append(molprops_sym_CH3)
                molprops_combo.append(molprops_CH1)

            #  combine the molprops dataframes into one dataframe
            molprops_df = pd.concat(molprops_combo)

        # order the ppm values in the two dataframes
        hsqc_CH3CH1_df = hsqc_CH3CH1_df.sort_values(by=["f1_ppm"], ascending=False)
        molprops_df = molprops_df.sort_values(by=["ppm"], ascending=False)

        # transfer the information from the molprops_df to the hsqc dataframe

        hsqc_idx = hsqc_CH3CH1_df.index
        molprops_idx = molprops_df.index

        hsqc.loc[hsqc_idx, "CH3"] = molprops_df.loc[molprops_idx, "CH3"].values
        hsqc.loc[hsqc_idx, "numProtons"] = molprops_df.loc[
            molprops_idx, "numProtons"
        ].values
        hsqc.loc[hsqc_idx, "CH1"] = molprops_df.loc[molprops_idx, "CH1"].values
        hsqc.loc[hsqc_idx, "CH3CH1"] = molprops_df.loc[molprops_idx, "CH3CH1"].values
        hsqc.loc[hsqc_idx, "atom_idx"] = molprops_df.loc[
            molprops_idx, "atom_idx"
        ].values
        hsqc.loc[hsqc_idx, "atomNumber"] = molprops_df.loc[
            molprops_idx, "atomNumber"
        ].values
        hsqc.loc[hsqc_idx, "sym_atom_idx"] = molprops_df.loc[
            molprops_idx, "sym_atom_idx"
        ].values
        hsqc.loc[hsqc_idx, "sym_atomNumber"] = molprops_df.loc[
            molprops_idx, "sym_atomNumber"
        ].values
        hsqc.loc[hsqc_idx, "x"] = molprops_df.loc[molprops_idx, "x"].values
        hsqc.loc[hsqc_idx, "y"] = molprops_df.loc[molprops_idx, "y"].values

        print("\nhsqc\n")
        print(hsqc[["f1_ppm", "atom_idx", "atomNumber", "CH3", "CH2", "CH1", "CH3CH1"]])

        print("\nmolprops_df\n")
        print(
            molprops_df[
                ["ppm", "atom_idx", "atomNumber", "CH3", "CH2", "CH1", "CH3CH1"]
            ]
        )

        molprops_CH2_df = molprops_df[molprops_df.CH2]
        print("\nmol.molprops_df CH2\n")
        print(
            molprops_CH2_df[
                ["ppm", "atom_idx", "atomNumber", "CH3", "CH2", "CH1", "CH3CH1"]
            ]
        )

        # for idx, row in molprops_df.iterrows():
        #     hsqc.loc[hsqc.f1_ppm == row["ppm"], "CH3"] = True
        #     hsqc.loc[hsqc.f1_ppm == row["ppm"], "numProtons"] = 3

        # # check the symmetry of the molecule first
        # ch3_mol_symmetry = False
        # ch1_mol_symmetry = False

        # ch3_mol_symmetry = len(mol.sym_molprops_df[mol.sym_molprops_df.CH3]) < len(
        #     mol.molprops_df[mol.molprops_df.CH3]
        # )
        # ch1_mol_symmetry = len(mol.sym_molprops_df[mol.sym_molprops_df.CH1]) < len(
        #     mol.molprops_df[mol.molprops_df.CH1]
        # )

        # print("ch3_mol_symmetry", ch3_mol_symmetry)
        # print("ch1_mol_symmetry", ch1_mol_symmetry)

        # print("len(c13)", len(c13))
        # print("len(hsqc)", len(hsqc))
        # print("len(mol.molprops_df)", len(mol.molprops_df))
        # print("len(mol.sym_molprops_df)", len(mol.sym_molprops_df))

        # len_c13 = len(c13)
        # len_c13_CH2 = len(c13[c13.CH2])
        # len_c13_CH3CH = len(c13[c13.CH3CH1])
        # len_c13_CH0 = len(c13[c13.CH0])

        # print("len_c13", len_c13)
        # print("len_c13_CH2", len_c13_CH2)
        # print("len_c13_CH3CH", len_c13_CH3CH)
        # print("len_c13_CH0", len_c13_CH0)

        # print(c13)

        # len_molprops_df = len(mol.molprops_df)
        # len_sym_molprops_df = len(mol.sym_molprops_df)

        # len_molprops_df_CH3 = len(mol.molprops_df[mol.molprops_df.CH3])
        # len_sym_molprops_df_CH3 = len(mol.sym_molprops_df[mol.sym_molprops_df.CH3])

        # len_molprops_df_CH1 = len(mol.molprops_df[mol.molprops_df.CH1])
        # len_sym_molprops_df_VH1 = len(mol.sym_molprops_df[mol.sym_molprops_df.CH1])

        # len_molprops_df_CH2 = len(mol.molprops_df[mol.molprops_df.CH2])
        # len_sym_molprops_df_CH2 = len(mol.sym_molprops_df[mol.sym_molprops_df.CH2])

        # len_molprops_df_ch0 = len(mol.molprops_df[mol.molprops_df.CH0])
        # len_sym_molprops_df_ch0 = len(mol.sym_molprops_df[mol.sym_molprops_df.CH0])

        # len_molprops_df_CH3CH = len(mol.molprops_df[mol.molprops_df.CH3CH1])
        # len_sym_molprops_df_CH3CH = len(mol.sym_molprops_df[mol.sym_molprops_df.CH3CH1])

        # if len_c13 == len_molprops_df:
        #     print("c13 and mol.molprops_df are the same length")
        #     molprops_df = mol.molprops_df

        # elif len(c13) == len_sym_molprops_df:
        #     print("c13 and mol.sym_molprops_df are the same length")
        #     molprops_df = mol.sym_molprops_df

        # else:
        #     molprops_combo = []
        #     if len_c13_CH2 == len_molprops_df_CH2:
        #         molprops_combo.append(mol.molprops_df[mol.molprops_df.CH2])
        #     elif len_c13_CH2 == len_sym_molprops_df_CH2:
        #         molprops_combo.append(mol.sym_molprops_df[mol.sym_molprops_df.CH2])

        #     if len_c13_CH3CH == len_molprops_df_CH3CH:
        #         molprops_combo.append(mol.molprops_df[mol.molprops_df.CH3CH1])
        #     elif len_c13_CH3CH == len_sym_molprops_df_CH3CH:
        #         molprops_combo.append(mol.sym_molprops_df[mol.sym_molprops_df.CH3])

        #     if len_c13_CH0 == len_molprops_df_ch0:
        #         molprops_combo.append(mol.molprops_df[mol.molprops_df.CH0])
        #     elif len_c13_CH0 == len_sym_molprops_df_ch0:
        #         molprops_combo.append(mol.sym_molprops_df[mol.sym_molprops_df.CH0])

        #     #  combine the molprops dataframes into one dataframe
        #     molprops_df = pd.concat(molprops_combo)

        # print("\nmolprops_df\n", molprops_df)

        # use the fact carbon atoms above 67 ppm mst be CH1
        # split the hsqc dataframe into two dataframes based on the ppm values and CH3CH1 values being True
        # hsqc_CH3CH1_df = hsqc[hsqc.CH3CH1]

        # hsqc_CH3CH1_gt_67_df = hsqc_CH3CH1_df[hsqc_CH3CH1_df.f1_ppm > 67]
        # hsqc_CH3CH1_lt_67_df = hsqc_CH3CH1_df[hsqc_CH3CH1_df.f1_ppm <= 67]

        # # hsqc_CH3CH1_gt_67_df.loc[idx, "numProtons"] = 1

        # # # set the CH1 values to true in the hsqc dataframe and set numProtons to 1

        # # # split the mol dataframe into two dataframes based on the ppm values and CH3CH1 values being True
        # # mol_sym_CH3CH1_df = mol.sym_molprops_df[mol.sym_molprops_df.CH3CH1]
        # # mol_sym_CH3CH1_gt_67_df = mol_sym_CH3CH1_df[mol_sym_CH3CH1_df.ppm > 67]
        # # mol_sym_CH3CH1_lt_67_df = mol_sym_CH3CH1_df[mol_sym_CH3CH1_df.ppm <= 67]

        # mol_CH3CH1_df = molprops_df[molprops_df.CH3CH1]
        # mol_CH3CH1_gt_67_df = mol_CH3CH1_df[mol_CH3CH1_df.ppm > 67]
        # mol_CH3CH1_lt_67_df = mol_CH3CH1_df[mol_CH3CH1_df.ppm <= 67]

        # len_molprops_gt_67_CH1 = len(mol_CH3CH1_gt_67_df[mol_CH3CH1_gt_67_df.CH1])
        # len_molprops_lt_67_CH1 = len(mol_CH3CH1_lt_67_df[mol_CH3CH1_lt_67_df.CH1])

        # len_hsqc_gt_67_CH1 = len(hsqc_CH3CH1_gt_67_df[hsqc_CH3CH1_gt_67_df.CH1])
        # len_hsqc_lt_67_CH1 = len(hsqc_CH3CH1_lt_67_df[hsqc_CH3CH1_lt_67_df.CH1])

        # # if length of mol_CH3CH1_gt_67_df == length of hsqc_CH3CH1_gt_67_df then we can assign the CH1 values
        # if len_molprops_gt_67_CH1 == len_hsqc_gt_67_CH1:
        #     idx = hsqc_CH3CH1_gt_67_df.index
        #     hsqc.loc[idx, "CH1"] = True
        #     hsqc.loc[idx, "numProtons"] = 1

        # else:  # we have a mis-match in ppm values so include all CH3CH1 values in CH3CH1_lt_67_df
        #     hsqc_CH3CH1_lt_67_df = hsqc_CH3CH1_df
        #     mol_CH3CH1_lt_67_df = mol_CH3CH1_df

        #     len_hsqc_

        # # check the lengths again
        # if len(mol_CH3CH1_lt_67_df) == len(hsqc_CH3CH1_lt_67_df):
        #     # just use molCH3CH1_lt_67_df to assign the CH3 and CH1 values
        #     mol_CH3CH1_df = mol_CH3CH1_lt_67_df
        # elif len(mol_sym_CH3CH1_lt_67_df) == len(hsqc_CH3CH1_lt_67_df):
        #     mol_CH3CH1_df = mol_sym_CH3CH1_lt_67_df

        # # # if length of mol_CH3CH1_lt_67_df == length of hsqc_CH3CH1_lt_67_df then we can assign the CH3 values
        # if len(mol_CH3CH1_lt_67_df[mol_CH3CH1_lt_67_df.CH3]) == len(
        #     hsqc_CH3CH1_lt_67_df
        # ):
        #     idx = hsqc_CH3CH1_lt_67_df.index
        #     hsqc.loc[idx, "CH3"] = True
        #     hsqc.loc[idx, "numProtons"] = 3
        # else:
        # # we have a mixture of CH3 and CH1 values, therefore we need to assign the CH3 and CH1 values based on the closest of the ppm values
        # # eliminate the CH3 values from the hsqc dataframe
        # mol_ch3_df = mol.molprops_df[mol.molprops_df.CH3]

        # mol_ch3ch1_df = mol.molprops_df[mol.molprops_df.CH3CH1]

        # mol_ch1_df = mol.molprops_df[mol.molprops_df.CH1]

        # # if no CH3 groups in expected molecule assign all CH3CH1 groups to CH1
        # if mol_ch3_df.shape[0] == 0 and mol_ch1_df.shape[0] > 0:
        #     hsqc.loc[hsqc.CH3CH1.index, "CH1"] = True

        #     hsqc.loc[hsqc.CH1.index, "numProtons"] = 1

        # # if no CH1 groups in expected molecule assign all CH3CH1 groups to CH3
        # elif mol_ch3_df.shape[0] > 0 and mol_ch1_df.shape[0] == 0:
        #     hsqc.loc[hsqc.CH3CH1.index, "CH3"] = True
        #     hsqc.loc[hsqc.CH3.index, "numProtons"] = 3

        # #     # sort mol_CH3CH1_lt_67_df in descending order based on ppm
        # mol_CH3CH1_lt_67_df = mol_CH3CH1_lt_67_df.sort_values(
        #     by=["ppm"], ascending=False
        # ).reset_index(drop=True)

        # # sort hsqc_CH3CH1_lt_67_df in descending order based on f1_ppm
        # hsqc_CH3CH1_lt_67_df = hsqc_CH3CH1_lt_67_df.sort_values(
        #     by=["f1_ppm"], ascending=False
        # )

        # # if the length  of the two dataframes are the not the same then we may a problem raise an error
        # if len(mol_CH3CH1_lt_67_df) != len(hsqc_CH3CH1_lt_67_df):
        #     # check if we use the fullmol.molprops_df are the same length as the hsqc dataframe
        #     mol_CH3CH1_lt_67_df = mol.molprops_df[mol.molprops_df.CH3CH1]
        #     mol_CH3CH1_lt_67_df = mol_CH3CH1_lt_67_df[mol_CH3CH1_lt_67_df.ppm <= 67]

        #     print("len(mol_CH3CH1_lt_67_df) = ", len(mol_CH3CH1_lt_67_df))
        #     print("len(hsqc_CH3CH1_lt_67_df) = ", len(hsqc_CH3CH1_lt_67_df))

        #     print("\nmol.molprops_df\n")
        #     print(mol.molprops_df[["ppm", "CH3", "CH2", "CH1", "CH0", "CH3CH1"]])

        #     print("\nmol_CH3CH1_lt_67_df\n")
        #     print(mol_CH3CH1_lt_67_df[["ppm", "CH3", "CH1", "CH3CH1"]])

        #     print("\nhsqc_CH3CH1_lt_67_df\n")
        #     print(
        #         hsqc_CH3CH1_lt_67_df[
        #             [
        #                 "f1_ppm",
        #                 "f2_ppm",
        #                 "atom_idx",
        #                 "CH3",
        #                 "CH1",
        #                 "CH3CH1",
        #                 "numProtons",
        #             ]
        #         ]
        #     )

        #     if len(mol_CH3CH1_lt_67_df) != len(hsqc_CH3CH1_lt_67_df):
        #         print(
        #             "611 ::: Length of mol_CH3CH1_lt_67_df is not the same as hsqc_CH3CH1_lt_67_df"
        #         )
        #         raise ValueError

        #         # print("Length of mol_CH3CH1_lt_67_df is not the same as hsqc_CH3CH1_lt_67_df")
        #         # raise ValueError

        # CH3 = mol_CH3CH1_lt_67_df.CH3.values
        # CH1 = mol_CH3CH1_lt_67_df.CH1.values

        # # assign the CH3 values to the hsqc dataframe
        # for idx, CH3_val in zip(hsqc_CH3CH1_lt_67_df.index, CH3):
        #     hsqc.loc[idx, "CH3"] = CH3_val

        # # assign the CH1 values to the hsqc dataframe
        # for idx, CH1_val in zip(hsqc_CH3CH1_lt_67_df.index, CH1):
        #     hsqc.loc[idx, "CH1"] = CH1_val

        # # set the numProtons values to 3 for CH3 values and 1 for CH1 values
        # hsqc.loc[hsqc[hsqc.CH3].index, "numProtons"] = 3
        # hsqc.loc[hsqc[hsqc.CH1].index, "numProtons"] = 1
        # # hsqc[hsqc.CH1]["numProtons"] = 1
        # # hsqc_CH3CH1_lt_67_df["CH3"] = CH3

    c13["x"] = c13["x"].astype(float)
    c13["y"] = c13["y"].astype(float)
    h1["x"] = h1["x"].astype(float)
    h1["y"] = h1["y"].astype(float)
    hsqc["x"] = hsqc["x"].astype(float)
    hsqc["y"] = hsqc["y"].astype(float)
    hmbc["x"] = hmbc["x"].astype(float)
    hmbc["y"] = hmbc["y"].astype(float)
    cosy["x"] = cosy["x"].astype(float)
    cosy["y"] = cosy["y"].astype(float)
    ddept["x"] = ddept["x"].astype(float)
    ddept["y"] = ddept["y"].astype(float)
    hsqc_clipcosy["x"] = hsqc_clipcosy["x"].astype(float)
    hsqc_clipcosy["y"] = hsqc_clipcosy["y"].astype(float)

    # print("hsqc", len(hsqc))
    # print("mol.sym_molprops_df", len(mol.sym_molprops_df))
    # print("mol.molprops_df", len(mol.molprops_df))

    # if len(hsqc) > len(mol.sym_molprops_df):
    #     assign_rows(
    #         hsqc,
    #         mol.molprops_df,
    #         "f1_ppm",
    #         "ppm",
    #         "CH2",
    #         "atom_idx",
    #         "atomNumber",
    #         "x",
    #         "y",
    #         "sym_atom_idx",
    #         "sym_atomNumber",
    #     )
    #     assign_rows(
    #         hsqc,
    #         mol.molprops_df,
    #         "f1_ppm",
    #         "ppm",
    #         "CH1",
    #         "atom_idx",
    #         "atomNumber",
    #         "x",
    #         "y",
    #         "sym_atom_idx",
    #         "sym_atomNumber",
    #     )
    #     assign_rows(
    #         hsqc,
    #         mol.molprops_df,
    #         "f1_ppm",
    #         "ppm",
    #         "CH3",
    #         "atom_idx",
    #         "atomNumber",
    #         "x",
    #         "y",
    #         "sym_atom_idx",
    #         "sym_atomNumber",
    #     )
    # else:
    #     assign_rows(
    #         hsqc,
    #         mol.sym_molprops_df,
    #         "f1_ppm",
    #         "ppm",
    #         "CH2",
    #         "atom_idx",
    #         "atomNumber",
    #         "x",
    #         "y",
    #         "sym_atom_idx",
    #         "sym_atomNumber",
    #     )
    #     assign_rows(
    #         hsqc,
    #         mol.sym_molprops_df,
    #         "f1_ppm",
    #         "ppm",
    #         "CH1",
    #         "atom_idx",
    #         "atomNumber",
    #         "x",
    #         "y",
    #         "sym_atom_idx",
    #         "sym_atomNumber",
    #     )
    #     assign_rows(
    #         hsqc,
    #         mol.sym_molprops_df,
    #         "f1_ppm",
    #         "ppm",
    #         "CH3",
    #         "atom_idx",
    #         "atomNumber",
    #         "x",
    #         "y",
    #         "sym_atom_idx",
    #         "sym_atomNumber",
    #     )

    h1 = hsqc.copy()
    h1["ppm"] = h1["f2_ppm"]

    # set correct dtypes for columns x and y in c13 and h1 dataframes to float
    # update c13 dataframe with the atom_idx and atomNumber values from hsqc
    for idx, row in c13.iterrows():
        hsqc_rows = hsqc[hsqc["f1_ppm"] == row["ppm"]]
        if len(hsqc_rows) > 0:
            c13.loc[idx, "CH3"] = hsqc_rows["CH3"].values[0]
            c13.loc[idx, "CH2"] = hsqc_rows["CH2"].values[0]
            c13.loc[idx, "CH1"] = hsqc_rows["CH1"].values[0]
            c13.loc[idx, "atom_idx"] = hsqc_rows["atom_idx"].values[0]
            c13.loc[idx, "atomNumber"] = hsqc_rows["atomNumber"].values[0]
            c13.loc[idx, "x"] = hsqc_rows["x"].values[0]
            c13.loc[idx, "y"] = hsqc_rows["y"].values[0]
            c13.loc[idx, "sym_atom_idx"] = hsqc_rows["sym_atom_idx"].values[0]
            c13.loc[idx, "sym_atomNumber"] = hsqc_rows["sym_atomNumber"].values[0]
            c13.loc[idx, "numProtons"] = hsqc_rows["numProtons"].values[0]

    # finally assgign CH0 values from mol.molprops_df to c13_1d based on the ppm values and the atom_idx being -1 in the c13_1d dataframe
    c13.loc[c13.atom_idx == -1, "CH0"] = True
    c13.loc[c13.atom_idx == -1, "numProtons"] = 0

    ch0_rows = (
        mol.molprops_df[mol.molprops_df.CH0]
        .sort_values(by=["ppm"], ascending=False)
        .reset_index(drop=True)
    )

    # copy over the atom_idx and atomNumber values from ch0_rows to c13_1d
    c13_idx = c13[c13.CH0].index

    print("mol.molprops_df")
    print(
        mol.molprops_df[
            ["ppm", "atom_idx", "atomNumber", "numProtons", "CH3", "CH2", "CH1", "CH0"]
        ]
    )

    print("c13")
    print(
        c13[["ppm", "atom_idx", "atomNumber", "numProtons", "CH3", "CH2", "CH1", "CH0"]]
    )

    print("hsqc")
    print(
        hsqc[
            [
                "f1_ppm",
                "f2_ppm",
                "atom_idx",
                "atomNumber",
                "numProtons",
                "CH3",
                "CH2",
                "CH1",
                "CH0",
            ]
        ]
    )

    print("c13[c13.CH0]")
    print(c13[c13.CH0])

    c13.loc[c13_idx, "atom_idx"] = ch0_rows["atom_idx"].values
    c13.loc[c13_idx, "atomNumber"] = ch0_rows["atomNumber"].values
    c13.loc[c13_idx, "x"] = ch0_rows["x"].values
    c13.loc[c13_idx, "y"] = ch0_rows["y"].values
    c13.loc[c13_idx, "sym_atom_idx"] = ch0_rows["sym_atom_idx"].values
    c13.loc[c13_idx, "sym_atomNumber"] = ch0_rows["sym_atomNumber"].values

    # change type of atom_idx and atomNumber to int
    c13["atom_idx"] = c13["atom_idx"].astype(int)
    c13["atomNumber"] = c13["atomNumber"].astype(int)

    # use hsqc to update c13 dataframe with the atom_idx and atomNumber values
    # assign_rows(c13, hsqc, "ppm", "f1_ppm", "CH3", "atom_idx", "atomNumber", "x", "y", "sym_atom_idx", "sym_atomNumber")

    return h1, c13, hsqc


def assign_jcouplings(df1, df2):
    # check if h1_1d is the same length as h1 and proceed
    if len(df1) != len(df2):
        return df1

    # sort the dataframes by ppm in descending order
    if ("ppm" not in df1.columns) and ("f2_ppm" not in df1.columns):
        return df1

    if ("ppm" not in df2.columns) and ("f2_ppm" not in df2.columns):
        return df1

    if "ppm" in df1.columns:
        df1 = df1.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
    elif "f2_ppm" in df1.columns:
        df1 = df1.sort_values(by=["f2_ppm"], ascending=False).reset_index(drop=True)

    if "ppm" in df2.columns:
        df2 = df2.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
    elif "f2_ppm" in df2.columns:
        df2 = df2.sort_values(by=["f2_ppm"], ascending=False).reset_index(drop=True)

    df1["jCouplingVals"] = df2["jCouplingVals"].values
    df1["jCouplingClass"] = df2["jCouplingClass"].values

    return df1


def assign_jcouplings_to_c13(c13, hsqc):

    for hsqc_idx, row in hsqc.iterrows():
        # find the row in c13 with the same atom_idx
        c13_rows = c13[c13.atom_idx == row["atom_idx"]]
        c13_idx = c13_rows.index
        c13.loc[c13_idx, "jCouplingVals"] = row["jCouplingVals"]
        c13.loc[c13_idx, "jCouplingClass"] = row["jCouplingClass"]

    return c13


def create_expected_molecule_from_molfilestring(
    molfilestr: str,
    carbonAtomsInfo: pd.DataFrame,
    df: pd.DataFrame,
    prediction_from_nmrshiftdb: bool,
) -> expectedmolecule.expectedMolecule:

    """
    Create an RDKit molecule object from a molfile string and a pandas dataframe
    """

    mol = expectedmolecule.expectedMolecule(
        molfilestr,
        is_smiles=False,
        carbonAtomsInfo=carbonAtomsInfo,
        mnova_c13predictions=df,
        predict_from_nmrshiftdb=prediction_from_nmrshiftdb,
    )

    svg_str, new_xy3 = create_svg_string(
        mol.mol,
        molWidth=svgDimensions.MOLWIDTH,
        molHeight=svgDimensions.MOLHEIGHT,
        svgWidth=svgDimensions.SVGWIDTH,
        svgHeight=svgDimensions.SVGHEIGHT,
    )

    mol.svg_str = svg_str
    mol.xy3 = new_xy3
    mol.molfile = molfilestr
    # update xy3 values in mol2
    expectedmolecule.copy_xy3_to_molecule(mol, new_xy3)

    return mol


def expand_c13_dataframe(c13: pd.DataFrame, mol: expectedmolecule.expectedMolecule):
    """
    Expand the c13 dataframe to include the number of protons attached to the carbon atom
    """

    if len(c13) == len(mol.molprops_df):
        return c13

    c13_sym = c13[c13["sym_atom_idx"] != ""].copy()

    # swap atom_idx and sym_atom_idx values
    c13_sym["atom_idx"], c13_sym["sym_atom_idx"] = (
        c13_sym["sym_atom_idx"],
        c13_sym["atom_idx"],
    )
    c13_sym["atomNumber"], c13_sym["sym_atomNumber"] = (
        c13_sym["sym_atomNumber"],
        c13_sym["atomNumber"],
    )

    # append c13_sym to c13
    # c13_all = c13.append(c13_sym, ignore_index=True)
    c13_all = pd.concat([c13, c13_sym], ignore_index=True)

    # set atomNumber and atom_idx type to int
    c13_all["atomNumber"] = c13_all["atomNumber"].astype(int)
    c13_all["atom_idx"] = c13_all["atom_idx"].astype(int)

    # set sym_atomNumber and sym_atom_idx type to str
    c13_all["sym_atomNumber"] = c13_all["sym_atomNumber"].astype(str)
    c13_all["sym_atom_idx"] = c13_all["sym_atom_idx"].astype(str)

    # update the x and y coord in c13_all from mol.molprops_df using the atom_idx
    for idx, row in c13_all.iterrows():
        atom_idx = row["atom_idx"]
        mol_row = mol.molprops_df[mol.molprops_df.atom_idx == atom_idx]
        if len(mol_row) == 0:
            continue
        c13_all.at[idx, "x"] = mol_row["x"].values[0]
        c13_all.at[idx, "y"] = mol_row["y"].values[0]

    # sort the dataframe by ppm in descending order descending order, reset the index and start the index at 1
    c13_all = c13_all.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
    c13_all.index += 1

    return c13_all


def expand_h1_dataframe(h1, mol):
    # expand h1 dataframe by adding new rows where sym_atom_idx is not ""
    h1_sym = h1[h1["sym_atom_idx"] != ""].copy()

    # swap atom_idx and sym_atom_idx values
    h1_sym["atom_idx"], h1_sym["sym_atom_idx"] = (
        h1_sym["sym_atom_idx"],
        h1_sym["atom_idx"],
    )
    h1_sym["atomNumber"], h1_sym["sym_atomNumber"] = (
        h1_sym["sym_atomNumber"],
        h1_sym["atomNumber"],
    )

    # append h1_sym to h1
    # h1_all = h1.append(h1_sym, ignore_index=True)
    h1_all = pd.concat([h1, h1_sym], ignore_index=True)

    # convert atom_idx and atomNumber to int
    h1_all["atom_idx"] = h1_all["atom_idx"].astype(int)
    h1_all["atomNumber"] = h1_all["atomNumber"].astype(int)
    h1_all["sym_atom_idx"] = h1_all["sym_atom_idx"].astype(str)
    h1_all["sym_atomNumber"] = h1_all["sym_atomNumber"].astype(str)

    # update the x and y coord in h1_all from mol.molprops_df using the atom_idx
    for idx, row in h1_all.iterrows():
        atom_idx = row["atom_idx"]
        mol_row = mol.molprops_df[mol.molprops_df.atom_idx == atom_idx]
        h1_all.at[idx, "x"] = mol_row["x"].values[0]
        h1_all.at[idx, "y"] = mol_row["y"].values[0]

    # sort the dataframe by ppm in descending order, reset the index and start the index at 1
    h1_all = h1_all.sort_values(by=["ppm"], ascending=False).reset_index(drop=True)
    h1_all.index += 1

    return h1_all


def create_network_graph_nodes(c13: pd.DataFrame, h1: pd.DataFrame):
    """
    Create a network graph from the c13 and h1 dataframes
    """
    c13_sorted = c13.sort_values(by=["atom_idx"], ascending=True).reset_index(drop=True)
    c13_sorted["f1_ppm"] = c13_sorted["ppm"]
    c13_sorted["iupacLabel"] = ""
    G2 = create_network_graph(c13_sorted, h1)

    return G2


def add_all_cosy_edges_to_graph(
    G2, cosy: pd.DataFrame, clipcosy: pd.DataFrame, h1_all: pd.DataFrame
):
    """
    Add edges to the graph G2 from the cosy dataframe
    """
    cosy_combinations = create_cosy_combinations(cosy, h1_all)

    clipcosy_combinations = create_clipcosy_combinations(clipcosy, h1_all)

    kept_cosy_combinations = remove_duplicate_cosy_combinations(cosy_combinations)
    kept_clipcosy_combinations = remove_duplicate_cosy_combinations(
        clipcosy_combinations
    )

    G2 = add_cosy_edges_to_graph(G2, kept_cosy_combinations)
    G2 = add_cosy_edges_to_graph(G2, kept_clipcosy_combinations)

    return G2


def add_all_hmbc_edges_to_graph(
    G2, hmbc: pd.DataFrame, h1_all: pd.DataFrame, c13_sorted: pd.DataFrame
):
    """
    Add edges to the graph G2 from the hmbc dataframe
    """
    G2 = add_hmbc_edges_to_graph(G2, hmbc, h1_all, c13_sorted)

    return G2


def complete_jcouplings(h1_1d, h1, c13, hsqc):

    h1 = assign_jcouplings(h1, h1_1d)
    hsqc = assign_jcouplings(hsqc, h1)
    c13 = assign_jcouplings_to_c13(c13, hsqc)

    # update ppm columns
    if "ppm" not in h1.columns and "f2_ppm" in h1.columns:
        h1["ppm"] = h1["f2_ppm"]

    if "f2_ppm" not in h1.columns and "ppm" in h1.columns:
        h1["f2_ppm"] = h1["ppm"]

    if "ppm" not in c13.columns and "f1_ppm" in c13.columns:
        c13["ppm"] = c13["f1_ppm"]

    if "f1_ppm" not in c13.columns and "ppm" in c13.columns:
        c13["f1_ppm"] = c13["ppm"]

    return h1, c13, hsqc


def html_from_predictions(nmrproblem):
    """
    Create an html page from the nmrproblem object
    """

    molfilestr = nmrproblem.dataframes["molfile"].loc[0, "molfile"]
    print(
        "nmrprpblem.prediction_from_nmrshiftdb2",
        nmrproblem.prediction_from_nmrshiftdb2(),
    )
    mol = create_expected_molecule_from_molfilestring(
        molfilestr,
        nmrproblem.dataframes["carbonAtomsInfo"],
        nmrproblem.dataframes["c13predictions"],
        nmrproblem.prediction_from_nmrshiftdb2(),
    )

    print(
        mol.molprops_df[
            ["ppm", "atom_idx", "atomNumber", "numProtons", "CH3", "CH2", "CH1", "CH0"]
        ]
    )
    # nmrAssignments = nmrproblem.dataframes["nmrAssignments"].copy()
    cosy = nmrproblem.dataframes["COSY"].copy()
    hmbc = nmrproblem.dataframes["HMBC"].copy()
    hsqc = nmrproblem.dataframes["HSQC"].copy()
    h1_1d = nmrproblem.dataframes["H1_1D"].copy()
    c13_1d = nmrproblem.dataframes["C13_1D"].copy()
    dept = nmrproblem.dataframes["DEPT135"].copy()
    noesy = nmrproblem.dataframes["NOESY"].copy()
    h1_pureshift = nmrproblem.dataframes["H1_pureshift"].copy()
    ddept_ch3_only = nmrproblem.dataframes["DDEPT_CH3_ONLY"].copy()
    hsqc_clipcosy = nmrproblem.dataframes["HSQC_CLIPCOSY"].copy()

    # check if mol has CH2 in dataframe then does hsqc have some negative intensities
    if (mol.num_CH2_carbon_atoms > 0) and (hsqc[hsqc.intensity < 0].shape[0] == 0):

        if len(dept) == 0:
            return (
                None,
                None,
                "hsqc needs to be multiplicity edited experiment as no DEPT data found",
            )

        # set intensity of CH2 values in hsqc to negative by using dept values
        dept_CH2 = dept[dept.intensity < 0].copy()

        if len(dept_CH2) == 0:
            return (
                None,
                None,
                "hsqc needs to be multiplicity edited experiment as no DEPT CH2 data found",
            )

        # sort dept_CH2 by ppm in descending order
        dept_CH2 = dept_CH2.sort_values(by=["ppm"], ascending=False).reset_index(
            drop=True
        )

        # sort hsqc by f1_ppm in descending order
        hsqc = hsqc.sort_values(by=["f1_ppm"], ascending=False).reset_index(drop=True)

        # loop through dept_CH2, find closest ppm value in hsqc and set the intensity to negative
        tol = CARBONSEPARATION
        for v1 in dept_CH2["ppm"]:

            print("v1", v1)

            hsqc["prob"] = hsqc.apply(
                lambda x: stats.norm(v1, tol).pdf(x["f1_ppm"]), axis=1
            )
            # set the intensity to -1 for all prob values greater than 0.0
            hsqc.loc[hsqc["prob"] > 0.0, "intensity"] = -1.0
            hsqc["prob"] = 0.0

    print("\n===================")
    print("hsqc")
    print("===================\n")

    print(hsqc)

    if (
        len(c13_1d) == 0
        or len(c13_1d) > mol.num_carbon_atoms
        or len(c13_1d) < mol.num_sym_carbon_atoms
    ):
        h1, c13, hsqc = complete_dataframes_hsqc_only(
            mol,
            hsqc,
            hmbc=hmbc,
            ddept=ddept_ch3_only,
            cosy=cosy,
            hsqc_clipcosy=hsqc_clipcosy,
        )
    elif len(c13_1d) == mol.num_carbon_atoms or len(c13_1d) == mol.num_sym_carbon_atoms:
        h1, c13, hsqc = complete_dataframes_c13_hsqc(
            mol,
            hsqc,
            c13_1d,
            hmbc=hmbc,
            ddept=ddept_ch3_only,
            cosy=cosy,
            hsqc_clipcosy=hsqc_clipcosy,
        )
    else:
        return None, None, "error in c13_1d dataframe"

    h1, c13, hsqc = complete_jcouplings(h1_1d, h1, c13, hsqc)

    print("c13")
    print(
        c13[
            [
                "ppm",
                "atom_idx",
                "sym_atom_idx",
                "atomNumber",
                "numProtons",
                "CH3",
                "CH2",
                "CH1",
                "CH0",
            ]
        ]
    )

    c13_all = expand_c13_dataframe(c13, mol)
    print("c13_all")
    print(
        c13_all[
            [
                "ppm",
                "atom_idx",
                "sym_atom_idx",
                "atomNumber",
                "sym_atomNumber",
                "numProtons",
                "CH3",
                "CH2",
                "CH1",
                "CH0",
            ]
        ]
    )

    # c13_all = c13.copy()

    print("c13_all")
    print(
        c13_all[
            [
                "ppm",
                "atom_idx",
                "sym_atom_idx",
                "atomNumber",
                "numProtons",
                "CH3",
                "CH2",
                "CH1",
                "CH0",
                "x",
                "y",
            ]
        ]
    )

    # add ppm values from mol.molprops_df to c13_all using atom_idx
    for idx, row in c13_all.iterrows():
        atom_idx = row["atom_idx"]
        mol_row = mol.molprops_df[mol.molprops_df.atom_idx == atom_idx]
        if len(mol_row) == 0:
            continue
        c13_all.at[idx, "ppm_calculated"] = mol_row["ppm"].values[0]

    h1_all = expand_h1_dataframe(h1, mol)
    print("h1_all")
    print(
        h1_all[
            [
                "ppm",
                "atom_idx",
                "sym_atom_idx",
                "atomNumber",
                "numProtons",
                "CH3",
                "CH2",
                "CH1",
                "CH0",
            ]
        ]
    )

    # h1_all = h1.copy()

    G2 = create_network_graph_nodes(c13_all, h1_all)

    G2 = add_all_cosy_edges_to_graph(G2, cosy, hsqc_clipcosy, h1_all)

    G2 = add_all_hmbc_edges_to_graph(G2, hmbc, h1_all, c13_all)

    jsonGraphData = json_graph.node_link_data(G2)

    svgWidth = svgDimensions.SVGWIDTH
    svgHeight = svgDimensions.SVGHEIGHT
    molWidth = svgDimensions.MOLWIDTH
    molHeight = svgDimensions.MOLHEIGHT

    jinjadata = {
        "svg_container": mol.svg_str,
        "graph_edges": jsonGraphData["links"],
        "graph_nodes": jsonGraphData["nodes"],
        "translateX": int((svgWidth - molWidth) / 2),
        "translateY": int((svgHeight - molHeight) / 2),
        "title": "dummy_title",
    }

    return jinjadata, G2, "ok"


if __name__ == "__main__":

    import platform
    from pathlib import Path
    from html_from_assignments import create_htmlpage_from_graph

    print("Running buildNMRdataframes.py")

    if platform.system() == "Windows":
        fn_json = Path(
            r"C:\Users\vsmw51\Dropbox\PC\Downloads\exampleProblems\dimethylbenzene\dimethylbenzene_assignments_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/predictC13test/predictC13test_predictions_mresnova.json"
        )
        fn_json = Path(
            r"C:\Users\vsmw51\Dropbox\PC\Downloads\exampleProblems\menthol_eeh_predict\menthol_eeh_predict_assignments_mresnova.json"
        )
        fn_json = Path(
            r"C:\Users\vsmw51\Dropbox\PC\Downloads\exampleProblems\EVB_330b_annotated\EVB_330b_annotated_assignments_mresnova.json"
        )
        fn_json = Path(
            r"C:\Users\vsmw51\Dropbox\PC\Downloads\exampleProblems\menthol_eeh_predict_and_assign\menthol_eeh_predict_and_assign_assignments_mresnova.json"
        )

    elif platform.system() == "Darwin":
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/dimethylbenzene/dimethylbenzene_assignments_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/predictC13test/predictC13test_predictions_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/menthol_eeh_predict\menthol_eeh_predict_assignments_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2024/python/awh/exampleProblems/EVB_330b_annotated/EVB_330b_annotated_assignments_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2022/python/python/rdkit37/projects/simpleNMR/exampleProblems/EVB_330b_annotated/EVB_330b_annotated_assignments_mresnova.json"
        )
        fn_json = Path(
            r"/Users/vsmw51/OneDrive - Durham University/projects/programming/2022/python/python/rdkit37/projects/simpleNMR/exampleProblems/santonin/santonin_assignments_mresnova.json"
        )

    if fn_json.exists():
        problemdata_json = NMRProblem.from_mnova_json_file(fn_json)
        for k, v in problemdata_json.dataframes.items():
            print(k)
    else:
        print(f"{fn_json.name} not found")

    jinjadata, G2, ok = html_from_predictions(problemdata_json)

    from html_from_assignments import create_htmlpage_from_graph

    create_htmlpage_from_graph(
        G2, "test.html", r"templates/d3molplotmnova_template.html", jinjadata
    )
