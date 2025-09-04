
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment
import networkx as nx
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from typing import Tuple, Union, List, Dict, Any

from flask import render_template


from rdkit.Chem import Draw

# from excelheaders import excel_orig_df_columns, excel_df_columns
# from html_from_assignments import tidyup_ppm_values
from html_from_assignments import NMRProblem
import expectedmolecule

def warning_dialog(return_message, title_message, qstarted=True):
    print(return_message)


class NMRsolution:
    def __init__(self, problemdata_json: NMRProblem):

        self.exact_ppm_values = problemdata_json.exact_ppm_values

        self.qtstarted = False

        self.nmrsolution_failed = False
        self.nmrsolution_error_message = "ok"
        self.nmrsolution_error_code = 200

        self.numtimes_HSQC_CH2_set_to_FALSE = 0

        self.problemdata_json = problemdata_json

        self.expts_available = set(problemdata_json.dataframes["chosenSpectra"]["expt"])

        print("self.expts_available", self.expts_available)

        # remove "SKIP" from expts_available
        if "SKIP" in self.expts_available:
            self.expts_available.remove("SKIP")

        if "HSQC" not in self.expts_available:
            self.nmrsolution_failed = True
            self.nmrsolution_error_message = "<p>HSQC experiment not present</p>"
            self.nmrsolution_error_code = 401
            return

        self.hsqc_df = problemdata_json.dataframes["HSQC"]
        self.hmbc_df = problemdata_json.dataframes["HMBC"]
        self.cosy_df = problemdata_json.dataframes["COSY"]
        self.dept_df = problemdata_json.dataframes["DEPT135"]
        self.ddept_ch3_only_df = problemdata_json.dataframes["DDEPT_CH3_ONLY"]
        self.hsqc_clipcosy_df = problemdata_json.dataframes["HSQC_CLIPCOSY"]
        self.h1_df = problemdata_json.dataframes["H1_1D"]
        self.pureshift_df = problemdata_json.dataframes["H1_pureshift"]
        self.c13_df = problemdata_json.dataframes["C13_1D"]

        print("self.c13_df\n", self.c13_df)

        # if hsqc_df is empty then return
        if self.hsqc_df.empty:
            self.nmrsolution_failed = True
            self.nmrsolution_error_message = "<p>HSQC experiment contains no data</p>"
            self.nmrsolution_error_code = 401
            return

        self.molstr = problemdata_json.dataframes["molfile"].loc[0, "molfile"]

        self.smilesstr = problemdata_json.dataframes["smiles"].loc[0, "smiles"]
        self.carbonAtomsInfo = problemdata_json.dataframes["carbonAtomsInfo"]
        self.c13predictions = problemdata_json.dataframes["c13predictions"]

        self.expected_molecule = self.create_expected_molecule_from_molfilestring(
            self.molstr,
            self.carbonAtomsInfo,
            self.c13predictions,
            self.problemdata_json.prediction_from_nmrshiftdb2(),
            self.problemdata_json.JEOL_prediction_used()
        )

        self.solution_error_message = self.expected_molecule.nmrshiftdb_failed_message
        self.solution_error_code = self.expected_molecule.nmrshiftdb_failed_code

        # check if we need to use dept to assign CH2 in hsqc or has the user annotated the CH2 in the hsqc dataframe
        # print("nmrsolution line 108")
        # print("self.hsqc_df before check_hsqc_assign_CH2_from_DEPT\n", self.hsqc_df)

        print("self.hsqc_df.columns\n", self.hsqc_df.columns)

        if (
            (self.expected_molecule.num_CH2_carbon_atoms > 0)
            and (self.hsqc_df[self.hsqc_df.intensity < 0].shape[0] == 0)
            and (self.hsqc_df[self.hsqc_df.Annotation == "CH2"].shape[0] == 0)
        ):

            self.check_hsqc_assign_CH2_from_DEPT135()

        elif (
            (self.expected_molecule.num_CH2_carbon_atoms > 0)
            and (self.hsqc_df[self.hsqc_df.intensity < 0].shape[0] == 0)
            and (self.hsqc_df[self.hsqc_df.Annotation == "CH2"].shape[0] > 0)
        ):

            # multiply the intensity and the integral by -1 in the hsqc_df where CH2 is present in the annotation column
            self.hsqc_df.loc[self.hsqc_df.Annotation == "CH2", "intensity"] = -1.0
            self.hsqc_df.loc[self.hsqc_df.Annotation == "CH2", "integral"] = -1.0

        # print("self.hsqc_df after check_hsqc_assign_CH2_from_DEPT\n", self.hsqc_df)
        # print(self.hsqc_df)

    def initiate_molgraph(self, json_data, G2):
        # Create a new NetworkX graph
        molgraph = nx.Graph()

        # if rubteralone remove - 1 from the atom numbers
        if (json_data["MNOVAcalcMethod"]["data"]["0"] == "NMRSHIFTDB2 Predict") or (json_data["MNOVAcalcMethod"]["data"]["0"] == "JEOL Predict"):
            nodes_offset = 0
        else:
            nodes_offset = 1

        # Add nodes for each atom in the molecule
        for atom in self.expected_molecule.GetAtoms():
            molgraph.add_node(
                atom.GetIdx(),
                symbol=atom.GetSymbol(),
                atom_number=atom.GetAtomicNum(),
                numProtons=atom.GetNumImplicitHs(),
                id0=atom.GetIdx(),
                x=self.expected_molecule.xy3_allatoms[atom.GetIdx()][0],
                y=self.expected_molecule.xy3_allatoms[atom.GetIdx()][1],
            )

        # Add edges for each bond in the molecule
        for bond in self.expected_molecule.GetBonds():
            molgraph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        # add the ppm values from the c13predictions in the json data

        # set the 'ppm' attribute for each node in the molgraph to 10000 to begin with
        for node in molgraph.nodes():
            molgraph.nodes[node]["ppm"] = 10000
            # molgraph.nodes[node]["numProtons"] = -1

        # c13predictions = json_data["oldjsondata"]["c13predictions"]["data"]

        for node in molgraph.nodes():

            G2_node = node + nodes_offset

            if G2_node in G2.nodes():
                molgraph.nodes[node]["ppm"] = G2.nodes[G2_node]["ppm"]
                # molgraph.nodes[node]["numProtons"] = G2.nodes[node]["numProtons"]
                # check if the atmomNumber can be converted to an integer if not dont try
                if isinstance(G2.nodes[G2_node]["atomNumber"], (float, int)):
                    molgraph.nodes[node]["atomNumber"] = int(
                        G2.nodes[G2_node]["atomNumber"]
                    )
                elif isinstance(G2.nodes[G2_node]["atomNumber"], str):
                    try:
                        molgraph.nodes[node]["atomNumber"] = int(
                            G2.nodes[G2_node]["atomNumber"]
                        )
                    except:
                        molgraph.nodes[node]["atomNumber"] = G2.nodes[G2_node][
                            "atomNumber"
                        ]
                molgraph.nodes[node]["id"] = int(node)

        self.molgraph = molgraph

    def assign_CH3_CH2_CH1_overall(self):

        self.assign_CH3_CH2_CH1_in_HSQC_using_Assignments()

        # check if everything has been assigned, no negative values left in numProtons

        if self.hsqc[self.hsqc.numProtons == -1].empty:
            # return True, "ok"
            return "ok", 200

        #  attempt ro remove assign CH3 CH2 CH1 from hsqc using different methods

        if self.hsqc[self.hsqc.numProtons == -1].shape[0] > 0:
            self.assign_CH3_CH1_in_HSQC_using_H1()

        if self.hsqc[self.hsqc.numProtons == -1].shape[0] > 0:
            self.assign_CH3_CH1_in_HSQC_using_expected_molecule()

        if self.hsqc[self.hsqc.numProtons == -1].shape[0] > 0:
            print("solution failed")
            self.nmrsolution_failed = True

        if self.nmrsolution_failed:
            print(self.nmrsolution_error_message, 400)
            return self.nmrsolution_error_message, 400
        else:
            print("ok", 200)
            return "ok", 200

    def initialise_prior_to_carbon_assignment(self):

        self.c13 = self.c13.sort_values("ppm", ascending=False, ignore_index=True)
        c13 = self.c13

        self.all_molprops_df = self.expected_molecule.molprops_df.sort_values(
            "ppm", ascending=False, ignore_index=True
        ).copy()
        self.sym_molprops_df = self.expected_molecule.sym_molprops_df.sort_values(
            "ppm", ascending=False, ignore_index=True
        ).copy()

        error_msg = ""
        if len(c13) == len(self.all_molprops_df):
            error_msg = "<p>c13 and all_molprops_df are the same length</p>"
            self.molprops_df = self.all_molprops_df

        elif len(c13) == len(self.sym_molprops_df):
            error_msg = error_msg + "<p>c13 and sym_molprops_df are the same length</p>"
            self.molprops_df = self.sym_molprops_df

        elif (len(c13) < len(self.all_molprops_df)) and (
            len(c13) > len(self.sym_molprops_df)
        ):
            self.molprops_df = self.all_molprops_df

        elif len(c13) <= len(self.all_molprops_df[~self.all_molprops_df.CH0]) or len(
            c13
        ) <= len(self.sym_molprops_df[~self.sym_molprops_df.CH0]):
            self.molprops_df = self.all_molprops_df

        elif len(c13) < len(self.all_molprops_df):
            self.molprops_df = self.all_molprops_df

        elif len(c13) > len(self.all_molprops_df):
            # No simple solution right now so create tables to show the differences
            # and return it as a html page

            c13_df_list = []
            for CHn_str in ["CH0", "CH1", "CH2", "CH3"]:
                df1_C13 = c13[c13[CHn_str]][
                    ["ppm", "numProtons", "CH0", "CH1", "CH2", "CH3"]
                ].copy()
                # replace True with 1 and False with 0
                df1_C13 = df1_C13.replace({True: 1, False: 0})
                # set ppm value to a string with 2 decimal places
                df1_C13["ppm"] = df1_C13["ppm"].map("{:.2f}".format)
                c13_df_list.append(df1_C13.values.tolist())

                self.molprops_df1 = self.all_molprops_df[self.all_molprops_df[CHn_str]][
                    ["ppm", "numProtons", "CH0", "CH1", "CH2", "CH3"]
                ].copy()
                # replace True with 1 and False with 0
                self.molprops_df1 = self.molprops_df1.replace({True: 1, False: 0})

                c13_df_list.append(self.molprops_df1.values.tolist())

            tab_headings = ["CH0", "CH1", "CH2", "CH3"]

            rtn_html = render_template(
                "error_table.html",
                colHeadings=["ppm", "numProtons", "CH0", "CH1", "CH2", "CH3"],
                df_data=c13_df_list,
                tab_headings=tab_headings,
            )

            return rtn_html, 400

        return "ok", 200

    def update_assignments_expt_dataframes(self):
        """
        Updates experimental NMR dataframes by synchronizing values from the C13 reference dataframe.

        Copies C13-related values into various experiment dataframes,
        augments the H1 dataframe with atom and coordinate information from HSQC, assigns J-coupling information,
        and resets specific columns for further analysis.
        """
        self.hmbc = copy_over_values_c13_all_to_hetero2D(self.hmbc, self.c13)
        self.hsqc = copy_over_values_c13_all_to_hetero2D(self.hsqc, self.c13)
        self.cosy = copy_over_values_c13_all_to_homo2D(self.cosy, self.c13)
        self.hsqc_clipcosy = copy_over_values_c13_all_to_hetero2D(
            self.hsqc_clipcosy, self.c13
        )
        self.ddept_ch3_only = copy_over_values_c13_all_to_hetero2D(
            self.ddept_ch3_only, self.c13
        )

        # add atom_idx, atomNumber, x and y columns to h1 dataframe using hsqc dataframe

        for idx, row in self.h1.iterrows():
            f2_ppm = row["ppm"]
            hsqc_row = self.hsqc[self.hsqc["f2_ppm"] == f2_ppm]
            if len(hsqc_row) == 0:
                continue
            self.h1.at[idx, "atom_idx"] = hsqc_row["f2_atom_idx"].values[0]
            self.h1.at[idx, "sym_atom_idx"] = hsqc_row["f2_sym_atom_idx"].values[0]
            self.h1.at[idx, "atomNumber"] = hsqc_row["f2_atomNumber"].values[0]
            self.h1.at[idx, "sym_atomNumber"] = hsqc_row["f2_sym_atomNumber"].values[0]
            self.h1.at[idx, "x"] = hsqc_row["f2_x"].values[0]
            self.h1.at[idx, "y"] = hsqc_row["f2_y"].values[0]

        self.hsqc = assign_jcouplings(self.hsqc, self.h1)
        self.c13 = assign_jcouplings_to_c13(self.c13, self.hsqc)

        self.h1["f1_ppm"] = self.h1["ppm"]
        self.c13["f1_ppm"] = self.c13["ppm"]
        self.c13["iupacLabel"] = ""

    def attempt_assignment_CH3_CH2_CH1_to_C13_table(self):
        """
        Attempts to assign CH3, CH2, CH1, and CH0 groups from the molecular properties DataFrame
        (`molprops_df`) to the C13 table (`c13`) based on the number of protons and chemical shift (ppm)
        values. The function iterates over each group type, matches rows between the two tables using
        optimal assignment (Hungarian algorithm) or nearest-neighbor matching, and updates the C13 table
        with corresponding atom indices and properties. If the assignment cannot be completed due to
        mismatched row counts, an error table is rendered for debugging.
        Returns:
            tuple: ("ok", 200) if assignment is successful, or (rendered_html, 400) if an error occurs.
        """

        c13 = self.c13

        for CHn, nProtons in zip(["CH3", "CH2", "CH1", "CH0"], [3, 2, 1, 0]):

            df_CHn = self.molprops_df[self.molprops_df[CHn]]
            c13_CHn = c13[c13["numProtons"] == nProtons]

            if c13_CHn.empty:
                print(f"{CHn} c13_CHn is empty")
                continue

            if len(c13_CHn) == len(df_CHn):
                # order dataframes based on ppm column

                # Compute the pairwise distance matrix
                distance_matrix = np.abs(
                    df_CHn["ppm"].values[:, np.newaxis] - c13_CHn["ppm"].values
                )

                # Find optimal assignment
                row_ind, col_ind = linear_sum_assignment(distance_matrix)

                # Extract matched pairs of indices and total minimum distance
                matches = [
                    (df_CHn.index[i], c13_CHn.index[j])
                    for i, j in zip(row_ind, col_ind)
                ]
                total_distance = distance_matrix[row_ind, col_ind].sum()

                for df_idx, c13_idx in matches:
                    c13.at[c13_idx, "atom_idx"] = df_CHn.at[df_idx, "atom_idx"]
                    c13.at[c13_idx, "sym_atom_idx"] = df_CHn.at[df_idx, "sym_atom_idx"]
                    c13.at[c13_idx, "atomNumber"] = df_CHn.at[df_idx, "atomNumber"]
                    c13.at[c13_idx, "sym_atomNumber"] = df_CHn.at[
                        df_idx, "sym_atomNumber"
                    ]
                    c13.at[c13_idx, "x"] = df_CHn.at[df_idx, "x"]
                    c13.at[c13_idx, "y"] = df_CHn.at[df_idx, "y"]
                    c13.at[c13_idx, "ppm_calculated"] = df_CHn.at[df_idx, "ppm"]

            elif len(c13_CHn) < len(df_CHn):
                # attempt to match the c13 rows individually to the molprops_df rows
                # based on the ppm values

                # calculate the number of unique ppm values in df_CHn
                nunique_ppm = df_CHn.ppm.nunique()

                if len(c13_CHn) == nunique_ppm:
                    # drop duplicates from df_CHn based on ppm
                    df_CHn = df_CHn.drop_duplicates(subset=["ppm"])

                for idx, row in c13_CHn.iterrows():
                    ppm = row["ppm"]

                    # find closet ppm value in df_CHn to ppm
                    mol_row = df_CHn.iloc[(df_CHn["ppm"] - ppm).abs().argsort()[:1]]
                    if len(mol_row) == 0:
                        continue
                    c13.at[idx, "atom_idx"] = mol_row["atom_idx"].values[0]
                    c13.at[idx, "sym_atom_idx"] = mol_row["sym_atom_idx"].values[0]
                    c13.at[idx, "atomNumber"] = mol_row["atomNumber"].values[0]
                    c13.at[idx, "sym_atomNumber"] = mol_row["sym_atomNumber"].values[0]
                    c13.at[idx, "x"] = mol_row["x"].values[0]
                    c13.at[idx, "y"] = mol_row["y"].values[0]

                    c13.at[idx, "ppm_calculated"] = mol_row["ppm"].values[0]

                    # remove the row from df_CHn
                    df_CHn = df_CHn.drop(mol_row.index)

            else:

                c13_df_list = []
                for CHn_str in ["CH0", "CH1", "CH2", "CH3"]:
                    df1_C13 = c13[c13[CHn_str]][
                        ["ppm", "numProtons", "CH0", "CH1", "CH2", "CH3"]
                    ].copy()
                    # replace True with 1 and False with 0
                    df1_C13 = df1_C13.replace({True: 1, False: 0})
                    # set ppm value to a string with 2 decimal places
                    df1_C13["ppm"] = df1_C13["ppm"].map("{:.2f}".format)
                    c13_df_list.append(df1_C13.values.tolist())

                    self.molprops_df1 = self.molprops_df[self.molprops_df[CHn_str]][
                        ["ppm", "numProtons", "CH0", "CH1", "CH2", "CH3"]
                    ].copy()
                    # replace True with 1 and False with 0
                    self.molprops_df1 = self.molprops_df1.replace({True: 1, False: 0})

                    c13_df_list.append(self.molprops_df1.values.tolist())

                tab_headings = ["CH0", "CH1", "CH2", "CH3"]

                print("c13_df_list\n", c13_df_list)

                rtn_html = render_template(
                    "error_table.html",
                    colHeadings=["ppm", "numProtons", "CH0", "CH1", "CH2", "CH3"],
                    df_data=c13_df_list,
                    tab_headings=tab_headings,
                )
                return rtn_html, 400

        return "ok", 200

    def check_hsqc_assign_CH2_from_DEPT135(self):
        """
        Assigns CH2 group intensities in the HSQC dataframe using DEPT135 experiment data.

        This function checks for the presence of CH2 carbons and updates the HSQC dataframe
        by setting the intensity of CH2 signals to negative values based on DEPT135 results.
        It handles both diastereotopic and non-diastereotopic CH2 cases and sets error messages
        if the data is inconsistent or missing.
        """
        mol = self.expected_molecule

        print("check_hsqc_assign_CH2_from_DEPT135")
        print("self.hsqc_df.shape[0]", self.hsqc_df.shape[0])
        print("self.dept_df.shape[0]", self.dept_df.shape[0])
        print("mol.num_carbon_atoms_with_protons", mol.num_carbon_atoms_with_protons)

        print(
            "self.hsqc_df.shape[0] == self.dept_df.shape[0]",
            self.hsqc_df.shape[0] == self.dept_df.shape[0],
        )

        # check if mol has CH2 in dataframe then does hsqc have some negative intensities
        if (
            (mol.num_CH2_carbon_atoms > 0)
            and (self.hsqc_df[self.hsqc_df.intensity < 0].shape[0] == 0)
            and
            # the number of entries in the hsqc_df with intensity > 0 is greater than the number of CH3CH1 entries in the mol
            (
                self.hsqc_df[self.hsqc_df.intensity > 0].shape[0]
                > (mol.num_CH1_carbon_atoms + mol.num_CH3_carbon_atoms)
            )
        ):

            if len(self.dept_df) == 0:
                self.nmrsolution_failed = True
                self.nmrsolution_error_message = "<p> hsqc needs to be multiplicity edited experiment as no DEPT data found</p>"
                self.nmrsolution_error_code = 401
                return

            # set intensity of CH2 values in hsqc to negative by using dept values
            dept_CH2 = self.dept_df[self.dept_df.intensity < 0].copy()
            dept_CH3CH1 = self.dept_df[self.dept_df.intensity > 0].copy()

            if len(dept_CH2) == 0:
                self.nmrsolution_failed = True
                self.nmrsolution_error_message = "<p> hsqc needs to be multiplicity edited experiment as no DEPT CH2 data found</p>"
                self.nmrsolution_error_code = 401
                return

            if len(dept_CH2) != mol.num_CH2_carbon_atoms:
                self.nmrsolution_failed = True
                self.nmrsolution_error_message = "<p> hsqc needs to be multiplicity edited experiment as DEPT CH2 data does not match expected number of CH2 carbon atoms</p>"
                self.nmrsolution_error_code = 401
                return

            if (len(self.dept_df) != mol.num_carbon_atoms_with_protons) and (len(self.dept_df) != mol.num_sym_carbon_atoms_with_protons):
                self.nmrsolution_failed = True
                self.nmrsolution_error_message = "<p> Please check the number of carbon peaks in the DEPT data as the number does not equal to the number of protonated carbons in the given molecule</p>"
                self.nmrsolution_error_code = 401
                return

            # sort dept_CH2 by ppm in descending order
            dept_CH2 = dept_CH2.sort_values(by=["ppm"], ascending=False).reset_index(
                drop=True
            )

            # sort hsqc by f1_ppm in descending order
            self.hsqc_df = self.hsqc_df.sort_values(
                by=["f1_ppm"], ascending=False
            ).reset_index(drop=True)

            hsqc_ppm_values = self.hsqc_df.f1_ppm.tolist()

            # check if the number of hsqc values is equal to the number of carbons in molprops_df with protons attached
            if len(hsqc_ppm_values) > mol.num_carbon_atoms_with_protons:
                # we have  diasterotopic CH2 carbons ie one proton ppm value for each carbon
                pass

            elif len(hsqc_ppm_values) == len(
                mol.sym_molprops_df[~mol.sym_molprops_df.CH0]
            ):
                # no diasterotopic CH2 carbons
                pass

            # check if the number of hsqc rows equals the number of dept rows
            # if yes no distasterotopic CH2 carbons, therefore each CH2 carbon has only one proton group

            if self.hsqc_df.shape[0] == self.dept_df.shape[0]:
                print("no diasterotopic CH2 carbons")
                # loop through dept_CH2, find closest ppm value in hsqc and set the intensity to negative
                # tol = 0.0025  # CARBONSEPARATION
                tol = self.problemdata_json.carbonSeparation
                for v1 in dept_CH2["ppm"]:

                    hsqc = self.hsqc_df[self.hsqc_df.intensity > 0]

                    # find the index of closest ppm value in hsqc
                    idx = hsqc["f1_ppm"].sub(v1).abs().idxmin()
                    # set the intensity of the idx loc in self.hsqc_df to -1
                    self.hsqc_df.loc[idx, "intensity"] = -1.0

            elif self.hsqc_df.shape[0] > self.dept_df.shape[0]:
                print("diasterotopic CH2 carbons")
                # loop through dept_CH2, find closest ppm value in hsqc and set the intensity to negative
                # tol = 0.0025  # CARBONSEPARATION
                tol = self.problemdata_json.carbonSeparation
                for v1 in dept_CH2["ppm"]:

                    self.hsqc_df["prob"] = self.hsqc_df.apply(
                        lambda x: stats.norm(v1, tol).pdf(x["f1_ppm"]), axis=1
                    )
                    # set the intensity to -1 for all prob values greater than 0.0
                    self.hsqc_df.loc[self.hsqc_df["prob"] > 0.0, "intensity"] = -1.0
                    self.hsqc_df["prob"] = 0.0

    def find_and_group_CH2s(self, df1):
        """
        Groups CH2 resonances in the provided DataFrame based on their chemical shift proximity.

        This function identifies CH2 signals and clusters them into groups where each group contains
        resonances close to each other in ppm, using a probability distribution threshold.
        It returns the indices and values of these grouped CH2 resonances.

        Args:
            df1 (pd.DataFrame): DataFrame containing NMR resonance data
                                with a 'CH2' boolean column and 'f1_ppm' values.

        Returns:
            tuple: (unique_idxs, unique_ch2s) where unique_idxs is a list of lists of indices for each group,
                                              and unique_ch2s is a list of lists of grouped ppm values.
        """
        # get a list of the CH2 resonances
        ch2_vals = df1[df1.CH2].f1_ppm.tolist()
        ch2_idx_vals = df1[df1.CH2].index.tolist()

        unique_ch2s = []
        unique_idxs = []

        # start with first CH2 value in the  list
        # create a probability distribution around it and obtain the probability of all the other values to
        # to see if they are close to the first value.
        # all values that have a +ve probability are close to the first value
        # all values with a zero probability are not.
        # add the +ve to a saved list of lists "similar_CH2s"
        # then remove them from the original list and repeat until original list length is zero

        while len(ch2_vals):
            # choose first from the list
            p0 = ch2_vals[0]
            print("p0", p0, type(p0))
            # find list of hmbc values that are similar to the first in the list
            similar_ch2s = [
                p
                for p in ch2_vals
                if stats.norm.pdf(
                    p, loc=p0, scale=self.problemdata_json.carbonSeparation
                )
                > 0
            ]
            similar_idxs = [
                i
                for i, p in zip(ch2_idx_vals, ch2_vals)
                if stats.norm.pdf(
                    p, loc=p0, scale=self.problemdata_json.carbonSeparation
                )
                > 0
            ]
            # if the length of the list is > 2 then we need to keep only the two closest values
            if len(similar_ch2s) > 2:
                # find the two closest values to the first value
                closest = np.argsort(np.abs(np.array(similar_ch2s) - p0))[:2]
                similar_ch2s = [similar_ch2s[i] for i in closest]
                similar_idxs = [similar_idxs[i] for i in closest]
            # save the list
            unique_ch2s.append(similar_ch2s)
            unique_idxs.append(similar_idxs)

            # keep only hmbc values that were not similar
            ch2_idx_vals = [
                i
                for i, p in zip(ch2_idx_vals, ch2_vals)
                if stats.norm.pdf(
                    p, loc=p0, scale=self.problemdata_json.carbonSeparation
                )
                == 0
            ]
            ch2_vals = [
                p
                for p in ch2_vals
                if stats.norm.pdf(
                    p, loc=p0, scale=self.problemdata_json.carbonSeparation
                )
                == 0
            ]

        return unique_idxs, unique_ch2s

    def create_svg_string(
        self, mol, molWidth=1000, molHeight=600, svgWidth=1200, svgHeight=700
    ):
        """
        Generates an SVG string representation of a molecule and computes normalized coordinates for carbon atoms.

        This function draws the given RDKit molecule as an SVG, centers it within the specified dimensions,
        and returns both the SVG string and a dictionary of normalized (x, y) coordinates for each carbon atom.

        Args:
            mol: RDKit molecule object to be drawn.
            molWidth (int, optional): Width of the molecule drawing area in pixels. Defaults to 1000.
            molHeight (int, optional): Height of the molecule drawing area in pixels. Defaults to 600.
            svgWidth (int, optional): Total SVG width in pixels. Defaults to 1200.
            svgHeight (int, optional): Total SVG height in pixels. Defaults to 700.

        Returns:
            tuple: (svg_string, new_xy3) where svg_string is the SVG representation of the molecule,
                   and new_xy3 is a dict mapping atom indices to normalized (x, y) coordinates.
        """

        translateWidth = int((svgWidth - molWidth) / 2)
        translateHeight = int((svgHeight - molHeight) / 2)

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

    def create_expected_molecule_from_molfilestring(
        self,
        molfilestr: str,
        carbonAtomsInfo: pd.DataFrame,
        df: pd.DataFrame,
        prediction_from_nmrshiftdb: bool,
        JEOL_predict: bool
    ) -> expectedmolecule.expectedMolecule:

        """
        Creates an expectedMolecule object from a molfile string and associated carbon atom information.

        This function initializes an expectedMolecule using the provided molfile string,
        carbon atom data, and optional NMR prediction data.

        Args:
            molfilestr (str): The molfile string representing the molecule.
            carbonAtomsInfo (pd.DataFrame): DataFrame containing information about carbon atoms.
            df (pd.DataFrame): DataFrame with C13 prediction data.
            prediction_from_nmrshiftdb (bool): Whether to use NMRShiftDB predictions.

        Returns:
            expectedmolecule.expectedMolecule: The constructed expectedMolecule object.
        """

        mol = expectedmolecule.expectedMolecule(
            molfilestr,
            is_smiles=False,
            carbonAtomsInfo=carbonAtomsInfo,
            mnova_c13predictions=df,
            predict_from_nmrshiftdb=prediction_from_nmrshiftdb,
            JEOL_predict=JEOL_predict
        )

        return mol

    def init_h1_and_pureshift_from_c13_hsqc_hmbc_cosy(
        self
    ) ->  Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Initializes the H1 and pureshift DataFrames from C13, HSQC, HMBC, and COSY experiment data.

        This function creates the H1 and pureshift DataFrames using information from the HSQC DataFrame,
        ensuring the correct columns and ordering for downstream NMR analysis. If the data is already present,
        it returns the existing DataFrames.

        Returns:
            tuple: (h1_df, pureshift_df) where both are pandas DataFrames containing H1 and pureshift data.
        """
        if self.PURESHIFT_data_present and self.H1_data_present:
            return self.h1_df, self.pureshift_df

        if (not self.h1_df.empty) and (not self.pureshift_df.empty):
            return self.h1_df, self.pureshift_df

        self.h1_df = pd.DataFrame(
            columns=[
                "ppm",
                "numProtons",
                "integral",
                "jCouplingVals",
                "jCouplingClass",
                "intensity",
                "range",
                "signaltype",
            ]
        )
        self.h1_df["ppm"] = self.hsqc_df["f2_ppm"].copy()
        self.h1_df["signaltype"] = self.hsqc_df["signaltype"].copy()

        # keep only rows where type is 'Compound' or  "" (empty)
        self.h1_df = self.h1_df[
            (self.h1_df["signaltype"] == "Compound") | (self.h1_df["signaltype"] == 0)
        ]

        self.h1_df["numProtons"] = -1
        self.h1_df["integral"] = -1
        self.h1_df["jCouplingVals"] = ""
        self.h1_df["jCouplingClass"] = ""
        self.h1_df["intensity"] = self.hsqc_df["intensity"].copy()
        self.h1_df["range"] = 0.0

        # order h1_df by ppm in place highest to lowest, reset index to 1,2,3,4,5...
        self.h1_df = self.h1_df.sort_values("ppm", ascending=False, ignore_index=True)
        self.h1_df.index = self.h1_df.index + 1

        # create pureshift_df from h1_df
        self.pureshift_df = self.h1_df.copy()

        return self.h1_df, self.pureshift_df

    def init_pureshift_from_h1(self) -> pd.DataFrame:
        """
        Initializes the pureshift DataFrame from the H1 DataFrame if not already present.

        This function copies the H1 DataFrame to the pureshift DataFrame if pureshift data is missing,
        and returns the resulting pureshift DataFrame.

        Returns:
            pd.DataFrame: The pureshift DataFrame.
        """

        print("init_pureshift_from_h1")

        if self.PURESHIFT_data_present:
            return self.pureshift_df

        if not self.h1_df.empty:
            self.pureshift_df = self.h1_df.copy()
        return self.pureshift_df

    def init_h1_from_pureshift(self) -> pd.DataFrame:
        """
        Initializes the H1 DataFrame from the pureshift DataFrame if H1 data is missing.

        This function creates the H1 DataFrame using information from the pureshift DataFrame,
        ensuring the correct columns and ordering for downstream NMR analysis. If H1 data is already present,
        it returns the existing H1 DataFrame.

        Returns:
            pd.DataFrame: The H1 DataFrame.
        """

        if self.H1_data_present:
            return self.h1_df
        if self.PURESHIFT_data_missing:
            return self.h1_df

        h1_df = pd.DataFrame(
            columns=[
                "ppm",
                "numProtons",
                "integral",
                "jCouplingVals",
                "jCouplingClass",
                "intensity",
                "range",
                "signaltype",
            ]
        )
        h1_df["ppm"] = self.pureshift_df["ppm"].copy()
        h1_df["signaltype"] = self.pureshift_df["signaltype"].copy()

        # keep only rows where type is 'Compound' or  "" (empty)

        h1_df = h1_df[(h1_df["signaltype"] == "Compound") | (h1_df["signaltype"] == 0)]

        h1_df["numProtons"] = -1
        h1_df["integral"] = -1
        h1_df["jCouplingVals"] = ""
        h1_df["jCouplingClass"] = ""
        h1_df["intensity"] = 1.0
        h1_df["range"] = 0.0

        # order h1_df by ppm in place highest to lowest, reset index to 1,2,3,4,5...
        h1_df = h1_df.sort_values("ppm", ascending=False, ignore_index=True)
        h1_df.index = h1_df.index + 1

        return h1_df

    def add_CH0_CH1_CH2_CH3_CH3CH1_to_C13_df(self) -> None:
        """
        Labels C13 DataFrame rows with CH0, CH1, CH2, CH3, and CH3CH1 group information using HSQC data.

        This function updates the C13 DataFrame to indicate the type of each carbon atom (e.g., CH0, CH1, CH2, CH3, CH3CH1, quaternary)
        based on matching chemical shifts with the HSQC DataFrame, and sets the number of protons accordingly.
        """

        # label c13 values as CH2 or not using the hsqc_df_
        print("add_CH0_CH1_CH2_CH3_CH3CH1_to_C13_df")
        self.c13_df["CH0"] = False
        self.c13_df["CH1"] = False
        self.c13_df["CH2"] = False
        self.c13_df["CH3"] = False
        self.c13_df["CH3CH1"] = False
        self.c13_df["quaternary"] = False
        self.c13["numProtons"] = -1

        for idx, hsqc_row in self.hsqc.iterrows():
            hsqc_ppm = hsqc_row["f1_ppm"]
            c13_rows = self.c13_df[self.c13_df.ppm == hsqc_ppm]

            if hsqc_ppm in self.c13_df.ppm.values:
                for column in ["numProtons", "CH1", "CH2", "CH3", "CH3CH1"]:
                    self.c13_df.loc[c13_rows.index, column] = hsqc_row[column].values[0]

        # label all rows in c13_df as CH0 where numProtons is -1
        self.c13_df.loc[self.c13_df.numProtons == -1, "CH0"] = True
        self.c13_df.loc[self.c13_df.numProtons == -1, "numProtons"] = 0

        self.c13

    def add_CH0_CH1_CH2_CH3_CH3CH1_to_C13(self) -> None:
        """
        Labels C13 DataFrame rows with CH0, CH1, CH2, CH3, and CH3CH1 group information using HSQC data.

        This function updates the C13 DataFrame to indicate the type of each carbon atom (e.g., CH0, CH1, CH2, CH3, CH3CH1, quaternary)
        based on matching chemical shifts with the HSQC DataFrame, and sets the number of protons accordingly.
        """

        # label c13 values as CH2 or not using the hsqc_df_
        self.c13["CH0"] = False
        self.c13["CH1"] = False
        self.c13["CH2"] = False
        self.c13["CH3"] = False
        self.c13["CH3CH1"] = False
        self.c13["quaternary"] = False
        self.c13["numProtons"] = -1

        for idx, ppm in zip(self.c13.index, self.c13.ppm):
            if ppm in self.hsqc.f1_ppm.values:
                self.c13.loc[idx, "CH2"] = self.hsqc_df[self.hsqc.f1_ppm == ppm][
                    "CH2"
                ].values[0]
                self.c13.loc[idx, "numProtons"] = 2

        # label c13 values as quaternary if not present in hsqc_df
        for idx, ppm in zip(self.c13.index, self.c13.ppm):
            if ppm not in self.hsqc.f1_ppm.values:
                self.c13.loc[idx, "quaternary"] = True
                self.c13.loc[idx, "CH0"] = True
                self.c13.loc[idx, "numProtons"] = 0

        # label c13 values as CH3CH if quaternary and CH2 both False
        for idx, ppm in zip(self.c13.index, self.c13.ppm):
            if not self.c13.loc[idx, "quaternary"] and not self.c13.loc[idx, "CH2"]:
                self.c13.loc[idx, "CH3CH1"] = True
                # self.c13.loc[idx, "numProtons"] = 1

    def init_CH_CH3_HSQC_from_DDEPT_CH3_only(self) -> None:
        """
        Updates HSQC DataFrames with CH3 and CH1 group assignments using DDEPT CH3-only experiment data.

        This function sets the CH3 and CH1 columns in the HSQC DataFrames based on the presence of CH3 signals
        in the DDEPT CH3-only experiment, ensuring correct group assignments for downstream analysis.
        """
        if self.ddept_ch3_only.empty:
            return

        if "CH1" not in self.hsqc_df.columns:
            self.hsqc_df["CH1"] = False
        if "CH3" not in self.hsqc_df.columns:
            self.hsqc_df["CH3"] = False
        if "CH1" not in self.hsqc.columns:
            self.hsqc["CH1"] = False
        if "CH3" not in self.hsqc.columns:
            self.hsqc["CH3"] = False

        for idx, ppm in zip(self.ddept_ch3_only.index, self.ddept_ch3_only.f1_ppm):
            if ppm in self.hsqc.f1_ppm.values:
                self.hsqc.loc[self.hsqc.f1_ppm == ppm, "CH3"] = True
                self.hsqc.loc[self.hsqc.f1_ppm == ppm, "CH1"] = False

        # update hsqc_df with CH3 and CH1 values from hsqc
        self.hsqc_df["CH3"] = self.hsqc["CH3"]
        self.hsqc_df["CH1"] = self.hsqc["CH1"]

    def init_c13_from_hsqc_and_hmbc(self) -> Tuple[pd.DataFrame, bool]:
        """
        Initializes the C13 DataFrame using HSQC and HMBC experiment data.

        This function generates a C13 DataFrame by combining and processing chemical shift values from the HSQC and HMBC DataFrames,
        assigning group labels (CH0, CH1, CH2, CH3, CH3CH1, quaternary) and the number of protons for each carbon resonance.
        If C13 data is already present, it returns the existing DataFrame and False.

        Returns:
            tuple: (c13_df, created) where c13_df is the resulting C13 DataFrame and created is a boolean indicating if new data was created.
        """
        print("init_c13_from_hsqc_and_hmbc")

        if self.C13_data_present:
            return self.c13_df, False

        # replace values for f2_ppm in hmbc starting from f2_ppm hsqc
        if not self.exact_ppm_values:
            self.hmbc_df = self.tidyup_ppm_values(
                self.hmbc_df,
                sorted(self.hsqc.f2_ppm.unique().tolist(), reverse=True),
                "f2_ppm",
                ppm_tolerance=self.problemdata_json.protonSeparation,
            )
        else:
            print("Exact ppm values only, no tidy up required")

        # if any hmbc f2_ppm_probs == 0 then drop the row
        self.hmbc_df.drop(
            self.hmbc_df[self.hmbc_df.f2_ppm_prob == 0].index, inplace=True
        )

        # find all f1_ppm HMBC idx resonances that are not showing up in the HSQC f2_ppm
        iii = []
        for i in self.hmbc_df.index:
            prob_vals = []
            for c in self.hsqc.f1_ppm.unique():
                prob_vals.append(
                    stats.norm.pdf(
                        self.hmbc_df.loc[i, "f1_ppm"],
                        loc=c,
                        scale=self.problemdata_json.carbonSeparation,
                    )
                )
            if np.array(prob_vals).sum() == 0:
                iii.append(i)
            else:
                pass

        # keep only the unique hmbc resonances not in f1_ppm HSQC
        # get a list of the hmbc resonances
        hmbcs = self.hmbc_df.loc[iii, "f1_ppm"].tolist()
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
                p
                for p in hmbcs
                if stats.norm.pdf(
                    p, loc=p0, scale=self.problemdata_json.carbonSeparation
                )
                > 0
            ]
            # save the list
            unique_hmbc.append(similar_hmbcs)

            # keep only hmbc values that were not similar
            hmbcs = [
                p
                for p in hmbcs
                if stats.norm.pdf(
                    p, loc=p0, scale=self.problemdata_json.carbonSeparation
                )
                == 0
            ]

        # create an array of mean values for hmbc f1_ppm not found in hsqc f1_ppm
        mean_unique_hmbc_vals = [np.mean(h) for h in unique_hmbc]

        # tidyup f1_ppm values in hmbc that are not in f1_ppm HSQC
        if not self.exact_ppm_values:

            hmbc_1 = self.tidyup_ppm_values(
                self.hmbc_df.loc[iii],
                mean_unique_hmbc_vals,
                "f1_ppm",
                ppm_tolerance=self.problemdata_json.carbonSeparation,
            )

            # tidyup f1_ppm values in hmbc that are in f1_ppm HSQC
            hmbc_2 = self.tidyup_ppm_values(
                self.hmbc_df.drop(iii),
                self.hsqc.f1_ppm.unique(),
                "f1_ppm",
                ppm_tolerance=self.problemdata_json.carbonSeparation,
            )
        else:
            print("Exact ppm values only, no tidy up required")

        # rejoin two parts of HMBC data
        self.hmbc_df = pd.concat([hmbc_1, hmbc_2])
        self.hmbc_df.sort_index(inplace=True)

        # add f2p_ppm column to HSQC and HMBC tables
        # f2p_ppm is C13 one to one relation between f1_ppm and f2_ppm in HSQC
        self.hmbc_df["f2p_ppm"] = 0.0
        for idx, f2ppmHSQC, f1ppmHSQC in zip(
            self.hsqc.index, self.hsqc.f2_ppm, self.hsqc.f1_ppm
        ):
            self.hmbc_df.loc[self.hmbc_df.f2_ppm == f2ppmHSQC, "f2p_ppm"] = f1ppmHSQC

        # return list of C13 values
        c13_list = sorted(
            set(self.hmbc_df.f1_ppm).union(
                set(self.hmbc_df.f2p_ppm), set(self.hsqc.f1_ppm)
            ),
            reverse=True,
        )

        # create a dataframe of C13 values
        self.c13_df = pd.DataFrame(c13_list, columns=["ppm"])
        self.c13_df["signaltype"] = "Compound"
        self.c13_df.index = self.c13_df.index + 1

        # label c13 values as CH2 or not using the hsqc
        self.c13_df["CH2"] = False
        if "numProtons" not in self.c13_df.columns:
            self.c13_df["numProtons"] = -1

        CH2_hsqc_df = self.hsqc[self.hsqc.CH2]
        for idx, ppm in zip(self.c13_df.index, self.c13_df.ppm):
            if ppm in CH2_hsqc_df.f1_ppm.values:
                self.c13_df.loc[idx, "CH2"] = True

        #  label c13 values as quaternary if not present in hsqc
        self.c13_df["quaternary"] = False
        self.c13_df["CH0"] = False
        for idx, ppm in zip(self.c13_df.index, self.c13_df.ppm):
            if ppm not in self.hsqc.f1_ppm.values:
                self.c13_df.loc[idx, "quaternary"] = True
                self.c13_df.loc[idx, "CH0"] = True

        # label c13 values as CH3CH if quaternary and CH2 both False
        self.c13_df["CH3CH1"] = False
        CH3CH1_hsqc_df = self.hsqc[self.hsqc.CH3CH1]
        for idx, ppm in zip(self.c13_df.index, self.c13_df.ppm):
            if ppm in CH3CH1_hsqc_df.f1_ppm.values:
                self.c13_df.loc[idx, "CH3CH1"] = True

        # label c13 values as CH3 if any CH3 values in hsqc
        self.c13_df["CH3"] = False
        self.c13_df["CH1"] = False
        if "numProtons" not in self.c13_df.columns:
            self.c13_df["numProtons"] = -1
        CH3_hsqc_df = self.hsqc[self.hsqc.CH3]
        if not CH3_hsqc_df.empty:
            for idx, ppm in zip(self.c13_df.index, self.c13_df.ppm):
                if ppm in CH3_hsqc_df.f1_ppm.values:
                    self.c13_df.loc[idx, "CH3"] = True

        # label c13 values as CH1 if any CH1 values in hsqc
        CH1_hsqc_df = self.hsqc[self.hsqc.CH1]
        if not CH1_hsqc_df.empty:
            for idx, ppm in zip(self.c13_df.index, self.c13_df.ppm):
                if ppm in CH1_hsqc_df.f1_ppm.values:
                    self.c13_df.loc[idx, "CH1"] = True

        self.c13_df.loc[self.c13_df.CH3, "numProtons"] = 3
        self.c13_df.loc[self.c13_df.CH2, "numProtons"] = 2
        self.c13_df.loc[self.c13_df.CH1, "numProtons"] = 1
        self.c13_df.loc[self.c13_df.CH0, "numProtons"] = 0

        return self.c13_df, True

    def check_what_experiments_are_present(self) -> None:
        """
        Checks which NMR experiment data types are present in the available experiments.

        This function sets boolean flags indicating the presence or absence of each expected NMR experiment
        (e.g., HSQC, H1_pureshift, C13_1D, H1_1D, COSY, HMBC, NOESY, HSQC_CLIPCOSY, DDEPT_CH3_ONLY, DEPT135)
        based on the contents of the expts_available attribute.
        """
        # do basic checks on excels sheets and decide how to proceed
        if {"HSQC"}.issubset(self.expts_available):
            self.HSQC_data_present = True
            self.HSQC_data_missing = False
        else:
            self.HSQC_data_present = False
            self.HSQC_data_missing = True
        if {"H1_pureshift"}.issubset(self.expts_available):
            self.PURESHIFT_data_present = True
            self.PURESHIFT_data_missing = False
        else:
            self.PURESHIFT_data_present = False
            self.PURESHIFT_data_missing = True

        if {"C13_1D"}.issubset(self.expts_available):
            self.C13_data_present = True
            self.C13_data_missing = False
        else:
            self.C13_data_present = False
            self.C13_data_missing = True

        if {"H1_1D"}.issubset(self.expts_available):
            self.H1_data_present = True
            self.H1_data_missing = False
        else:
            self.H1_data_present = False
            self.H1_data_missing = True

        if {"COSY"}.issubset(self.expts_available):
            self.COSY_data_present = True
            self.COSY_data_missing = False
        else:
            self.COSY_data_present = False
            self.COSY_data_missing = True

        if {"HMBC"}.issubset(self.expts_available):
            self.HMBC_data_present = True
            self.HMBC_data_missing = False
        else:
            self.HMBC_data_present = False
            self.HMBC_data_missing = True

        if {"NOESY"}.issubset(self.expts_available):
            self.NOESY_data_present = True
            self.NOESY_data_missing = False
        else:
            self.NOESY_data_present = False
            self.NOESY_data_missing = True

        if {"HSQC_CLIPCOSY"}.issubset(self.expts_available):
            self.HSQC_CLIPCOSY_data_present = True
            self.HSQC_CLIPCOSY_data_missing = False
        else:
            self.HSQC_CLIPCOSY_data_present = False
            self.HSQC_CLIPCOSY_data_missing = True

        if {"DDEPT_CH3_ONLY"}.issubset(self.expts_available):
            self.DDEPT_CH3_ONLY_data_present = True
            self.DDEPT_CH3_ONLY_data_missing = False
        else:
            self.DDEPT_CH3_ONLY_data_present = False
            self.DDEPT_CH3_ONLY_data_missing = True

        if {"DEPT135"}.issubset(self.expts_available):
            self.DEPT_data_present = True
            self.DEPT_data_missing = False
        else:
            self.DEPT_data_present = False
            self.DEPT_data_missing = True

    def check_h1_pureshift_not_empty(self) -> Tuple[bool, str]:
        """
        Checks if the H1 or pureshift DataFrames are not empty and initializes the H1 DataFrame accordingly.

        This function verifies the presence of H1 or pureshift data, attempts to extract valid H1 signals,
        and returns a tuple indicating success and a status message.

        Returns:
            tuple: (success, message) where success is a boolean indicating if H1 data is present,
                   and message is a string describing the result.
        """
        if not self.pureshift_df.empty:
            self.h1 = self.pureshift_df[
                (self.pureshift_df.signaltype == "Compound")
                | (self.pureshift_df.signaltype == 0)
            ][["ppm"]].copy()
            if self.h1.empty:
                # maybe user forgot to set type to compound so search for "" in Type column
                self.h1 = self.pureshift_df[self.pureshift_df.signaltype == ""][
                    ["ppm"]
                ].copy()
                self.h1["signaltype"] = "Compound"
            if self.h1.empty:
                warning_dialog(
                    "H1 Excel Sheet is Empty", "Error in Excel Sheet", self.qtstarted
                )
                return False, "h1 is empty"
        elif not self.h1_df.empty:
            self.h1 = self.h1_df[["ppm"]].copy()
        else:
            warning_dialog(
                "H1_df or pureshift_df excel sheet is empty",
                "Error in Excel Sheet",
                self.qtstarted,
            )
            return False, "No pureshift_df or h1_df"

        return True, "h1 is not empty"

    def check_c13_not_empty(self) -> Tuple[bool, str]:
        """
        Checks if the C13 DataFrame is not empty and initializes the C13 DataFrame accordingly.

        This function verifies the presence of C13 data, attempts to extract valid C13 signals,
        and returns a tuple indicating success and a status message.

        Returns:
            tuple: (success, message) where success is a boolean indicating if C13 data is present,
                   and message is a string describing the result.
        """
        if self.c13_df.empty:
            warning_dialog(
                "C13 Excel Sheet is Empty", "Error in Excel File", self.qtstarted
            )
            return False, "c13 is empty"
        else:
            self.c13 = self.c13_df[
                (self.c13_df.signaltype == "Compound") | (self.c13_df.signaltype == 0)
            ][["ppm"]].copy()
            if self.c13.empty:
                # maybe user forgot to set type to compound so search for "" in Type column
                self.c13 = self.c13_df[self.c13_df.signaltype == ""][["ppm"]].copy()
                self.c13["signaltype"] = "Compound"
            if self.c13.empty:
                warning_dialog(
                    "C13 Excel Sheet is Empty", "Error in Excel File", self.qtstarted
                )
                return False, "c13 is empty"

        if "numProtons" not in self.c13.columns:
            self.c13["numProtons"] = -1
        if "attached_protons" not in self.c13.columns:
            self.c13["attached_protons"] = -1
        if "CH0" not in self.c13.columns:
            self.c13["CH0"] = False
        if "CH1" not in self.c13.columns:
            self.c13["CH1"] = False
        if "CH2" not in self.c13.columns:
            self.c13["CH2"] = False
        if "CH3" not in self.c13.columns:
            self.c13["CH3"] = False
        if "CH3CH1" not in self.c13.columns:
            self.c13["CH3CH1"] = False
        if "quaternary" not in self.c13.columns:
            self.c13["quaternary"] = False

        for CHn in ["CH0", "CH1", "CH2", "CH3", "CH3CH1", "quaternary", "numProtons"]:
            if CHn in self.c13_df.columns:
                self.c13[CHn] = self.c13_df[
                    (self.c13_df.signaltype == "Compound")
                    | (self.c13_df.signaltype == 0)
                ][CHn].copy()

        # add two to numProtons for CH2
        self.c13.loc[self.c13.CH2, "numProtons"] = 2
        # add one to numProtons for CH3CH
        self.c13.loc[self.c13.CH1, "numProtons"] = 1
        # add zero to numProtons for quaternary
        self.c13.loc[self.c13.quaternary, "numProtons"] = 0
        self.c13.loc[self.c13.CH0, "numProtons"] = 0
        # add three to numProtons for CH3
        self.c13.loc[self.c13.CH3, "numProtons"] = 3

        self.c13.loc[:, "attached_protons"] = self.c13.numProtons

        return True, "c13 is not empty"

    def check_hsqc_df_not_empty(self) -> Tuple[bool, str]:
        """
        Checks if the HSQC DataFrame is not empty and initializes the HSQC DataFrame accordingly.

        This function verifies the presence of HSQC data, attempts to extract valid HSQC signals,
        and returns a tuple indicating success and a status message.

        Returns:
            tuple: (success, message) where success is a boolean indicating if HSQC data is present,
                   and message is a string describing the result.
        """

        print("self.hsqc_df:\n", self.hsqc_df)
        if not self.hsqc_df.empty:
            # copy rows where signaltype is Compound or 0
            self.hsqc = self.hsqc_df[
                (self.hsqc_df.signaltype == "Compound") | (self.hsqc_df.signaltype == 0)
            ][["f2_ppm", "f1_ppm", "intensity", "signaltype", "Annotation"]].copy()

        if self.hsqc.empty:
            # maybe user forgot to set type to compound so search for "" in Type column
            self.hsqc = self.hsqc_df[self.hsqc_df.signaltype == ""][
                ["f2_ppm", "f1_ppm", "intensity", "signaltype", "Annotation"]
            ].copy()
            self.hsqc["signaltype"] = "Compound"

        if self.hsqc.empty:
            warning_dialog(
                "HSQC excel Sheet is Empty", "Error in Excel File", self.qtstarted
            )
            return False, "hsqc is empty"
        else:
            self.hsqc["f1_i"] = -1
            self.hsqc["f2_i"] = -1
            self.hsqc["f2p_i"] = -1
            self.hsqc["f1C_i"] = -1
            self.hsqc["f2H_i"] = -1
            self.hsqc["f2Cp_i"] = -1
            self.hsqc["f2p_ppm"] = -1

            self.hsqc["CH2"] = False
            self.numtimes_HSQC_CH2_set_to_FALSE += 1
            # print(
            #     "line 1043, check_hsqc_df_not_empty, self.numtimes_HSQC_CH2_set_to_FALSE"
            # )
            self.hsqc["CH3CH1"] = False
            self.hsqc["CH3"] = False
            self.hsqc["CH1"] = False
            self.hsqc["CH0"] = False
            self.hsqc["quaternary"] = False

            self.hsqc["numProtons"] = -1
            self.hsqc["attached_protons"] = -1
            self.hsqc["f2_integral"] = -1
            return True, "hsqc is not empty"

    def tidyup_ppm_values(
        self, df: pd.DataFrame, true_values: list, column_name: str, ppm_tolerance=0.005
    ) -> pd.DataFrame:
        """
        Adjusts chemical shift values in a DataFrame to their nearest true values within a specified tolerance.

        This function creates new columns for the original values and their probability of matching the adjusted value,
        then replaces each value in the specified column with its nearest value from a provided list.

        Args:
            df (pd.DataFrame): The DataFrame containing chemical shift values to adjust.
            true_values (list): List of reference values to match against.
            column_name (str): Name of the column in df to adjust.
            ppm_tolerance (float, optional): Standard deviation for probability calculation. Defaults to 0.005.

        Returns:
            pd.DataFrame: The DataFrame with adjusted values and probability columns.
        """

        # make a copy of the column_name adding a suffix orig
        df[f"{column_name}_orig"] = df[column_name]

        # make a probability column to see how far replacement is from original
        df[f"{column_name}_prob"] = 0

        # create dataframe with ppm values and their nearest true value
        # dfnew = df.assign(nearest_true_value = df[column_name].apply(lambda x: find_nearest(true_values, x)))
        df[column_name] = df[column_name].apply(
            lambda x: self.find_nearest(true_values, x)
        )

        # calculate probabilities
        for idx in df.index:
            df.loc[idx, f"{column_name}_prob"] = stats.norm.pdf(
                df.loc[idx, column_name],
                loc=df.loc[idx, f"{column_name}_orig"],
                scale=ppm_tolerance,
            )

        return df

    def add_CH2_CH3CH_to_hsqc_dataframes(self) -> None:
        """
        Updates the HSQC DataFrames with the number of protons and related columns for CH3, CH2, and CH1 groups.

        This function sets the 'numProtons', 'integral', and 'attached_protons' columns in both the main and
        DataFrame versions of HSQC for each group type, ensuring correct assignment for downstream NMR analysis.
        """

        print("add_CH2_CH3CH_to_hsqc_dataframes")
        self.hsqc["numProtons"] = -1
        self.hsqc["integral"] = -1
        self.hsqc["attached_protons"] = -1

        self.hsqc_df["numProtons"] = -1
        self.hsqc_df["integral"] = -1
        self.hsqc_df["attached_protons"] = -1

        for CHn, nHs in zip(["CH3", "CH2", "CH1"], [3, 2, 1]):

            for col_id in ["numProtons", "integral", "attached_protons"]:
                self.hsqc.loc[self.hsqc[CHn], col_id] = nHs
                self.hsqc_df.loc[self.hsqc_df[CHn], col_id] = nHs

    def transfer_hsqc_info_to_h1(self) -> None:
        """
        Transfers group and proton information from the HSQC DataFrame to the H1 DataFrame.

        This function synchronizes CH3CH1, CH3, CH2, CH1, and numProtons columns from HSQC to H1,
        updates attached protons and integrals, and handles special cases for diastereotopic protons.
        """
        # transfer information from the HSQC to the H1
        # check if CH3CH1, CH3, CH2, and CH1 columns are in H1 if not add and set to False

        h1 = self.h1
        hsqc = self.hsqc

        for col in ["CH3CH1", "CH3", "CH2", "CH1"]:
            if col not in h1.columns:
                h1[col] = False

        # copy information from HSQC to H1
        for index, row in hsqc.iterrows():
            h1.loc[h1["ppm"] == row["f2_ppm"], "CH3CH1"] = row["CH3CH1"]
            h1.loc[h1["ppm"] == row["f2_ppm"], "CH3"] = row["CH3"]
            h1.loc[h1["ppm"] == row["f2_ppm"], "CH2"] = row["CH2"]
            h1.loc[h1["ppm"] == row["f2_ppm"], "CH1"] = row["CH1"]

        # transfer numProtons from hsqc to h1 if numProtons in h1 == -1
        if h1[h1.numProtons == -1].shape[0] > 0:
            for idx, row in h1[h1.numProtons == -1].iterrows():
                h1.loc[idx, "numProtons"] = hsqc[hsqc.f2_ppm == row.ppm][
                    "numProtons"
                ].values[0]

        h1["attached_protons"] = h1["numProtons"]

        # check if h1 integrals all negative then set to intensity
        if h1[h1.integral < 0].shape[0] > 0:
            h1["integral"] = h1["numProtons"]

        # set any integral that equals 2 to -2 in h1
        h1.loc[h1.integral == 2, "integral"] = -2

        # if we have diasteremeric protons then set integral to -1
        # this can be tested if we count the frequency of the values in column f1p_i in h1 and if an equal 2 set the integraal to -1
        unique_f1p_i = h1.f1p_i.unique()
        for f1p_i in unique_f1p_i:
            if h1[h1.f1p_i == f1p_i].shape[0] == 2:
                h1.loc[h1.f1p_i == f1p_i, "integral"] = -1

    def transfer_hsqc_info_to_c13(self) -> None:
        """
        Transfers group and proton information from the HSQC DataFrame to the C13 DataFrame.

        This function synchronizes numProtons, CH3, CH2, CH1, and CH3CH1 columns from HSQC to C13,
        and updates CH0 and quaternary flags for carbons with zero protons.

        Returns:
            None
        """
        # transfer information from hsqc to c13

        c13: pd.DataFrame = self.c13
        hsqc: pd.DataFrame = self.hsqc

        for idx, row in hsqc.iterrows():
            # search c13 by ppm
            c13.loc[c13.ppm == row.f1_ppm, "numProtons"] = row.numProtons
            c13.loc[c13.ppm == row.f1_ppm, "CH3"] = row.CH3
            c13.loc[c13.ppm == row.f1_ppm, "CH2"] = row.CH2
            c13.loc[c13.ppm == row.f1_ppm, "CH1"] = row.CH1
            c13.loc[c13.ppm == row.f1_ppm, "CH3CH1"] = row.CH3CH1

        # any leftover numprotons in c13 == -1 set to 0
        c13.loc[c13.numProtons == -1, "numProtons"] = 0
        c13.loc[c13.numProtons == 0, "CH0"] = True
        c13.loc[c13.numProtons == 0, "quaternary"] = True

        c13["attached_protons"] = c13["numProtons"]

    def assign_CH3_CH1_in_HSQC_using_expected_molecule(self) -> None:
        """
        Assigns CH3 and CH1 group information in the HSQC DataFrame using the expected molecule's properties.

        This function matches and labels CH3 and CH1 groups in the HSQC data based on the expected
        molecule's group counts and chemical shifts, handling various scenarios for group assignment
        and updating related columns accordingly.
        """
        print("c13_from_hsqc and H1_data_missing")

        expected_molecule = self.expected_molecule
        hsqc = self.hsqc

        CH3_mol_df = expected_molecule.molprops_df[expected_molecule.molprops_df.CH3]

        CH3CH1_mol_df = expected_molecule.molprops_df[
            expected_molecule.molprops_df.CH3CH1
        ]

        CH1_mol_df = expected_molecule.molprops_df[expected_molecule.molprops_df.CH1]

        print(f"CH1_mol_df.shape[0]: {CH1_mol_df.shape[0]}")

        CH3CH1_sym_mol_df = expected_molecule.sym_molprops_df[
            expected_molecule.sym_molprops_df.CH3CH1
        ]
        CH3_sym_mol_df = expected_molecule.sym_molprops_df[
            expected_molecule.sym_molprops_df.CH3
        ]
        CH1_sym_mol_df = expected_molecule.sym_molprops_df[
            expected_molecule.sym_molprops_df.CH1
        ]
        CH1_sym_mol_gt_67_df = CH1_sym_mol_df[CH1_sym_mol_df.ppm >= 67]
        CH1_sym_mol_lt_67_df = CH1_sym_mol_df[CH1_sym_mol_df.ppm < 67]

        CH3CH1_hsqc_df = hsqc[hsqc.CH3CH1].copy()

        num_CH3CH1_hsqc = CH3CH1_hsqc_df.shape[0]
        print(f"num_CH3CH1_hsqc: {num_CH3CH1_hsqc}")

        print(f"CH3CH1_hsqc_df.shape[0]: {CH3CH1_hsqc_df.shape[0]}")
        print(f"CH3CH1_mol_df.shape[0]: {CH3CH1_mol_df.shape[0]}")
        print(f"CH3CH1_sym_mol_df.shape[0]: {CH3CH1_sym_mol_df.shape[0]}")

        if num_CH3CH1_hsqc == 0:
            #  skip if no CH3CH1 groups in HSQC
            return
        
        if num_CH3CH1_hsqc > CH3CH1_mol_df.shape[0]:
            #  if more CH3CH1 groups in HSQC than expected molecule then we have a problem
            print(
                "more CH3CH1 groups in HSQC than expected molecule, this is a problem"
            )
            self.nmrsolution_failed = True
            self.nmrsolution_error_message = "CH3CH1 HSQC > CH3CH1 expected molecule"
            self.nmrsolution_error_code = 401
            return

        #  check if CH3 groups already assigned and if so assign CH1 groups
        if (CH3_mol_df.shape[0] > 0) and hsqc.CH3.sum() > 0:
            print("CH3 groups already assigned")
            #  set the CH1 to True where hsqc.numprotons < 0
            hsqc.loc[hsqc.numProtons < 0, "CH1"] = True
            hsqc.loc[hsqc.numProtons < 0, "numProtons"] = 1
            return

        # if no CH3 groups in expected molecule assign all CH3CH1 groups to CH1
        if CH3_mol_df.shape[0] == 0:

            hsqc.loc[CH3CH1_hsqc_df.index, "CH1"] = True
            hsqc.loc[CH3CH1_hsqc_df.index, "numProtons"] = 1
            print("no CH3 groups in expected molecule setting all CH3CH1 groups to CH1")
            return

        # if no CH1 groups in expected molecule assign all CH3CH1 groups to CH3
        elif CH1_mol_df.shape[0] == 0:
            print("no CH1 groups in expected molecule")
            hsqc.loc[CH3CH1_hsqc_df.index, "CH3"] = True
            hsqc.loc[CH3CH1_hsqc_df.index, "numProtons"] = 3
            print("no CH1 groups in expected molecule setting all CH3CH1 groups to CH3")
            return
        
        #  if not all CH3CH1 peaks are picked

        #  check if we think we have symmetry in the expected molecule
        if CH3CH1_sym_mol_df.shape[0] < CH3CH1_mol_df.shape[0] and CH3CH1_hsqc_df.shape[0] < CH3CH1_sym_mol_df.shape[0]:

            print("not all HSQC peaks picked yet")
            # loop through the CH3CH1_hsqc_df and try to match with CH3CH1_mol_df 
            CH3CH1_sym_mol_df = CH3CH1_sym_mol_df[~CH3CH1_sym_mol_df.picked]
            for idx, row in CH3CH1_hsqc_df.iterrows():
                # find the closest match in the CH3CH1_mol_df
                closest_match = CH3CH1_sym_mol_df.iloc[
                    (CH3CH1_sym_mol_df["ppm"] - row.f1_ppm).abs().argsort()[:1]
                ]
                print(f"closest_match: \n{closest_match}")
                print(f"closest_matc.columns: {closest_match.columns} ")
                print(f"closest_match.totalNumHs: {closest_match.totalNumHs.values[0]}")
                if not closest_match.empty:
                    # update the hsqc dataframe with the closest match
                    hsqc.loc[row.name, "numProtons"] = closest_match.totalNumHs.values[0]
                    hsqc.loc[row.name, "CH3"] = closest_match.CH3.values[0]
                    hsqc.loc[row.name, "CH1"] = closest_match.CH1.values[0]

                    CH3CH1_sym_mol_df.loc[closest_match.index, "picked"] = True
                    CH3CH1_sym_mol_df = CH3CH1_sym_mol_df[~CH3CH1_sym_mol_df.picked]
            
        # elif CH3CH1_hsqc_df.shape[0] < CH3CH1_mol_df.shape[0]:

        #     # loop through the CH3CH1_hsqc_df and try to match with CH3CH1_mol_df 
        #     CH3CH1_mol_df = CH3CH1_mol_df[~CH3CH1_mol_df.picked]
        #     for idx, row in CH3CH1_hsqc_df.iterrows():
        #         # find the closest match in the CH3CH1_mol_df
        #         closest_match = CH3CH1_mol_df.iloc[
        #             (CH3CH1_mol_df["ppm"] - row.f1_ppm).abs().argsort()[:1]
        #         ]
        #         if not closest_match.empty:
        #             # update the hsqc dataframe with the closest match
        #             hsqc.loc[row.name, "numProtons"] = closest_match.totalNumHs.values[0]
        #             hsqc.loc[row.name, "CH3"] = closest_match.CH3.values[0]
        #             hsqc.loc[row.name, "CH1"] = closest_match.CH1.values[0]
        #             # mark the closest match as picked in the CH3CH1_mol_df

        #             CH3CH1_mol_df.loc[closest_match.index, "picked"] = True
        #             CH3CH1_mol_df = CH3CH1_mol_df[~CH3CH1_mol_df.picked]

  


        else:

            print("we have both CH3 and CH1 groups in the expected molecule")
            # we have both CH3 and CH1 groups in the expected molecule

            # split the CH3CH1 based on ppm value above and below 67 ppm
            # CH3 groups expected to be below 67 ppm

            CH3CH1_hsqc_df_lessthan_67 = CH3CH1_hsqc_df[CH3CH1_hsqc_df.f1_ppm < 67]
            CH3CH1_hsqc_df_morethan_67 = CH3CH1_hsqc_df[CH3CH1_hsqc_df.f1_ppm >= 67]

            print(
                f"CH3CH1_hsqc_df_lessthan_67.shape[0]: {CH3CH1_hsqc_df_lessthan_67.shape[0]}"
            )
            print(
                f"CH3CH1_hsqc_df_morethan_67.shape[0]: {CH3CH1_hsqc_df_morethan_67.shape[0]}"
            )

            #  label the CH3CH1 groups above 67 ppm as CH1
            hsqc.loc[CH3CH1_hsqc_df_morethan_67.index, "CH1"] = True
            hsqc.loc[CH3CH1_hsqc_df_morethan_67.index, "numProtons"] = 1

            # split the expected molecule CH3CH1 groups based on 67 ppm
            # CH3 groups expected to be below 67 ppm
            CH3CH1_mol_df_lessthan_67 = CH3CH1_mol_df[CH3CH1_mol_df.ppm < 67]
            CH3CH1_mol_df_morethan_67 = CH3CH1_mol_df[CH3CH1_mol_df.ppm >= 67]

            print(
                f"CH3CH1_mol_df_lessthan_67.shape[0]: {CH3CH1_mol_df_lessthan_67.shape[0]}"
            )
            print(
                f"CH3CH1_mol_df_morethan_67.shape[0]: {CH3CH1_mol_df_morethan_67.shape[0]}"
            )

            CH3_mol_df_lessthan_67 = CH3CH1_mol_df_lessthan_67[
                CH3CH1_mol_df_lessthan_67.CH3
            ]
            CH1_mol_df_lessthan_67 = CH3CH1_mol_df_lessthan_67[
                CH3CH1_mol_df_lessthan_67.CH1
            ]

            # print(f"CH3_mol_df_lessthan_67.shape[0]: {CH3_mol_df_lessthan_67.shape[0]}")
            # print(f"CH1_mol_df_lessthan_67.shape[0]: {CH1_mol_df_lessthan_67.shape[0]}")
            # print(
            #     f"CH3CH1_hsqc_df_lessthan_67.shape[0]: {CH3CH1_hsqc_df_lessthan_67.shape[0]}"
            # )

            # if there are no CH1 groups below 67 ppm in the expected molecule then we can set the hsqc CH3CH1 groups to CH3
            if CH1_mol_df_lessthan_67.shape[0] == 0:
                print("no CH1 groups below 67 ppm in expected molecule")
                hsqc.loc[CH3CH1_hsqc_df_lessthan_67.index, "CH3"] = True
                hsqc.loc[CH3CH1_hsqc_df_lessthan_67.index, "numProtons"] = 3

            # check the numbers
            elif CH3CH1_hsqc_df.shape[0] == CH3CH1_mol_df.shape[0]:
                print("CH3CH1_hsqc_df.shape[0] == CH3CH1_mol_df.shape[0]")
                for idx, row in CH3_mol_df.iterrows():
                    # find the closest match in the hsqc dataframe
                    # and update the numProtons column to 3
                    closest_match = CH3CH1_hsqc_df.iloc[
                        (CH3CH1_hsqc_df["f1_ppm"] - row.ppm).abs().argsort()[:1]
                    ]
                    hsqc.loc[closest_match.index, "numProtons"] = 3
                    hsqc.loc[closest_match.index, "CH3"] = True
                    CH3CH1_hsqc_df.drop(closest_match.index, inplace=True)

                for idx, row in CH1_mol_df.iterrows():
                    # find the closest match in the hsqc dataframe
                    # and update the numProtons column to 1
                    closest_match = CH3CH1_hsqc_df.iloc[
                        (CH3CH1_hsqc_df["f1_ppm"] - row.ppm).abs().argsort()[:1]
                    ]
                    hsqc.loc[closest_match.index, "numProtons"] = 1
                    hsqc.loc[closest_match.index, "CH1"] = True
                    CH3CH1_hsqc_df.drop(closest_match.index, inplace=True)

            elif CH3CH1_hsqc_df.shape[0] == CH3CH1_sym_mol_df.shape[0]:
                print("CH3CH1_hsqc_df.shape[0] == CH3CH1_sym_mol_df.shape[0]")
                for idx, row in CH3_sym_mol_df.iterrows():
                    # find the closest match in the hsqc dataframe
                    # and update the numProtons column to 3
                    closest_match = CH3CH1_hsqc_df.iloc[
                        (CH3CH1_hsqc_df["f1_ppm"] - row.ppm).abs().argsort()[:1]
                    ]
                    hsqc.loc[closest_match.index, "numProtons"] = 3
                    hsqc.loc[closest_match.index, "CH3"] = True
                    CH3CH1_hsqc_df.drop(closest_match.index, inplace=True)

                for idx, row in CH1_sym_mol_df.iterrows():
                    # find the closest match in the hsqc dataframe
                    # and update the numProtons column to 1
                    closest_match = CH3CH1_hsqc_df.iloc[
                        (CH3CH1_hsqc_df["f1_ppm"] - row.ppm).abs().argsort()[:1]
                    ]
                    hsqc.loc[closest_match.index, "numProtons"] = 1
                    hsqc.loc[closest_match.index, "CH1"] = True
                    CH3CH1_hsqc_df.drop(closest_match.index, inplace=True)

            elif (
                (CH3CH1_mol_df_morethan_67.shape[0]
                == CH3CH1_hsqc_df_morethan_67.shape[0]) and (CH3CH1_mol_df_morethan_67.shape[0]>0)
            ):
                print("then match up CH1 and CH3s base on ranking ppm values")
                # then match up CH1 and CH3s base on ranking ppm values

                # sort the CH3CH1_hsqc_df_lessthan_67 by f1_ppm
                CH3CH1_hsqc_df_lessthan_67 = CH3CH1_hsqc_df_lessthan_67.sort_values(
                    by=["f1_ppm"]
                )
                CH3CH1_mol_df_lessthan_67 = CH3CH1_mol_df_lessthan_67.sort_values(
                    by=["ppm"]
                )
                CH3CH1_sym_mol_df_lessthan_67 = CH3CH1_sym_mol_df[
                    CH3CH1_sym_mol_df.ppm < 67
                ].sort_values(by=["ppm"])

                mol_idx = CH3CH1_mol_df_lessthan_67.index
                mol_sym_idx = CH3CH1_sym_mol_df_lessthan_67.index
                hsqc_idx = CH3CH1_hsqc_df_lessthan_67.index

                if len(hsqc_idx) == len(mol_idx):
                    hsqc.loc[hsqc_idx, "numProtons"] = CH3CH1_mol_df_lessthan_67.loc[
                        mol_idx, "numProtons"
                    ].values
                    hsqc.loc[hsqc.numProtons == 3, "CH3"] = True
                    hsqc.loc[hsqc.numProtons == 1, "CH1"] = True
                elif len(hsqc_idx) == len(mol_sym_idx):
                    hsqc.loc[hsqc_idx, "numProtons"] = CH3CH1_mol_df_lessthan_67.loc[
                        mol_sym_idx, "numProtons"
                    ].values
                    hsqc.loc[hsqc.numProtons == 3, "CH3"] = True
                    hsqc.loc[hsqc.numProtons == 1, "CH1"] = True
                else:
                    print("No solution found")
                    self.nmrsolution_failed = True
                    self.nmrsolution_error_message = "Line 2020 nmrsolution.py CH3CH1 groups in HSQC do not match expected molecule"
                    self.nmrsolution_error_code = 401

            # elif num_CH3CH1_hsqc < CH3CH1_mol_df.shape[0]:
            #     # attempt to match up CH3 and CH1 based on ppm values
            elif CH3CH1_hsqc_df.shape[0] < CH3CH1_mol_df.shape[0]:

                # loop through the CH3CH1_hsqc_df and try to match with CH3CH1_mol_df 
                CH3CH1_mol_df = CH3CH1_mol_df[~CH3CH1_mol_df.picked]
                for idx, row in CH3CH1_hsqc_df.iterrows():
                    # find the closest match in the CH3CH1_mol_df
                    closest_match = CH3CH1_mol_df.iloc[
                        (CH3CH1_mol_df["ppm"] - row.f1_ppm).abs().argsort()[:1]
                    ]
                    if not closest_match.empty:
                        # update the hsqc dataframe with the closest match
                        hsqc.loc[row.name, "numProtons"] = closest_match.totalNumHs.values[0]
                        hsqc.loc[row.name, "CH3"] = closest_match.CH3.values[0]
                        hsqc.loc[row.name, "CH1"] = closest_match.CH1.values[0]
                        # mark the closest match as picked in the CH3CH1_mol_df

                        CH3CH1_mol_df.loc[closest_match.index, "picked"] = True
                        CH3CH1_mol_df = CH3CH1_mol_df[~CH3CH1_mol_df.picked]
            else:
                print("No solution found")
                self.nmrsolution_failed = True
                self.nmrsolution_error_message = "Line 2042 nmrsolution.py CH3CH1 groups in HSQC do not match expected molecule"
                self.nmrsolution_error_code = 401

        print("self.error_code", self.nmrsolution_error_code)

        if not self.nmrsolution_failed:
            # updated integral and f2_integral columns based on numProtons
            hsqc["integral"] = hsqc["numProtons"]
            hsqc["f2_integral"] = hsqc["numProtons"]
            hsqc["attached_protons"] = hsqc["numProtons"]
            # set the integral to negative if CH2 columns is True
            hsqc.loc[hsqc.CH2, "integral"] = (
                np.abs(hsqc.loc[hsqc.CH2, "integral"].values) * -1
            )
            hsqc.loc[hsqc.CH2, "f2_integral"] = (
                np.abs(hsqc.loc[hsqc.CH2, "f2_integral"].values) * -1
            )

    def assign_CH3_CH1_in_HSQC_using_H1(self) -> None:
        """
        Assigns CH3 and CH1 group information in the HSQC DataFrame using the H1 DataFrame.

        This function matches and labels CH3 and CH1 groups in the HSQC data based on the number of CH3 groups in the H1 data,
        updating the number of protons and group flags accordingly.

        Returns:
            None
        """
        h1 = self.h1
        hsqc = self.hsqc

        num_CH3_in_H1 = h1[h1.numProtons == 3].shape[0]

        if (num_CH3_in_H1 == self.expected_molecule.num_CH3_carbon_atoms) or (
            num_CH3_in_H1 == self.expected_molecule.num_sym_CH3_carbon_atoms
        ):
            print("number of CH3 in H1 matches expected number")
            for index, row in h1[h1.numProtons == 3].iterrows():
                hsqc_rows = hsqc[hsqc.f2_ppm == row.ppm]
                if hsqc_rows.shape[0] == 1:
                    idx = hsqc_rows.index[0]
                    hsqc.loc[idx, "CH3"] = True
                    hsqc.loc[idx, "numProtons"] = 3
                else:
                    print("something went wrong")
                    self.nmrsolution_failed = True
                    self.nmrsolution_error_message = "Error attempting to assign HSQC CH3 protons using H1 data source"
                    self.nmrsolution_error_code = 401

            # now assign the rest of the -1 numProtons in HSQC to 1 and update CH1 flag in HSQC
            hsqc.loc[hsqc.numProtons == -1, "numProtons"] = 1
            hsqc.loc[hsqc.numProtons == 1, "CH1"] = True

        hsqc["attached_protons"] = hsqc["numProtons"]
        hsqc["f2_integral"] = hsqc["numProtons"]
        hsqc["integral"] = hsqc["numProtons"]
        #  set the integral to negative if CH2 columns is True
        hsqc.loc[hsqc.CH2, "integral"] = (
            np.abs(hsqc.loc[hsqc.CH2, "integral"].values) * -1
        )
        hsqc.loc[hsqc.CH2, "f2_integral"] = (
            np.abs(hsqc.loc[hsqc.CH2, "f2_integral"].values) * -1
        )

    def assign_CH3_CH2_CH1_in_HSQC_using_Assignments(self) -> None:
        """
        Assigns CH3, CH2, and CH1 group information in the HSQC DataFrame using the Annotation column.

        This function labels CH3, CH2, and CH1 groups in the HSQC data based on the Annotation values and intensity,
        sets the number of protons, and updates related columns for integrals and attached protons.
        """

        print(
            "====================================\nassign_CH3_CH2_CH1_in_HSQC_using_Assignments\n===================================="
        )
        # print(self.hsqc["Annotation"])
        # print(self.hsqc["Annotation"].str.contains(""))

        hsqc = self.hsqc

        hsqc["numProtons"] = -1

        #  check if all Assignments in the column are "" then just return

        hsqc["CH3"] = False
        hsqc["CH2"] = False
        hsqc["CH1"] = False
        hsqc["CH0"] = False
        hsqc["quaternary"] = False
        hsqc["CH3CH1"] = False

        # set column CH2 to True if intensity < 0
        hsqc.loc[hsqc.intensity < 0, "CH2"] = True

        hsqc.loc[hsqc["Annotation"] == "CH3", "CH3"] = True
        hsqc.loc[hsqc["Annotation"] == "CH2", "CH2"] = True
        hsqc.loc[hsqc["Annotation"] == "CH1", "CH1"] = True

        # set CH3CH1 to true for all other rows
        hsqc.loc[hsqc.CH2 == False, "CH3CH1"] = True

        hsqc.loc[hsqc.CH3, "numProtons"] = 3
        hsqc.loc[hsqc.CH2, "numProtons"] = 2
        hsqc.loc[hsqc.CH1, "numProtons"] = 1

        # set the integral, attached_protons and f2_integral
        hsqc["f2_integral"] = hsqc["numProtons"]
        hsqc["integral"] = hsqc["numProtons"]
        hsqc["attached_protons"] = hsqc["numProtons"]
        hsqc["f2_integral"] = hsqc["numProtons"]
        hsqc["integral"] = hsqc["numProtons"]

        #  set the integral to negative if CH2 columns is True
        hsqc.loc[hsqc.CH2, "integral"] = (
            np.abs(hsqc.loc[hsqc.CH2, "integral"].values) * -1
        )
        hsqc.loc[hsqc.CH2, "f2_integral"] = (
            np.abs(hsqc.loc[hsqc.CH2, "f2_integral"].values) * -1
        )

    def assign_CH3_CH2_CH1_in_HSQC_using_DoubleDept(self) -> None:
        """
        Assigns CH3, CH2, and CH1 group information in the HSQC DataFrame using Double DEPT experiment data.

        This function labels CH3, CH2, and CH1 groups in the HSQC data based on matches to the Double DEPT CH3-only experiment,
        sets the number of protons, and updates related columns for integrals and attached protons.
        """

        hsqc = self.hsqc
        hsqc_CH3 = self.hsqc.copy()

        ddept_ch3_only_df = self.ddept_ch3_only_df

        hsqc["numProtons"] = -1

        # set CH3CH1 for all other rows
        hsqc.loc[hsqc.CH2 == False, "CH3CH1"] = True

        # if ddept_ch3_only not empty set CH3 based on f1_ppm values closest to hsqc f1_ppm values
        if not ddept_ch3_only_df.empty:
            ddept_ch3_only_df = ddept_ch3_only_df.assign(
                CH3=lambda x: False, CH2=lambda x: False, CH1=lambda x: False
            )
            for idx, ppm in zip(ddept_ch3_only_df.index, ddept_ch3_only_df.f1_ppm):
                # find the closest match in the hsqc dataframe
                closest_match = hsqc_CH3.iloc[
                    (hsqc_CH3["f1_ppm"] - ppm).abs().argsort()[:1]
                ]
                print(f"{ppm} closest_match: {closest_match['f1_ppm']}")
                ddept_ch3_only_df.loc[idx, "CH3"] = True
                hsqc.loc[closest_match.index, "CH3"] = True
                # remove row from hsqc_CH3
                hsqc_CH3.drop(closest_match.index, inplace=True)

            # set CH1 based on CH3 values and CH3CH1 values
            CH3CH1_hsqc_df = hsqc[hsqc.CH3CH1]
            CH1_hsqc_df = CH3CH1_hsqc_df[hsqc.CH3 == False]
            hsqc.loc[CH1_hsqc_df.index, "CH1"] = True

        # set the numprotons for CH3 CH2 and CH1 of self.hsqc
        hsqc.loc[hsqc.CH3, "numProtons"] = 3
        hsqc.loc[hsqc.CH2, "numProtons"] = 2
        hsqc.loc[hsqc.CH1, "numProtons"] = 1

        # print("\nhsqc after setting CH3 based on ddept_ch3_only_df\n")
        # print(hsqc[["f1_ppm", "f2_ppm", "numProtons", "CH3", "CH2", "CH1", "CH3CH1"]])

        # set the integral, attached_protons and f2_integral
        hsqc["f2_integral"] = hsqc["numProtons"]
        hsqc["integral"] = hsqc["numProtons"]
        hsqc["attached_protons"] = hsqc["numProtons"]
        hsqc["f2_integral"] = hsqc["numProtons"]
        hsqc["integral"] = hsqc["numProtons"]

        #  set the integral to negative if CH2 columns is True
        hsqc.loc[hsqc.CH2, "integral"] = (
            np.abs(hsqc.loc[hsqc.CH2, "integral"].values) * -1
        )
        hsqc.loc[hsqc.CH2, "f2_integral"] = (
            np.abs(hsqc.loc[hsqc.CH2, "f2_integral"].values) * -1
        )

    def init_class_from_json(self) -> Union[Tuple[bool, str], Tuple[bool, str, int]]:
        """
        Initializes the class instance from a JSON configuration and sets up all experiment dataframes.

        This function loads and processes all relevant NMR experiment data, sets up molecule and experiment attributes,
        and returns a tuple indicating success and a status message or error code if initialization fails.

        Returns:
            tuple: (success, message) where success is a boolean indicating if initialization succeeded,
                   and message is a string or error code describing the result.
        """

        # finally create rdkit molecule from smiles string
        self.dbe = self.expected_molecule.dbe
        self.elements = self.expected_molecule.elements
        self.png = self.expected_molecule.png
        # reconstruct moleculeAtomsStr from elements
        self.moleculeAtomsStr = CalcMolFormula(self.expected_molecule.mol)

        # check what experiments are present
        self.check_what_experiments_are_present()

        print("self.HSQC_data_present:", self.HSQC_data_present)

        # check if hsqc is present
        if not self.HSQC_data_present:
            self.nmrsolution_failed = True
            self.nmrsolution_error_message = "HSQC data missing"
            self.nmrsolution_error_code = 401
            return False, self.nmrsolution_error_message, self.nmrsolution_error_code

        # create short views of the dataframes
        hsqc_not_empty, return_message = self.check_hsqc_df_not_empty()

        if self.hsqc.empty:
            self.nmrsolution_failed = True
            self.nmrsolution_error_message = "HSQC data missing"
            self.nmrsolution_error_code = 401
            return False, return_message, self.nmrsolution_error_code

        # find CH2 groups in hsqc_df
        self.hsqc.loc[self.hsqc.intensity < 0, "CH2"] = True
        unique_idxs, unique_ch2s = self.find_and_group_CH2s(self.hsqc)

        # replace hsqc CH2 values in f1_ppm of HSQC
        for idx, ch2 in zip(unique_idxs, unique_ch2s):
            self.hsqc.loc[idx, "f1_ppm"] = np.mean(ch2)

        # count how many unique CH2 groups are in the hsqc dataframe
        # create missing dataframes
        self.c13_df, self.c13_from_hsqc = self.init_c13_from_hsqc_and_hmbc()

        self.h1_df = self.init_h1_from_pureshift()

        self.pureshift_df = self.init_pureshift_from_h1()

        (
            self.h1_df,
            self.pureshift_df,
        ) = self.init_h1_and_pureshift_from_c13_hsqc_hmbc_cosy()

        self.c13_from_hsqc = True

        c13_not_empty, return_message = self.check_c13_not_empty()
        if not c13_not_empty:
            return False, return_message

        h1_pureshift_not_empty, return_message = self.check_h1_pureshift_not_empty()
        if not h1_pureshift_not_empty:
            return False, return_message

        # check if integral column was missing from h1_df, if so replace integral values by using hsqc
        self.numCarbonGroups = self.c13.shape[0]
        self.numProtonGroups = self.h1.shape[0]

        self.symmetric_molecule = False
        if self.elements["C"] > 0 and self.elements["C"] > self.numCarbonGroups:
            self.symmetric_molecule = True

        self.c13["attached_protons"] = 0
        self.c13["ppmH1s"] = None

        self.hsqc["f1_i"] = 0
        self.hsqc["f2_i"] = 0
        self.hsqc["f2p_i"] = 0
        self.hsqc["f1C_i"] = 0
        self.hsqc["f2H_i"] = 0
        self.hsqc["f2Cp_i"] = 0
        self.hsqc["f2p_ppm"] = 0

        self.hmbc = self.hmbc_df[["f1_ppm", "f2_ppm", "intensity"]].copy()
        self.hmbc["f2p_ppm"] = 0
        self.hmbc["f1_i"] = 0
        self.hmbc["f2_i"] = 0
        self.hmbc["f2p_i"] = 0

        self.pureshift = self.pureshift_df.copy()

        self.cosy = self.cosy_df[["f1_ppm", "f2_ppm", "intensity"]].copy()
        self.cosy = self.cosy.assign(f1_i=lambda x: 0)
        self.cosy = self.cosy.assign(f1p_i=lambda x: 0)
        self.cosy = self.cosy.assign(f2_i=lambda x: 0)
        self.cosy = self.cosy.assign(f2p_i=lambda x: 0)

        self.cosy = self.cosy.assign(f1p_ppm=lambda x: np.nan)
        self.cosy = self.cosy.assign(f2p_ppm=lambda x: np.nan)

        # f1H_i	f2H_i	f1Cp_i	f2Cp_i
        self.cosy["f1H_i"] = ""
        self.cosy["f2H_i"] = ""
        self.cosy["f1Cp_i"] = ""
        self.cosy["f2Cp_i"] = ""

        self.c13["max_bonds"] = 4

        self.h1["integral"] = self.h1_df["integral"]
        self.h1["numProtons"] = self.h1_df["numProtons"]
        self.h1["jCouplingClass"] = self.h1_df["jCouplingClass"]
        self.h1["jCouplingVals"] = self.h1_df["jCouplingVals"]
        # self.h1["range"] = self.h1_df["range"]
        # tidy up chemical shift values by replacing cosy, hsqc and hmbc picked peaks with values from c13ppm and h1ppm dataframes

        # HMBC
        if not self.exact_ppm_values:
            self.hmbc = self.tidyup_ppm_values(
                self.hmbc,
                self.c13.ppm.tolist(),
                "f1_ppm",
                ppm_tolerance=self.problemdata_json.carbonSeparation,
            )
            self.hmbc = self.tidyup_ppm_values(
                self.hmbc,
                self.h1.ppm.tolist(),
                "f2_ppm",
                ppm_tolerance=self.problemdata_json.protonSeparation,
            )

            self.hmbc.drop(self.hmbc[self.hmbc.f1_ppm_prob == 0].index, inplace=True)
            self.hmbc.drop(self.hmbc[self.hmbc.f2_ppm_prob == 0].index, inplace=True)

        else:
            print("Exact ppm values only, no tidy up required")

        # HSQC
        if not self.exact_ppm_values:
            self.hsqc = self.tidyup_ppm_values(
                self.hsqc,
                self.c13.ppm.tolist(),
                "f1_ppm",
                self.problemdata_json.carbonSeparation,
            )
            self.hsqc = self.tidyup_ppm_values(
                self.hsqc,
                self.h1.ppm.tolist(),
                "f2_ppm",
                self.problemdata_json.protonSeparation,
            )

            # tidy up cosy H1 shifts
            self.cosy = self.tidyup_ppm_values(
                self.cosy,
                self.h1.ppm.tolist(),
                "f1_ppm",
                self.problemdata_json.protonSeparation,
            )
            self.cosy = self.tidyup_ppm_values(
                self.cosy,
                self.h1.ppm.tolist(),
                "f2_ppm",
                self.problemdata_json.protonSeparation,
            )

            # check if any probability equals zero and remove the row
            # because it is likely that proton is not connected directly to carbon
            self.cosy.drop(self.cosy[self.cosy.f1_ppm_prob == 0].index, inplace=True)
            self.cosy.drop(self.cosy[self.cosy.f2_ppm_prob == 0].index, inplace=True)
        
        else:
            print("Exact ppm values only, no tidy up required")

        # add index columns to h1
        self.h1["label"] = ["H" + str(i) for i in self.h1.index]
        self.h1["f1H_i"] = ["H" + str(i) for i in self.h1.index]
        self.h1["f2H_i"] = ["H" + str(i) for i in self.h1.index]
        self.h1["f1_i"] = self.h1.index
        self.h1["f2_i"] = self.h1.index

        # add lookup dicts for dataframe h1
        self.H1ppmH1label = dict(zip(self.h1.ppm, self.h1.label))
        self.H1labelH1ppm = dict(zip(self.h1.label, self.h1.ppm))

        self.H1indexH1label = dict(zip(self.h1.index, self.h1.label))
        self.H1labelH1index = dict(zip(self.h1.label, self.h1.index))

        self.H1ppmH1index = dict(zip(self.h1.ppm, self.h1.index))
        self.H1indexH1ppm = dict(zip(self.h1.index, self.h1.ppm))

        # add index columns to c13
        self.c13["label"] = ["C" + str(i) for i in self.c13.index]
        self.c13["f2C_i"] = ["C" + str(i) for i in self.c13.index]
        self.c13["f1C_i"] = ["C" + str(i) for i in self.c13.index]
        self.c13["f1_i"] = self.c13.index
        self.c13["f2_i"] = self.c13.index

        # add lookup dicts for dataframe c13
        self.C13ppmC13label = dict(zip(self.c13.ppm, self.c13.label))
        self.C13labelC13ppm = dict(zip(self.c13.label, self.c13.ppm))

        self.C13indexC13label = dict(zip(self.c13.index, self.c13.label))
        self.C13labelC13index = dict(zip(self.c13.label, self.c13.index))

        self.C13ppmC13index = dict(zip(self.c13.ppm, self.c13.index))
        self.C13indexC13ppm = dict(zip(self.h1.index, self.c13.ppm))

        print("self.C13indexC13ppm", self.C13indexC13ppm)

        # open and read excel file into datframe and then process it
        # with open(self.yamlFiles[0], 'r') as fp:
        #    info = yaml.safe_load(fp)
        #    self.init_variables_from_dict(info)

        # add index columns to hsqc
        for i in self.hsqc.index:
            print(i)
            self.hsqc.loc[i, "f2_i"] = self.H1ppmH1index[self.hsqc.loc[i, "f2_ppm"]]
            self.hsqc.loc[i, "f2H_i"] = self.H1ppmH1label[self.hsqc.loc[i, "f2_ppm"]]
            self.hsqc.loc[i, "f1_i"] = self.C13ppmC13index[self.hsqc.loc[i, "f1_ppm"]]
            self.hsqc.loc[i, "f1C_i"] = self.C13ppmC13label[self.hsqc.loc[i, "f1_ppm"]]
            self.hsqc.loc[i, "f2Cp_i"] = self.C13ppmC13label[self.hsqc.loc[i, "f1_ppm"]]

        self.hsqc["f2p_i"] = self.hsqc["f1_i"]
        self.hsqc["f2p_ppm"] = self.hsqc["f1_ppm"]

        # add lookup dicts for hsqc
        self.hsqcH1ppmC13index = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2p_i))
        self.hsqcH1ppmC13label = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2Cp_i))
        self.hsqcH1ppmC13ppm = dict(zip(self.hsqc.f2_ppm, self.hsqc.f2p_ppm))

        self.hsqcH1indexC13index = dict(zip(self.hsqc.f2_i, self.hsqc.f2p_i))
        self.hsqcH1indexC13ppm = dict(zip(self.hsqc.f2_i, self.hsqc.f2p_ppm))
        self.hsqcH1indexC13label = dict(zip(self.hsqc.f2_i, self.hsqc.f2Cp_i))

        self.hsqcH1labelC13label = dict(zip(self.hsqc.f2H_i, self.hsqc.f1C_i))
        self.hsqcH1labelC13index = dict(zip(self.hsqc.f2H_i, self.hsqc.f1_i))
        self.hsqcH1labelC13ppm = dict(zip(self.hsqc.f2H_i, self.hsqc.f1_ppm))

        if not self.exact_ppm_values:
            # add index columns to hsqc_clipcosy
            self.hsqc_clipcosy = self.tidyup_hsqc_clipcosy(
                self.hsqc_clipcosy_df, self.c13, self.h1
            )
        else:
            print("Exact ppm values only, no tidy up required")
            self.hsqc_clipcosy = self.hsqc_clipcosy_df.copy()

        self.hsqc_clipcosy = self.process_hsqc_clipcosy(self.hsqc_clipcosy, self.hsqc)
        # loguru.logger.debug(f"hsqc_clipcosy\n{self.hsqc_clipcosy}")
        # loguru.logger.debug(f"hsqc\n{self.hsqc}")

        # process ddept_ch3_only
        self.ddept_ch3_only = self.ddept_ch3_only_df.copy()

        if not self.exact_ppm_values:
            self.ddept_ch3_only = self.tidyup_ppm_values(
                self.ddept_ch3_only,
                self.c13["ppm"],
                "f1_ppm",
                ppm_tolerance=self.problemdata_json.carbonSeparation,
            )
            self.ddept_ch3_only = self.tidyup_ppm_values(
                self.ddept_ch3_only,
                self.h1["ppm"],
                "f2_ppm",
                ppm_tolerance=self.problemdata_json.protonSeparation,
            )
        else:
            print("Exact ppm values only, no tidy up required")

        self.ddept_ch3_only = self.process_ddept_ch3_only(
            self.ddept_ch3_only, self.hsqc
        )
        # loguru.logger.debug(f"ddept_ch3_only\n{self.ddept_ch3_only}")

        # add index columns to cosy
        # fill in cosy dataframe
        for hppm in self.h1.ppm:
            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1_i"
            ] = self.H1ppmH1index.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2_i"
            ] = self.H1ppmH1index.get(hppm)

            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1p_i"
            ] = self.hsqcH1ppmC13index.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2p_i"
            ] = self.hsqcH1ppmC13index.get(hppm)

            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1p_ppm"
            ] = self.hsqcH1ppmC13ppm.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2p_ppm"
            ] = self.hsqcH1ppmC13ppm.get(hppm)

            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1H_i"
            ] = self.H1ppmH1label.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2H_i"
            ] = self.H1ppmH1label.get(hppm)

            self.cosy.loc[
                self.cosy[self.cosy.f1_ppm == hppm].index, "f1Cp_i"
            ] = self.hsqcH1ppmC13label.get(hppm)
            self.cosy.loc[
                self.cosy[self.cosy.f2_ppm == hppm].index, "f2Cp_i"
            ] = self.hsqcH1ppmC13label.get(hppm)

        # combine cosy and hsqc_clipcosy dataframes
        # loguru.logger.debug(f"\ncosy\n{self.cosy}")
        self.cosy = self.process_and_merge_cosy_dfs(self.hsqc_clipcosy, self.cosy)
        # loguru.logger.debug(f"\ncosy\n{self.cosy}")

        # loguru.logger.debug("Add CH3 to HSQC from DDEPT_CH3_ONLY")
        # loguru.logger.debug(f"\nhsqc\n{self.hsqc}")
        # self.init_CH_CH3_HSQC_from_DDEPT_CH3_only()
        # loguru.logger.debug(f"\nhsqc\n{self.hsqc}")

        # add index columns to hmbc
        self.hmbc["f1C_i"] = ""
        self.hmbc["f2H_i"] = ""
        self.hmbc["f2Cp_i"] = ""

        # fill in hmbc dataframe
        for i in self.hmbc.index:
            self.hmbc.loc[i, "f2p_ppm"] = self.hsqcH1ppmC13ppm.get(
                self.hmbc.loc[i, "f2_ppm"]
            )
            self.hmbc.loc[i, "f2p_i"] = self.hsqcH1ppmC13index.get(
                self.hmbc.loc[i, "f2_ppm"]
            )

            self.hmbc.loc[i, "f2_i"] = self.H1ppmH1index.get(self.hmbc.loc[i, "f2_ppm"])
            self.hmbc.loc[i, "f1_i"] = self.C13ppmC13index.get(
                self.hmbc.loc[i, "f1_ppm"]
            )

            self.hmbc.loc[i, "f1C_i"] = self.C13ppmC13label.get(
                self.hmbc.loc[i, "f1_ppm"]
            )
            self.hmbc.loc[i, "f2H_i"] = self.H1ppmH1label.get(
                self.hmbc.loc[i, "f2_ppm"]
            )
            self.hmbc.loc[i, "f2Cp_i"] = self.hsqcH1ppmC13label.get(
                self.hmbc.loc[i, "f2_ppm"]
            )

        # add f1p_ppm and f1p_i columns to self.h1 based on hsqc dataframe
        self.h1["f1p_ppm"] = -1e6
        self.h1["f1p_i"] = -1
        print("Assigning f1p_ppm and f1p_i to h1 dataframe based on hsqc dataframe")
        print(f"self.h1.:\n {self.h1}")
        print(f"self.hsqc:\n {self.hsqc}")
        for i in self.h1.index:
            f2_ppm = self.h1.loc[i, "ppm"]
            if self.hsqc[self.hsqc.f2_ppm == f2_ppm].empty:
                print(f"f2_ppm {f2_ppm} not found in hsqc dataframe")
                self.h1.loc[i, "f1p_ppm"] = -1e6
                self.h1.loc[i, "f1p_i"] = -1
                break

            self.h1.loc[i, "f1p_ppm"] = self.hsqc[
                self.hsqc.f2_ppm == f2_ppm
            ].f1_ppm.values[0]
            self.h1.loc[i, "f1p_i"] = self.hsqc[self.hsqc.f2_ppm == f2_ppm].f1_i.values[
                0
            ]

        # if "f2_integral" not in self.hsqc.columns:
        #     self.define_hsqc_f2integral()
        # self.check_and_ammend_CH3groups_greater_than_3()
        # self.check_and_ammend_CH2groups_greater_than_2()
        # self.check_and_ammend_CH1groups_greater_than_1()

        # check if h1 contains any -1e6 values in f1p_ppm and if so return false and error message
        if (self.h1["f1p_ppm"] == -1e6).any():
            self.nmrsolution_failed = True
            self.nmrsolution_error_message = f"Misalignment between 1D-proton and HSQC values. Please check the the 1-D proton and hsqc spectra around {f2_ppm:.4} ppm for alignment."
            self.nmrsolution_error_code = 401
            return False, self.nmrsolution_error_message

        return True, "ok"

    def process_and_merge_cosy_dfs(self, df: pd.DataFrame, cosy: pd.DataFrame) -> pd.DataFrame:
        """
        Merges two COSY DataFrames, aligning columns and removing duplicate cross-peaks.

        This function drops columns from the first DataFrame that are not present in the second,
        concatenates both DataFrames, and removes duplicate rows based on the 'f1_ppm' and 'f2_ppm' columns.

        Args:
            df (pd.DataFrame): The first COSY DataFrame to merge.
            cosy (pd.DataFrame): The second COSY DataFrame to merge.

        Returns:
            pd.DataFrame: The merged COSY DataFrame with duplicates removed.
        """
        # Remove columns in df that are not in the cosy dataframe
        df = df.drop(columns=set(df.columns).difference(set(cosy.columns)))

        # Merge the two dataframes
        cosy_all = pd.concat([cosy, df], axis=0, ignore_index=True)

        # Remove duplicate rows
        cosy_all = cosy_all.drop_duplicates(subset=["f1_ppm", "f2_ppm"], keep="first")

        return cosy_all

    def process_ddept_ch3_only(self, df: pd.DataFrame, hsqc: pd.DataFrame) -> pd.DataFrame:
        """
        Adds index and group information to the DDEPT CH3-only DataFrame using corresponding indices from the HSQC DataFrame.

        This function populates the DDEPT CH3-only DataFrame with index columns and group assignments by matching chemical shifts to the HSQC DataFrame.

        Args:
            df (pd.DataFrame): The DDEPT CH3-only DataFrame to process.
            hsqc (pd.DataFrame): The HSQC DataFrame used for index and group assignment.

        Returns:
            pd.DataFrame: The processed DDEPT CH3-only DataFrame with added indices and group columns.
        """
        if df.empty:
            return df
        df["f1_i"] = 0
        df["f2_i"] = 0
        df["f2p_i"] = 0
        df["f2p_ppm"] = 0

        df["f2H_i"] = ""

        df["f1C_i"] = ""
        df["f2Cp_i"] = ""
        df["integral"] = 3
        df["numProtons"] = 3
        df["CH3"] = True
        df["CH3CH1"] = True
        df["CH2"] = False

        for f2_ppm, f1_ppm, f2p_ppm, f2p_i, f1_i, f2_i in zip(
            hsqc.f2_ppm, hsqc.f1_ppm, hsqc.f2p_ppm, hsqc.f2p_i, hsqc.f1_i, hsqc.f2_i
        ):
            df.loc[df.f1_ppm == f1_ppm, "f1_i"] = f1_i
            df.loc[df.f2_ppm == f2_ppm, "f2_i"] = f2_i
            df.loc[df.f2_ppm == f2_ppm, "f2p_ppm"] = f2p_ppm
            df.loc[df.f2_ppm == f2_ppm, "f2p_i"] = f2p_i
            df.loc[df.f2_ppm == f2_ppm, "f2Cp_i"] = f"C{f2p_i}"
            df.loc[df.f2_ppm == f2_ppm, "f2H_i"] = f"H{f2_i}"

            df.loc[df.f2_ppm == f2_ppm, "f1H_i"] = f"H{f2_i}"
            df.loc[df.f1_ppm == f1_ppm, "f1C_i"] = f"C{f1_i}"

        return df

    def tidyup_hsqc_clipcosy(
        self,
        hsqc_clipcosy_df: pd.DataFrame,
        c13: pd.DataFrame,
        h1: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Cleans and aligns the HSQC-CLIP-COSY DataFrame by filtering for negative intensity crosspeaks and adjusting ppm values.

        This function filters the input DataFrame for negative intensity values, then aligns the chemical shift columns to the nearest values in the C13 and H1 DataFrames using specified tolerances.

        Args:
            hsqc_clipcosy_df (pd.DataFrame): The HSQC-CLIP-COSY DataFrame to tidy.
            c13 (pd.DataFrame): The C13 DataFrame for reference ppm values.
            h1 (pd.DataFrame): The H1 DataFrame for reference ppm values.

        Returns:
            pd.DataFrame: The tidied HSQC-CLIP-COSY DataFrame with aligned ppm values.
        """

        # Filter the dataframe for intensity < 0
        hsqc_clipcosy = hsqc_clipcosy_df[hsqc_clipcosy_df.intensity < 0].copy()


        if not self.exact_ppm_values:
            # Tidy up ppm values
            hsqc_clipcosy = self.tidyup_ppm_values(
                hsqc_clipcosy,
                c13["ppm"],
                "f1_ppm",
                ppm_tolerance=self.problemdata_json.carbonSeparation,
            )
            hsqc_clipcosy = self.tidyup_ppm_values(
                hsqc_clipcosy,
                h1["ppm"],
                "f2_ppm",
                ppm_tolerance=self.problemdata_json.protonSeparation,
            )
        else:
            print("Exact ppm values only, no tidy up required")

        return hsqc_clipcosy

    def process_hsqc_clipcosy(self, df, hsqc):
        if df.empty:
            return df
        ## adds indices to the hsqc_clipcosy dataframe using corresponding indices from hsqc dataframe
        df["f1_i"] = 0
        df["f1p_i"] = 0
        df["f2_i"] = 0
        df["f2p_i"] = 0
        df["f2p_ppm"] = 0

        df["f1H_i"] = ""
        df["f2H_i"] = ""

        df["f1Cp_i"] = ""
        df["f2Cp_i"] = ""
        # df["signaltype"] = "compound"

        # copy f1_ppm column to f1p_ppm
        df["f1p_ppm"] = df["f1_ppm"]

        for f2_ppm, f1_ppm, f2p_ppm, f2p_i, f1_i, f2_i in zip(
            hsqc.f2_ppm, hsqc.f1_ppm, hsqc.f2p_ppm, hsqc.f2p_i, hsqc.f1_i, hsqc.f2_i
        ):
            df.loc[df.f1_ppm == f1_ppm, "f1_ppm"] = f2_ppm
            df.loc[df.f2_ppm == f2_ppm, "f2_i"] = f2_i
            df.loc[df.f2_ppm == f2_ppm, "f2p_i"] = f2p_i
            df.loc[df.f2_ppm == f2_ppm, "f2p_ppm"] = f2p_ppm
            df.loc[df.f2_ppm == f2_ppm, "f2Cp_i"] = f"C{f2p_i}"
            df.loc[df.f2_ppm == f2_ppm, "f2H_i"] = f"H{f2_i}"

            df.loc[df.f1_ppm == f2_ppm, "f1_i"] = f2_i
            df.loc[df.f1_ppm == f2_ppm, "f1p_i"] = f2p_i
            df.loc[df.f1_ppm == f2_ppm, "f1H_i"] = f"H{f2_i}"
            df.loc[df.f1_ppm == f2_ppm, "f1Cp_i"] = f"C{f2p_i}"

        return df

    def find_nearest(self, array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]


def copy_over_values_c13_all_to_hetero2D(df, c13_all):
    """
    Copy over values from c13_all to hmbc dataframe
    """
    if df.empty:
        columns = list(df.columns)
        columns.append("f1_atom_idx")
        columns.append("f1_sym_atom_idx")
        columns.append("f1_atomNumber")
        columns.append("f1_sym_atomNumber")
        columns.append("f1_x")
        columns.append("f1_y")
        columns.append("f2_atom_idx")
        columns.append("f2_sym_atom_idx")
        columns.append("f2_atomNumber")
        columns.append("f2_sym_atomNumber")
        columns.append("f2_x")
        columns.append("f2_y")
        df = pd.DataFrame(columns=columns)
        return df

    df["f1_atom_idx"] = None
    df["f1_sym_atom_idx"] = None
    df["f1_atomNumber"] = None
    df["f1_sym_atomNumber"] = None
    df["f1_x"] = None
    df["f1_y"] = None
    df["f2_atom_idx"] = None
    df["f2_sym_atom_idx"] = None
    df["f2_atomNumber"] = None
    df["f2_sym_atomNumber"] = None
    df["f2_x"] = None
    df["f2_y"] = None


    for idx, row in c13_all.iterrows():
        c13_ppm = row["ppm"]
        f1_idx = df[df["f1_ppm"] == c13_ppm].index
        for ii in f1_idx:
            df.loc[ii, "f1_atom_idx"] = row["atom_idx"]
            df.at[ii, "f1_sym_atom_idx"] = row["sym_atom_idx"]
            df.at[ii, "f1_atomNumber"] = row["atomNumber"]
            df.at[ii, "f1_x"] = row["x"]
            df.at[ii, "f1_y"] = row["y"]

        f2p_idx = (df[df["f2p_ppm"] == c13_ppm]).index.tolist()
        for ii in f2p_idx:
            df.at[ii, "f2_atom_idx"] = row["atom_idx"]
            df.at[ii, "f2_sym_atom_idx"] = row["sym_atom_idx"]
            df.at[ii, "f2_atomNumber"] = row["atomNumber"]
            df.at[ii, "f2_x"] = row["x"]
            df.at[ii, "f2_y"] = row["y"]


    return df


def copy_over_values_c13_all_to_homo2D(df, c13_all):
    """
    Copy over values from c13_all to hmbc dataframe
    """
    if df.empty:
        columns = list(df.columns)
        columns.append("f1_atom_idx")
        columns.append("f1_asym_tom_idx")
        columns.append("f1_atomNumber")
        columns.append("f1_sym_atomNumber")
        columns.append("f1_x")
        columns.append("f1_y")
        columns.append("f2_atom_idx")
        columns.append("f2_sym_atom_idx")
        columns.append("f2_atomNumber")
        columns.append("f2_sym_atomNumber")
        columns.append("f2_x")
        columns.append("f2_y")
        df = pd.DataFrame(columns=columns)
        return df
    for idx, row in c13_all.iterrows():
        c13_ppm = row["ppm"]
        f1_idx = df[df["f1p_ppm"] == c13_ppm].index
        for ii in f1_idx:
            df.at[ii, "f1_atom_idx"] = row["atom_idx"]
            df.at[ii, "f1_sym_atom_idx"] = row["sym_atom_idx"]
            df.at[ii, "f1_atomNumber"] = row["atomNumber"]
            df.at[ii, "f1_sym_atomNumber"] = row["sym_atomNumber"]
            df.at[ii, "f1_x"] = row["x"]
            df.at[ii, "f1_y"] = row["y"]

        f2p_idx = df[df["f2p_ppm"] == c13_ppm].index
        for ii in f2p_idx:
            df.at[ii, "f2_atom_idx"] = row["atom_idx"]
            df.at[ii, "f2_sym_atom_idx"] = row["sym_atom_idx"]
            df.at[ii, "f2_atomNumber"] = row["atomNumber"]
            df.at[ii, "f2_sym_atomNumber"] = row["sym_atomNumber"]
            df.at[ii, "f2_x"] = row["x"]
            df.at[ii, "f2_y"] = row["y"]

    return df


def create_network_graph(c13, h1):

    # create nodes of graph
    G2 = nx.Graph()

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
        idx = int(i)
        # iupacLabel = iupacLabel.replace("\'", "\\'")
        G2.add_node(idx)
        G2.nodes[idx]["symbol"] = "C"
        # check if the atomNumber is a numeric value string

        if isinstance(atomNumber, (int, float)):
            G2.nodes[idx]["atomNumber"] = int(atomNumber)
        elif isinstance(atomNumber, str):
            if atomNumber.isnumeric():
                G2.nodes[idx]["atomNumber"] = int(atomNumber)
            else:
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

    if "jCouplingVals" not in df2.columns:
        df1["jCouplingVals"] = ""
        df1["jCouplingClass"] = ""
    else:
        df1["jCouplingVals"] = df2["jCouplingVals"].values
        df1["jCouplingClass"] = df2["jCouplingClass"].values

    return df1


def assign_jcouplings_to_c13(c13, hsqc):

    c13["jCouplingVals"] = ""
    c13["jCouplingClass"] = ""

    if ("jCouplingVals" not in hsqc.columns) or ("jCouplingClass" not in hsqc.columns):
        return c13
    for hsqc_idx, row in hsqc.iterrows():
        # find the row in c13 with the same atom_idx
        c13_rows = c13[c13.atom_idx == row["f1_atom_idx"]]
        c13_idx = c13_rows.index
        c13.loc[c13_idx, "jCouplingVals"] = row["jCouplingVals"]
        c13.loc[c13_idx, "jCouplingClass"] = row["jCouplingClass"]

    return c13


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


def add_cosy_edges_to_graph(G2: nx.Graph, keep: dict) -> nx.Graph:

    # add cosy edges to graph
    for k, v in keep.items():
        for vv in v:
            node1 = int(vv[0])
            node2 = int(vv[1])
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
            node1 = int(vv[0])
            node2 = int(vv[1])
            if G2.has_edge(node1, node2):
                # set cosy attribute to True
                G2.get_edge_data(node1, node2)["cosy"] = True
            else:
                G2.add_edge(node1, node2)
                # set cosy attribute to True
                G2.get_edge_data(node1, node2)["cosy"] = True

    return G2


def add_all_hmbc_edges_to_graph(
    G2, hmbc: pd.DataFrame, h1_all: pd.DataFrame, c13_sorted: pd.DataFrame
):
    """
    Add edges to the graph G2 from the hmbc dataframe
    """
    G2 = add_hmbc_edges_to_graph(G2, hmbc, h1_all, c13_sorted)

    return G2


def add_hmbc_edges_to_graph(
    G2: nx.Graph, hmbc: pd.DataFrame, h1: pd.DataFrame, c13: pd.DataFrame
) -> nx.Graph:
    # add HMBC edges to G2
    print("Adding HMBC edges 1")
    count = 0
    for i, row in hmbc.iterrows():
        # find the rows in h1 that have the same ppm vaslue as the f1_ppm value in the row
        df1 = h1[h1.ppm == row.f2_ppm]
        df2 = c13[c13.f1_ppm == row.f1_ppm]

        for j, row1 in df1.iterrows():
            for k, row2 in df2.iterrows():
                e1 = int(row1.atom_idx)
                e2 = int(row2.atom_idx)
                ppm1 = row1.f1_ppm
                ppm2 = row2.f1_ppm
                if (e1 == e2) or (ppm1 == ppm2):
                    continue
                if not G2.has_edge(e1, e2):
                    G2.add_edge(e1, e2, weight=1)
                    count += 1
                G2.get_edge_data(e1, e2)["hmbc"] = True

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
