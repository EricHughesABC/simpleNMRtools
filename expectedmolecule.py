# create nmrmol which is a subclass of rdkit.Chem.rdchem.Mol
import os
import platform
import subprocess
from io import StringIO
import numpy as np
import pandas as pd


import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


# from decorators import count_calls


# import java
import javaUtils

from functools import lru_cache

from globals import svgDimensions

XYDIM = 800


@lru_cache(maxsize=None)
def calc_c13_chemical_shifts_using_nmrshift2D(molstr: str) -> pd.DataFrame:

    mol_df = javaUtils.calculate_nmrpredictions_nmrshiftdb(molstr)

    if len(mol_df) > 0:
        return mol_df
    else:
        return False


class expectedMolecule:
    def copy_xy3_to_molecule(self, new_xy3, molprops_df):
        # copy new_xy3 to molprops_df x and y columns
        for idx, row in molprops_df.iterrows():
            if idx in new_xy3:
                molprops_df.at[idx, "x"] = new_xy3[idx][0]
                molprops_df.at[idx, "y"] = new_xy3[idx][1]
            else:
                print(f"idx {idx} not found in new_xy3")

    def __init__(
        self,
        smiles_str,
        carbonAtomsInfo=None,
        is_smiles=True,
        mnova_c13predictions=None,
        predict_from_nmrshiftdb=False,
        JEOL_predict=False
    ):
        self.nmrshiftdb_failed = False
        self.nmrshiftdb_failed_message = ""
        self.nmrshiftdb_failed_code = 200

        if is_smiles:
            self.smiles_str = smiles_str
            self.mol = Chem.MolFromSmiles(smiles_str)

        else:
            self.mol_str = smiles_str
            self.mol = Chem.MolFromMolBlock(smiles_str)
            # create smiles string from rdkit mol

            self.smiles_str = Chem.MolToSmiles(self.mol)

        # self.java_available = java.JAVA_AVAILABLE
        # self.java_command = java.JAVA_COMMAND
        # self.java_params = java.JAVA_PARAMS

        # fix coordinates of molecule before creating png

        # self.mol = Chem.AddHs(self.mol)
        # AllChem.EmbedMolecule(self.mol, randomSeed=3)
        # self.mol = Chem.RemoveHs(self.mol)

        # self.mol.Compute2DCoords()

        self.svg_str, self.xy3, self.xy3_allatoms = self.create_svg_string(
            molWidth=svgDimensions.MOLWIDTH,
            molHeight=svgDimensions.MOLHEIGHT,
            svgWidth=svgDimensions.SVGWIDTH,
            svgHeight=svgDimensions.SVGHEIGHT,
        )
        # print("self.svg_str\n", self.svg_str)

        # copy_xy3_to_molecule(self, self.xy3)

        self.png = Draw.MolToImage(self.mol, size=(XYDIM, XYDIM))

        self.dbe = self.calc_dbe()
        self.elements = self.init_elements_dict()
        self.num_carbon_atoms = self.elements.get("C", 0)
        self.num_hydrogen_atoms = self.elements.get("H", 0)

        molprops = [
            [
                atom.GetIdx(),
                atom.GetIdx(),
                atom.GetNumImplicitHs(),
                atom.GetTotalNumHs(),
                atom.GetDegree(),
                atom.GetHybridization(),
                atom.GetIsAromatic(),
            ]
            for atom in self.mol.GetAtoms()
            if atom.GetSymbol() == "C"
        ]

        self.molprops_df = pd.DataFrame(
            data=molprops,
            columns=[
                "idx",
                "atom_idx",
                "implicitHs",
                "totalNumHs",
                "degree",
                "hybridization",
                "aromatic",
            ],
        )

        #  add 1 to atom_idx

        # self.molprops_df["atom_idx"] = self.molprops_df["atom_idx"] + 1
        if JEOL_predict:    
            self.molprops_df["atom_idx"] = self.molprops_df["atom_idx"]
        else:
            self.molprops_df["atom_idx"] = self.molprops_df["atom_idx"] + 1

        self.molprops_df["numProtons"] = self.molprops_df["totalNumHs"].astype(int)
        self.molprops_df = self.molprops_df.set_index(["idx"])
        self.molprops_df["idx"] = self.molprops_df.index

        # define quaternary carbon atoms column for totalNumHs == 0
        self.molprops_df["quaternary"] = False
        self.molprops_df["CH0"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 0, "quaternary"] = True
        self.molprops_df["CH0"] = self.molprops_df["quaternary"]

        # define CH1 column for totalNumHs == 1
        self.molprops_df["CH1"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 1, "CH1"] = True

        self.molprops_df["CH2"] = False
        # set CH2 to True if carbon has 2 protons attached
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 2, "CH2"] = True

        # define CH3 column for totalNumHs == 3
        self.molprops_df["CH3"] = False
        self.molprops_df.loc[self.molprops_df["totalNumHs"] == 3, "CH3"] = True

        # define CH3CH1 column where CH3 or CH1 are True
        self.molprops_df["CH3CH1"] = False
        self.molprops_df["CH3CH1"] = self.molprops_df["CH3"] | self.molprops_df["CH1"]

        self.molprops_df["picked"] = False

        # decide where to calculate C13 NMR chemical shifts from
        # if mnova_c13predictions is not None then use mnova_c13predictions
        # else use nmrshift2D to calculate C13 NMR chemical shifts

        if (mnova_c13predictions is not None) and (predict_from_nmrshiftdb == False):

            print("=====================================================")
            print("Using mnova_c13predictions to calculate C13 NMR chemical shifts")
            print("=====================================================")

            # self.c13ppm = {}
            # self.c13ppm['mean'] = mnova_c13predictions["ppm"].values

            # self.molprops_df["atom_idx"] = self.molprops_df.index

            print("n=====================================================")
            print("mnova_c13predictions", mnova_c13predictions.columns)
            print("=====================================================")

            self.molprops_df = pd.merge(
                mnova_c13predictions, self.molprops_df, on="atom_idx", how="left"
            )
            self.molprops_df["numProtons"] = self.molprops_df["numProtons_x"]
            # reset the index to atom index
            self.molprops_df.set_index("idx", inplace=True)

        else:
            print("Using nmrshift2D to calculate C13 NMR chemical shifts")
            self.c13ppm = self.calculated_c13_chemical_shifts()
            try:
                print("self.c13ppm", len(self.c13ppm), len(self.molprops_df))
                if len(self.c13ppm) == len(self.molprops_df):
                    self.molprops_df["ppm"] = self.c13ppm["mean"]

                    self.molprops_df["atom_idx"] = self.molprops_df.index
                    # this is wrong we need to use the carbonAtomInfo dataframe!!!!!
                    self.molprops_df["atomNumber"] = carbonAtomsInfo[
                        "atomNumber"
                    ].values
                    self.nmrshiftdb_failed = False

                else:
                    self.nmrshiftdb_failed_message = "<p> An error occurred during the calculation of the chemical shifts using nmrshiftDB  </p><p>- number of atoms in c13ppm does not match number of atoms in th molecule </p>"
                    self.molprops_df["ppm"] = -100000.0
                    self.molprops_df["atom_idx"] = self.molprops_df.index
                    self.molprops_df["atomNumber"] = carbonAtomsInfo[
                        "atomNumber"
                    ].values
                    self.nmrshiftdb_failed = True
                    self.nmrshiftdb_failed_code = 401

            except:
                print(
                    "An error occurred during the calculation of the chemical shifts using nmrshiftDB"
                )
                self.nmrshiftdb_failed_message = "<p> An error occurred during the calculation of the chemical shifts using nmrshiftDB </p>"
                self.molprops_df["ppm"] = -100000.0
                self.molprops_df["atom_idx"] = self.molprops_df.index
                self.molprops_df["atomNumber"] = carbonAtomsInfo["atomNumber"].values
                self.nmrshiftdb_failed = True
                self.nmrshiftdb_failed_code = 401

        # check if there are rings in the molecule
        self.molprops_df["ring_idx"] = -1
        self.molprops_df["ring_size"] = 0
        ring_atoms = self.GetRingSystems()

        # print(f"self.molprops_df.index {self.molprops_df.index}")

        # create sets of carbon atoms in rings store them in molprops_df column aromatic_rings
        for ring_idx, ring in enumerate(ring_atoms):
            carbon_atoms_in_ring = [
                i for i in ring if self.mol.GetAtomWithIdx(i).GetSymbol() == "C"
            ]
            # print(f"ring_idx {ring_idx} carbon_atoms_in_ring {carbon_atoms_in_ring}")
            self.molprops_df.loc[carbon_atoms_in_ring, "ring_idx"] = ring_idx

        # set the ring size for each carbon atom
        for ring in ring_atoms:
            for atom in ring:
                # check if atom is carbon
                if self.mol.GetAtomWithIdx(atom).GetSymbol() == "C":
                    self.molprops_df.loc[atom, "ring_size"] = len(ring)

        # add ring info to molprops_df
        self.add_ring_info_to_dataframe()

        self.has_symmetry = (
            len(self.mol.GetSubstructMatches(self.mol, uniquify=False, maxMatches=3))
            > 1
        )

        # idx_list, xxx, yyy = self.calc_carbon_xy_positions_png(self.mol)
        self.molprops_df["x"] = 0.0
        self.molprops_df["y"] = 0.0
        # self.molprops_df.loc[idx_list, "x"] = xxx
        # self.molprops_df.loc[idx_list, "y"] = yyy

        self.copy_xy3_to_molecule(self.xy3, self.molprops_df)

        self.molprops_df["symmetry_idx1"] = -1
        self.molprops_df["symmetry_idx2"] = -1
        self.molprops_df["sym_atom_idx"] = ""
        self.molprops_df["sym_atomNumber"] = ""

        print(
            "\nfind the symmetry atoms in each ring and non ring group and assign the symmetry indices to each\n"
        )

        # find the symmetry atoms in each ring and non ring group and assign the symmetry indices to each
        for ring_idx in self.molprops_df.ring_idx.unique():
            ring_df = self.molprops_df[self.molprops_df.ring_idx == ring_idx]
            for ppm in ring_df.ppm.unique():
                ring_df_ppm = ring_df[ring_df.ppm == ppm]
                if len(ring_df_ppm) == 2:
                    self.molprops_df.loc[
                        ring_df_ppm.index[0], "symmetry_idx1"
                    ] = ring_df_ppm.index[1]
                    self.molprops_df.loc[
                        ring_df_ppm.index[1], "symmetry_idx1"
                    ] = ring_df_ppm.index[0]
                elif len(ring_df_ppm) == 4 and ring_idx > -1:
                    self.molprops_df.loc[
                        ring_df_ppm.index[0], "symmetry_idx1"
                    ] = ring_df_ppm.index[3]
                    self.molprops_df.loc[
                        ring_df_ppm.index[1], "symmetry_idx1"
                    ] = ring_df_ppm.index[2]
                    self.molprops_df.loc[
                        ring_df_ppm.index[2], "symmetry_idx1"
                    ] = ring_df_ppm.index[1]
                    self.molprops_df.loc[
                        ring_df_ppm.index[3], "symmetry_idx1"
                    ] = ring_df_ppm.index[0]

        # update symmetry indices sym_atom_idx and sym_atomNumber

        self.molprops_df["sym_atom_idx"] = ""
        self.molprops_df["sym_atomNumber"] = ""

        for idx, row in self.molprops_df.iterrows():
            df = self.molprops_df[self.molprops_df.ppm == row["ppm"]]

            if len(df) > 1:
                # remove row from df that has the same index as the current row
                df = df.drop(idx)

                for idx2, row2 in df.iterrows():
                    if self.molprops_df.at[idx2, "sym_atom_idx"] == "":
                        self.molprops_df.at[idx2, "sym_atom_idx"] = f"{row['atom_idx']}"
                        self.molprops_df.at[
                            idx2, "sym_atomNumber"
                        ] = f"{row['atomNumber']}"
                    else:
                        self.molprops_df.at[
                            idx2, "sym_atom_idx"
                        ] = f"{self.molprops_df.at[idx2, 'sym_atom_idx']}, {row['atom_idx']}"
                        self.molprops_df.at[
                            idx2, "sym_atomNumber"
                        ] = f"{self.molprops_df.at[idx2, 'sym_atomNumber']}, {row['atomNumber']}"

        # calculate number of carbons without protons attached
        self.num_quaternary_carbons = self.molprops_df[
            self.molprops_df.quaternary
        ].shape[0]
        self.num_CH0_carbon_atoms = self.molprops_df[self.molprops_df.quaternary].shape[
            0
        ]

        # calculate number of carbon with two protons attached
        self.num_CH2_carbon_atoms = self.molprops_df[self.molprops_df.CH2].shape[0]

        # calculate number of carbon with three protons attached
        self.num_CH3_carbon_atoms = self.molprops_df[self.molprops_df.CH3].shape[0]

        # calculate number of carbon with one proton  attached
        self.num_CH1_carbon_atoms = self.molprops_df[self.molprops_df.CH1].shape[0]

        # calculate number of carbons with protons attached
        self.num_carbon_atoms_with_protons = (
            self.num_CH2_carbon_atoms
            + self.num_CH3_carbon_atoms
            + self.num_CH1_carbon_atoms
        )

        self.num_hsqc_carbon_atoms = self.num_carbon_atoms_with_protons

        # check if there are aromatic rings that map symmetrically to each other
        self.mapped_symmetric_aromatic_rings = self.map_symmetric_aromatic_rings(
            self.molprops_df
        )

        # create reduced  molprops_df based on symmetry of NMR ppm values
        self.sym_molprops_df = self.molprops_df.drop_duplicates(
            subset=["ppm", "aromatic", "numProtons"], inplace=False
        )

        # self.copy_xy3_to_molecule(self.xy3, self.molprops_df)

        # molecule has hose-code symmetry if there are less rows in sym_molprops_df than in molprops_df
        self.has_hose_code_symmetry = (
            self.sym_molprops_df.shape[0] < self.molprops_df.shape[0]
        )

        # calculate number of aromatic carbon atoms in sym_molprops_df
        self.num_sym_aromatic_carbon_atoms = self.sym_molprops_df[
            self.sym_molprops_df.aromatic
        ].shape[0]
        self.num_sym_CH0_carbon_atoms = self.sym_molprops_df[
            self.sym_molprops_df.CH0
        ].shape[0]

        # calculate the number of quaternary carbon atoms in sym_molprops_df
        self.num_sym_quaternary_carbons = self.sym_molprops_df[
            self.sym_molprops_df.quaternary
        ].shape[0]

        # calculate the number of carbon atoms with two protons attached in sym_molprops_df
        self.num_sym_CH2_carbon_atoms = self.sym_molprops_df[
            self.sym_molprops_df.CH2
        ].shape[0]

        # calculate the number of carbon atoms with three protons attached in sym_molprops_df
        self.num_sym_CH3_carbon_atoms = self.sym_molprops_df[
            self.sym_molprops_df.CH3
        ].shape[0]

        # calculate the number of carbon atoms with one proton attached in sym_molprops_df
        self.num_sym_CH1_carbon_atoms = self.sym_molprops_df[
            self.sym_molprops_df.CH1
        ].shape[0]

        # calculate the number of carbon atoms with protons attached in sym_molprops_df
        self.num_sym_carbon_atoms_with_protons = (
            self.num_sym_CH2_carbon_atoms
            + self.num_sym_CH3_carbon_atoms
            + self.num_sym_CH1_carbon_atoms
        )

        self.num_sym_hsqc_carbon_atoms = self.num_sym_carbon_atoms_with_protons

        # calculate number of carbon atoms in sym_molprops_df
        self.num_sym_carbon_atoms = self.sym_molprops_df.shape[0]

        # calculate number of aromatic rings
        self.num_aromatic_rings = rdkit.Chem.rdMolDescriptors.CalcNumAromaticRings(
            self.mol
        )

        # calculate number of all rings
        self.num_rings = rdkit.Chem.rdMolDescriptors.CalcNumRings(self.mol)

        # create a list of pairs of symmetry atoms if symmetry_idx1 is not -1
        self.symmetry_pairs = []
        for idx in self.molprops_df.index:
            if self.molprops_df.loc[idx, "symmetry_idx1"] != -1:
                pair = {idx, self.molprops_df.loc[idx, "symmetry_idx1"]}
                # sort pair so that the lower index is first
                # pair.sort()
                # add pair to list if it is not already in the list
                if pair not in self.symmetry_pairs:
                    self.symmetry_pairs.append(pair)

    def map_symmetric_aromatic_rings(self, df):
        # if there are no aromatic rings then return an empty list
        if df.ring_idx.max() <= 0:
            return []

        # if there are more than one aromatic ring then map the aromatic rings to the symmetry pairs
        # first get the symmetry pairs
        symmetry_pairs = []
        for x, y in np.array(np.triu_indices(df.ring_idx.max() + 1, k=1)).T:
            s1 = df.query("ring_idx == @x")["ppm"].to_list()
            s2 = df.query("ring_idx == @y")["ppm"].to_list()
            s1.sort()
            s2.sort()
            if s1 == s2:
                symmetry_pairs.append((x, y))

        return symmetry_pairs

    def GetRingSystems(self, includeSpiro=False):
        ri = self.mol.GetRingInfo()
        systems = []
        for ring in ri.AtomRings():
            ringAts = set(ring)
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system))
                if nInCommon and (includeSpiro or nInCommon > 1):
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
            nSystems.append(ringAts)
            systems = nSystems
        return systems

    def add_ring_info_to_dataframe(self):
        ring_info = self.mol.GetRingInfo()

        for i, ring in enumerate(ring_info.AtomRings()):
            if f"ring_idx{i+1}" not in self.molprops_df.columns:
                self.molprops_df[f"ring_idx{i}"] = -1
                self.molprops_df[f"ring_size{i}"] = 0
            for atom_idx in ring:
                # check if atom is a carbon atom
                if atom_idx in self.molprops_df.index:

                    self.molprops_df.loc[atom_idx, f"ring_idx{i}"] = i
                    self.molprops_df.loc[atom_idx, f"ring_size{i}"] = len(ring)
                else:
                    print(f"atom_idx {atom_idx} not in mol.molprops_df.index")

    def _repr_png_(self):
        return self.mol._repr_png_()

    def _repr_svg_(self):
        return self.mol._repr_svg_()

    def __repr__(self):
        return self.mol.__repr__()

    def __str__(self):
        return self.mol.__str__()

    def GetAtoms(self):
        return self.mol.GetAtoms()

    def GetBonds(self):
        return self.mol.GetBonds()

    def GetAtomWithIdx(self, idx):
        return self.mol.GetAtomWithIdx(idx)

    def GetBondWithIdx(self, idx):
        return self.mol.GetBondWithIdx(idx)

    def GetNumAtoms(self):
        return self.mol.GetNumAtoms()

    def GetNumBonds(self):
        return self.mol.GetNumBonds()

    def GetAromaticAtoms(self):
        return self.mol.GetAromaticAtoms()

    def GetBondBetweenAtoms(self, i, j):
        return self.mol.GetBondBetweenAtoms(i, j)

    def GetRingInfo(self):
        return self.mol.GetRingInfo()

    def GetConformer(self):
        return self.mol.GetConformer()

    def init_elements_dict(self):
        return (
            pd.DataFrame(
                [
                    [atom.GetIdx(), atom.GetSymbol()]
                    for atom in Chem.AddHs(self.mol).GetAtoms()
                ],
                columns=["atom_index", "atom_symbol"],
            )["atom_symbol"]
            .value_counts()
            .to_dict()
        )

    # calculate DBE for molecule
    def calc_dbe(self) -> int:
        elements = self.init_elements_dict()
        if "C" in elements:
            dbe_value = elements["C"]
        if "N" in elements:
            dbe_value += elements["N"] / 2
        for e in ["H", "F", "Cl", "Br"]:
            if e in elements:
                dbe_value -= elements[e] / 2

        return dbe_value + 1

    # def calculated_c13_chemical_shifts(self) -> pd.DataFrame:
    #     return self.calc_c13_chemical_shifts_using_nmrshift2D()

    def calculated_c13_chemical_shifts(self) -> pd.DataFrame:
        return calc_c13_chemical_shifts_using_nmrshift2D(self.mol_str)

    # return dictionary of dictionarys, first key is the number of protons, second key is carbon atom index, value is calculated C13 NMR chemical shift for molecule
    def c13_nmr_shifts(self) -> dict:
        c13_nmr_shifts = {
            k: {i: None for i in v} for k, v in self.proton_groups().items()
        }

        c13ppm_df = calc_c13_chemical_shifts_using_nmrshift2D(self.mol_str)
        if isinstance(c13ppm_df, pd.DataFrame):
            # reset index to atom index
            for v in c13_nmr_shifts.values():
                for k2, v2 in v.items():
                    # find row in c13ppm_df with atom index k2
                    v[k2] = c13ppm_df.loc[k2, "mean"]

        return c13_nmr_shifts

    # @cache
    # def calc_c13_chemical_shifts_using_nmrshift2D(self) -> pd.DataFrame:
    #     print("*****************************************")
    #     print("calc_c13_chemical_shifts_using_nmrshift2D")
    #     print("*****************************************")

    #     with open("mol.mol", "w") as fp:
    #         fp.write(Chem.MolToMolBlock(self.mol))

    #     ret = os.system(self.java_command)

    #     if ret == 1:
    #         print("NMRShift2D failed to calculate C13 chemical shifts")
    #         return False
    #     else:
    #         mol_df = pd.read_csv("mol.csv", index_col=0)
    #         mol_df.index = mol_df.index - 1
    #         return mol_df

    # return dictionary of lists key is the number of protons attached to carbon, value is list of carbon atom indices
    def proton_groups(self) -> dict:
        proton_groups = {
            atom.GetTotalNumHs(): []
            for atom in self.mol.GetAtoms()
            if atom.GetAtomicNum() == 6
        }
        for atom in self.mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                proton_groups[atom.GetTotalNumHs()].append(atom.GetIdx())
        return proton_groups

    def create_png(self):
        """Creates a png image from a smiles string via rdkit"""
        png = None

        # mol2 = Chem.AddHs(self.mol)
        # AllChem.EmbedMolecule(mol2, randomSeed=3)
        # rdkit_molecule = Chem.RemoveHs(mol2)

        # rdkit_molecule.Compute2DCoords()

        return Draw.MolToImage(rdkit_molecule, size=(XYDIM, XYDIM))

    # def create_png_from_smiles(smiles_str: str) -> PIL.Image.Image:
    #     """Creates a png image from a smiles string via rdkit"""
    #     png = None
    #     # rdkit_molecule = Chem.MolFromSmiles(smiles_str)
    #     rdkit_molecule = expectedmolecule.expectedMolecule(smiles_str)

    #     mol2 = Chem.AddHs(rdkit_molecule)
    #     AllChem.EmbedMolecule(mol2, randomSeed=3)
    #     rdkit_molecule = Chem.RemoveHs(mol2)

    #     rdkit_molecule.Compute2DCoords()

    #     return Draw.MolToImage(rdkit_molecule, size=(XYDIM, XYDIM))

    def calc_carbon_xy_positions_png(self, rdkit_molecule: Chem.Mol) -> list:
        """Returns the xy3 positions of the carbon atoms in the molecule"""

        d2d = Draw.rdMolDraw2D.MolDraw2DSVG(XYDIM, XYDIM)
        d2d.DrawMolecule(rdkit_molecule)
        d2d.FinishDrawing()
        idx_list = []
        xxx = []
        yyy = []
        for atom in rdkit_molecule.GetAtoms():
            if atom.GetSymbol() == "C":
                idx = atom.GetIdx()
                point = d2d.GetDrawCoords(idx)
                idx_list.append(idx)
                xxx.append(point.x / XYDIM)
                yyy.append(point.y / XYDIM)
        return idx_list, xxx, yyy

    def calc_allatom_xy_positions_png(self) -> list:
        """Returns the xy3 positions of all heavy atoms in the molecule"""

        d2d = Draw.rdMolDraw2D.MolDraw2DSVG(XYDIM, XYDIM)
        d2d.DrawMolecule(self.mol)
        d2d.FinishDrawing()
        idx_list = []
        xxx = []
        yyy = []
        for atom in self.mol.GetAtoms():

            idx = atom.GetIdx()
            point = d2d.GetDrawCoords(idx)
            idx_list.append(idx)
            xxx.append(point.x / XYDIM)
            yyy.append(point.y / XYDIM)
        return idx_list, xxx, yyy

    # def create_svg_string(self,
    #     mol, molWidth=1000, molHeight=600, svgWidth=1200, svgHeight=700
    # ):

    #     translateWidth = int((svgWidth - molWidth) / 2)
    #     translateHeight = int((svgHeight - molHeight) / 2)

    #     d2d = Draw.rdMolDraw2D.MolDraw2DSVG(molWidth, molHeight)
    #     d2d.DrawMolecule(mol)
    #     d2d.TagAtoms(mol)
    #     d2d.FinishDrawing()
    #     sss = d2d.GetDrawingText().replace(
    #         f"width='{molWidth}px' height='{molHeight}px'",
    #         f"width={molWidth} height={molHeight}",
    #     )
    #     sss = sss.replace("fill:#FFFFFF", "fill:none").replace(
    #         "<svg", '<svg class="center"'
    #     )

    #     sss = sss.replace(
    #         f"<!-- END OF HEADER -->",
    #         f"<!-- END OF HEADER -->\n<g transform='translate({translateWidth}, {translateHeight})'>",
    #     )
    #     sss = sss.replace("</svg>", "</g>\n</svg>")
    #     sss = sss.replace(
    #         f"width={molWidth} height={molHeight} viewBox='0 0 {molWidth} {molHeight}'",
    #         f"width={svgWidth} height={svgHeight} viewBox='0 0 {svgWidth} {svgHeight}'",
    #     )

    #     idx_list = []
    #     xxx = []
    #     yyy = []

    #     new_xy3 = {}
    #     for atom in mol.GetAtoms():
    #         if atom.GetSymbol() == "C":
    #             idx = atom.GetIdx()
    #             point = d2d.GetDrawCoords(idx)
    #             idx_list.append(idx)
    #             xxx.append(point.x / molWidth)
    #             yyy.append(point.y / molHeight)
    #             new_xy3[idx] = (point.x / molWidth, point.y / molHeight)

    #     return sss, new_xy3

    def create_svg_string(
        self, molWidth=1000, molHeight=600, svgWidth=1200, svgHeight=700
    ):

        translateWidth = int((svgWidth - molWidth) / 2)
        translateHeight = int((svgHeight - molHeight) / 2)

        # self.mol = Chem.AddHs(self.mol)
        # AllChem.EmbedMolecule(self.mol)
        # self.mol = Chem.RemoveHs(self.mol)

        # d2d = Draw.rdMolDraw2D.MolDraw2DSVG(molWidth, molHeight)
        # d2d.drawOptions().minFontSize = 15
        # d2d.DrawMolecule(self.mol)
        # d2d.TagAtoms(self.mol)
        # d2d.FinishDrawing()

        # fixed_coords = {}
        # for atom in self.mol.GetAtoms():
        #     idx = atom.GetIdx()
        #     point = d2d.GetDrawCoords(idx)
        #     fixed_coords[idx] = point

        # AllChem.Compute2DCoords(self.mol)

        d2d = Draw.rdMolDraw2D.MolDraw2DSVG(molWidth, molHeight)
        d2d.drawOptions().minFontSize = 20
        # d2d.drawOptions().multipleBondOffset = 1
        d2d.DrawMolecule(self.mol)
        d2d.TagAtoms(self.mol)
        d2d.FinishDrawing()

        sss = d2d.GetDrawingText()
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
            f"width={svgWidth} height={svgHeight} viewBox='0 0 {svgWidth} {svgHeight}'",
        )

        idx_list = []
        xxx = []
        yyy = []

        new_xy3 = {}
        new_xy3_all = {}
        # for atom in self.mol.GetAtoms():
        #     if atom.GetSymbol() == "C":
        #         idx = atom.GetIdx()
        #         point = d2d.GetDrawCoords(idx)
        #         idx_list.append(idx)
        #         xxx.append(point.x / molWidth)
        #         yyy.append(point.y / molHeight)
        #         new_xy3[idx] = (point.x / molWidth, point.y / molHeight)

        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            point = d2d.GetDrawCoords(idx)
            idx_list.append(idx)
            xxx.append(point.x / molWidth)
            yyy.append(point.y / molHeight)

            new_xy3_all[idx] = (point.x / molWidth, point.y / molHeight)

            if atom.GetSymbol() == "C":
                new_xy3[idx] = (point.x / molWidth, point.y / molHeight)

        return sss, new_xy3, new_xy3_all


# if __name__ == "__main__":

#     from matplotlib import pyplot as plt

#     # mol = expectedMolecule("COc1cc(/C=C/C=O)cc(OC)c1O")
#     # mol = expectedMolecule("CO[C@H]1[C@@H](O)[C@@H](C)O[C@H](O[C@@H]2[C@@H](O)[C@@H](O)[C@@H](C)O[C@H]2OC2=CC=CC3=C(O)C4=C5C(OC(=O)C6=C5C(OC4=O)=CC=C6C)=C23)[C@@H]1O")
#     # mol = expectedMolecule("CCCC1CC(N(C1)C)C(=O)NC(C2C(C(C(C(O2)SC)OP(=O)(O)O)O)O)C(C)Cl")
#     # mol = expectedMolecule("CC(=O)OC(C)(C)\C=C\C(=O)C(C)(O)C1C(O)CC2(C)C3CC=C4C(C=C(O)C(=O)C4(C)C)C3(C)C(=O)CC12C")
#     # mol = expectedMolecule("CNC(=O)OC1=CC=C2N(C)C3N(C)CCC3(C)C2=C1")
#     # mol = expectedMolecule("CC1=CC(=CC(=C1C2=CC=CC(=C2)COC3=CC4=C(C=C3)[C@@H](CO4)CC(=O)O)C)OCCCS(=O)(=O)C")
#     # mol = expectedMolecule("c1cc2ccc1CC2")
#     # mol = expectedMolecule("c1cc2cc(c1)CC2")

#     molstr = """
# RDKit          2D

#  20 22  0  0  0  0  0  0  0  0999 V2000
#     2.6401   -3.1891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#     1.2135   -2.7256    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
#     0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#     1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#     0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.5000    2.5981    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.7500    3.8971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.5000    5.1962    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.7500    6.4952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#     0.7500    3.8971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.2135   -2.7256    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.4271   -1.8439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.4271   -3.6073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.9635   -5.0339    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.4635   -5.0339    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
#     0.4182   -6.2474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.0000   -3.6073    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   1  2  1  0
#   2  3  1  0
#   3  4  2  0
#   4  5  1  0
#   5  6  2  0
#   6  7  1  0
#   7  8  1  0
#   8  9  1  0
#   9 10  1  0
#   8 11  2  0
#   6 12  1  0
#  12 13  2  0
#  13 14  1  0
#  14 15  1  0
#  14 16  1  0
#  16 17  1  0
#  17 18  1  0
#  18 19  1  0
#  18 20  1  0
#  20  2  1  0
#  13  3  1  0
#  20 14  1  0
# M  END
# """

#     mol = expectedMolecule(molstr, is_smiles=False)

#     print(mol.molprops_df)

#     plt.imshow(
#         mol.png,
#         aspect="auto",
#         extent=[0, 1, 1, 0],
#     )

#     xxx = mol.molprops_df["x"]
#     yyy = mol.molprops_df["y"]
#     plt.scatter(xxx, yyy, c="red", s=50)

#     plt.xlim(-0.1, 1.1)
#     plt.ylim(1.1, -0.1)
#     plt.show()
