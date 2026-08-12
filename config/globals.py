"""
Global configuration constants for the application.

Usage Examples:
--------------

# Import specific constants
from globals import SVG_DIMENSIONS, CARBONSEPARATION, NODE_COLOR_MAP, NMREXPERIMENTS

# Use SVG dimensions
canvas_width = SVG_DIMENSIONS.svg_width
molecule_height = SVG_DIMENSIONS.mol_height

# Use chemical constants
spacing = CARBONSEPARATION * 2
proton_gap = PROTONSEPARATION

# Use color mapping
atom_color = NODE_COLOR_MAP[1]  # Light green
carbon_color = NODE_COLOR_MAP[0]  # Orange

# Use NMR experiments
if "HSQC" in NMREXPERIMENTS:
    experiment_type = NMREXPERIMENTS["HSQC"]

# Alternative: Import entire module
import globals as config
width = config.SVG_DIMENSIONS.svg_width
color = config.NODE_COLOR_MAP[-1]

# Alternative: Import with aliases
from globals import SVG_DIMENSIONS as dims, NODE_COLOR_MAP as colors
size = dims.mol_width
node_color = colors[2]
"""

import collections
import types

# globals.py
import sys

class _ConstModule(types.ModuleType):
    def __setattr__(self, name, value):
        if name in self.__dict__:
            raise AttributeError(f"Cannot reassign constant '{name}'")
        super().__setattr__(name, value)

XYDIM = 800
CH1 = "CH1"
CH3 = "CH3"



# SVG dimensions configuration
SVGDimensions = collections.namedtuple(
    "SVGDimensions", ["mol_width", "mol_height", "svg_width", "svg_height"]
)

SVG_DIMENSIONS = SVGDimensions(
    mol_width=1000, mol_height=600, svg_width=1000, svg_height=600
)

# Global constant for XY dimension used in visualizations of molecule PNG image
XYDIM = 800

# Chemical separation constants
CARBONSEPARATION = 0.0025
PROTONSEPARATION = 0.005

IDX = "idx"
IMPLICITHS = "implicitHs"
TOTALNUMHS = "totalNumHs"
DEGREE = "degree"
HYBRIDIZATION = "hybridization"
AROMATIC = "aromatic"

CH1 = "CH1"
CH2 = "CH2"
CH3 = "CH3"
CH0 = "CH0"
CH3CH1 = "CH3CH1"  
CH3plusCH1 = "CH3+CH1"
QUATERNARY_CARBON = "quaternary"  # Quaternary carbon (no attached hydrogens)
NUMPROTONS = "numProtons"  # Number of attached protons (0, 1, 2, 3)
ATTACHED_PROTONS = "attached_protons"  # Number of attached protons (0, 1, 2, 3)
INTEGRAL = "integral"  # Integral value of the signal
F2_INTEGRAL = "f2_integral"  # Integral value of the signal in F2 dimension
INTENSITY = "intensity"  # Intensity of the signal
PPM = "ppm"  # Chemical shift in parts per million (ppm)
PPM_CALCULATED = "ppm_calculated"
F1_PPM = "f1_ppm"  # Chemical shift in F1 dimension (ppm)
F2_PPM = "f2_ppm"  # Chemical shift in F2 dimension (ppm)
F2P_PPM = "f2p_ppm"  # Chemical shift in F2 dimension (ppm) for HSQC one-to-one mapping
F1P_PPM = "f1p_ppm"  # Chemical shift in F1 dimension (ppm) for HMBC one-to-one mapping
J_COUPLING_CLASS = "jCouplingClass"  # J-coupling class (e.g., doublet, triplet)
J_COUPLING_VALS = "jCouplingVals"  # J-coupling
RANGE = "range"  # Range of the signal (start and end ppm)
SIGNALTYPE = "signaltype"  # Type of signal (e.g., multiplet, singlet)
COMPOUND = "Compound"  # Compound signal type


ATOMNUMBER = "atomNumber"  # Atom number in the molecule
SYM_ATOMNUMBER = "sym_atomNumber"  # Symmetric atom number in the molecule
ATOMIDX = "atom_idx"  # Atom index in the molecule
SYM_ATOMIDX = "sym_atom_idx"  # Symmetric atom index in the molecule

F1_ATOMIDX = "f1_atom_idx"  # Atom index in the F1 dimension
F1_SYM_ATOMIDX = "f1_sym_atom_idx"  # Symmetric atom index in the F1 dimension
F1_ATOMNUMBER = "f1_atomNumber"  # Atom number in the F1 dimension
F1_SYM_ATOMNUMBER = "f1_sym_atomNumber"  # Symmetric atom number in the F1 dimension
F2_ATOMIDX = "f2_atom_idx"  # Atom index in the F2 dimension
F2_SYM_ATOMIDX = "f2_sym_atom_idx"  # Symmetric atom index in the F2 dimension
F2_ATOMNUMBER = "f2_atomNumber"  # Atom number in the F2 dimension
F2_SYM_ATOMNUMBER = "f2_sym_atomNumber"  # Symmetric atom number in the F2 dimension


X = "x"  # X-coordinate of the atom in the molecule
Y = "y"  # Y-coordinate of the atom in the molecule
F1_X = "f1_x"  # X-coordinate of the atom in the F1 dimension
F1_Y = "f1_y"  # Y-coordinate of the atom in the F1 dimension
F2_X = "f2_x"  # X-coordinate of the atom in the F2 dimension
F2_Y = "f2_y"  # Y-coordinate of the atom in the F2 dimension

IUPAC_LABEL = "iupacLabel"  # IUPAC label for the atom
H1_PPM = "H1_ppm"  # Proton chemical shift in parts per million (ppm)
# HTML CODES
BADREQUEST = 400
NOTFOUND = 404
GOODREQUEST = 200

DATATYPE = "datatype"
MULTIPLETS = "multiplets"
NMRSPECTRUM = "nmrspectrum"
PEAKS = "peaks"
INTEGRALS = "integrals"
COUNT = "count"
DATA = "data"


DELTA1 = "delta1"
DELTA2 = "delta2"
TYPE = "type"
ANNOTATION = "annotation"





# CH3/CH1 split heuristic boundary (ppm)
# CH1 groups (e.g. C-O, C-N) typically resonate at or above this value;
# CH3 groups typically fall below it.  Reliable for >99% of encountered
# molecules.  Adjust here if a specific compound class requires a different
# boundary — do not hard-code 67 elsewhere.
CH1_CH3_PPM_BOUNDARY: int = 67

# NMR experiment types
NMREXPERIMENTS = types.MappingProxyType(
    {
        "H1_1D": "H1_1D",
        "C13_1D": "C13_1D",
        "DEPT135": "DEPT135",
        "HSQC": "HSQC",
        "HMBC": "HMBC",
        "COSY": "COSY",
        "NOESY": "NOESY",
        "H1_pureshift": "H1_pureshift",
        "HSQC_CLIPCOSY": "HSQC_CLIPCOSY",
        "DDEPTCH3ONLY": "DDEPTCH3ONLY",
    }
)

# Node color mapping for visualization
NODE_COLOR_MAP = types.MappingProxyType(
    {
        0: "#FFA500",  # Orange
        1: "#98FB98",  # Light green
        2: "#FFFF00",  # Yellow
        3: "#00FFFF",  # Cyan
        -1: "#ADD8E6",  # Light blue (non-carbon atoms)
        -2: "#D3D3D3",  # Light grey
    }
)

# Define a global constant for column renaming
NMREXPERIMENTS_COLUMN_RENAME_MAP = types.MappingProxyType(
    {
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
    }
)

# Simulated annealing parameters
SIMULATED_ANNEALING_RUNS = 50
SIMULATED_ANNEALING_COOLING_RATE = 0.995

COSY = "cosy"

sys.modules[__name__].__class__ = _ConstModule

# Example usage (remove in production)
if __name__ == "__main__":
    print(
        f"Molecule dimensions: {SVG_DIMENSIONS.mol_width}x{SVG_DIMENSIONS.mol_height}"
    )
    print(f"SVG dimensions: {SVG_DIMENSIONS.svg_width}x{SVG_DIMENSIONS.svg_height}")
    print(f"Available NMR experiments: {len(NMREXPERIMENTS)}")


