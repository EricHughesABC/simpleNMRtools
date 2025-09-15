# import collections
# import types

# global svgDimensions

# SVGDIMENSIONS = collections.namedtuple(
#     "SVGDIMENSIONS", "MOLWIDTH MOLHEIGHT SVGWIDTH SVGHEIGHT"
# )
# svgDimensions = SVGDIMENSIONS(1000, 600, 1000, 600)


# global CARBONSEPARATION
# global PROTONSEPARATION

# CARBONSEPARATION = 0.0025
# PROTONSEPARATION = 0.005

# global NMREXPERIMENTS

# NMREXPERIMENTS = (
#     "H1_1D",
#     "C13_1D",
#     "DEPT135",
#     "HSQC",
#     "HMBC",
#     "COSY",
#     "NOESY",
#     "H1_pureshift",
#     "HSQC_CLIPCOSY",
#     "DDEPT_CH3_ONLY",
# )



# NODECOLORMAP = types.MappingProxyType({
#     0: "#FFA500",
#     1: "#98FB98",
#     2: "yellow", 
#     3: "#00FFFF",
#     -1: "lightblue",  # For non-carbon atoms
#     -2: "lightgrey",
# })

# print(svgDimensions.MOLWIDTH)
# print(svgDimensions.MOLHEIGHT)
# print(svgDimensions.SVGWIDTH)
# print(svgDimensions.SVGHEIGHT)
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

# SVG dimensions configuration
SVGDimensions = collections.namedtuple(
    "SVGDimensions", ["mol_width", "mol_height", "svg_width", "svg_height"]
)

SVG_DIMENSIONS = SVGDimensions(
    mol_width=1000,
    mol_height=600, 
    svg_width=1000,
    svg_height=600
)

# Chemical separation constants
CARBONSEPARATION = 0.0025
PROTONSEPARATION = 0.005

# NMR experiment types
NMREXPERIMENTS = types.MappingProxyType({
    "H1_1D": "H1_1D",
    "C13_1D": "C13_1D", 
    "DEPT135": "DEPT135",
    "HSQC": "HSQC",
    "HMBC": "HMBC",
    "COSY": "COSY",
    "NOESY": "NOESY",
    "H1_pureshift": "H1_pureshift",
    "HSQC_CLIPCOSY": "HSQC_CLIPCOSY",
    "DDEPT_CH3_ONLY": "DDEPT_CH3_ONLY",
})

# Node color mapping for visualization
NODE_COLOR_MAP = types.MappingProxyType({
    0: "#FFA500",    # Orange
    1: "#98FB98",    # Light green
    2: "#FFFF00",    # Yellow
    3: "#00FFFF",    # Cyan
    -1: "#ADD8E6",   # Light blue (non-carbon atoms)
    -2: "#D3D3D3",   # Light grey
})

# Example usage (remove in production)
if __name__ == "__main__":
    print(f"Molecule dimensions: {SVG_DIMENSIONS.mol_width}x{SVG_DIMENSIONS.mol_height}")
    print(f"SVG dimensions: {SVG_DIMENSIONS.svg_width}x{SVG_DIMENSIONS.svg_height}")
    print(f"Available NMR experiments: {len(NMREXPERIMENTS)}")