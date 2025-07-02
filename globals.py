import collections

global svgDimensions

SVGDIMENSIONS = collections.namedtuple(
    "SVGDIMENSIONS", "MOLWIDTH MOLHEIGHT SVGWIDTH SVGHEIGHT"
)
svgDimensions = SVGDIMENSIONS(1000, 600, 1000, 600)


global CARBONSEPARATION
global PROTONSEPARATION

CARBONSEPARATION = 0.00075
PROTONSEPARATION = 0.0025

global NMREXPERIMENTS

NMREXPERIMENTS = (
    "H1_1D",
    "C13_1D",
    "DEPT135",
    "HSQC",
    "HMBC",
    "COSY",
    "NOESY",
    "H1_pureshift",
    "HSQC_CLIPCOSY",
    "DDEPT_CH3_ONLY",
)

print(svgDimensions.MOLWIDTH)
print(svgDimensions.MOLHEIGHT)
print(svgDimensions.SVGWIDTH)
print(svgDimensions.SVGHEIGHT)
