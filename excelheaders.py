"""
NMR DataFrame Column Definitions

This module defines the column structures for pandas DataFrames used in NMR data processing.
Contains both processed and original column mappings for various NMR experiment types.
"""

# Processed DataFrame Columns
# These are the standardized column names used in the processed data pipeline
EXCEL_DF_COLUMNS = {
    "nmrAssignments": [
        "atom_idx",
        "atomNumber", 
        "numProtons",
        "f1_ppm",
        "iupacLabel",
        "f2_ppm1",
        "f2_ppm2",
        "f2_ppm3",
    ],
    
    "c13predictions": [
        "atomNumber",
        "atom_idx", 
        "numProtons",
        "ppm"
    ],
    
    "h1": [
        "ppm",
        "integral",
        "jCouplingClass",
        "jCouplingVals",
        "range",
        "label",
        "f1H_i",
        "f2H_i",
        "f1_i",
        "f2_i",
    ],
    
    "c13": [
        "ppm",
        "attached_protons",
        "ppmH1s",
        "max_bonds",
        "label",
        "f2C_i",
        "f1C_i",
        "f1_i",
        "f2_i",
    ],
    
    "DEPT135": [
        "ppm",
        "attached_protons",
        "ppmH1s",
        "max_bonds",
        "label",
        "f2C_i",
        "f1C_i",
        "f1_i",
        "f2_i",
    ],
    
    "pureshift": [
        "ppm",
        "intensity",
        "Width",
        "Area",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    # 2D NMR Experiments - Common column structure
    "cosy": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f1_i",
        "f1p_i",
        "f2_i",
        "f2p_i",
        "f1p_ppm",
        "f2p_ppm",
        "f1H_i",
        "f2H_i",
        "f1Cp_i",
        "f2Cp_i",
    ],
    
    "noesy": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f1_i",
        "f1p_i",
        "f2_i",
        "f2p_i",
        "f1p_ppm",
        "f2p_ppm",
        "f1H_i",
        "f2H_i",
        "f1Cp_i",
        "f2Cp_i",
    ],
    
    "hsqc": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f1_i",
        "f2_i",
        "f2p_i",
        "f1C_i",
        "f2H_i",
        "f2Cp_i",
        "f2p_ppm",
    ],
    
    "hsqc_clip_cosy": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f1_i",
        "f2_i",
        "f2p_i",
        "f1C_i",
        "f2H_i",
        "f2Cp_i",
        "f2p_ppm",
    ],
    
    "ddept_ch3_only": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f1_i",
        "f2_i",
        "f2p_i",
        "f1C_i",
        "f2H_i",
        "f2Cp_i",
        "f2p_ppm",
    ],
    
    "hmbc": [
        "f1_ppm",
        "f2_ppm",
        "intensity",
        "f2p_ppm",
        "f1_i",
        "f2_i",
        "f2p_i",
        "f1C_i",
        "f2H_i",
        "f2Cp_i",
    ],
}


# Original DataFrame Columns  
# These represent the raw column names from Excel files before processing
EXCEL_ORIG_DF_COLUMNS = {
    "nmrAssignments": [
        "atom_idx",
        "atomNumber",
        "numProtons",
        "f1_ppm",
        "iupacLabel",
        "f2_ppm1",
        "f2_ppm2",
        "f2_ppm3",
    ],
    
    "c13predictions": [
        "atomNumber",
        "atom_idx",
        "numProtons", 
        "ppm"
    ],
    
    "molecule": [
        "molecule",
        "smiles"
    ],
    
    # 1D NMR Experiments
    "H1_1D": [
        "Name",
        "Shift",
        "Range",
        "H's",
        "Integral",
        "Class",
        "J's",
        "Method"
    ],
    
    "C13_1D": [
        "ppm",
        "Intensity",
        "Width",
        "Area",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    "DEPT135": [
        "ppm",
        "Intensity",
        "Width",
        "Area",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    "H1_pureshift": [
        "ppm",
        "Intensity",
        "Width",
        "Area",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    # 2D NMR Experiments - Original Excel format
    "COSY": [
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    "HSQC": [
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    "HSQC_CLIPCOSY": [
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    "DDEPT_CH3_ONLY": [
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    "HMBC": [
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
    
    "NOESY": [
        "f2 (ppm)",
        "f1 (ppm)",
        "Intensity",
        "Width f2",
        "Width f1",
        "Volume",
        "Type",
        "Flags",
        "Impurity/Compound",
        "Annotation",
    ],
}


# # Helper functions for working with column definitions
# def get_processed_columns(experiment_type: str) -> list[str]:
#     """
#     Get the processed column names for a specific NMR experiment type.
    
#     Args:
#         experiment_type: The NMR experiment type key
        
#     Returns:
#         List of column names for the processed DataFrame
        
#     Raises:
#         KeyError: If experiment_type is not found in EXCEL_DF_COLUMNS
#     """
#     return EXCEL_DF_COLUMNS[experiment_type]


# def get_original_columns(experiment_type: str) -> list[str]:
#     """
#     Get the original Excel column names for a specific NMR experiment type.
    
#     Args:
#         experiment_type: The NMR experiment type key
        
#     Returns:
#         List of column names from the original Excel format
        
#     Raises:
#         KeyError: If experiment_type is not found in EXCEL_ORIG_DF_COLUMNS
#     """
#     return EXCEL_ORIG_DF_COLUMNS[experiment_type]


# def get_available_experiments() -> tuple[list[str], list[str]]:
#     """
#     Get lists of available experiment types for both processed and original formats.
    
#     Returns:
#         Tuple containing (processed_experiments, original_experiments)
#     """
#     processed_experiments = list(EXCEL_DF_COLUMNS.keys())
#     original_experiments = list(EXCEL_ORIG_DF_COLUMNS.keys())
#     return processed_experiments, original_experiments

# excel_df_columns = {

#     "nmrAssignments": [
#         "atom_idx",
#         "atomNumber",
#         "numProtons",
#         "f1_ppm",
#         "iupacLabel",
#         "f2_ppm1",
#         "f2_ppm2",
#         "f2_ppm3",
#     ],
#     "c13predictions": ["atomNumber", "atom_idx", "numProtons", "ppm"],
#     "h1": [
#         "ppm",
#         "integral",
#         "jCouplingClass",
#         "jCouplingVals",
#         "range",
#         "label",
#         "f1H_i",
#         "f2H_i",
#         "f1_i",
#         "f2_i",
#     ],
#     "c13": [
#         "ppm",
#         "attached_protons",
#         "ppmH1s",
#         "max_bonds",
#         "label",
#         "f2C_i",
#         "f1C_i",
#         "f1_i",
#         "f2_i",
#     ],
#     "DEPT135": [
#         "ppm",
#         "attached_protons",
#         "ppmH1s",
#         "max_bonds",
#         "label",
#         "f2C_i",
#         "f1C_i",
#         "f1_i",
#         "f2_i",
#     ],
#     "pureshift": [
#         "ppm",
#         "intensity",
#         "Width",
#         "Area",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "cosy": [
#         "f1_ppm",
#         "f2_ppm",
#         "intensity",
#         "f1_i",
#         "f1p_i",
#         "f2_i",
#         "f2p_i",
#         "f1p_ppm",
#         "f2p_ppm",
#         "f1H_i",
#         "f2H_i",
#         "f1Cp_i",
#         "f2Cp_i",
#     ],
#     "noesy": [
#         "f1_ppm",
#         "f2_ppm",
#         "intensity",
#         "f1_i",
#         "f1p_i",
#         "f2_i",
#         "f2p_i",
#         "f1p_ppm",
#         "f2p_ppm",
#         "f1H_i",
#         "f2H_i",
#         "f1Cp_i",
#         "f2Cp_i",
#     ],
#     "hsqc": [
#         "f1_ppm",
#         "f2_ppm",
#         "intensity",
#         "f1_i",
#         "f2_i",
#         "f2p_i",
#         "f1C_i",
#         "f2H_i",
#         "f2Cp_i",
#         "f2p_ppm",
#     ],
#     "hsqc_clip_cosy": [
#         "f1_ppm",
#         "f2_ppm",
#         "intensity",
#         "f1_i",
#         "f2_i",
#         "f2p_i",
#         "f1C_i",
#         "f2H_i",
#         "f2Cp_i",
#         "f2p_ppm",
#     ],
#     "ddept_ch3_only": [
#         "f1_ppm",
#         "f2_ppm",
#         "intensity",
#         "f1_i",
#         "f2_i",
#         "f2p_i",
#         "f1C_i",
#         "f2H_i",
#         "f2Cp_i",
#         "f2p_ppm",
#     ],
#     "hmbc": [
#         "f1_ppm",
#         "f2_ppm",
#         "intensity",
#         "f2p_ppm",
#         "f1_i",
#         "f2_i",
#         "f2p_i",
#         "f1C_i",
#         "f2H_i",
#         "f2Cp_i",
#     ],
# }

# excel_orig_df_columns = {

#         "nmrAssignments": [
#         "atom_idx",
#         "atomNumber",
#         "numProtons",
#         "f1_ppm",
#         "iupacLabel",
#         "f2_ppm1",
#         "f2_ppm2",
#         "f2_ppm3",
#     ],
#     "c13predictions": ["atomNumber", "atom_idx", "numProtons", "ppm"],
#     "molecule": ["molecule", "smiles"],
#     "H1_1D": ["Name", "Shift", "Range", "H's", "Integral", "Class", "J's", "Method"],
#     "C13_1D": [
#         "ppm",
#         "Intensity",
#         "Width",
#         "Area",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "DEPT135": [
#         "ppm",
#         "Intensity",
#         "Width",
#         "Area",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "H1_pureshift": [
#         "ppm",
#         "Intensity",
#         "Width",
#         "Area",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "COSY": [
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "HSQC": [
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "HSQC_CLIPCOSY": [
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "DDEPT_CH3_ONLY": [
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "HMBC": [
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
#     "NOESY": [
#         "f2 (ppm)",
#         "f1 (ppm)",
#         "Intensity",
#         "Width f2",
#         "Width f1",
#         "Volume",
#         "Type",
#         "Flags",
#         "Impurity/Compound",
#         "Annotation",
#     ],
# }
