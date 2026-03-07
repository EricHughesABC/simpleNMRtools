# simpleNMR JSON File Structure

## Overview

The simpleNMR workflow involves two JSON files:

**Input file** (`*_mresnova.json`) — assembled by the MNova plugin and posted to the simpleNMR server. It carries all the structural and spectral information needed to run the NMR assignment.

**Output file** (`*_from_simplemnova.json`) — returned by the simpleNMR server once the user accepts the result. This file is read back by the MNova plugin to display the final assignments on the molecule.

---

## Part 1 — Input File (`*_mresnova.json`)

The input file contains two categories of data:

- **Metadata** — molecular structure, prediction values, and configuration
- **Spectra data** — peaks, integrals, and multiplets picked from each NMR experiment

### 1.1 General metadata structure

All metadata blocks share a common envelope:

```json
"block_name": {
    "datatype": "block_name",   // must match the key
    "count": 1,                  // number of items in data
    "data": {
        "0": ...                 // items indexed as string keys "0", "1", ...
    }
}
```

---

### Metadata blocks

#### smiles

The SMILES string for the molecule. Only one is allowed per file.

```json
"smiles": {
    "datatype": "smiles",
    "count": 1,
    "data": {
        "0": "Brc1ccccc1O"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### molfile

The MDL molfile for the molecule. Only one is allowed per file. The atom order defined here determines the `atom_idx` values used throughout the rest of the file.

```json
"molfile": {
    "datatype": "molfile",
    "count": 1,
    "data": {
        "0": "\n  Mrv2011 ...\nM  END\n"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### hostname

A unique identifier for the machine submitting the job. Used for licensing and logging.

```json
"hostname": {
    "datatype": "hostname",
    "count": 1,
    "data": {
        "0": "5PAN0-H2J3HEPH-WQWA0-N2EHQRNY"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### workingDirectory

Full path to the directory containing the data files.

```json
"workingDirectory": {
    "datatype": "workingDirectory",
    "count": 1,
    "data": {
        "0": "Y:/downloads/Eric/student_projects/2025/simpleNMR_examples/2-bromophenol"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### workingFilename

The dataset name without file extension.

```json
"workingFilename": {
    "datatype": "workingFilename",
    "count": 1,
    "data": {
        "0": "2-bromophenol"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### MNOVAcalcMethod

Identifies which method was used to generate the ¹³C predictions. Common values are `"MNOVA Predict"`, `"NMRSHIFTDB2 Predict"`, `"JEOL Predict"`, and `"MNova Manually Assigned data"`.

```json
"MNOVAcalcMethod": {
    "datatype": "MNOVAcalcMethod",
    "count": 1,
    "data": {
        "0": "MNova Manually Assigned data"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### carbonCalcPositionsMethod

Kept for backward compatibility; not actively used by the program.

```json
"carbonCalcPositionsMethod": {
    "datatype": "carbonCalcPositionsMethod",
    "count": 1,
    "data": {
        "0": "Calculated Positions"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### allAtomsInfo

Information about all heavy atoms in the molecule (carbons and heteroatoms). The `atom_idx` and `id` fields must always be equal and correspond to the atom index in the molfile.

```json
"allAtomsInfo": {
    "datatype": "allAtomsInfo",
    "count": 2,
    "data": {
        "0": {
            "atom_idx": 0,       // molfile atom index (0-based)
            "id": 0,             // must equal atom_idx
            "atomNumber": "1",   // label displayed on the molecule
            "symbol": "C",
            "numProtons": 0      // number of directly attached protons
        },
        "1": {
            "atom_idx": 1,
            "id": 1,
            "atomNumber": "2",
            "symbol": "C",
            "numProtons": 1
        }
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### carbonAtomsInfo

The same structure as `allAtomsInfo` but restricted to carbon atoms only.

```json
"carbonAtomsInfo": {
    "datatype": "carbonAtomsInfo",
    "count": 6,
    "data": {
        "0": {
            "atom_idx": 0,
            "id": 0,
            "atomNumber": "1",
            "symbol": "C",
            "numProtons": 1
        }
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### c13predictions

Predicted ¹³C chemical shifts for each carbon. When NMRShiftDB2 is used for prediction, this block is sent empty (`count: 0`, `data: {}`) because the server calls NMRShiftDB2 directly after receiving the file. When MNova Predict or JEOL Jason is used, the predictions are populated here before posting.

Empty (NMRShiftDB2 path):
```json
"c13predictions": {
    "datatype": "c13predictions",
    "count": 0,
    "data": {}
}
```

Populated (MNova / JEOL path):
```json
"c13predictions": {
    "datatype": "c13predictions",
    "count": 6,
    "data": {
        "0": {
            "atomNumber": 1,     // label displayed on molecule
            "atom_idx": 0,       // molfile atom index
            "numProtons": 1,     // protons attached to this carbon
            "ppm": 129.302       // predicted 13C chemical shift (ppm)
        },
        "1": {
            "atomNumber": 2,
            "atom_idx": 1,
            "numProtons": 1,
            "ppm": 121.77
        }
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### nmrAssignments

The experimentally assigned ¹³C and ¹H chemical shifts for each carbon, derived from HSQC (and HSQC-CLIPCOSY) cross-peaks. Quaternary carbons have no HSQC correlation, so their `f2_ppm` values are `"undefined"`. CH₂ carbons may have two distinct proton shifts (`f2_ppm1` and `f2_ppm2`). Keys start at `"1"` (not `"0"`).

```json
"nmrAssignments": {
    "datatype": "nmrAssignments",
    "count": 6,
    "data": {
        "1": {
            "atom_idx": 1,                   // molfile atom index
            "atomNumber": "1",               // label displayed on molecule
            "numProtons": 1,                 // 0 = quat C, 1 = CH, 2 = CH2, 3 = CH3
            "f1_ppm": 129.1695876106564,     // assigned 13C chemical shift (ppm)
            "iupacLabel": "",
            "f2_ppm1": 7.2485803027607805,   // first (or only) 1H shift; "undefined" if quat
            "f2_ppm2": "undefined",          // second 1H shift (CH2 only)
            "f2_ppm3": "undefined"           // third 1H shift (CH3 only)
        },
        "4": {
            "atom_idx": 4,
            "atomNumber": "5",
            "numProtons": 0,                 // quaternary carbon — no HSQC peak
            "f1_ppm": 152.27407251687237,
            "iupacLabel": "",
            "f2_ppm1": "undefined",
            "f2_ppm2": "undefined",
            "f2_ppm3": "undefined"
        }
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### chosenSpectra

A list of all spectra selected by the user in MNova, with their experiment type labels appended. The experiment type label (last token) matches the entries in `exptIdentifiers`. A value of `SKIP` means the spectrum was excluded from the analysis.

```json
"chosenSpectra": {
    "datatype": "chosenSpectra",
    "count": 3,
    "data": {
        "0": "[13C, 1H] Unknown gns_noah3-BSScc.eeh 04114305.15001.ser_5 HMBC",
        "1": "[13C, 1H] Unknown gns_noah3-BSScc.eeh 04114305.15002.ser_6 HSQC",
        "2": "[13C, 1H] Unknown gns_noah3-BSScc.eeh 04114305.15003.ser_7 HSQC_CLIPCOSY"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### spectraWithPeaks

Lists the spectra for which peaks or integrals have been picked. The format is the same as `chosenSpectra` but without the experiment type label suffix.

```json
"spectraWithPeaks": {
    "datatype": "spectraWithPeaks",
    "count": 3,
    "data": {
        "0": "[13C, 1H] Unknown gns_noah3-BSScc.eeh 04114305.15001.ser_5",
        "1": "[13C, 1H] Unknown gns_noah3-BSScc.eeh 04114305.15002.ser_6",
        "2": "[13C, 1H] Unknown gns_noah3-BSScc.eeh 04114305.15003.ser_7"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### exptIdentifiers

The experiment type label for each spectrum in `chosenSpectra`, in the same order. Valid labels are: `SKIP`, `H1_1D`, `C13_1D`, `H1_pureshift`, `DEPT135`, `COSY`, `HSQC`, `HMBC`, `HSQC_CLIPCOSY`, `DDEPT_CH3_ONLY`.

```json
"exptIdentifiers": {
    "datatype": "exptIdentifiers",
    "count": 3,
    "data": {
        "0": "HMBC",
        "1": "HSQC",
        "2": "HSQC_CLIPCOSY"
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### simulatedAnnealing

Whether to use simulated annealing when optimising the carbon assignment. Set to `true` in normal use.

```json
"simulatedAnnealing": {
    "datatype": "simulatedAnnealing",
    "count": 1,
    "data": {
        "0": true
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

#### ml_consent

Whether the user has consented to their data being used for machine learning. Defaults to `false`.

```json
"ml_consent": {
    "datatype": "ml_consent",
    "count": 1,
    "data": {
        "0": false
    }
}
```

[↑ Back to top](#part-1--input-file-_mresnovajson)

---

### 1.2 Spectra data blocks

Each NMR experiment contributes a separate top-level block. The key is constructed from the spectrum filename and an index suffix (e.g. `HMBC_0`, `HSQC_0`, `C13_1D_0`). Multiple acquisitions of the same experiment type are indexed `_0`, `_1`, etc.

Required: at least one `HSQC` or `HSQC_CLIPCOSY` block.
Optional: `C13_1D`, `H1_1D`, `H1_pureshift`, `HMBC`, `COSY`, `DDEPT_CH3_ONLY`.

#### Spectrum attributes

Each spectrum block begins with a set of scalar attributes describing the acquisition:

```json
"HSQC_0": {
    "datatype": "nmrspectrum",
    "experimenttype": "HSQC",         // experiment type label
    "pulsesequence": "hsqcedetgpsp.3",
    "datafilename": "/path/to/data/3",
    "solvent": "CDCl3",
    "expt_fn": "3",                    // experiment folder number
    "specfrequency": [400.13, 100.62], // [1H freq, 13C freq] in MHz; single value for 1D
    "type": "2D",                      // "1D" or "2D"
    "temperature": 26.55,              // in °C
    "nucleus": ["1H", "13C"],          // single string for 1D spectra
    ...
}
```

#### peaks

Peak positions picked from the spectrum. For 2D spectra, `delta1` is the F2 (direct, usually ¹H) dimension and `delta2` is the F1 (indirect, usually ¹³C) dimension. For 1D spectra, `delta2` is 0.

```json
"peaks": {
    "datatype": "peaks",
    "count": 4,
    "data": {
        "0": {
            "intensity": 51278040.375,
            "delta1": 7.249,       // F2 chemical shift (ppm)
            "delta2": 129.17,      // F1 chemical shift (ppm); 0 for 1D
            "type": 0,             // 0 = genuine signal
            "annotation": ""       // optional: "CH", "CH2", "CH3" for edited HSQC
        }
    }
}
```

The `annotation` field on HSQC peaks can be used to indicate multiplicity (`CH`, `CH2`, `CH3`) when an edited HSQC is available, helping the assignment algorithm.

#### integrals

Integral regions. For 1D spectra, `rangeMin2`/`rangeMax2`/`delta2` are 0.

```json
"integrals": {
    "datatype": "integrals",
    "count": 4,
    "normValue": 1,
    "data": {
        "0": {
            "intensity": 23545327372.44,
            "rangeMin1": 5.97,    // F2 lower bound (ppm)
            "rangeMax1": 6.07,    // F2 upper bound (ppm)
            "rangeMin2": 131.82,  // F1 lower bound (ppm)
            "rangeMax2": 132.63,  // F1 upper bound (ppm)
            "delta1": 6.028,      // F2 centre (ppm)
            "delta2": 132.23,     // F1 centre (ppm)
            "annotation": "",
            "type": 0
        }
    }
}
```

#### multiplets

Multiplet data from 1D ¹H spectra. This block is only meaningful for 1D proton data; for all other spectra it is present as an empty placeholder for compatibility.

Empty placeholder (all non-1D spectra):
```json
"multiplets": {
    "datatype": "multiplets",
    "normValue": 1,
    "count": 0,
    "data": {}
}
```

Populated (1D ¹H spectrum, see Appendix A for full detail):
```json
"multiplets": {
    "datatype": "multiplets",
    "normValue": 17635.57,
    "count": 11,
    "data": {
        "0": {
            "delta1": 6.028,          // 1H chemical shift (ppm)
            "nH": 1,                  // number of protons (from integral)
            "realH": 1,               // actual number of protons
            "integralValue": 17635.57,
            "category": "dd",         // multiplicity code: s, d, t, q, dd, dt, tt, ...
            "type": 0,
            "jlistcount": 2,          // number of J coupling values
            "jvals": [0.808, 15.784]  // J coupling constants (Hz)
        }
    }
}
```

---

## Part 2 — Output File (`*_from_simplemnova.json`)

This file is returned by the simpleNMR server and loaded back into MNova to render the final assignment on the molecule. It is a flat JSON object (no `datatype`/`count`/`data` envelope).

### Top-level structure

| Key | Type | Description |
|-----|------|-------------|
| `nodes_orig` | array | Carbon node data as received from the server (original) |
| `nodes_now` | array | Current carbon node data after any user edits; initially identical to `nodes_orig` |
| `links` | array | HMBC and COSY correlations connecting carbon nodes |
| `catoms` | array | Carbon atom list used by the MNova display (current state) |
| `catoms_orig` | array | Carbon atom list as originally returned (pre-edit snapshot) |
| `best_results` | object | Summary statistics for the best assignment found |
| `molgraph` | object | NetworkX-format molecular graph of all heavy atoms, used for graph operations |
| `shortest_paths` | object | All-pairs shortest bond-path distance matrix across all heavy atoms |
| `svg` | string | RDKit-rendered SVG of the molecule for display in MNova |
| `smilesString` | string | SMILES of the molecule |
| `molfile` | string | MDL molfile of the molecule |
| `dataFrom` | string | Source identifier; always `"mnova"` |
| `oldjsondata` | object | A copy of the original spectra data from the input file |
| `workingDirectory` | string | Full path to the working directory |
| `workingFilename` | string | Dataset name without extension |
| `title` | string | Display title (placeholder value `"dummy_title"` if not set) |

---

### 2.1 Node objects (`nodes_orig`, `nodes_now`, `catoms`, `catoms_orig`)

All four arrays share the same node structure. Each element represents one carbon atom.

```json
{
    "atomNumber": 1,            // atom label displayed on molecule (integer)
    "id": 0,                    // 0-based index; matches atom_idx in input file
    "symbol": "C",
    "numProtons": 1,            // 0 = quat, 1 = CH, 2 = CH2, 3 = CH3
    "ppm": 129.1695876106564,   // assigned 13C chemical shift (ppm)
    "ppm_calculated": 129.302,  // predicted 13C chemical shift (ppm)
    "iupacLabel": "",           // IUPAC label if available
    "jCouplingClass": "",       // J coupling classification (if used)
    "jCouplingVals": "",        // J coupling values (if used)
    "visible": true,            // whether the node is shown on the molecule display
    "x": 0.7404,                // normalised x coordinate for display (0–1)
    "y": 0.8033,                // normalised y coordinate for display (0–1)
    "H1_ppm": [7.249]           // list of assigned 1H shifts; absent for quaternary carbons
}
```

Note: `nodes_orig` / `nodes_now` and `catoms_orig` / `catoms` appear redundant but serve different roles in the MNova plugin — the `catoms` arrays drive the structure diagram rendering while `nodes` arrays drive the NMR display graph.

---

### 2.2 Link objects (`links`)

Each link represents a spectral correlation between two carbon atoms.

```json
{
    "source": 4,       // id of the source carbon node
    "target": 2,       // id of the target carbon node
    "hmbc": true,      // present and true for HMBC correlations
    "weight": 1        // correlation weight used in the assignment scoring
}
```

```json
{
    "source": 2,
    "target": 1,
    "cosy": true       // present and true for COSY correlations (no "weight" field)
}
```

A link carries either `"hmbc": true` or `"cosy": true`, never both. HMBC links always include a `weight` field; COSY links do not.

---

### 2.3 molgraph

A NetworkX-compatible graph representation of the full molecular structure, including all heavy atoms (carbons and heteroatoms). This is used internally by the server for graph-based operations such as computing shortest paths and validating HMBC correlations against bond topology.

The structure follows the standard NetworkX node-link JSON format.

```json
"molgraph": {
    "directed": false,
    "multigraph": false,
    "graph": {},
    "nodes": [
        {
            "id": 0,              // node index (matches atom_idx / id in input file)
            "id0": 0,             // duplicate of id; kept for compatibility
            "symbol": "C",        // element symbol
            "atom_number": 6,     // atomic number
            "numProtons": 1,      // directly attached protons
            "atomNumber": 1,      // atom label displayed on the molecule (carbons only)
            "ppm": 129.17,        // assigned 13C shift; heteroatoms use placeholder 10000
            "x": 0.7404,          // normalised display x coordinate (0–1)
            "y": 0.8033           // normalised display y coordinate (0–1)
        },
        {
            "id": 6,
            "id0": 6,
            "symbol": "O",
            "atom_number": 8,
            "numProtons": 1,
            "ppm": 10000,         // placeholder — heteroatoms have no 13C assignment
            "x": 0.5880,
            "y": 0.0696
                                  // note: no atomNumber field for heteroatoms
        }
    ],
    "links": [
        {
            "source": 0,          // id of first atom
            "target": 1           // id of second atom (one entry per bond)
        }
    ]
}
```

Key differences from the `nodes_orig` / `catoms` arrays: `molgraph` includes all heavy atoms (not just carbons), carries atomic number (`atom_number`), and uses `10000` as a sentinel ppm value for heteroatoms.

---

### 2.4 shortest_paths

A precomputed all-pairs shortest-path distance matrix for the molecular graph, expressed as bond counts. Both keys are string representations of atom `id` values. The matrix covers all heavy atoms including heteroatoms.

This is used by the assignment algorithm to check whether HMBC correlations (which should span 2–3 bonds) are topologically plausible for a given assignment.

```json
"shortest_paths": {
    "0": { "0": 0, "1": 1, "2": 2, "3": 3, "4": 2, "5": 1, "6": 3, "7": 4 },
    "1": { "0": 1, "1": 0, "2": 1, "3": 2, "4": 3, "5": 2, "6": 4, "7": 3 },
    ...
}
```

Reading the matrix: `shortest_paths["0"]["3"]` gives the number of bonds on the shortest path between atom `id=0` and atom `id=3`. A value of `0` means the atom is querying itself. The matrix is symmetric.

---

### 2.5 svg

An SVG string of the molecule rendered by RDKit, used by the MNova plugin to display the structure diagram. The image is 1000 × 600 pixels.

```json
"svg": "<?xml version='1.0' encoding='iso-8859-1'?>\n<svg class=\"center\" version='1.1' ..."
```

The string contains a complete, self-contained SVG document and can be written directly to a `.svg` file or embedded in HTML. Atom coordinates within the SVG are independent of the normalised `x`/`y` coordinates stored in the node arrays.

---

### 2.6 best_results

Summary of the quality of the best assignment found.

```json
"best_results": {
    "best_weight": 0,                    // total assignment weight score
    "best_mae": 0.2119,                  // mean absolute error vs predicted shifts (ppm)
    "best_lae": 0.5223,                  // largest absolute error vs predicted shifts (ppm)
    "best_lae_atomNumber": 6             // atomNumber with the largest error
}
```

A low `best_mae` indicates good overall agreement between the assigned and predicted ¹³C shifts. The `best_lae_atomNumber` flags which carbon shows the greatest discrepancy, which is useful for reviewing the assignment.

---

## Appendix A — Full multiplet structure (1D ¹H spectra)

```json
"multiplets": {
    "datatype": "multiplets",
    "normValue": 17635.57,
    "count": 11,
    "data": {
        "0": {
            "delta1": 6.028,
            "nH": 1,
            "realH": 1,
            "integralValue": 17635.57,
            "category": "dd",
            "type": 0,
            "jlistcount": 2,
            "jvals": [0.808, 15.784]
        },
        "1": {
            "delta1": 5.479,
            "nH": 1,
            "realH": 1,
            "integralValue": 17993.81,
            "category": "tt",
            "type": 0,
            "jlistcount": 2,
            "jvals": [1.606, 3.891]
        }
    }
},
"integrals": {
    "datatype": "integrals",
    "count": 11,
    "normValue": 17635.57,
    "data": {
        "0": {
            "intensity": 17635.57,
            "rangeMin1": 5.985,
            "rangeMax1": 6.082,
            "rangeMin2": 0,
            "rangeMax2": 0,
            "delta1": 6.034,
            "delta2": 0,
            "type": 0
        },
        "1": {
            "intensity": 17993.81,
            "rangeMin1": 5.452,
            "rangeMax1": 5.530,
            "rangeMin2": 0,
            "rangeMax2": 0,
            "delta1": 5.491,
            "delta2": 0,
            "type": 0
        }
    }
}
```
