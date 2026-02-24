# JSON File Structure

## Overview

This file is sent to the simpleNMR server and holds all the information for the program to attempt to assign the data.

The file is made up of two types of JSON data structures:

    - Data structures that hold information that is required to run the analysis
      
    - Data structures that hold information about the NMR spectra and the peaks and integrals picked

## List of the components in the file

- [smiles](#11-smiles)

- [molfile](#12-molfile)

- [hostname](#13-hostname)

- [workingDirectory](#14-workingdirectory)

- [workingFilename](#15-workingfilename)

- [MNOVAcalcMethod](#16-mnovacalcmethod)

- [carbonCalcPositionsMethod](#17-carboncalcpositionsmethod)

- [allAtomsInfo](#18-allatomsinfo)

- [carbonAtomsInfo](#19-carbonatomsinfo)

- [c13predictions](#110-c13predictions)

- [nmrAssignments](#111-nmrassignments)

- [chosenSpectra](#112-chosenspectra)

- [spectraWithPeaks](#113-spectrawithpeaks)

- [exptIdentifiers](#114-exptidentifiers)

- [simulatedAnnealing](#115-simulatedannealing)

- [ml_consent](#116-ml_consent)

- ----------

- C13_1D_0,  C13_1D_1, ... (optional)

- H1_1D_0, H1_1D_1, ... (optional)

- H1_pureshift_0, H1_pureshift_1, (optional)

- HSQC_0, HSQC_1, ... (required)

- HMBC_0, HMBC_1, ... (optional)

- COSY_0, COSY_1, ... (optional)

- HSQC_CLIPCOSY_0,  HSQC_CLIPCOSY_1, ... (optional)

- DDEPT_CH3_ONLY_0, DDEPT_CH3_ONLY_1, ... (optional)



## Overview of the two data structure types

### 1 Meta data

```python
{
    "data_name": {
        "datatype": str,      # Identifier for the type of data
        "count": int,         # Number of data items
        "data": dict or {}    # Dictionary of data items indexed by string keys
    }
}
```
- #### 1.1 smiles
```python
    "smiles": {
        "datatype": "smiles",                        # identifier must match the key
        "count": 1,                                 # always 1, only one smile string allowed for now
        "data": {
            "0": "CC(=O)C=CC1C(C)=CCCC1(C)C"
        },  
    },
```
[↑ Back to list](#list-of-the-components-in-the-file)

- #### 1.2 molfile
```python
    "molfile": {
        "datatype": "molfile",
        "count": 1,                                # always 1, only one molfile string allowed for now
        "data": {
        "0": "\n     RDKit          2D\n\n 14 14  0  0  0  0  0  0  0  0999 V2000\n    4.7546    0.1821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1565    0.1821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4556    0.9321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8575    2.4321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4416    3.1821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8575    0.9321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.4416    0.1821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7407    0.9321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5226   -0.9670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0537    0.9321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.4057   -0.9670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7407    2.4321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1565    3.1821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7546   -1.3179    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1 10  1  0\n  1  3  1  0\n  1 14  2  0\n  2  6  1  0\n  2  3  2  0\n  4 13  1  0\n  4  6  1  0\n  4  5  2  0\n  5 12  1  0\n  6  7  1  0\n  7 11  1  0\n  7  9  1  0\n  7  8  1  0\n  8 12  1  0\nM  END"
      },
  },
```
[↑ Back to list](#list-of-the-components-in-the-file)
  - #### 1.3 hostname
```json
    "hostname": {
        "datatype": "hostname",
        "count": 1,
      
        "data": {
            "0": "0x624bb2cbbxxx"
        },
    },
```
[↑ Back to list](#list-of-the-components-in-the-file)

- #### 1.4 workingDirectory

```json
    "workingDirectory": {
        "datatype": "workingDirectory",
        "count": 1,
      
        "data": {
            "0": "/opt/topspin4.5.0/examdata/exam_CMCse_1"     # Full directory path where json file is
        },   
     },
```
[↑ Back to list](#list-of-the-components-in-the-file)
- #### 1.5 workingFilename

```json
    "workingFilename": {
        "datatype": "workingFilename"
        "count": 1,
        "data": {
            "0": "alpha_ionone"                          # file name of data set without extension
        },
    },
```
[↑ Back to list](#list-of-the-components-in-the-file)

  - #### 1.6 MNOVAcalcMethod

```json
    "MNOVAcalcMethod": {
        "datatype": "MNOVAcalcMethod",
        "count": 1,
      
        "data": {
            "0": "JEOL Predict",                     # "MNOVA Predict", "NMRSHIFTDB2 Predict", etc
        },
    },
```
[↑ Back to list](#list-of-the-components-in-the-file)

- #### 1.7 carbonCalcPositionsMethod
  Not really used by the program, but included for backward compatibility
  ```json
      "carbonCalcPositionsMethod": {
          "datatype": "carbonCalcPositionsMethod",
          "count": 1,
          "data": {
              "0": "Calculated Positions"
          }
      },
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)

- #### 1.8 allAtomsInfo
  
  All hetero atoms
  
  ```python
      "allAtomsInfo": {
          "datatype": "allAtomsInfo",
          "count": 2,                     
          "data": {
              "0": {
                  "atom_idx": 0,                    # atom_idx and id must have the same value
                  "id": 0,                          # corresponds to the molfile atom index
                  "atomNumber": "1",                # can be a number or string and is the value displayed
                                                    # on the molecule
                  "symbol": "C", 
                  "numProtons": 0
              },
              "1": {
                  "atom_idx": 1,
                  "id": 1,
                  "atomNumber": "2",
                  "symbol": "C",
                  "numProtons": 1
              },
  				 }
       }
  
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
- #### 1.9 carbonAtomsInfo
  
  Carbon atoms only
  
  ```python
      "carbonAtomsInfo": {
          "datatype": "carbonAtomsInfo",
          "count": 2, 
          "data": {
              "0": {
                  "atom_idx": 0,
                  "id": 0,
                  "atomNumber": "1",
                  "symbol": "C",
                  "numProtons": 0
              },
              "1": {
                  "atom_idx": 1,
                  "id": 1,
                  "atomNumber": "2",
                  "symbol": "C",
                  "numProtons": 1
              },
          }
      }
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
- #### 1.10 c13predictions
  
  This structure holds the carbon 13 ppm predicted values. When nmshiftDB2 is used to predict the values the data structure is empty as nmrshiftdb is called from the simpleNMR server after the json file has been posted.
  
  When Mestrenova (with prediction license) or JEOL Jason (prediction included) then this data structure is used to transfer the predictions across to the simpleNMR server
  
  ```python
      "c13predictions": {
          "datatype": "c13predictions",
          "count": 0,                            # 0 if nmrshiftDB2 is used
          "data": {}
      },
  ```
  
  ```python
      "c13predictions": {
        "count": 3,                               # prediction from Mestrenova or JEOL Jason
        "data": {
          "0": {
            "atomNumber": 1,                      # value displayed on the molecule
            "atom_idx": 0,                        # mol file atom index
            "numProtons": 0,                      # number of protons attached to carbon 0,1,2,3
            "ppm": 197.1999969482422              # carbon chemical shift ppm
          },
          "1": {
            "atomNumber": 2,
            "atom_idx": 1,
            "numProtons": 1,
            "ppm": 148.25
          },
          "2": {
            "atomNumber": 3,
            "atom_idx": 2,
            "numProtons": 1,
            "ppm": 132.4499969482422
          }
        }
      }
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
- #### 1.11 nmrAssignments
  
  ```python
      "nmrAssignments": {
          "datatype": "nmrAssignments",
          "count": 13,
          "data": {
              "1": {
                  "atom_idx": 1,                         # mol file atom index number
                  "atomNumber": "14",                    # carbon value displayed on molecule
                  "numProtons": 3,                       # number of protons attached to carbon
                  "f1_ppm": 26.88671562420045,           # 13C ppm value
                  "iupacLabel": "",
                  "f2_ppm1": 2.227323226450975,          # 1H ppm value 
                  "f2_ppm2": 2.227323226450975,
                  "f2_ppm3": 2.227323226450975
              },
              "2": {
                  "atom_idx": 2,
                  "atomNumber": "11",
                  "numProtons": 0,
                  "f1_ppm": 198.3511985110063,
                  "iupacLabel": "",
                  "f2_ppm1": "undefined",                # since quaternary carbon 1H ppm values undefined
                  "f2_ppm2": "undefined",
                  "f2_ppm3": "undefined"
              },
              "4": {
                  "atom_idx": 4,
                  "atomNumber": "10",
                  "numProtons": 1,
                  "f1_ppm": 132.2982230359739,
                  "iupacLabel": "",
                  "f2_ppm1": 6.027334964222518,         # CH group so only one 1H ppm value set
                  "f2_ppm2": "undefined",
                  "f2_ppm3": "undefined"
              },
              "11": {
                  "atom_idx": 11,
                  "atomNumber": "6",
                  "numProtons": 2,
                  "f1_ppm": 31.178625044027,
                  "iupacLabel": "",
                  "f2_ppm1": 1.4321982837216,          # CH group, 2 1H ppm values set, can be different
                  "f2_ppm2": 1.204580518086928,
                  "f2_ppm3": "undefined"
              },
          }
      }
            
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
- #### 1.12 chosenSpectra
  
  ```python
      "chosenSpectra": {
          "datatype": "chosenSpectra",
          "count": 5,
          "data": {
              "0": "1H 1D zg exam_CMCse_1.1.fid_1 SKIP",   #nucleus, type, pulsesequence, exptIdentifier
              "1": "13C 1D zgpg30 exam_CMCse_1.5.fid_2 C13_1D",
              "2": "[1H, 1H] COSY cosygpmfqf exam_CMCse_1.2.ser_3 COSY",
              "3": "[13C, 1H] HSQC-EDITED hsqcedetgpsp.3 exam_CMCse_1.3.ser_4 HSQC",
              "4": "[13C, 1H] HMBC hmbcetgpl3nd exam_CMCse_1.4.ser_5 HMBC"
          }
      },
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
- #### 1.13 spectraWithPeaks
  
  ```python
      "spectraWithPeaks": {
          "datatype": "spectraWithPeaks",
          "count": 5,
          "data": {
              "0": "1H 1D zg exam_CMCse_1.1.fid_1",    # nucleus, type, pulsesequence, 
              "1": "13C 1D zgpg30 exam_CMCse_1.5.fid_2",
              "2": "[1H, 1H] COSY cosygpmfqf exam_CMCse_1.2.ser_3",
              "3": "[13C, 1H] HSQC-EDITED hsqcedetgpsp.3 exam_CMCse_1.3.ser_4",
              "4": "[13C, 1H] HMBC hmbcetgpl3nd exam_CMCse_1.4.ser_5"
          }
      }
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
- #### 1.14 exptIdentifiers
  
  ```python
      "exptIdentifiers": {
          "count": 5,
          "datatype": "exptIdentifiers",
          "data": {
              "0": "SKIP",             # SKIP, H1_1D, C13_1D, H1_pureshift, DEPT135,  
              "1": "C13_1D",           # COSY, HSQC, HMBC, HSQC_CLIPCOSY, DDEPT_CH3_ONLY
              "2": "COSY",
              "3": "HSQC",
              "4": "HMBC"
          }
      },
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
- #### 1.15 simulatedAnnealing
  
  ```python
      "simulatedAnnealing": {
          "datatype": "simulatedAnnealing",
          "count": 1,
          "data": {
              "0": true
          }
      },
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
- #### 1.16 ml_consent
  
  ```python
      "ml_consent": {
          "datatype": "ml_consent",
          "count": 1,
          "data": {
              "0": false
          }
      },
  ```
  [↑ Back to list](#list-of-the-components-in-the-file)
  
  


### 2 Spectra data

```json

data_spectrum_name: {
		"spectrum_attributes": str, int, or list
		.
		.
		.
		"peaks": {},
		"integrals": {},
		"multiplets": {}
}
```

- #### 2.1 spectrum_attributes

  ```json
        "datatype": str,                 # Identifier for the type of data
        "pulsesequence": str,
        "experimenttype": str,           # HSQC, HMBC, COSY, C13_1D, 
        "datafilename": str,             # "/opt/topspin4.5.0/examdata/exam_CMCse_1/3"
        "solvent": str,                  # "CDCl3"
        "expt_fn": str                   # "3"
        "specfrequency": [float, float], # [400.13,100.62] for 2D [100.62] for 1D
        "type": str,                     # "2D" or "1D"
        "temperature": float,            # 26.55,
        "nucleus": [ str, str].          # ["1H","13C"] for 2D, "1H" or "13C" for 1D
  ```
  
- #### 2.2 peaks
  ```json
    "peaks": {
        "count": int,                # number of peaks in data
        "datatype": str,             # "peaks"
            
        "data": {
            "0": {
                "intensity": float,               # peak height,
                "delta2": float,                  # 0.0      
                "delta1": float,                  # 198.3  direct dimension if 1D
                "type": int                       # 0 means that it is a proper signal
                "annotation": str,                # if spectrum is HSQC then one can use annotation
                                                  # to tell if the peak is CH3, CH2 or CH1 to help
                                                  # simpleNMR out
                    
             },
             "1": {
                    "intensity": 51278040.375,
                    "delta2": 0.0,
                    "delta1": 148.92,
                    "annotation": "",
                    "type": 0
             },     
        },
    }
  ```
  
- #### 2.3 integrals

  ```json
      "integrals": {
          "count": int                               # number of integrals
          "normValue": int,                          # integral normalisation value not used
          "datatype": str,                           # "integrals"
        
          "data": {
             "0": {
                 "intensity": float,               # peak intensity
                 "rangeMin1": float,               # f2 min value
                 "rangeMax1": float,               # f2 max value
                 "rangeMin2": float,               # f1 min value
                 "rangeMax2": float,               # f1 max value
                 "delta1": float,                  # f2 centre of integral
                 "delta2": float,                  # f1 centre of integral
                 "annotation": str,                # not really used
                 "type": int                       # 0  means it is a real peak
                  },
             "1": {
                 "intensity": 23545327372.44,
                 "rangeMin1": 5.97,
                 "rangeMax1": 6.07,
                 "rangeMin2": 131.82,
                 "rangeMax2": 132.63,
                 "delta1": 6.028,
                 "delta2": 132.23,
                 "annotation": "",
                 "type": 0
             },
         }                                           # end of data
     }                                               # end of integrals
  ```

  

- #### 2.4 multiplets

The multiplets data is normally not used by simpleNMR, but is added to each spectrum dataset for compatibiliity. It is only used in a 1-D proton dataset when the simpleNMR program was being developed. It is kept for compatibility. Appendix A describes what information is recorded when it is used. for now, the minimum information is described below that is added to the file.

```json
      "multiplets": {
        "datatype": "multiplets",
        "normValue": 1,
        "count": 0,
        "data": {},
      },
```



## Appendix A

### Multiplet example of 1D proton dataset

```python
        			.
  						.
    					.
  			"multiplets": {
            "datatype": "multiplets",
            "normValue": 17635.57460361719,
            "count": 11,
            "data": {
                "0": {
                    "delta1": 6.02835345402265,
                    "nH": 1,
                    "realH": 1,
                    "integralValue": 17635.57460361719,
                    "category": "dd",
                    "type": 0,
                    "jlistcount": 2,
                    "jvals": [
                        0.8076171875,
                        15.78369140625
                    ]
                },
                "1": {
                    "delta1": 5.479175846577852,
                    "nH": 1,
                    "realH": 1,
                    "integralValue": 17993.810027956963,
                    "category": "tt",
                    "type": 0,
                    "jlistcount": 2,
                    "jvals": [
                        1.6058349609375,
                        3.89111328125
                    ]
                },
            },
            "integrals": {
            "datatype": "integrals",
            "count": 11,
            "normValue": 17635.57460361719,
            "data": {
                "0": {
                    "intensity": 17635.57460361719,
                    "rangeMin1": 5.985485831704805,
                    "rangeMin2": 0,
                    "rangeMax1": 6.08238015447256,
                    "rangeMax2": 0,
                    "delta1": 6.0339329930886825,
                    "delta2": 0,
                    "type": 0
                },
                "1": {
                    "intensity": 17993.810027956963,
                    "rangeMin1": 5.451617112141305,
                    "rangeMin2": 0,
                    "rangeMax1": 5.529512548091851,
                    "rangeMax2": 0,
                    "delta1": 5.490564830116578,
                    "delta2": 0,
                    "type": 0
                }
            }
```

