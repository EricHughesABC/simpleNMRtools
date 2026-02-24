# NMR Assignments JSON Structure Specification

## Document Overview

This document describes the structure of NMR (Nuclear Magnetic Resonance) assignment data files used for molecular structure analysis and spectral data processing. The format supports multiple experiment types, molecular structure representations, and comprehensive spectral data with peaks, integrals, and multiplets.

## Top-Level Structure

The JSON file is a dictionary containing multiple named sections. Each section follows a consistent pattern with metadata about the data type, count, and the actual data.

### Standard Section Format

```python
{
    "section_name": {
        "datatype": str,      # Identifier for the type of data
        "count": int,         # Number of data items
        "data": dict or {}    # Dictionary of data items indexed by string keys
    }
}
```

------

## Section Definitions

### 1. Molecular Structure Data

#### 1.1 SMILES

**Purpose**: Stores the molecular structure in SMILES notation

```python
{
    "smiles": {
        "datatype": "smiles",
        "count": 1,
        "data": {
            "0": str  # SMILES string, e.g., "CC(=O)/C=C/C1C(C)=CCCC1(C)C"
        }
    }
}
```

**Fields**:

- `datatype`: Always "smiles"
- `count`: Always 1 for single molecule
- `data["0"]`: SMILES string with stereochemistry notation

#### 1.2 MOL File

**Purpose**: Stores the molecular structure in MDL MOL format (V2000)

```python
{
    "molfile": {
        "datatype": "molfile",
        "count": 1,
        "data": {
            "0": str  # Complete MOL file content with \n line endings
        }
    }
}
```

**Fields**:

- `datatype`: Always "molfile"
- `count`: Always 1
- `data["0"]`: Multi-line string containing complete V2000 MOL file format
  - Header lines
  - Counts line
  - Atom block
  - Bond block
  - Properties block (optional)
  - `M  END` terminator
  - `$$$$` file terminator

------

### 2. System Metadata

#### 2.1 Hostname

```python
{
    "hostname": {
        "datatype": "hostname",
        "count": 1,
        "data": {
            "0": str  # Machine identifier, e.g., "0x624bb2cbbfd1"
        }
    }
}
```

#### 2.2 Working Directory

```python
{
    "workingDirectory": {
        "datatype": "workingDirectory",
        "count": 1,
        "data": {
            "0": str  # Full path, e.g., "/opt/topspin4.5.0/examdata/exam_CMCse_1"
        }
    }
}
```

#### 2.3 Working Filename

```python
{
    "workingFilename": {
        "datatype": "workingFilename",
        "count": 1,
        "data": {
            "0": str  # Base filename without path or extension
        }
    }
}
```

------

### 3. Atom Information

#### 3.1 All Atoms Info

**Purpose**: Complete information about all atoms in the molecule

```python
{
    "allAtomsInfo": {
        "datatype": "allAtomsInfo",
        "count": int,  # Total number of atoms
        "data": {
            "0": {
                "atom_idx": int,      # Index in atom array (0-based)
                "id": int,            # Unique identifier (same as atom_idx)
                "atomNumber": str,    # Atom number as string ("1", "2", etc.)
                "symbol": str,        # Element symbol ("C", "H", "O", "N", etc.)
                "numProtons": int     # Number of attached protons (0-3)
            },
            # ... additional atoms indexed "1", "2", etc.
        }
    }
}
```

**Field Descriptions**:

- `atom_idx`: Zero-based index matching the atom's position
- `id`: Unique identifier (typically equals atom_idx)
- `atomNumber`: String representation of atom number (1-based)
- `symbol`: Chemical element symbol
- `numProtons`: Number of hydrogen atoms attached (important for NMR analysis)

#### 3.2 Carbon Atoms Info

**Purpose**: Filtered information containing only carbon atoms

```python
{
    "carbonAtomsInfo": {
        "datatype": "carbonAtomsInfo",
        "count": int,  # Number of carbon atoms
        "data": {
            # Same structure as allAtomsInfo but filtered to symbol == "C"
            "0": {
                "atom_idx": int,
                "id": int,
                "atomNumber": str,
                "symbol": "C",  # Always "C"
                "numProtons": int
            }
        }
    }
}
```

------

### 4. NMR Assignment Data

#### 4.1 NMR Assignments

**Purpose**: User-created or computed assignments of peaks to atoms

```python
{
    "nmrAssignments": {
        "datatype": "nmrAssignments",
        "count": int,  # Number of assignments (can be 0)
        "data": {}     # Empty if no assignments exist
    }
}
```

**Note**: When populated, each entry would map spectral peaks to specific atoms.

#### 4.2 C13 Predictions

**Purpose**: Predicted carbon-13 chemical shifts

```python
{
    "c13predictions": {
        "datatype": "c13predictions",
        "count": int,
        "data": {}  # Structure not shown in sample (empty)
    }
}
```

------

### 5. NMR Spectrum Data

Each NMR experiment is stored as a separate section with a naming pattern: `{EXPERIMENT_TYPE}_{INDEX}`

Common experiment types:

- `HMBC_0` - Heteronuclear Multiple Bond Correlation
- `HSQC_0` - Heteronuclear Single Quantum Coherence
- `COSY_0` - Correlation Spectroscopy
- `NOESY_0` - Nuclear Overhauser Effect Spectroscopy
- `H1_1D_0` - Proton 1D spectrum
- `C13_1D_0` - Carbon-13 1D spectrum
- `DEPT135_0` - Distortionless Enhancement by Polarization Transfer

#### 5.1 NMR Spectrum Structure (2D Example: HMBC)

```python
{
    "HMBC_0": {
        "datatype": "nmrspectrum",
        "origin": str,           # e.g., "Bruker XWIN-NMR"
        "type": str,             # "1D" or "2D"
        "subtype": str,          # e.g., "13C1H, HMBC"
        "experimenttype": str,   # e.g., "2D-HMBC"
        "experiment": str,       # e.g., "HMBC"
        "class": str,            # Optional classification
        "spectype": str,         # Optional spectrum type
        "pulsesequence": str,    # e.g., "hmbcetgpl3nd"
        "intrument": str,        # Instrument name, e.g., "Avance"
        "probe": str,            # Probe identifier, e.g., "36"
        "datafilename": str,     # Full path to raw data
        "nucleus": str,          # String representation, e.g., "['1H', '13C']"
        "specfrequency": str,    # String representation, e.g., "[100.612769, 400.13]"
        "temperature": str,      # Temperature in Kelvin, e.g., "300"
        
        "peaks": {
            "datatype": "peaks",
            "count": int,
            "data": {
                "0": {
                    "intensity": float,    # Peak intensity/height
                    "type": int,           # Peak type (0 = normal)
                    "annotation": str,     # User annotation (often empty)
                    "delta1": float,       # Chemical shift dimension 1 (ppm)
                    "delta2": float        # Chemical shift dimension 2 (ppm)
                }
                # ... additional peaks
            }
        },
        
        "integrals": {
            "datatype": "integrals",
            "count": int,
            "normValue": float,  # Normalization value
            "data": {}           # Individual integral regions (if any)
        },
        
        "multiplets": {
            "datatype": "multiplets",
            "count": int,
            "normValue": float,
            "data": {}  # Multiplet analysis data (if any)
        },
        
        "filename": str  # Base identifier, e.g., "HMBC_0"
    }
}
```

#### 5.2 NMR Spectrum Structure (1D Example: C13_1D)

```python
{
    "C13_1D_0": {
        "datatype": "nmrspectrum",
        "origin": "Bruker XWIN-NMR",
        "type": "1D",
        "subtype": "13C",
        "experimenttype": "1D",
        "experiment": "C13_1D",
        "class": "",
        "spectype": "",
        "pulsesequence": "zgpg30",
        "intrument": "Avance",
        "probe": "36",
        "datafilename": str,
        "nucleus": "13C",           # Single nucleus (not a list)
        "specfrequency": float,     # Single float value
        "temperature": "300",
        
        "peaks": {
            "datatype": "peaks",
            "count": int,
            "data": {
                "0": {
                    "intensity": float,
                    "type": int,
                    "annotation": str,
                    "delta1": float,  # Chemical shift (ppm)
                    "delta2": int     # Always 0 for 1D spectra
                }
            }
        },
        
        "integrals": {
            "datatype": "integrals",
            "count": 0,
            "normValue": 1,
            "data": {}
        },
        
        "multiplets": {
            "datatype": "multiplets",
            "count": 0,
            "normValue": 1,
            "data": {}
        },
        
        "filename": str
    }
}
```

**Key Differences Between 1D and 2D Spectra**:

- `nucleus`: String (1D) vs. string representation of list (2D)
- `specfrequency`: Float (1D) vs. string representation of list (2D)
- `delta2`: Always 0 (1D) vs. second dimension chemical shift (2D)

------

### 6. Spectrum Selection and Organization

#### 6.1 Chosen Spectra

**Purpose**: List of spectra selected for analysis

```python
{
    "chosenSpectra": {
        "datatype": "chosenSpectra",
        "count": int,
        "data": {
            "0": str,  # Format: "[nuclei] dimension pulseseq filename experiment"
            "1": str,  # e.g., "[1H, 13C] 2D hmbcetgpl3nd HMBC_0 HMBC"
            # ... additional entries
        }
    }
}
```

#### 6.2 Experiment Identifiers

**Purpose**: List of experiment types present in the file

```python
{
    "exptIdentifiers": {
        "datatype": "exptIdentifiers",
        "count": int,
        "data": {
            "0": str,  # e.g., "H1_1D", "HMBC", "HSQC", "COSY", "C13_1D"
            # ... additional experiment types
        }
    }
}
```

#### 6.3 Spectra With Peaks

**Purpose**: List of spectra that contain peak data

```python
{
    "spectraWithPeaks": {
        "datatype": "spectraWithPeaks",
        "count": int,
        "data": {
            "0": str,  # Format: "nucleus dimension pulseseq datafile_index"
            # e.g., "1H 1D zg 1.fid_0"
            # e.g., "[1H, 13C] HMBC hmbcetgpl3nd 4.ser_0"
        }
    }
}
```

------

### 7. Calculation and Algorithm Parameters

#### 7.1 Carbon Calculation Method

```python
{
    "carbonCalcPositionsMethod": {
        "datatype": "carbonCalcPositionsMethod",
        "count": 1,
        "data": {
            "0": str  # e.g., "Calculated Positions"
        }
    }
}
```

#### 7.2 Prediction Method

```python
{
    "MNOVAcalcMethod": {
        "datatype": "MNOVAcalcMethod",
        "count": 1,
        "data": {
            "0": str  # e.g., "NMRSHIFTDB2 Predict"
        }
    }
}
```

#### 7.3 Simulated Annealing Parameters

```python
{
    "randomizeStart": {
        "datatype": "randomizeStart",
        "count": 1,
        "data": {"0": bool}  # Whether to randomize initial state
    },
    
    "startingTemperature": {
        "datatype": "startingTemperature",
        "count": 1,
        "data": {"0": float}  # e.g., 1000
    },
    
    "endingTemperature": {
        "datatype": "endingTemperature",
        "count": 1,
        "data": {"0": float}  # e.g., 0.1
    },
    
    "coolingRate": {
        "datatype": "coolingRate",
        "count": 1,
        "data": {"0": float}  # e.g., 0.999 (exponential decay)
    },
    
    "numberOfSteps": {
        "datatype": "numberOfSteps",
        "count": 1,
        "data": {"0": int}  # e.g., 10000
    },
    
    "ppmGroupSeparation": {
        "datatype": "ppmGroupSeparation",
        "count": 1,
        "data": {"0": float}  # e.g., 2 ppm
    },
    
    "simulatedAnnealing": {
        "datatype": "simulatedAnnealing",
        "count": 1,
        "data": {"0": bool}  # Whether SA is enabled
    }
}
```

#### 7.4 User Preferences

```python
{
    "ml_consent": {
        "datatype": "ml_consent",
        "count": 1,
        "data": {"0": bool}  # Machine learning data sharing consent
    }
}
```

------

## Python Implementation Examples

### Complete Data Structure Definition

```python
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field
from enum import Enum

class ExperimentType(Enum):
    """NMR experiment types"""
    H1_1D = "H1_1D"
    C13_1D = "C13_1D"
    DEPT135 = "DEPT135"
    HSQC = "HSQC"
    HMBC = "HMBC"
    COSY = "COSY"
    NOESY = "NOESY"
    TOCSY = "TOCSY"

@dataclass
class AtomInfo:
    """Information about a single atom"""
    atom_idx: int
    id: int
    atomNumber: str
    symbol: str
    numProtons: int
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "atom_idx": self.atom_idx,
            "id": self.id,
            "atomNumber": self.atomNumber,
            "symbol": self.symbol,
            "numProtons": self.numProtons
        }

@dataclass
class Peak:
    """NMR peak data"""
    intensity: float
    type: int
    annotation: str
    delta1: float
    delta2: float  # 0 for 1D spectra
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "intensity": self.intensity,
            "type": self.type,
            "annotation": self.annotation,
            "delta1": self.delta1,
            "delta2": self.delta2
        }

@dataclass
class PeakData:
    """Collection of peaks"""
    datatype: str = "peaks"
    count: int = 0
    data: Dict[str, Peak] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "datatype": self.datatype,
            "count": self.count,
            "data": {k: v.to_dict() for k, v in self.data.items()}
        }

@dataclass
class IntegralData:
    """Integral region data"""
    datatype: str = "integrals"
    count: int = 0
    normValue: float = 1.0
    data: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "datatype": self.datatype,
            "count": self.count,
            "normValue": self.normValue,
            "data": self.data
        }

@dataclass
class MultipletData:
    """Multiplet analysis data"""
    datatype: str = "multiplets"
    count: int = 0
    normValue: float = 1.0
    data: Dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "datatype": self.datatype,
            "count": self.count,
            "normValue": self.normValue,
            "data": self.data
        }

@dataclass
class NMRSpectrum:
    """Complete NMR spectrum data"""
    datatype: str = "nmrspectrum"
    origin: str = ""
    type: str = ""  # "1D" or "2D"
    subtype: str = ""
    experimenttype: str = ""
    experiment: str = ""
    class_name: str = ""  # 'class' is a reserved keyword
    spectype: str = ""
    pulsesequence: str = ""
    intrument: str = ""  # Note: original has typo "intrument"
    probe: str = ""
    datafilename: str = ""
    nucleus: str = ""  # String or string representation of list
    specfrequency: Union[str, float] = ""  # String for 2D, float for 1D
    temperature: str = ""
    peaks: PeakData = field(default_factory=PeakData)
    integrals: IntegralData = field(default_factory=IntegralData)
    multiplets: MultipletData = field(default_factory=MultipletData)
    filename: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        result = {
            "datatype": self.datatype,
            "origin": self.origin,
            "type": self.type,
            "subtype": self.subtype,
            "experimenttype": self.experimenttype,
            "experiment": self.experiment,
            "class": self.class_name,
            "spectype": self.spectype,
            "pulsesequence": self.pulsesequence,
            "intrument": self.intrument,
            "probe": self.probe,
            "datafilename": self.datafilename,
            "nucleus": self.nucleus,
            "specfrequency": self.specfrequency,
            "temperature": self.temperature,
            "peaks": self.peaks.to_dict(),
            "integrals": self.integrals.to_dict(),
            "multiplets": self.multiplets.to_dict(),
            "filename": self.filename
        }
        return result

@dataclass
class NMRAssignmentData:
    """Complete NMR assignment file structure"""
    # Molecular structure
    smiles: str = ""
    molfile: str = ""
    
    # System metadata
    hostname: str = ""
    workingDirectory: str = ""
    workingFilename: str = ""
    
    # Atom information
    allAtomsInfo: List[AtomInfo] = field(default_factory=list)
    carbonAtomsInfo: List[AtomInfo] = field(default_factory=list)
    
    # Assignments and predictions
    nmrAssignments: Dict[str, Any] = field(default_factory=dict)
    c13predictions: Dict[str, Any] = field(default_factory=dict)
    
    # NMR Spectra (stored by experiment key)
    spectra: Dict[str, NMRSpectrum] = field(default_factory=dict)
    
    # Spectrum organization
    chosenSpectra: List[str] = field(default_factory=list)
    exptIdentifiers: List[str] = field(default_factory=list)
    spectraWithPeaks: List[str] = field(default_factory=list)
    
    # Calculation parameters
    carbonCalcPositionsMethod: str = ""
    MNOVAcalcMethod: str = ""
    
    # Simulated annealing parameters
    randomizeStart: bool = False
    startingTemperature: float = 1000.0
    endingTemperature: float = 0.1
    coolingRate: float = 0.999
    numberOfSteps: int = 10000
    ppmGroupSeparation: float = 2.0
    simulatedAnnealing: bool = True
    
    # User preferences
    ml_consent: bool = False
    
    def to_json_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable dictionary matching file format"""
        
        result = {
            "smiles": {
                "datatype": "smiles",
                "count": 1,
                "data": {"0": self.smiles}
            },
            "molfile": {
                "datatype": "molfile",
                "count": 1,
                "data": {"0": self.molfile}
            },
            "hostname": {
                "datatype": "hostname",
                "count": 1,
                "data": {"0": self.hostname}
            },
            "workingDirectory": {
                "datatype": "workingDirectory",
                "count": 1,
                "data": {"0": self.workingDirectory}
            },
            "workingFilename": {
                "datatype": "workingFilename",
                "count": 1,
                "data": {"0": self.workingFilename}
            },
            "allAtomsInfo": {
                "datatype": "allAtomsInfo",
                "count": len(self.allAtomsInfo),
                "data": {
                    str(i): atom.to_dict() 
                    for i, atom in enumerate(self.allAtomsInfo)
                }
            },
            "carbonAtomsInfo": {
                "datatype": "carbonAtomsInfo",
                "count": len(self.carbonAtomsInfo),
                "data": {
                    str(i): atom.to_dict() 
                    for i, atom in enumerate(self.carbonAtomsInfo)
                }
            },
            "nmrAssignments": {
                "datatype": "nmrAssignments",
                "count": len(self.nmrAssignments),
                "data": self.nmrAssignments
            },
            "c13predictions": {
                "datatype": "c13predictions",
                "count": len(self.c13predictions),
                "data": self.c13predictions
            },
            "chosenSpectra": {
                "datatype": "chosenSpectra",
                "count": len(self.chosenSpectra),
                "data": {str(i): spec for i, spec in enumerate(self.chosenSpectra)}
            },
            "exptIdentifiers": {
                "datatype": "exptIdentifiers",
                "count": len(self.exptIdentifiers),
                "data": {str(i): ident for i, ident in enumerate(self.exptIdentifiers)}
            },
            "spectraWithPeaks": {
                "datatype": "spectraWithPeaks",
                "count": len(self.spectraWithPeaks),
                "data": {str(i): spec for i, spec in enumerate(self.spectraWithPeaks)}
            },
            "carbonCalcPositionsMethod": {
                "datatype": "carbonCalcPositionsMethod",
                "count": 1,
                "data": {"0": self.carbonCalcPositionsMethod}
            },
            "MNOVAcalcMethod": {
                "datatype": "MNOVAcalcMethod",
                "count": 1,
                "data": {"0": self.MNOVAcalcMethod}
            },
            "randomizeStart": {
                "datatype": "randomizeStart",
                "count": 1,
                "data": {"0": self.randomizeStart}
            },
            "startingTemperature": {
                "datatype": "startingTemperature",
                "count": 1,
                "data": {"0": self.startingTemperature}
            },
            "endingTemperature": {
                "datatype": "endingTemperature",
                "count": 1,
                "data": {"0": self.endingTemperature}
            },
            "coolingRate": {
                "datatype": "coolingRate",
                "count": 1,
                "data": {"0": self.coolingRate}
            },
            "numberOfSteps": {
                "datatype": "numberOfSteps",
                "count": 1,
                "data": {"0": self.numberOfSteps}
            },
            "ppmGroupSeparation": {
                "datatype": "ppmGroupSeparation",
                "count": 1,
                "data": {"0": self.ppmGroupSeparation}
            },
            "ml_consent": {
                "datatype": "ml_consent",
                "count": 1,
                "data": {"0": self.ml_consent}
            },
            "simulatedAnnealing": {
                "datatype": "simulatedAnnealing",
                "count": 1,
                "data": {"0": self.simulatedAnnealing}
            }
        }
        
        # Add spectrum data
        for spectrum_key, spectrum in self.spectra.items():
            result[spectrum_key] = spectrum.to_dict()
        
        return result
```

### Usage Example: Creating a New File

```python
import json
from typing import List

def create_nmr_assignment_file(
    smiles: str,
    molfile: str,
    atoms: List[AtomInfo],
    spectra: Dict[str, NMRSpectrum],
    output_file: str
) -> None:
    """
    Create a new NMR assignment JSON file
    
    Args:
        smiles: SMILES string of molecule
        molfile: MOL file content
        atoms: List of atom information
        spectra: Dictionary of NMR spectra (key = experiment name)
        output_file: Path to output JSON file
    """
    
    # Create data structure
    data = NMRAssignmentData(
        smiles=smiles,
        molfile=molfile,
        hostname="generated",
        workingDirectory="/path/to/data",
        workingFilename="output",
        allAtomsInfo=atoms,
        carbonAtomsInfo=[a for a in atoms if a.symbol == "C"],
        spectra=spectra
    )
    
    # Populate experiment identifiers
    data.exptIdentifiers = list(set(s.experiment for s in spectra.values()))
    
    # Convert to JSON and write
    with open(output_file, 'w') as f:
        json.dump(data.to_json_dict(), f, indent=4)

# Example: Create a simple 1D C13 spectrum
atoms = [
    AtomInfo(atom_idx=0, id=0, atomNumber="1", symbol="C", numProtons=3),
    AtomInfo(atom_idx=1, id=1, atomNumber="2", symbol="C", numProtons=2),
]

peaks = PeakData(
    count=2,
    data={
        "0": Peak(intensity=100.0, type=0, annotation="", delta1=20.5, delta2=0),
        "1": Peak(intensity=95.0, type=0, annotation="", delta1=30.2, delta2=0),
    }
)

spectrum = NMRSpectrum(
    origin="Generated",
    type="1D",
    subtype="13C",
    experimenttype="1D",
    experiment="C13_1D",
    pulsesequence="zgpg30",
    intrument="Simulated",
    probe="1",
    nucleus="13C",
    specfrequency=100.0,
    temperature="300",
    peaks=peaks,
    filename="C13_1D_0"
)

create_nmr_assignment_file(
    smiles="CC",
    molfile="...",  # Full MOL file content
    atoms=atoms,
    spectra={"C13_1D_0": spectrum},
    output_file="output.json"
)
```

### Usage Example: Reading an Existing File

```python
def parse_nmr_assignment_file(filepath: str) -> NMRAssignmentData:
    """Parse an NMR assignment JSON file"""
    
    with open(filepath, 'r') as f:
        raw_data = json.load(f)
    
    # Parse basic metadata
    data = NMRAssignmentData(
        smiles=raw_data["smiles"]["data"]["0"],
        molfile=raw_data["molfile"]["data"]["0"],
        hostname=raw_data["hostname"]["data"]["0"],
        workingDirectory=raw_data["workingDirectory"]["data"]["0"],
        workingFilename=raw_data["workingFilename"]["data"]["0"],
    )
    
    # Parse atoms
    for atom_dict in raw_data["allAtomsInfo"]["data"].values():
        atom = AtomInfo(**atom_dict)
        data.allAtomsInfo.append(atom)
    
    for atom_dict in raw_data["carbonAtomsInfo"]["data"].values():
        atom = AtomInfo(**atom_dict)
        data.carbonAtomsInfo.append(atom)
    
    # Parse spectra
    for key, value in raw_data.items():
        if isinstance(value, dict) and value.get("datatype") == "nmrspectrum":
            spectrum = parse_spectrum(value)
            data.spectra[key] = spectrum
    
    # Parse algorithm parameters
    data.startingTemperature = raw_data["startingTemperature"]["data"]["0"]
    data.endingTemperature = raw_data["endingTemperature"]["data"]["0"]
    data.coolingRate = raw_data["coolingRate"]["data"]["0"]
    data.numberOfSteps = raw_data["numberOfSteps"]["data"]["0"]
    data.simulatedAnnealing = raw_data["simulatedAnnealing"]["data"]["0"]
    data.ml_consent = raw_data["ml_consent"]["data"]["0"]
    
    return data

def parse_spectrum(spec_dict: Dict[str, Any]) -> NMRSpectrum:
    """Parse a single spectrum from dictionary"""
    
    # Parse peaks
    peaks = PeakData(
        count=spec_dict["peaks"]["count"],
        data={
            k: Peak(**v) for k, v in spec_dict["peaks"]["data"].items()
        }
    )
    
    # Create spectrum object
    spectrum = NMRSpectrum(
        origin=spec_dict["origin"],
        type=spec_dict["type"],
        subtype=spec_dict["subtype"],
        experimenttype=spec_dict["experimenttype"],
        experiment=spec_dict["experiment"],
        class_name=spec_dict["class"],
        spectype=spec_dict["spectype"],
        pulsesequence=spec_dict["pulsesequence"],
        intrument=spec_dict["intrument"],
        probe=spec_dict["probe"],
        datafilename=spec_dict["datafilename"],
        nucleus=spec_dict["nucleus"],
        specfrequency=spec_dict["specfrequency"],
        temperature=spec_dict["temperature"],
        peaks=peaks,
        filename=spec_dict["filename"]
    )
    
    return spectrum

# Usage
data = parse_nmr_assignment_file("exam_CMCse_1_assignments.json")
print(f"Molecule: {data.smiles}")
print(f"Number of atoms: {len(data.allAtomsInfo)}")
print(f"Number of carbons: {len(data.carbonAtomsInfo)}")
print(f"Spectra available: {list(data.spectra.keys())}")
```

------

## Data Validation Rules

### Required Fields

All sections must have:

- `datatype` - Must match section purpose
- `count` - Must equal number of items in data dict
- `data` - Dictionary with string keys "0", "1", "2", etc.

### Atom Information Rules

1. `atom_idx` must equal the dictionary key (as integer)
2. `id` typically equals `atom_idx`
3. `atomNumber` is 1-based (1, 2, 3...) while `atom_idx` is 0-based
4. `numProtons` must be 0-3 for carbon atoms
5. Carbon atoms in `carbonAtomsInfo` must be subset of `allAtomsInfo`

### NMR Spectrum Rules

1. For 1D spectra:
   - `nucleus` is a simple string (e.g., "13C", "1H")
   - `specfrequency` is a float
   - `delta2` in peaks must be 0
2. For 2D spectra:
   - `nucleus` is string representation of list (e.g., "['1H', '13C']")
   - `specfrequency` is string representation of list (e.g., "[400.13, 100.61]")
   - `delta2` in peaks contains actual chemical shift values
3. Peak dictionary keys must be sequential strings: "0", "1", "2", etc.

### Simulated Annealing Rules

1. `startingTemperature` > `endingTemperature` > 0
2. `coolingRate` must be 0 < rate < 1 (typically 0.99-0.999)
3. `numberOfSteps` must be positive integer

------

## Common Patterns and Best Practices

### 1. String Keys for All Data Dictionaries

All data dictionaries use string keys even for numeric indices:

```python
# Correct
{"0": value, "1": value, "2": value}

# Incorrect
{0: value, 1: value, 2: value}
```

### 2. Empty Data Sections

Sections with no data still include full structure:

```python
{
    "nmrAssignments": {
        "datatype": "nmrAssignments",
        "count": 0,
        "data": {}  # Empty dict, not None
    }
}
```

### 3. Type Conversions

Be careful with type conversions for spectrum metadata:

- Nucleus and frequency are strings for 2D, native types for 1D
- Temperature is always string
- Probe is always string (even if numeric)

### 4. MOL File Format

MOL files must include:

- Complete header (3 lines)
- Counts line
- All atom records
- All bond records
- `M  END` line
- `$$$$` terminator
- Preserve all whitespace and newlines

------

## File Size Considerations

- Typical file size: 10KB - 1MB depending on number of peaks
- Large files may have:
  - 50-100+ peaks per 2D spectrum
  - Multiple experiment types (5-10 different spectra)
  - Detailed multiplet analysis
  - Comprehensive integral data

For large datasets, consider:

- Streaming JSON parsing (ijson library)
- Lazy loading of spectra
- Database storage for collections of files

------

## Version Compatibility Notes

This format specification is based on:

- Bruker XWIN-NMR output format
- MDL MOL V2000 format
- Standard NMR experiment nomenclature

When implementing from other sources:

- Ensure consistent use of ppm units for chemical shifts
- Maintain proper dimension ordering (F2, F1 for 2D)
- Preserve original peak intensities (don't auto-normalize)
- Include all required metadata fields even if empty