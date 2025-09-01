# PROTAC Conformer Generator (PCG)

Supporting code for "Data-driven generation of conformational ensembles and ternary complexes for PROTAC and other Chimera systems" (DOI pending).

## Overview

The PROTAC Conformer Generator (PCG) is a CSD Python API workflow for generating conformational ensembles of PROTAC (Proteolysis Targeting Chimera) molecules with optional protein-protein interaction scoring. The algorithm samples the linker region while keeping the warhead and E3 binding group conformations fixed to provided observed (experimental or modeled) structures.

Minimal required inputs:

 - Full PROTAC 3D structure (mol2/sdf; can be built from SMILES via the CSD API if needed)
 - E3 binding group 3D substructure
 - Warhead 3D substructure

Optional:

 - Protein–ligand complex structures for ternary complex modeling

## Algorithm Description

1. Read a 3D model of the full PROTAC.
2. Read 3D models of the warhead and E3 binding group.
3. Perform maximum common substructure (MCS) searches to match these fragments onto the full molecule.
4. Drive matching acyclic torsions in the PROTAC to those in the provided fragment conformations.
5. Lock ring conformations and matched rotatable torsions.
6. Use the CSD Conformer Generator to explore linker conformational space and build an ensemble.
7. (Optional) If protein complexes are supplied, overlay ligands, build ternary models, and score them.

## Ternary Complex Modeling (Optional)

If protein–ligand complexes for the E3 ligase and protein of interest are provided:

 - Overlay their bound ligand fragments onto each PROTAC conformer.
 - Transform the corresponding protein coordinates to assemble ternary models.
 - Score each model:
   - Clash Score: counts heavy atom pairs < 2.0 Å (lower is better).
   - Surface Score: counts heavy atom pairs between 2.0–4.0 Å (higher suggests more putative interface contacts).

## Dependencies & Requirements

 - CSD Python API with a valid "CSD-Discovery" or "CSD-Core + CSD-Conformer-Generator" license.

## Installation

Ensure a working CSD Python API environment. No separate package installation is required beyond dependencies already provided with the API.

## Usage

### Basic

```bash
python protacs_conformer_generator.py protac.mol2 e3_binding_group.mol2 warhead.mol2
```

### Advanced (multithreading, RMSD reference, ternary modeling)

```bash
python protacs_conformer_generator.py protac_model.mol2 e3_binding_group.mol2 warhead.mol2 --number_of_threads 16 --max_conformations 100 --find_closest_to protac_reference.mol2 --proteins_to_align e3_complex.pdb poi_complex.pdb
```

Show all options:

```bash
python protacs_conformer_generator.py -h
```

## Example case

An example folder (example_6HAX) is included containing:
- protac.mol2
- e3_binding_group.mol2
- warhead.mol2
- e3_complex.pdb 
- poi_complex.pdb 

To run the example:
1. Copy protacs_conformer_generator.py and all the files in example_6HAX to the same directory
2. Open a powershell, command-promt, or temrminal and navigate to the directory 
3. Run:
```bash
python protacs_conformer_generator.py protac.mol2 e3_binder.mol2 warhead.mol2 --proteins_to_align e3_complex.pdb poi_complex.pdb --max_conformations 10 --number_of_threads 8
```

### Key Parameters

 - `--number_of_threads`: Parallel threads (default: 1)
 - `--max_conformations`: Maximum conformers to generate
 - `--find_closest_to`: Reference structure for RMSD (defaults to input PROTAC if omitted)
 - `--proteins_to_align`: Two protein–ligand complex files (following same order as E3 binding group and warhead)
 - `--output_file`: Conformers output (default: protacs_conformers.mol2)
 - `--output_ranking`: CSV with scoring/results (default: protacs_full_results.csv)

### Filtering

 - `--clashscore_range` min max
 - `--surfacescore_range` min max

Only models within specified ranges are written.

### Scoring Parameter Overrides
Changing the parameters has not been extensively tested; use at your own risk.

**Clash Score:**

 - `--cs_min_d` (default 1.0)
 - `--cs_max_d` (default 2.0)
 - `--cs_max_height` (default 10.0)

**Surface Score:**

 - `--ss_min_d` (default 2.0)
 - `--ss_max_d` (default 4.0)
 - `--ss_max_height` (default 10.0)


## Output Files

 - **protacs_conformers.mol2**: Generated conformer ensemble
 - **protacs_full_results.csv**: Scores (RMSD, clash, surface, rankings)
 - **closest.mol2**: Conformer closest to reference (if requested)
 - **full_model_#.mol2**: Individual ternary complex models (if proteins supplied)
 - **protac_best_conformer_[metric].mol2**: Best per scoring metric

## Citation

Please cite the associated scientific paper (add DOI when available).

## License

Provided by the Cambridge Crystallographic Data Centre (CCDC).  See license notice in the script header for detailed terms and conditions.

## Support

This script is provided as-is and is not formally supported by CCDC at this time. For questions or issues, please refer to the CSD Python API documentation or contact the authors.

## Authors

- **Jason Cole** (cole@ccdc.cam.ac.uk) - Original implementation
- **Kepa K. Burusco-Goni** (kburusco@ccdc.cam.ac.uk) - Parallel implementation and bug fixes
- **Fabio Montisci** (fmontisci@ccdc.cam.ac.uk) - Documentation and code cleanup
