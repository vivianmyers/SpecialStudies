# Protein data-curation pipeline

## Requirements

- Python 3.10 +
- UCSF ChimeraX 1.12 (from `/Applications` on macOS)

## Run the pipeline

```bash
python3 run_pipeline.py cif/7HVP.cif
```

After the outputs are created, the user may open the
`.cxs` session in the ChimeraX GUI.

The prototype uses these fixed defaults:

- Output directory: `pipeline_output/`
- ChimeraX recommended H-bond tolerances
- Intra-residue H-bonds included
- Solvent handling

## Outputs

Each input receives its own directory under `pipeline_output/`:

```text
pipeline_output/7HVP/
├── 7HVP_cleaned.cif
├── 7HVP_curated.cxs
├── 7HVP_atoms.csv
├── 7HVP_covalent_bonds.csv
├── 7HVP_hydrogen_bonds.csv
└── 7HVP_curation.txt
```

The cleaned mmCIF and atom table contain the explicit hydrogen atoms added by
ChimeraX. Hydrogen bonds are noncovalent pseudobonds rather than mmCIF covalent
bonds: the ChimeraX session preserves them for visualization, and the hydrogen-
bond CSV exports their donor, hydrogen, acceptor, distance, angle, and criteria
for later analysis.

The plain-text report records all selections, warnings, removed contents,
processing choices, and output counts, including the number of ChimeraX H-bond
pseudobonds. Missing residues or atoms are reported but are not repaired in
this prototype. Alternate locations are preserved and reported.

## Pipeline order

```text
PDB mmCIF -> inspect -> select model/chains/components/solvent handling -> remove
unselected material -> optionally delete all solvent -> add hydrogens -> create
ChimeraX H-bond pseudobonds -> export and report
```

## ChimeraX commands

The pipeline launches ChimeraX in noninteractive mode to inspect
the input and to perform the selected curation. The launcher
has this form:

```bash
"<ChimeraX executable>" --nogui --exit \
  --script "<project>/chimerax_pipeline_worker.py --job <temporary-job.json>"
```

## Citations

- Meng EC, Goddard TD, Pettersen EF, et al. UCSF ChimeraX: Tools for structure
  building and analysis. *Protein Science*. 2023;32(11):e4792.
- Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids
  Research*. 2000;28:235-242. https://doi.org/10.1093/nar/28.1.235
