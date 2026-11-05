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

The pipeline uses these fixed defaults:

- Output directory: `pipeline_output/`

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
pseudobonds. Missing residues or atoms are reported but are not repaired.

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

### Adding hydrogens

After the selected model has been curated and any excluded chains,
components, or solvent molecules have been removed, the pipeline adds
hydrogen atoms with ChimeraX's standard protonation behavior. For model `#1`,
the exact ChimeraX command is:

```chimerax
addh #1
```

ChimeraX uses its default settings to determine protonation states and hydrogen placement from the
retained structure. Hydrogenation occurs before H-bond detection.

The added hydrogens become part of the curated atomic model. They are saved in
the `*_cleaned.cif` structure and `*_curated.cxs` session and are included in
the `*_atoms.csv` and `*_covalent_bonds.csv` tables. The curation report records
the atom count before hydrogenation and the final exported atom count.

### H-bonds

```chimerax
hbonds #1 interModel false interSubmodel false relax true intraRes true makePseudobonds true
```

Hydrogen bonds are identified after ChimeraX adds explicit hydrogen atoms. The pipeline uses ChimeraX's donor and acceptor typing and its recommended relaxed distance and angle tolerances. Only bonds within the selected model are considered; bonds between models or submodels are excluded. Intra-residue hydrogen bonds are included by default.

ChimeraX stores the detected H-bonds as noncovalent pseudobonds in the
`*_curated.cxs` session so they can be displayed in GUI. These are
not covalent bonds, so they are not written to the mmCIF bond table. The pipeline writes each
donor--acceptor pair to `*_hydrogen_bonds.csv`.

For each pair, the export records the donor, acceptor, and associated donor
hydrogen; the donor--acceptor and hydrogen--acceptor distances; and the
donor--hydrogen--acceptor angle. 

## Citations

- Meng EC, Goddard TD, Pettersen EF, et al. UCSF ChimeraX: Tools for structure
  building and analysis. *Protein Science*. 2023;32(11):e4792.
- Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids
  Research*. 2000;28:235-242. https://doi.org/10.1093/nar/28.1.235
