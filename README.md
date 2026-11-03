# Protein data-curation pipeline

## Requirements

- Python 3.10 +
- UCSF ChimeraX 1.12 (detected automatically from `/Applications` on macOS)

## Run the pipeline

```bash
python3 run_pipeline.py cif/7HVP.cif
```

After the outputs are created, the pipeline asks whether to open the curated
`.cxs` session in the ChimeraX GUI. The default is `no`.

The prototype uses these fixed defaults:

- Output directory: `pipeline_output/`
- ChimeraX recommended H-bond tolerances
- Intra-residue H-bonds included
- Solvent handling: keep all by default, with explicit delete-all choice
- ChimeraX detected automatically, or from the `CHIMERAX` environment variable

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

The pipeline launches ChimeraX in noninteractive mode twice: once to inspect
the input and once to perform the selected curation. The launcher invocation
has this form:

```bash
"<ChimeraX executable>" --nogui --exit \
  --script "<project>/chimerax_pipeline_worker.py --job <temporary-job.json>"
```

Angle brackets below mark values filled in by the pipeline; they are not typed
literally. Paths are quoted when sent to ChimeraX.

### 1. Inspect the input

```chimerax
open "<absolute-input-path>"
```

The chain, component, water, atom, residue, and warning tables are then built
through the ChimeraX Python API. No additional ChimeraX text command is used
for inspection.

### 2. Reopen and curate the selected model

```chimerax
open "<absolute-input-path>"
```

Unselected models are closed through the ChimeraX Python API. Unselected
polymer residues and non-polymer components are removed with the structure
Python API, not with a ChimeraX text command.

If the user chooses to delete all solvent, the pipeline sends exactly:

```chimerax
delete solvent
```

If the user keeps solvent, that command is not sent.

### 3. Add hydrogens

```chimerax
addh #<selected-ChimeraX-model-ID>
```

For the usual model ID `1`, the literal command is:

```chimerax
addh #1
```

### 4. Calculate hydrogen bonds and create pseudobonds

The current fixed settings send:

```chimerax
hbonds #<selected-ChimeraX-model-ID> interModel false interSubmodel false relax true intraRes true makePseudobonds true
```

For model ID `1`, the literal command is:

```chimerax
hbonds #1 interModel false interSubmodel false relax true intraRes true makePseudobonds true
```

The H-bond CSV is subsequently produced through ChimeraX's `find_hbonds()`
Python API using the same model restrictions, relaxed distance and angle
tolerances, and inclusion of intra-residue H-bonds.

### 5. Set the session to atom view

Before saving the session, the pipeline hides the ribbon/cartoon and shows all
atoms and their bonds:

```chimerax
hide #<selected-ChimeraX-model-ID> cartoons
show #<selected-ChimeraX-model-ID> atoms
```

For model ID `1`, the literal commands are:

```chimerax
hide #1 cartoons
show #1 atoms
```

The H-bond pseudobonds created in the previous step remain displayed.

### 6. Save the curated structure and ChimeraX session

```chimerax
save "<temporary-output-path>/<input-stem>_cleaned.cif" format mmcif
save "<temporary-output-path>/<input-stem>_curated.cxs"
```

The pipeline moves these successfully completed temporary files into
`pipeline_output/<input-stem>/`. Atom, covalent-bond, H-bond, and text-report
exports are written through the ChimeraX/Python APIs rather than additional
ChimeraX text commands.

### 7. Optionally open the completed session in the GUI

If the user answers `yes` to the final GUI prompt, the pipeline launches a new
ChimeraX process without `--nogui` or `--exit`:

```bash
"<ChimeraX executable>" "<output-path>/<input-stem>_curated.cxs"
```

The GUI process is detached from the pipeline so the terminal command can
finish while ChimeraX remains open. If the user answers `no`, no GUI process
is launched.

## Citations

- Meng EC, Goddard TD, Pettersen EF, et al. UCSF ChimeraX: Tools for structure
  building and analysis. *Protein Science*. 2023;32(11):e4792.
- Berman HM, Westbrook J, Feng Z, et al. The Protein Data Bank. *Nucleic Acids
  Research*. 2000;28:235-242. https://doi.org/10.1093/nar/28.1.235
