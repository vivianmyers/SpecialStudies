import csv
from pathlib import Path
from math import sqrt

from chimerax.atomic import AtomicStructure # pyright: ignore[reportMissingImports]


output_path = Path.home() / "Documents/Classes/SpecialStudies" / "covalent_bonds.csv"

structures = session.models.list(type=AtomicStructure) # pyright: ignore[reportUndefinedVariable]

if not structures:
    raise RuntimeError("No atomic structure is currently open.")

columns = [
    "model",
    "atom1_serial",
    "atom1_chain",
    "atom1_residue_name",
    "atom1_residue_number",
    "atom1_insertion_code",
    "atom1_name",
    "atom1_element",
    "atom2_serial",
    "atom2_chain",
    "atom2_residue_name",
    "atom2_residue_number",
    "atom2_insertion_code",
    "atom2_name",
    "atom2_element",
    "bond_length_A",
    "same_residue",
    "contains_hydrogen",
]

rows = []

for structure in structures:
    model_id = ".".join(str(value) for value in structure.id)

    for bond in structure.bonds:
        atom1, atom2 = bond.atoms
        residue1 = atom1.residue
        residue2 = atom2.residue

        dx = float(atom1.coord[0] - atom2.coord[0])
        dy = float(atom1.coord[1] - atom2.coord[1])
        dz = float(atom1.coord[2] - atom2.coord[2])
        length = sqrt(dx * dx + dy * dy + dz * dz)

        same_residue = (
            residue1.chain_id == residue2.chain_id
            and residue1.number == residue2.number
            and residue1.name == residue2.name
        )

        contains_hydrogen = (
            atom1.element.name == "H" or atom2.element.name == "H"
        )

        rows.append([
            model_id,
            getattr(atom1, "serial_number", ""),
            residue1.chain_id,
            residue1.name,
            residue1.number,
            getattr(residue1, "insertion_code", ""),
            atom1.name,
            atom1.element.name,
            getattr(atom2, "serial_number", ""),
            residue2.chain_id,
            residue2.name,
            residue2.number,
            getattr(residue2, "insertion_code", ""),
            atom2.name,
            atom2.element.name,
            f"{length:.4f}",
            same_residue,
            contains_hydrogen,
        ])

with output_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(columns)
    writer.writerows(rows)

session.logger.info( # pyright: ignore[reportUndefinedVariable]
    f"Exported {len(rows)} covalent bonds to {output_path}"
)