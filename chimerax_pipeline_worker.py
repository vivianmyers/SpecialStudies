"""ChimeraX worker for hydrogenation and bond-table export.

This script must run inside ChimeraX. Use run_pipeline.py as the public CLI.
"""

from __future__ import annotations

import argparse
import csv
import json
from math import acos, degrees, sqrt
from pathlib import Path

import numpy as np

from chimerax.atomic import AtomicStructure
from chimerax.core.commands import run
from chimerax.hbonds import find_hbonds, rec_angle_slop, rec_dist_slop


WATER_NAMES = {"HOH", "WAT", "DOD"}

ATOM_COLUMNS = [
    "atom_id",
    "model",
    "serial",
    "chain",
    "residue_name",
    "residue_number",
    "insertion_code",
    "atom_name",
    "element",
    "x",
    "y",
    "z",
    "is_hydrogen",
]

COVALENT_COLUMNS = [
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

HYDROGEN_COLUMNS = [
    "model",
    "donor_serial",
    "donor_chain",
    "donor_residue_name",
    "donor_residue_number",
    "donor_insertion_code",
    "donor_atom_name",
    "donor_element",
    "hydrogen_serial",
    "hydrogen_chain",
    "hydrogen_residue_name",
    "hydrogen_residue_number",
    "hydrogen_insertion_code",
    "hydrogen_atom_name",
    "hydrogen_element",
    "acceptor_serial",
    "acceptor_chain",
    "acceptor_residue_name",
    "acceptor_residue_number",
    "acceptor_insertion_code",
    "acceptor_atom_name",
    "acceptor_element",
    "donor_acceptor_distance_A",
    "hydrogen_acceptor_distance_A",
    "donor_hydrogen_acceptor_angle_deg",
    "same_residue",
    "same_chain",
    "involves_water",
    "criteria",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job", required=True, type=Path)
    return parser.parse_args()


def command_path(path: Path) -> str:
    """Quote a filesystem path for use in a ChimeraX text command."""
    return json.dumps(str(path))


def model_label(structure: AtomicStructure) -> str:
    return ".".join(str(part) for part in structure.id)


def atom_id(structure: AtomicStructure, atom) -> str:
    return f"{model_label(structure)}:{atom.serial_number}"


def distance(coord1, coord2) -> float:
    return float(np.linalg.norm(np.asarray(coord1) - np.asarray(coord2)))


def residue_sort_value(value):
    try:
        return (0, int(value))
    except (TypeError, ValueError):
        return (1, str(value))


def atom_fields(prefix, atom) -> dict[str, object]:
    if atom is None:
        return {
            f"{prefix}_serial": "",
            f"{prefix}_chain": "",
            f"{prefix}_residue_name": "",
            f"{prefix}_residue_number": "",
            f"{prefix}_insertion_code": "",
            f"{prefix}_atom_name": "",
            f"{prefix}_element": "",
        }

    residue = atom.residue
    return {
        f"{prefix}_serial": getattr(atom, "serial_number", ""),
        f"{prefix}_chain": residue.chain_id,
        f"{prefix}_residue_name": residue.name,
        f"{prefix}_residue_number": residue.number,
        f"{prefix}_insertion_code": getattr(residue, "insertion_code", ""),
        f"{prefix}_atom_name": atom.name,
        f"{prefix}_element": atom.element.name,
    }


def covalent_atom_fields(prefix, atom) -> dict[str, object]:
    residue = atom.residue
    return {
        f"{prefix}_serial": getattr(atom, "serial_number", ""),
        f"{prefix}_chain": residue.chain_id,
        f"{prefix}_residue_name": residue.name,
        f"{prefix}_residue_number": residue.number,
        f"{prefix}_insertion_code": getattr(residue, "insertion_code", ""),
        f"{prefix}_name": atom.name,
        f"{prefix}_element": atom.element.name,
    }


def donor_hydrogen(donor, acceptor):
    hydrogens = [
        neighbor for neighbor in donor.neighbors if neighbor.element.number == 1
    ]
    if not hydrogens:
        return None
    return min(
        hydrogens,
        key=lambda hydrogen: distance(hydrogen.scene_coord, acceptor.scene_coord),
    )


def dha_angle(donor, hydrogen, acceptor):
    if hydrogen is None:
        return None

    h_coord = np.asarray(hydrogen.scene_coord, dtype=float)
    dh = np.asarray(donor.scene_coord, dtype=float) - h_coord
    ha = np.asarray(acceptor.scene_coord, dtype=float) - h_coord
    denominator = np.linalg.norm(dh) * np.linalg.norm(ha)
    if denominator == 0:
        return None

    cosine = float(np.dot(dh, ha) / denominator)
    cosine = max(-1.0, min(1.0, cosine))
    return degrees(acos(cosine))


def export_atoms(structures, output_path: Path) -> int:
    rows: list[dict[str, object]] = []
    for structure in structures:
        for atom in structure.atoms:
            residue = atom.residue
            x, y, z = atom.scene_coord
            rows.append(
                {
                    "atom_id": atom_id(structure, atom),
                    "model": model_label(structure),
                    "serial": atom.serial_number,
                    "chain": residue.chain_id,
                    "residue_name": residue.name,
                    "residue_number": residue.number,
                    "insertion_code": getattr(residue, "insertion_code", ""),
                    "atom_name": atom.name,
                    "element": atom.element.name,
                    "x": f"{float(x):.4f}",
                    "y": f"{float(y):.4f}",
                    "z": f"{float(z):.4f}",
                    "is_hydrogen": atom.element.number == 1,
                }
            )

    rows.sort(key=lambda row: (str(row["model"]), int(row["serial"])))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=ATOM_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def export_covalent_bonds(structures, output_path: Path) -> int:
    rows: list[dict[str, object]] = []
    for structure in structures:
        for bond in structure.bonds:
            atom1, atom2 = bond.atoms
            delta = atom1.scene_coord - atom2.scene_coord
            length = sqrt(sum(float(component) ** 2 for component in delta))

            row: dict[str, object] = {"model": model_label(structure)}
            row.update(covalent_atom_fields("atom1", atom1))
            row.update(covalent_atom_fields("atom2", atom2))
            row.update(
                {
                    "bond_length_A": f"{length:.4f}",
                    "same_residue": atom1.residue is atom2.residue,
                    "contains_hydrogen": (
                        atom1.element.number == 1 or atom2.element.number == 1
                    ),
                }
            )
            rows.append(row)

    rows.sort(
        key=lambda row: (
            str(row["model"]),
            str(row["atom1_chain"]),
            residue_sort_value(row["atom1_residue_number"]),
            str(row["atom1_name"]),
            str(row["atom2_chain"]),
            residue_sort_value(row["atom2_residue_number"]),
            str(row["atom2_name"]),
        )
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=COVALENT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def export_hydrogen_bonds(
    session,
    structures,
    output_path: Path,
    *,
    relaxed_criteria: bool,
    include_intraresidue: bool,
) -> tuple[int, int]:
    distance_slop = rec_dist_slop if relaxed_criteria else 0.0
    angle_slop = rec_angle_slop if relaxed_criteria else 0.0
    rows: list[dict[str, object]] = []

    for structure in structures:
        pairs = find_hbonds(
            session,
            [structure],
            inter_model=False,
            intra_model=True,
            inter_submodel=False,
            dist_slop=distance_slop,
            angle_slop=angle_slop,
            status=True,
        )
        if not include_intraresidue:
            pairs = [
                (donor, acceptor)
                for donor, acceptor in pairs
                if donor.residue is not acceptor.residue
            ]

        for donor, acceptor in pairs:
            hydrogen = donor_hydrogen(donor, acceptor)
            da_distance = distance(donor.scene_coord, acceptor.scene_coord)
            ha_distance = (
                distance(hydrogen.scene_coord, acceptor.scene_coord)
                if hydrogen is not None
                else None
            )
            angle = dha_angle(donor, hydrogen, acceptor)

            row: dict[str, object] = {"model": model_label(structure)}
            row.update(atom_fields("donor", donor))
            row.update(atom_fields("hydrogen", hydrogen))
            row.update(atom_fields("acceptor", acceptor))
            row.update(
                {
                    "donor_acceptor_distance_A": f"{da_distance:.4f}",
                    "hydrogen_acceptor_distance_A": (
                        f"{ha_distance:.4f}" if ha_distance is not None else ""
                    ),
                    "donor_hydrogen_acceptor_angle_deg": (
                        f"{angle:.2f}" if angle is not None else ""
                    ),
                    "same_residue": donor.residue is acceptor.residue,
                    "same_chain": donor.residue.chain_id == acceptor.residue.chain_id,
                    "involves_water": (
                        donor.residue.name in WATER_NAMES
                        or acceptor.residue.name in WATER_NAMES
                    ),
                    "criteria": "relaxed" if relaxed_criteria else "strict",
                }
            )
            rows.append(row)

    rows.sort(
        key=lambda row: (
            str(row["model"]),
            str(row["donor_chain"]),
            residue_sort_value(row["donor_residue_number"]),
            str(row["donor_atom_name"]),
            str(row["acceptor_chain"]),
            residue_sort_value(row["acceptor_residue_number"]),
            str(row["acceptor_atom_name"]),
        )
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=HYDROGEN_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    without_explicit_hydrogen = sum(not row["hydrogen_atom_name"] for row in rows)
    return len(rows), without_explicit_hydrogen


def main(session) -> None:
    args = parse_args()
    job = json.loads(args.job.read_text(encoding="utf-8"))
    input_path = Path(job["input_path"])
    structure_output = Path(job["structure_output"])
    atom_output = Path(job["atom_output"])
    covalent_output = Path(job["covalent_output"])
    hydrogen_output = Path(job["hydrogen_output"])

    run(session, f"open {command_path(input_path)}")
    structures = session.models.list(type=AtomicStructure)
    if not structures:
        raise RuntimeError(f"No atomic structure was read from {input_path}")

    run(session, "addh")
    structures = session.models.list(type=AtomicStructure)

    structure_output.parent.mkdir(parents=True, exist_ok=True)
    run(session, f"save {command_path(structure_output)} format mmcif")

    atom_count = export_atoms(structures, atom_output)
    covalent_count = export_covalent_bonds(structures, covalent_output)
    hydrogen_count, without_hydrogen = export_hydrogen_bonds(
        session,
        structures,
        hydrogen_output,
        relaxed_criteria=bool(job["relaxed_hbonds"]),
        include_intraresidue=bool(job["include_intraresidue"]),
    )
    session.logger.info(
        f"Pipeline complete for {input_path}: {atom_count} atoms, "
        f"{covalent_count} covalent bonds, "
        f"{hydrogen_count} hydrogen bonds, {without_hydrogen} hydrogen bonds "
        "without an explicit donor hydrogen."
    )


main(session)
