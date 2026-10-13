import csv
from math import acos, degrees
from pathlib import Path

import numpy as np

from chimerax.atomic import AtomicStructure
from chimerax.hbonds import find_hbonds, rec_angle_slop, rec_dist_slop


# ------------------------------ Configuration ------------------------------

OUTPUT_PATH = Path.home() / "Documents/Classes/SpecialStudies" / "hydrogen_bonds.csv"

# True matches the normal ChimeraX "hbonds" command default by using the
# recommended distance and angle tolerances.  False applies strict criteria.
RELAX_CRITERIA = True

# ChimeraX includes intra-residue H-bonds by default.  Set this to False if
# only interactions between different residues are wanted.
INCLUDE_INTRARESIDUE = True


# ------------------------------- Utilities ---------------------------------

def model_label(structure):
    return "#" + ".".join(str(part) for part in structure.id)


def distance(coord1, coord2):
    return float(np.linalg.norm(np.asarray(coord1) - np.asarray(coord2)))


def donor_hydrogen(donor, acceptor):
    """Return the donor-bound H nearest the acceptor, as ChimeraX reports it."""
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
    """Return the D-H-A angle in degrees; 180 degrees is ideally linear."""
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


def atom_fields(prefix, atom):
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


# -------------------------------- Export ------------------------------------

structure = session.models.list(type=AtomicStructure)[0]

distance_slop = rec_dist_slop if RELAX_CRITERIA else 0.0
angle_slop = rec_angle_slop if RELAX_CRITERIA else 0.0

hbond_pairs = find_hbonds(
    session,
    [structure],
    inter_model=False,
    intra_model=True,
    inter_submodel=False,
    dist_slop=distance_slop,
    angle_slop=angle_slop,
    status=True,
)

if not INCLUDE_INTRARESIDUE:
    hbond_pairs = [
        (donor, acceptor)
        for donor, acceptor in hbond_pairs
        if donor.residue is not acceptor.residue
    ]

fieldnames = [
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

water_names = {"HOH", "WAT", "DOD"}
rows = []

for donor, acceptor in hbond_pairs:
    hydrogen = donor_hydrogen(donor, acceptor)
    da_distance = distance(donor.scene_coord, acceptor.scene_coord)
    ha_distance = (
        distance(hydrogen.scene_coord, acceptor.scene_coord)
        if hydrogen is not None
        else None
    )
    angle = dha_angle(donor, hydrogen, acceptor)

    row = {"model": model_label(structure)}
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
                donor.residue.name in water_names
                or acceptor.residue.name in water_names
            ),
            "criteria": "relaxed" if RELAX_CRITERIA else "strict",
        }
    )
    rows.append(row)

rows.sort(
    key=lambda row: (
        str(row["donor_chain"]),
        int(row["donor_residue_number"]),
        str(row["donor_atom_name"]),
        str(row["acceptor_chain"]),
        int(row["acceptor_residue_number"]),
        str(row["acceptor_atom_name"]),
    )
)

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as output_file:
    writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

without_explicit_h = sum(not row["hydrogen_atom_name"] for row in rows)
session.logger.info(
    f"Exported {len(rows)} hydrogen bonds for {model_label(structure)} to "
    f"{OUTPUT_PATH}. Interactions without an explicit donor hydrogen: "
    f"{without_explicit_h}."
)
