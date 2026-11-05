"""ChimeraX pipeline worker script"""


from __future__ import annotations

import argparse
import csv
import json
from math import acos, degrees, sqrt
from pathlib import Path

import numpy as np

from chimerax.atomic import AtomicStructure, Residue
from chimerax.core.commands import run
from chimerax.hbonds import find_hbonds, rec_angle_slop, rec_dist_slop


WATER_NAMES = {"HOH", "WAT", "DOD"}

ATOM_COLUMNS = [
    "atom_id",
    "stable_atom_id",
    "model",
    "serial",
    "chain",
    "residue_name",
    "residue_number",
    "insertion_code",
    "atom_name",
    "element",
    "alt_loc",
    "occupancy",
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


def display_chain_id(chain_id: str) -> str:
    """Give blank chain IDs a selectable, printable representation."""
    return chain_id.strip() or "_"


def residue_id(residue) -> str:
    insertion_code = str(getattr(residue, "insertion_code", "")).strip()
    position = f"{residue.number}{insertion_code}"
    return f"{display_chain_id(residue.chain_id)}:{position}:{residue.name}"


def stable_atom_id(structure: AtomicStructure, atom) -> str:
    residue = atom.residue
    alt_loc = str(getattr(atom, "alt_loc", "")).strip() or "."
    return ":".join(
        [
            model_label(structure),
            residue_id(residue),
            atom.name,
            alt_loc,
        ]
    )


def is_water(residue) -> bool:
    categories = {atom.structure_category for atom in residue.atoms}
    return residue.name in Residue.water_res_names or "solvent" in categories


def component_class(residue) -> str:
    categories = {atom.structure_category for atom in residue.atoms}
    if "ions" in categories:
        return "ion"
    if "ligand" in categories:
        return "ligand"
    return "non-polymer"


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
        f"{prefix}_chain": display_chain_id(residue.chain_id),
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
        f"{prefix}_chain": display_chain_id(residue.chain_id),
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
                    "stable_atom_id": stable_atom_id(structure, atom),
                    "model": model_label(structure),
                    "serial": atom.serial_number,
                    "chain": display_chain_id(residue.chain_id),
                    "residue_name": residue.name,
                    "residue_number": residue.number,
                    "insertion_code": getattr(residue, "insertion_code", ""),
                    "atom_name": atom.name,
                    "element": atom.element.name,
                    "alt_loc": str(getattr(atom, "alt_loc", "")).strip(),
                    "occupancy": f"{float(atom.occupancy):.3f}",
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


def chain_class(polymer_types: set[int]) -> str:
    if polymer_types == {Residue.PT_AMINO}:
        return "protein"
    if polymer_types == {Residue.PT_NUCLEIC}:
        return "nucleic acid"
    if polymer_types:
        return "mixed polymer"
    return "non-polymer"


def inspect_atomic_structure(structure: AtomicStructure) -> dict[str, object]:
    residues = list(structure.residues)
    chain_groups: dict[str, list] = {}
    for residue in residues:
        chain_id = display_chain_id(residue.chain_id)
        chain_groups.setdefault(chain_id, []).append(residue)

    chain_descriptions = {
        display_chain_id(chain.chain_id): str(getattr(chain, "description", "") or "")
        for chain in structure.chains
    }
    chains: list[dict[str, object]] = []
    for chain_id, chain_residues in sorted(chain_groups.items()):
        polymer = [
            residue
            for residue in chain_residues
            if residue.polymer_type != Residue.PT_NONE
        ]
        nonpolymer = [
            residue
            for residue in chain_residues
            if residue.polymer_type == Residue.PT_NONE and not is_water(residue)
        ]
        waters = [residue for residue in chain_residues if is_water(residue)]
        chains.append(
            {
                "chain_id": chain_id,
                "class": chain_class({residue.polymer_type for residue in polymer}),
                "polymer_residues": len(polymer),
                "nonpolymer_residues": len(nonpolymer),
                "water_residues": len(waters),
                "atoms": sum(len(residue.atoms) for residue in chain_residues),
                "description": chain_descriptions.get(chain_id, ""),
            }
        )

    components: list[dict[str, object]] = []
    water_residues = []
    for residue in residues:
        if residue.polymer_type != Residue.PT_NONE:
            continue
        if is_water(residue):
            water_residues.append(residue)
            continue
        components.append(
            {
                "id": residue_id(residue),
                "chain_id": display_chain_id(residue.chain_id),
                "residue_name": residue.name,
                "residue_number": residue.number,
                "insertion_code": str(
                    getattr(residue, "insertion_code", "")
                ).strip(),
                "class": component_class(residue),
                "atoms": len(residue.atoms),
            }
        )
    components.sort(
        key=lambda item: (
            str(item["chain_id"]),
            residue_sort_value(item["residue_number"]),
            str(item["residue_name"]),
        )
    )

    warnings: list[str] = []
    alternate_residues = {
        residue_id(atom.residue)
        for atom in structure.atoms
        if getattr(atom, "num_alt_locs", 0) > 1
    }
    if alternate_residues:
        preview = ", ".join(sorted(alternate_residues)[:8])
        suffix = " ..." if len(alternate_residues) > 8 else ""
        warnings.append(
            f"Alternate locations occur in {len(alternate_residues)} residue(s): "
            f"{preview}{suffix}"
        )

    missing_heavy = []
    for residue in residues:
        if residue.polymer_type == Residue.PT_NONE:
            continue
        try:
            if residue.is_missing_heavy_template_atoms():
                missing_heavy.append(residue_id(residue))
        except RuntimeError:
            continue
    if missing_heavy:
        preview = ", ".join(missing_heavy[:8])
        suffix = " ..." if len(missing_heavy) > 8 else ""
        warnings.append(
            f"Missing heavy atoms detected in {len(missing_heavy)} residue(s): "
            f"{preview}{suffix}"
        )

    for chain in structure.chains:
        missing_count = max(0, chain.num_residues - chain.num_existing_residues)
        if missing_count:
            warnings.append(
                f"Chain {display_chain_id(chain.chain_id)} has "
                f"{missing_count} unresolved residue position(s)."
            )
    if len(structure.bonds) == 0:
        warnings.append("No covalent bonds were found or inferred on import.")
    if getattr(structure, "num_coordsets", 1) > 1:
        warnings.append(
            f"The model contains {structure.num_coordsets} coordinate sets; "
            "the active coordinate set will be used for tabular exports."
        )

    return {
        "model_id": model_label(structure),
        "name": structure.name,
        "atoms": len(structure.atoms),
        "residues": len(residues),
        "covalent_bonds": len(structure.bonds),
        "chains": chains,
        "components": components,
        "waters": {
            "residues": len(water_residues),
            "atoms": sum(len(residue.atoms) for residue in water_residues),
        },
        "warnings": warnings,
    }


def inspection_manifest(input_path: Path, structures) -> dict[str, object]:
    models = [inspect_atomic_structure(structure) for structure in structures]
    if len(models) > 1:
        for model in models:
            model["warnings"].insert(
                0,
                f"The file contains {len(models)} atomic models; select one explicitly.",
            )
    return {
        "schema_version": 1,
        "input_path": str(input_path),
        "models": models,
    }


def format_curation_report(report: dict[str, object]) -> str:
    """Render the exported audit information as a readable text report."""
    lines = [
        "Protein curation report",
        "=======================",
        "",
        f"Input file: {report['input_path']}",
        f"Selected model: {report['selected_model']} "
        f"(ChimeraX #{report['selected_model_id']})",
        "Selected chains: " + ", ".join(report["selected_chains"]),
        "Kept components: "
        + (", ".join(report["kept_components"]) or "none"),
        f"Solvent handling: {str(report['solvent_handling']).replace('_', ' ')}",
        f"Protonation: {report['protonation']}",
        f"Hydrogen-bond criteria: {report['hbond_criteria']}",
        "Include intra-residue hydrogen bonds: "
        + ("yes" if report["include_intraresidue_hbonds"] else "no"),
        "",
        "Warnings",
        "--------",
    ]

    warnings = report["warnings"]
    if warnings:
        lines.extend(f"- {warning}" for warning in warnings)
    else:
        lines.append("- none")

    removed = report["removed"]
    removed_components = removed.get("components", [])
    lines.extend(
        [
            "",
            "Removed",
            "-------",
            f"- Polymer residues: {removed.get('polymer_residues', 0)}",
            "- Components: "
            + (", ".join(removed_components) if removed_components else "none"),
            f"- Solvent residues: {removed.get('solvent_residues', 0)}",
            f"- Solvent atoms: {removed.get('solvent_atoms', 0)}",
            "",
            "Counts",
            "------",
        ]
    )
    for name, value in report["counts"].items():
        lines.append(f"- {name.replace('_', ' ').capitalize()}: {value}")

    lines.extend(["", "Output files", "------------"])
    for name, value in report["outputs"].items():
        lines.append(f"- {name.replace('_', ' ').capitalize()}: {value}")
    lines.append("")
    return "\n".join(lines)


def curate_model(structure: AtomicStructure, job: dict) -> dict[str, object]:
    selected_chains = set(job["selected_chains"])
    kept_components = set(job["kept_components"])
    removed = {
        "polymer_residues": 0,
        "components": [],
        "solvent_residues": 0,
        "solvent_atoms": 0,
    }

    residues_to_remove = []
    for residue in list(structure.residues):
        if residue.polymer_type != Residue.PT_NONE:
            keep = display_chain_id(residue.chain_id) in selected_chains
            if not keep:
                removed["polymer_residues"] += 1
        elif is_water(residue):
            # Solvent is handled separately with ChimeraX's native
            # ``delete solvent`` command when the user explicitly requests it.
            keep = True
        else:
            keep = residue_id(residue) in kept_components
            if not keep:
                removed["components"].append(residue_id(residue))
        if not keep:
            residues_to_remove.append(residue)

    for residue in residues_to_remove:
        structure.delete_residue(residue)
    if len(structure.atoms) == 0:
        raise RuntimeError("Curation removed every atom from the selected model.")
    return removed


def main(session) -> None:
    args = parse_args()
    job = json.loads(args.job.read_text(encoding="utf-8"))
    input_path = Path(job["input_path"])
    run(session, f"open {command_path(input_path)}")
    structures = session.models.list(type=AtomicStructure)
    if not structures:
        raise RuntimeError(f"No atomic structure was read from {input_path}")

    action = job.get("action", "curate")
    if action == "inspect":
        output_path = Path(job["inspection_output"])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(inspection_manifest(input_path, structures), indent=2),
            encoding="utf-8",
        )
        session.logger.info(f"Wrote structure inspection to {output_path}")
        return
    if action != "curate":
        raise RuntimeError(f"Unsupported worker action: {action}")

    model_number = int(job["model_number"])
    if model_number < 1 or model_number > len(structures):
        raise RuntimeError(
            f"Model {model_number} is unavailable; found {len(structures)} model(s)."
        )
    structure = structures[model_number - 1]
    other_structures = [item for item in structures if item is not structure]
    if other_structures:
        session.models.close(other_structures)

    structure_output = Path(job["structure_output"])
    session_output = Path(job["session_output"])
    atom_output = Path(job["atom_output"])
    covalent_output = Path(job["covalent_output"])
    hydrogen_output = Path(job["hydrogen_output"])
    report_output = Path(job["report_output"])

    before_atoms = len(structure.atoms)
    before_residues = len(structure.residues)
    removed = curate_model(structure, job)
    if job.get("delete_solvent", False):
        atoms_before_solvent_deletion = len(structure.atoms)
        residues_before_solvent_deletion = len(structure.residues)
        run(session, "delete solvent")
        removed["solvent_residues"] = (
            residues_before_solvent_deletion - len(structure.residues)
        )
        removed["solvent_atoms"] = (
            atoms_before_solvent_deletion - len(structure.atoms)
        )
    retained_before_hydrogenation = len(structure.atoms)

    run(session, f"addh #{model_label(structure)}")

    relax_value = "true" if job["relaxed_hbonds"] else "false"
    intraresidue_value = "true" if job["include_intraresidue"] else "false"
    run(
        session,
        f"hbonds #{model_label(structure)} interModel false "
        f"interSubmodel false relax {relax_value} "
        f"intraRes {intraresidue_value} makePseudobonds true",
    )
    hbond_group = structure.pseudobond_group("hydrogen bonds", create_type=None)
    hbond_pseudobond_count = (
        len(hbond_group.pseudobonds) if hbond_group is not None else 0
    )

    run(session, f"hide #{model_label(structure)} cartoons")
    run(session, f"show #{model_label(structure)} atoms")

    structure_output.parent.mkdir(parents=True, exist_ok=True)
    run(session, f"save {command_path(structure_output)} format mmcif")
    run(session, f"save {command_path(session_output)}")

    selected_structures = [structure]
    atom_count = export_atoms(selected_structures, atom_output)
    covalent_count = export_covalent_bonds(selected_structures, covalent_output)
    hydrogen_count, without_hydrogen = export_hydrogen_bonds(
        session,
        selected_structures,
        hydrogen_output,
        relaxed_criteria=bool(job["relaxed_hbonds"]),
        include_intraresidue=bool(job["include_intraresidue"]),
    )

    inspection = job.get("inspection", {})
    inspected_models = inspection.get("models", [])
    model_warnings = (
        inspected_models[model_number - 1].get("warnings", [])
        if model_number <= len(inspected_models)
        else []
    )
    report = {
        "schema_version": 1,
        "input_path": str(input_path),
        "selected_model": model_number,
        "selected_model_id": model_label(structure),
        "selected_chains": job["selected_chains"],
        "kept_components": job["kept_components"],
        "solvent_handling": (
            "delete_all" if job.get("delete_solvent", False) else "keep_all"
        ),
        "protonation": job.get("protonation", "ChimeraX addh defaults"),
        "hbond_criteria": "relaxed" if job["relaxed_hbonds"] else "strict",
        "include_intraresidue_hbonds": bool(job["include_intraresidue"]),
        "warnings": model_warnings,
        "removed": removed,
        "counts": {
            "input_atoms": before_atoms,
            "input_residues": before_residues,
            "atoms_before_hydrogenation": retained_before_hydrogenation,
            "exported_atoms": atom_count,
            "exported_covalent_bonds": covalent_count,
            "exported_hydrogen_bonds": hydrogen_count,
            "chimerax_hbond_pseudobonds": hbond_pseudobond_count,
            "hydrogen_bonds_without_explicit_hydrogen": without_hydrogen,
        },
        "outputs": job.get(
            "final_outputs",
            {
                "cleaned_structure": str(structure_output),
                "chimerax_session": str(session_output),
                "atoms": str(atom_output),
                "covalent_bonds": str(covalent_output),
                "hydrogen_bonds": str(hydrogen_output),
                "curation_report": str(report_output),
            },
        ),
    }
    report_output.parent.mkdir(parents=True, exist_ok=True)
    report_output.write_text(format_curation_report(report), encoding="utf-8")
    session.logger.info(
        f"Pipeline complete for {input_path}: {atom_count} atoms, "
        f"{covalent_count} covalent bonds, "
        f"{hydrogen_count} hydrogen bonds, {without_hydrogen} hydrogen bonds "
        "without an explicit donor hydrogen."
    )


main(session)
