#!/usr/bin/env python3
"""Inspect and curate PDB mmCIF structures with UCSF ChimeraX."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Any


SUPPORTED_SUFFIXES = {".cif", ".mmcif"}
OUTPUT_ROOT = Path("pipeline_output")
COMPONENT_DISTANCE_A = 4.0


def csv_values(value: str | None) -> list[str]:
    """Split a comma-separated selection, preserving item order."""
    if not value:
        return []
    return list(dict.fromkeys(item.strip() for item in value.split(",") if item.strip()))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Interactively inspect and curate one PDB mmCIF structure with ChimeraX."
        )
    )
    parser.add_argument("input", type=Path, help="A .cif or .mmcif file.")
    return parser.parse_args(sys.argv[1:] if argv is None else argv)


def input_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    return Path(name).stem


def validate_input(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise ValueError(f"Input file does not exist: {path}")
    suffix_target = Path(resolved.stem) if resolved.suffix.lower() == ".gz" else resolved
    if suffix_target.suffix.lower() not in SUPPORTED_SUFFIXES:
        allowed = ", ".join(sorted(SUPPORTED_SUFFIXES))
        raise ValueError(f"Unsupported input type for {path}; expected: {allowed}")
    return resolved


def find_chimerax() -> Path:
    candidates: list[Path] = []
    environment_path = os.environ.get("CHIMERAX")
    if environment_path:
        candidates.append(Path(environment_path).expanduser())
    candidates.extend(
        [
            Path("/Applications/ChimeraX.app/Contents/MacOS/ChimeraX"),
            Path.home() / "Applications/ChimeraX.app/Contents/MacOS/ChimeraX",
        ]
    )
    candidates.extend(
        sorted(
            Path("/Applications").glob("ChimeraX*.app/Contents/MacOS/ChimeraX"),
            reverse=True,
        )
    )
    candidates.extend(
        sorted(
            (Path.home() / "Applications").glob(
                "ChimeraX*.app/Contents/MacOS/ChimeraX"
            ),
            reverse=True,
        )
    )
    for command_name in ("ChimeraX", "chimerax"):
        executable = shutil.which(command_name)
        if executable:
            candidates.append(Path(executable))
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise RuntimeError(
        "Could not find ChimeraX. Install it in /Applications or set the "
        "CHIMERAX environment variable to its executable."
    )


def run_chimerax(
    chimerax: Path, worker: Path, job_path: Path, *, show_output: bool = False
) -> None:
    script_argument = shlex.join([str(worker), "--job", str(job_path)])
    completed = subprocess.run(
        [str(chimerax), "--nogui", "--exit", "--script", script_argument],
        check=False,
        stdout=None if show_output else subprocess.PIPE,
        stderr=None if show_output else subprocess.STDOUT,
        text=not show_output,
    )
    if completed.returncode != 0:
        if not show_output and completed.stdout:
            print(completed.stdout, file=sys.stderr, end="")
        raise RuntimeError(f"ChimeraX failed with exit code {completed.returncode}.")


def inspect_structure(
    input_path: Path, chimerax: Path, worker: Path
) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="curation_inspect_") as temporary_name:
        temporary_dir = Path(temporary_name)
        output_path = temporary_dir / "inspection.json"
        job_path = temporary_dir / "job.json"
        job_path.write_text(
            json.dumps(
                {
                    "action": "inspect",
                    "input_path": str(input_path),
                    "inspection_output": str(output_path),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        run_chimerax(chimerax, worker, job_path)
        if not output_path.is_file():
            raise RuntimeError("ChimeraX did not create the inspection manifest.")
        return json.loads(output_path.read_text(encoding="utf-8"))


def model_data(inspection: dict[str, Any], model_number: int) -> dict[str, Any]:
    models = inspection.get("models", [])
    if model_number < 1 or model_number > len(models):
        raise ValueError(
            f"Model {model_number} is unavailable; choose 1 through {len(models)}."
        )
    return models[model_number - 1]


def render_table(headers: list[str], rows: list[list[object]]) -> str:
    text_rows = [[str(value) for value in row] for row in rows]
    widths = [len(header) for header in headers]
    for row in text_rows:
        for index, value in enumerate(row):
            widths[index] = max(widths[index], len(value))
    lines = ["  ".join(h.ljust(widths[i]) for i, h in enumerate(headers))]
    lines.append("  ".join("-" * width for width in widths))
    lines.extend(
        "  ".join(value.ljust(widths[i]) for i, value in enumerate(row))
        for row in text_rows
    )
    return "\n".join(lines)


def print_inspection(
    inspection: dict[str, Any], model_number: int, selected_chains: list[str] | None = None
) -> None:
    selected = set(selected_chains or [])
    model = model_data(inspection, model_number)
    print(f"\nInput: {inspection['input_path']}")
    print(
        f"Model {model_number} of {len(inspection['models'])} "
        f"(ChimeraX model #{model['model_id']})"
    )
    chain_rows = [
        [
            chain["chain_id"],
            chain["class"],
            chain["polymer_residues"],
            chain["nonpolymer_residues"],
            chain["water_residues"],
            chain["atoms"],
            chain.get("description", ""),
        ]
        for chain in model["chains"]
    ]
    print("\nChains")
    print(
        render_table(
            [
                "CHAIN",
                "CLASS",
                "POLYMER",
                "OTHER",
                "WATER",
                "ATOMS",
                "DESCRIPTION",
            ],
            chain_rows,
        )
    )

    component_rows = []
    for component in model["components"]:
        distances = component.get("distances_to_chains", {})
        relevant = [distances[c] for c in selected if c in distances]
        distance = min(relevant) if relevant else component.get("nearest_distance_A")
        component_rows.append(
            [
                component["id"],
                component["class"],
                component["atoms"],
                "" if distance is None else f"{distance:.2f}",
                ",".join(component.get("covalently_attached_chains", [])) or "no",
            ]
        )
    if component_rows:
        print("\nNon-polymer components")
        print(
            render_table(
                ["ID", "CLASS", "ATOMS", "MIN DIST (A)", "COVALENT TO"],
                component_rows,
            )
        )
    waters = model["waters"]
    print(f"\nWaters: {waters['residues']} residues, {waters['atoms']} atoms")
    warnings = model.get("warnings", [])
    if warnings:
        print("\nWarnings")
        for warning in warnings:
            print(f"- {warning}")


def prompt_value(prompt: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{prompt}{suffix}: ").strip()
    return value or default


def validate_chain_ids(model: dict[str, Any], chains: list[str]) -> None:
    available = {chain["chain_id"] for chain in model["chains"]}
    unknown = [chain for chain in chains if chain not in available]
    if unknown:
        raise ValueError(
            "Unknown chain ID(s): "
            + ", ".join(unknown)
            + ". Available chains: "
            + ", ".join(sorted(available))
        )
    polymer_chains = {
        chain["chain_id"]
        for chain in model["chains"]
        if chain["polymer_residues"] > 0
    }
    without_polymer = [chain for chain in chains if chain not in polymer_chains]
    if without_polymer:
        raise ValueError(
            "Chain selection must contain polymer chains; these have no polymer residues: "
            + ", ".join(without_polymer)
        )


def suggested_components(
    model: dict[str, Any], chains: list[str], distance_cutoff: float
) -> list[str]:
    selected = set(chains)
    suggestions: list[str] = []
    for component in model["components"]:
        covalent = set(component.get("covalently_attached_chains", [])) & selected
        distances = component.get("distances_to_chains", {})
        near = any(
            chain in distances and distances[chain] <= distance_cutoff
            for chain in selected
        )
        if covalent or near:
            suggestions.append(component["id"])
    return suggestions


def resolve_model_selection(inspection: dict[str, Any]) -> int:
    model_count = len(inspection["models"])
    if model_count == 1:
        return 1
    value = prompt_value("Select model", "1")
    try:
        model_number = int(value)
    except ValueError as error:
        raise ValueError("Model selection must be a number.") from error
    model_data(inspection, model_number)
    return model_number


def resolve_selection(
    inspection: dict[str, Any], model_number: int
) -> tuple[list[str], list[str]]:
    model = model_data(inspection, model_number)
    defaults = [
        chain["chain_id"]
        for chain in model["chains"]
        if chain["polymer_residues"] > 0
    ]
    chains = csv_values(prompt_value("Select polymer chains", ",".join(defaults)))
    if not chains:
        raise ValueError("At least one polymer chain must be selected.")
    validate_chain_ids(model, chains)

    component_defaults = suggested_components(
        model, chains, COMPONENT_DISTANCE_A
    )
    if model["components"]:
        print_inspection(inspection, model_number, chains)
        print(
            f"\nComponents within {COMPONENT_DISTANCE_A:.1f} A or covalently "
            "attached are suggested by default."
        )
        components = csv_values(
            prompt_value("Components to retain", ",".join(component_defaults))
        )
    else:
        components = []

    available_components = {item["id"] for item in model["components"]}
    unknown_components = [item for item in components if item not in available_components]
    if unknown_components:
        raise ValueError(
            "Unknown component ID(s): " + ", ".join(unknown_components)
        )

    return chains, components


def output_paths(output_root: Path, stem: str) -> dict[str, Path]:
    result_dir = output_root / stem
    return {
        "directory": result_dir,
        "structure": result_dir / f"{stem}_cleaned.cif",
        "session": result_dir / f"{stem}_curated.cxs",
        "atoms": result_dir / f"{stem}_atoms.csv",
        "covalent": result_dir / f"{stem}_covalent_bonds.csv",
        "hydrogen": result_dir / f"{stem}_hydrogen_bonds.csv",
        "report": result_dir / f"{stem}_curation.json",
    }


def curate_structure(
    *,
    input_path: Path,
    chimerax: Path,
    worker: Path,
    inspection: dict[str, Any],
    model_number: int,
    chains: list[str],
    components: list[str],
    overwrite: bool,
) -> dict[str, Path]:
    stem = input_stem(input_path)
    destinations = output_paths(OUTPUT_ROOT.resolve(), stem)
    final_keys = (
        "structure",
        "session",
        "atoms",
        "covalent",
        "hydrogen",
        "report",
    )
    existing = [destinations[key] for key in final_keys if destinations[key].exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Output already exists: "
            + ", ".join(str(path) for path in existing)
        )

    with tempfile.TemporaryDirectory(prefix="curation_run_") as temporary_name:
        temporary_dir = Path(temporary_name)
        temporary_outputs = {
            key: temporary_dir / destinations[key].name for key in final_keys
        }
        job_path = temporary_dir / "job.json"
        job = {
            "action": "curate",
            "input_path": str(input_path),
            "model_number": model_number,
            "selected_chains": chains,
            "kept_components": components,
            "structure_output": str(temporary_outputs["structure"]),
            "session_output": str(temporary_outputs["session"]),
            "atom_output": str(temporary_outputs["atoms"]),
            "covalent_output": str(temporary_outputs["covalent"]),
            "hydrogen_output": str(temporary_outputs["hydrogen"]),
            "report_output": str(temporary_outputs["report"]),
            "relaxed_hbonds": True,
            "include_intraresidue": True,
            "inspection": inspection,
            "protonation": "ChimeraX addh defaults",
            "final_outputs": {
                "cleaned_structure": str(destinations["structure"]),
                "chimerax_session": str(destinations["session"]),
                "atoms": str(destinations["atoms"]),
                "covalent_bonds": str(destinations["covalent"]),
                "hydrogen_bonds": str(destinations["hydrogen"]),
                "curation_report": str(destinations["report"]),
            },
        }
        job_path.write_text(json.dumps(job, indent=2), encoding="utf-8")
        run_chimerax(chimerax, worker, job_path, show_output=True)
        missing = [
            str(temporary_outputs[key])
            for key in final_keys
            if not temporary_outputs[key].is_file()
        ]
        if missing:
            raise RuntimeError(
                "ChimeraX exited without creating expected output(s): "
                + ", ".join(missing)
            )
        destinations["directory"].mkdir(parents=True, exist_ok=True)
        for key in final_keys:
            os.replace(temporary_outputs[key], destinations[key])
    return destinations


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if not sys.stdin.isatty():
            raise RuntimeError("This prototype requires an interactive terminal.")
        input_path = validate_input(args.input)
        chimerax = find_chimerax()
        worker = Path(__file__).with_name("chimerax_pipeline_worker.py").resolve()
        if not worker.is_file():
            raise RuntimeError(f"Missing ChimeraX worker script: {worker}")
        inspection = inspect_structure(input_path, chimerax, worker)

        model_number = resolve_model_selection(inspection)
        print_inspection(inspection, model_number)
        chains, components = resolve_selection(inspection, model_number)
        print("\nCuration selection")
        print(f"- Model: {model_number}")
        print(f"- Polymer chains: {', '.join(chains)}")
        print(f"- Components: {', '.join(components) if components else 'none'}")

        destinations = output_paths(OUTPUT_ROOT.resolve(), input_stem(input_path))
        final_keys = (
            "structure",
            "session",
            "atoms",
            "covalent",
            "hydrogen",
            "report",
        )
        existing = [destinations[key] for key in final_keys if destinations[key].exists()]
        overwrite = False
        if existing:
            print("\nExisting output will be replaced:")
            for path in existing:
                print(f"- {path}")
            replace = prompt_value("Replace existing output?", "no").lower()
            if replace not in {"y", "yes"}:
                print("Cancelled; no output was changed.")
                return 0
            overwrite = True

        answer = prompt_value("Continue?", "yes").lower()
        if answer not in {"y", "yes"}:
            print("Cancelled; no output was changed.")
            return 0

        destinations = curate_structure(
            input_path=input_path,
            chimerax=chimerax,
            worker=worker,
            inspection=inspection,
            model_number=model_number,
            chains=chains,
            components=components,
            overwrite=overwrite,
        )
    except (FileExistsError, OSError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print("\nCreated curated dataset:")
    for key in ("structure", "session", "atoms", "covalent", "hydrogen", "report"):
        print(f"- {destinations[key]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
