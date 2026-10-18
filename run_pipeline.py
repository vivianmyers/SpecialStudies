#!/usr/bin/env python3
"""Command-line launcher for the ChimeraX atom and bond export pipeline."""

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


SUPPORTED_SUFFIXES = {".cif", ".mmcif"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Add hydrogens to mmCIF structures with ChimeraX and export "
            "atom, covalent-bond, and hydrogen-bond tables."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="One or more .cif or .mmcif structure files.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("pipeline_output"),
        help="Output root directory (default: ./pipeline_output).",
    )
    parser.add_argument(
        "--strict-hbonds",
        action="store_true",
        help="Use strict ChimeraX hydrogen-bond criteria instead of recommended tolerances.",
    )
    parser.add_argument(
        "--exclude-intraresidue-hbonds",
        action="store_true",
        help="Exclude hydrogen bonds whose donor and acceptor are in the same residue.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing output files for the supplied structures.",
    )
    parser.add_argument(
        "--chimerax",
        type=Path,
        help="Path to the ChimeraX executable; normally detected automatically on macOS.",
    )
    return parser.parse_args()


def input_stem(path: Path) -> str:
    """Return a stable structure name, including for files ending in .gz."""
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    return Path(name).stem


def validate_input(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise ValueError(f"Input file does not exist: {path}")

    suffix_target = resolved
    if resolved.suffix.lower() == ".gz":
        suffix_target = Path(resolved.stem)
    if suffix_target.suffix.lower() not in SUPPORTED_SUFFIXES:
        allowed = ", ".join(sorted(SUPPORTED_SUFFIXES))
        raise ValueError(f"Unsupported input type for {path}; expected one of: {allowed}")
    return resolved


def find_chimerax(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit.expanduser())

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
        "Could not find ChimeraX. Supply its executable with --chimerax, for "
        "example --chimerax /Applications/ChimeraX-1.12.app/Contents/MacOS/ChimeraX"
    )


def output_paths(output_root: Path, stem: str) -> dict[str, Path]:
    result_dir = output_root / stem
    return {
        "directory": result_dir,
        "structure": result_dir / f"{stem}_H.cif",
        "atoms": result_dir / f"{stem}_atoms.csv",
        "covalent": result_dir / f"{stem}_covalent_bonds.csv",
        "hydrogen": result_dir / f"{stem}_hydrogen_bonds.csv",
    }


def process_structure(
    *,
    input_path: Path,
    destinations: dict[str, Path],
    chimerax: Path,
    worker: Path,
    relaxed_hbonds: bool,
    include_intraresidue: bool,
    overwrite: bool,
) -> None:
    final_files = [
        destinations["structure"],
        destinations["atoms"],
        destinations["covalent"],
        destinations["hydrogen"],
    ]
    existing = [path for path in final_files if path.exists()]
    if existing and not overwrite:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(
            f"Output already exists: {names}. Use --overwrite to replace it."
        )

    with tempfile.TemporaryDirectory(prefix="bond_pipeline_") as temporary_name:
        temporary_dir = Path(temporary_name)
        temporary_outputs = {
            "structure": temporary_dir / destinations["structure"].name,
            "atoms": temporary_dir / destinations["atoms"].name,
            "covalent": temporary_dir / destinations["covalent"].name,
            "hydrogen": temporary_dir / destinations["hydrogen"].name,
        }
        job_path = temporary_dir / "job.json"
        job = {
            "input_path": str(input_path),
            "structure_output": str(temporary_outputs["structure"]),
            "atom_output": str(temporary_outputs["atoms"]),
            "covalent_output": str(temporary_outputs["covalent"]),
            "hydrogen_output": str(temporary_outputs["hydrogen"]),
            "relaxed_hbonds": relaxed_hbonds,
            "include_intraresidue": include_intraresidue,
        }
        job_path.write_text(json.dumps(job, indent=2), encoding="utf-8")

        script_argument = shlex.join([str(worker), "--job", str(job_path)])
        command = [
            str(chimerax),
            "--nogui",
            "--exit",
            "--script",
            script_argument,
        ]
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            raise RuntimeError(
                f"ChimeraX failed for {input_path} with exit code "
                f"{completed.returncode}."
            )

        missing = [
            str(path) for path in temporary_outputs.values() if not path.is_file()
        ]
        if missing:
            raise RuntimeError(
                "ChimeraX exited without creating expected output(s): "
                + ", ".join(missing)
            )

        destinations["directory"].mkdir(parents=True, exist_ok=True)
        for key in ("structure", "atoms", "covalent", "hydrogen"):
            os.replace(temporary_outputs[key], destinations[key])


def main() -> int:
    args = parse_args()
    try:
        inputs = [validate_input(path) for path in args.inputs]
        stems = [input_stem(path) for path in inputs]
        duplicate_stems = sorted({stem for stem in stems if stems.count(stem) > 1})
        if duplicate_stems:
            raise ValueError(
                "Input filenames must have unique stems; duplicates: "
                + ", ".join(duplicate_stems)
            )

        chimerax = find_chimerax(args.chimerax)
    except (RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    output_root = args.output_dir.expanduser().resolve()
    worker = Path(__file__).with_name("chimerax_pipeline_worker.py").resolve()
    if not worker.is_file():
        print(f"error: missing ChimeraX worker script: {worker}", file=sys.stderr)
        return 2

    failures: list[tuple[Path, str]] = []
    for input_path, stem in zip(inputs, stems):
        destinations = output_paths(output_root, stem)
        print(f"Processing {input_path}", flush=True)
        try:
            process_structure(
                input_path=input_path,
                destinations=destinations,
                chimerax=chimerax,
                worker=worker,
                relaxed_hbonds=not args.strict_hbonds,
                include_intraresidue=not args.exclude_intraresidue_hbonds,
                overwrite=args.overwrite,
            )
        except (FileExistsError, OSError, RuntimeError) as error:
            failures.append((input_path, str(error)))
            print(f"FAILED: {error}", file=sys.stderr)
            continue

        print(f"Created {destinations['structure']}")
        print(f"Created {destinations['atoms']}")
        print(f"Created {destinations['covalent']}")
        print(f"Created {destinations['hydrogen']}")

    if failures:
        print(f"\n{len(failures)} input(s) failed:", file=sys.stderr)
        for input_path, message in failures:
            print(f"  {input_path}: {message}", file=sys.stderr)
        return 1

    print(f"\nCompleted {len(inputs)} structure(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
