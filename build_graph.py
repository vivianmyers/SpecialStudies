#!/usr/bin/env python3
"""Build a NetworkX graph from one pipeline result directory."""

import argparse
import csv
from pathlib import Path

import networkx as nx


def atom_id(model: str, serial: str) -> str:
    return f"{model}:{serial}"


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as input_file:
        yield from csv.DictReader(input_file)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path)
    args = parser.parse_args()

    result_dir = args.result_dir
    name = result_dir.name
    graph = nx.MultiDiGraph()

    for row in rows(result_dir / f"{name}_atoms.csv"):
        node = row.pop("atom_id")
        row["serial"] = int(row["serial"])
        row["residue_number"] = int(row["residue_number"])
        row["x"] = float(row["x"])
        row["y"] = float(row["y"])
        row["z"] = float(row["z"])
        row["is_hydrogen"] = row["is_hydrogen"] == "True"
        graph.add_node(node, **row)

    for row in rows(result_dir / f"{name}_covalent_bonds.csv"):
        atom1 = atom_id(row["model"], row["atom1_serial"])
        atom2 = atom_id(row["model"], row["atom2_serial"])
        attributes = {
            "edge_type": "covalent",
            "distance_A": float(row["bond_length_A"]),
            "same_residue": row["same_residue"] == "True",
            "contains_hydrogen": row["contains_hydrogen"] == "True",
        }
        graph.add_edge(atom1, atom2, **attributes)
        graph.add_edge(atom2, atom1, **attributes)

    for row in rows(result_dir / f"{name}_hydrogen_bonds.csv"):
        attributes = {
            "edge_type": "hydrogen_bond",
            "donor_acceptor_distance_A": float(
                row["donor_acceptor_distance_A"]
            ),
            "same_residue": row["same_residue"] == "True",
            "same_chain": row["same_chain"] == "True",
            "involves_water": row["involves_water"] == "True",
            "criteria": row["criteria"],
        }
        if row["hydrogen_serial"]:
            attributes.update(
                {
                    "hydrogen_id": atom_id(
                        row["model"], row["hydrogen_serial"]
                    ),
                    "hydrogen_acceptor_distance_A": float(
                        row["hydrogen_acceptor_distance_A"]
                    ),
                    "angle_deg": float(
                        row["donor_hydrogen_acceptor_angle_deg"]
                    ),
                }
            )
        graph.add_edge(
            atom_id(row["model"], row["donor_serial"]),
            atom_id(row["model"], row["acceptor_serial"]),
            **attributes,
        )

    output_path = result_dir / f"{name}_graph.graphml"
    nx.write_graphml(graph, output_path)
    print(f"Created {output_path}")
    print(f"{graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")


if __name__ == "__main__":
    main()
