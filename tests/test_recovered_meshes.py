#!/usr/bin/env python3
"""Regression and cell-complex invariants for the recovered Netgen meshes."""

import hashlib
import sys
from collections import Counter
from pathlib import Path


FIXTURES = Path(__file__).parent / "fixtures" / "netgen"
REGRESSIONS = {
    "furch.mesh": (1132, 5691, 846, "4f6a816b674e91ea"),
    "torus.mesh": (3742, 19874, 3154, "11a823165d0da2a6"),
}


def read_mesh(path, nodes_per_cell, nodes_per_facet):
    tokens = path.read_text(encoding="ascii").split()
    cursor = 0

    def section(width):
        nonlocal cursor
        count = int(tokens[cursor])
        cursor += 1
        rows = [tokens[cursor + i * width:cursor + (i + 1) * width]
                for i in range(count)]
        cursor += count * width
        return rows

    points = section(3)
    cells = section(nodes_per_cell + 1)
    facets = section(nodes_per_facet + 1)
    assert cursor == len(tokens), f"unparsed tokens in {path.name}"
    return points, cells, facets


def check(path, nodes_per_cell, nodes_per_facet):
    points, cells, boundary = read_mesh(path, nodes_per_cell, nodes_per_facet)
    cell_nodes = [tuple(map(int, row[1:])) for row in cells]
    boundary_nodes = [tuple(map(int, row[1:])) for row in boundary]

    assert all(len(set(cell)) == nodes_per_cell for cell in cell_nodes)
    assert all(len(set(facet)) == nodes_per_facet for facet in boundary_nodes)
    assert all(1 <= node <= len(points)
               for entity in cell_nodes + boundary_nodes for node in entity)
    assert len(set(map(frozenset, boundary_nodes))) == len(boundary_nodes)

    if nodes_per_cell == 4:
        local_facets = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    else:
        local_facets = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 2, 4, 6),
                        (1, 3, 5, 7), (0, 1, 4, 5), (2, 3, 6, 7))
    incidence = Counter(frozenset(cell[i] for i in face)
                        for cell in cell_nodes for face in local_facets)
    # Netgen's surface section contains exterior boundaries and labelled
    # material interfaces, so every declaration must be a volume facet but an
    # interface may legitimately have two incident cells.
    assert set(map(frozenset, boundary_nodes)) <= set(incidence)
    assert set(incidence.values()) <= {1, 2}, "mesh is non-manifold"
    return len(points), len(cells), len(boundary)


def main():
    for name, expected in REGRESSIONS.items():
        path = FIXTURES / name
        actual = check(path, 4, 3)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
        assert actual == expected[:3], (name, actual)
        assert digest == expected[3], (name, digest)
    assert check(FIXTURES / "unit_hex.mesh", 8, 4) == (8, 1, 6)
    print("recovered tetrahedral and hexahedral mesh invariants passed")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as error:
        print(f"FAILED: {error}", file=sys.stderr)
        raise
