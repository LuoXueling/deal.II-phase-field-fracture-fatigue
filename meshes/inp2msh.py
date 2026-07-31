#!/usr/bin/env python3
"""
Convert an Abaqus linear-tetrahedron (C3D4) .inp mesh into a gmsh 2.2 .msh
that deal.II's GridIn::read_msh can read, preserving boundary ids.

Why this exists
---------------
deal.II's GridIn::read_abaqus reads a fixed
    1 + GeometryInfo<dim>::vertices_per_cell
entries per element -- 8 nodes in 3D -- and its *SURFACE handler assumes six
quadrilateral faces per cell. It therefore cannot read C3D4 tetrahedra at all
(it fails inside read_ucd with a bogus vertex index). read_msh, by contrast,
handles tetrahedra and carries boundary ids through.

Boundary ids
------------
    *Surface, type=ELEMENT, name=Surf-N   ->   boundary id N

matching the convention the solver's *_boundary.txt files already use, so
those files need no change when switching a model from hexes to tets.

Usage
-----
    ./inp2msh.py model.inp model.msh

Then in the parameter file:

    subsection Project
      set Mesh from = ../meshes/model.msh
    end
    subsection Finite element system
      set Element type = tet
      set Physical dimension = 3
      set Polynomial degree = 1
      set Refine = false            # simplex meshes cannot be refined
    end
"""

import argparse
import re
import sys
from collections import defaultdict

# Abaqus C3D4 face -> local node indices (1-based).
# deal.II matches a boundary triangle to a cell face by node set, so the
# winding within each face does not affect the resulting boundary id.
TET_FACES = {1: (1, 2, 3), 2: (1, 4, 2), 3: (2, 4, 3), 4: (3, 4, 1)}

GMSH_TRIANGLE = 2
GMSH_TETRAHEDRON = 4


class ConversionError(Exception):
    pass


def parse_inp(path):
    """Return (nodes, elements, elsets, surfaces) from an Abaqus .inp."""
    nodes, elements = {}, {}
    elsets = defaultdict(list)
    surfaces = defaultdict(list)
    mode = None
    current = None

    try:
        handle = open(path)
    except OSError as exc:
        raise ConversionError("cannot open %s: %s" % (path, exc))

    with handle as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("**"):
                continue

            if line.startswith("*"):
                upper = line.upper()
                keyword = upper.split(",")[0].strip()
                name = None
                match = re.search(r"(?:ELSET|NSET|NAME)\s*=\s*([^,]+)", line, re.I)
                if match:
                    name = match.group(1).strip()

                if keyword == "*NODE":
                    mode = "node"
                elif keyword == "*ELEMENT":
                    etype = re.search(r"TYPE\s*=\s*([^,]+)", upper)
                    etype = etype.group(1).strip() if etype else "<unspecified>"
                    if etype != "C3D4":
                        raise ConversionError(
                            "line %d: element type %r is not supported; this "
                            "converter handles linear tetrahedra (C3D4) only.\n"
                            "Hexahedral meshes do not need converting -- point "
                            "'Mesh from' at the .inp and leave "
                            "'Element type = hex'." % (lineno, etype))
                    mode = "element"
                elif keyword == "*ELSET":
                    mode = "elset_generate" if "GENERATE" in upper else "elset"
                    current = name
                    elsets.setdefault(current, [])
                elif keyword == "*SURFACE":
                    mode, current = "surface", name
                    surfaces.setdefault(current, [])
                else:
                    mode = None
                continue

            fields = [f.strip() for f in line.split(",") if f.strip()]
            if mode == "node":
                if len(fields) < 4:
                    raise ConversionError(
                        "line %d: expected 'id, x, y, z'" % lineno)
                nodes[int(fields[0])] = tuple(float(v) for v in fields[1:4])
            elif mode == "element":
                if len(fields) < 5:
                    raise ConversionError(
                        "line %d: a C3D4 element needs 4 nodes, got %d"
                        % (lineno, len(fields) - 1))
                elements[int(fields[0])] = [int(v) for v in fields[1:5]]
            elif mode == "elset":
                elsets[current] += [int(v) for v in fields]
            elif mode == "elset_generate":
                first, last, step = (int(v) for v in fields[:3])
                elsets[current] += list(range(first, last + 1, step))
            elif mode == "surface":
                if len(fields) >= 2 and fields[1].upper().startswith("S"):
                    surfaces[current].append((fields[0], int(fields[1][1:])))

    if not nodes:
        raise ConversionError("no *Node block found in " + path)
    if not elements:
        raise ConversionError("no C3D4 *Element block found in " + path)
    return nodes, elements, elsets, surfaces


def boundary_triangles(elements, elsets, surfaces):
    """Expand *Surface definitions into (boundary_id, [n1, n2, n3])."""
    triangles = []
    for name, members in sorted(surfaces.items()):
        match = re.search(r"(\d+)", name or "")
        if not match:
            raise ConversionError(
                "surface %r has no number in its name; boundary ids are taken "
                "from names of the form 'Surf-N'. Rename it in Abaqus, or the "
                "boundary condition for it would be silently lost." % name)
        bid = int(match.group(1))
        for elset_name, face_no in members:
            if face_no not in TET_FACES:
                raise ConversionError(
                    "surface %r references face S%d; a tetrahedron only has "
                    "faces S1-S4." % (name, face_no))
            if elset_name not in elsets:
                raise ConversionError(
                    "surface %r references unknown element set %r"
                    % (name, elset_name))
            for eid in elsets[elset_name]:
                if eid not in elements:
                    raise ConversionError(
                        "surface %r references unknown element %d"
                        % (name, eid))
                conn = elements[eid]
                triangles.append(
                    (bid, [conn[i - 1] for i in TET_FACES[face_no]]))
    return triangles


def write_msh(path, nodes, elements, triangles):
    with open(path, "w") as out:
        out.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        out.write("$Nodes\n%d\n" % len(nodes))
        for nid in sorted(nodes):
            out.write("%d %.17g %.17g %.17g\n" % ((nid,) + nodes[nid]))
        out.write("$EndNodes\n")
        out.write("$Elements\n%d\n" % (len(triangles) + len(elements)))
        index = 0
        # Boundary triangles first; both tags carry the boundary id.
        for bid, tri in triangles:
            index += 1
            out.write("%d %d 2 %d %d %s\n"
                      % (index, GMSH_TRIANGLE, bid, bid,
                         " ".join(str(v) for v in tri)))
        for eid in sorted(elements):
            index += 1
            out.write("%d %d 2 0 0 %s\n"
                      % (index, GMSH_TETRAHEDRON,
                         " ".join(str(v) for v in elements[eid])))
        out.write("$EndElements\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert an Abaqus C3D4 tetrahedral mesh to gmsh 2.2 .msh "
                    "for deal.II, preserving *Surface boundary ids.",
        epilog="Boundary ids come from surface names: 'Surf-3' -> id 3.")
    parser.add_argument("input", help="Abaqus .inp file (C3D4 elements)")
    parser.add_argument("output", help="gmsh .msh file to write")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="suppress the summary")
    args = parser.parse_args()

    try:
        nodes, elements, elsets, surfaces = parse_inp(args.input)
        triangles = boundary_triangles(elements, elsets, surfaces)
    except ConversionError as exc:
        sys.exit("error: %s" % exc)

    write_msh(args.output, nodes, elements, triangles)

    if not args.quiet:
        ids = sorted({bid for bid, _ in triangles})
        print("wrote %s" % args.output)
        print("  %d nodes, %d tetrahedra, %d boundary triangles"
              % (len(nodes), len(elements), len(triangles)))
        print("  boundary ids: %s"
              % (", ".join(str(i) for i in ids) if ids else "none"))
        if not ids:
            print("  warning: no *Surface blocks found, so the mesh carries no "
                  "boundary ids and no boundary condition can be applied.")


if __name__ == "__main__":
    main()
