"""
symmetry_unique_carbons.py

Given a molecule, determine the minimum connected set of atoms containing
at least one representative of every NMR/graph-symmetry-equivalent carbon
class. Useful for picking a minimal, contiguous fragment that still shows
every chemically-unique carbon environment.

All atom indices used/returned are RDKit atom indices (0-indexed), which
equal (molfile 1-indexed atom number - 1) when the molecule is loaded
directly from a molfile with Chem.MolFromMolFile / Chem.MolFromMolBlock.
"""

from collections import defaultdict
import networkx as nx
from rdkit import Chem


def build_molgraph(mol):
    """Build a NetworkX graph of all heavy atoms in `mol`, using RDKit atom idx as node id."""
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            symbol=atom.GetSymbol(),
            atom_number=atom.GetAtomicNum(),
        )
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), weight=1)
    return G


def get_carbon_symmetry_classes(mol):
    """
    Return a list of symmetry-equivalence classes restricted to carbon atoms.

    Each class is a list of RDKit atom indices that are graph-symmetry
    equivalent (same CanonicalRankAtoms rank) AND all carbon. Singleton
    carbon classes (a carbon with no symmetry partner) are included as
    single-element lists, so callers don't need to special-case them.
    """
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    classes_by_rank = defaultdict(list)
    for idx, rank in enumerate(ranks):
        classes_by_rank[rank].append(idx)

    carbon_classes = []
    for atoms in classes_by_rank.values():
        if all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in atoms):
            carbon_classes.append(atoms)
    return carbon_classes


def greedy_group_steiner_tree(G, classes, weight="weight"):
    """
    Find a minimal connected subgraph of G containing at least one atom
    from each class in `classes` (a "Group Steiner Tree" heuristic).

    Parameters
    ----------
    G : nx.Graph
        Molecular graph (all heavy atoms, bonds as edges)
    classes : list of lists
        Each sublist is a set of interchangeable (symmetry-equivalent)
        atom indices - one class per unique carbon environment.
    weight : str
        Edge weight attribute name (bonds are typically unweighted, i.e. 1)

    Returns
    -------
    tree_nodes : set
        Minimal atom index set forming a connected subgraph, containing
        at least one representative of every class.
    chosen_reps : dict
        Maps class index (position in `classes`) -> the representative
        atom actually selected for that class.
    """
    classes = [list(c) for c in classes]
    remaining = list(range(len(classes)))

    first_class_idx = remaining.pop(0)
    seed = min(classes[first_class_idx])
    tree_nodes = {seed}
    chosen_reps = {first_class_idx: seed}

    while remaining:
        best = None  # (dist, class_idx, atom, path)
        lengths, paths = nx.multi_source_dijkstra(G, sources=tree_nodes, weight=weight)

        for ci in remaining:
            for atom in classes[ci]:
                if atom in tree_nodes:
                    dist, path = 0, [atom]
                else:
                    dist = lengths.get(atom, float("inf"))
                    path = paths.get(atom)
                if path is None:
                    continue
                if best is None or dist < best[0]:
                    best = (dist, ci, atom, path)

        if best is None:
            raise ValueError("Graph is disconnected w.r.t. remaining classes")

        _, ci, atom, path = best
        tree_nodes.update(path)
        chosen_reps[ci] = atom
        remaining.remove(ci)

    return tree_nodes, chosen_reps


def draw_minimal_atom_set(
    mol,
    tree_nodes,
    carbon_nodes,
    connector_nodes,
    output_path,
    size=(700, 700),
    carbon_color=(0.4, 0.85, 0.4),
    connector_color=(1.0, 0.55, 0.15),
    highlight_radius=0.5,
):
    """
    Render the molecule as a 2D PNG, highlighting:
      - carbon_nodes  (the chosen symmetry-class representatives) in green
      - connector_nodes (non-carbon bridge atoms needed for connectivity) in orange

    Parameters
    ----------
    mol : rdkit.Chem.Mol
    tree_nodes, carbon_nodes, connector_nodes : iterables of RDKit atom idx (0-indexed)
    output_path : str or Path
        Where to save the PNG.
    highlight_radius : float
        Radius (in molecule-drawing units) of the highlight circles. RDKit's
        default (~0.3) can look quite small/subtle on larger structures;
        increase this (e.g. 0.6-0.8) for more visible highlighting.
    """
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    mol = Chem.Mol(mol)  # don't mutate caller's molecule

    # Molfiles loaded straight from ChemDraw etc. usually already have 2D
    # coords; only recompute if none are present (e.g. a mol built purely
    # from a molblock with placeholder 0,0,0 coordinates).
    conf_ok = mol.GetNumConformers() > 0 and not all(
        mol.GetConformer().GetAtomPosition(i).x == 0
        and mol.GetConformer().GetAtomPosition(i).y == 0
        for i in range(mol.GetNumAtoms())
    )
    if not conf_ok:
        AllChem.Compute2DCoords(mol)

    highlight_atoms = list(carbon_nodes) + list(connector_nodes)
    highlight_colors = {i: carbon_color for i in carbon_nodes}
    highlight_colors.update({i: connector_color for i in connector_nodes})
    highlight_radii = {i: highlight_radius for i in highlight_atoms}

    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.addAtomIndices = False
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_colors,
        highlightAtomRadii=highlight_radii,
    )
    drawer.FinishDrawing()
    with open(output_path, "wb") as f:
        f.write(drawer.GetDrawingText())
    return output_path


def get_minimal_unique_carbon_set(mol, molgraph=None):
    """
    High-level entry point: given an RDKit mol, return the minimal connected
    atom set spanning all unique (symmetry-distinct) carbon environments.

    Returns
    -------
    tree_nodes : set of RDKit atom indices in the minimal connected subgraph
                 (includes carbon representatives + any non-carbon connectors)
    chosen_reps : dict of class_idx -> representative carbon atom idx
    carbon_classes : the symmetry classes used (list of lists of atom idx)
    """
    if molgraph is None:
        molgraph = build_molgraph(mol)
    carbon_classes = get_carbon_symmetry_classes(mol)
    tree_nodes, chosen_reps = greedy_group_steiner_tree(molgraph, carbon_classes)
    return tree_nodes, chosen_reps, carbon_classes


def main():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        description="Find the minimal connected atom set spanning every "
        "symmetry-unique carbon environment in a molfile."
    )
    parser.add_argument("molfile", type=Path, help="Path to a .mol / .sdf file")
    parser.add_argument(
        "--highlight-radius",
        type=float,
        default=0.5,
        help="Radius of the highlight circles in the output image "
        "(RDKit default look is ~0.3; larger molecules/denser structures "
        "often need a smaller value, simple molecules can use a larger "
        "one for visibility). Default: 0.5",
    )
    args = parser.parse_args()

    molfile_path: Path = args.molfile
    if not molfile_path.is_file():
        raise FileNotFoundError(f"No such file: {molfile_path}")

    mol = Chem.MolFromMolFile(str(molfile_path), sanitize=True)
    if mol is None:
        raise ValueError(f"RDKit could not parse molfile: {molfile_path}")

    molgraph = build_molgraph(mol)
    tree_nodes, chosen_reps, carbon_classes = get_minimal_unique_carbon_set(mol, molgraph)

    # Report everything using 1-indexed atom numbers (matches molfile numbering)
    tree_1indexed = sorted(n + 1 for n in tree_nodes)
    carbon_nodes = sorted(n + 1 for n in tree_nodes if molgraph.nodes[n]["symbol"] == "C")
    connector_nodes = sorted(
        n + 1 for n in tree_nodes if molgraph.nodes[n]["symbol"] != "C"
    )

    print(f"File: {molfile_path}")
    print(f"Total atoms in molecule: {mol.GetNumAtoms()}")
    print(f"Total carbon atoms: {sum(len(c) for c in carbon_classes)}")
    print(f"Number of unique carbon symmetry classes: {len(carbon_classes)}")
    print()
    print("Symmetry classes (1-indexed atom numbers):")
    for i, cls in enumerate(carbon_classes):
        cls_1indexed = sorted(a + 1 for a in cls)
        rep_1indexed = chosen_reps[i] + 1
        print(f"  Class {i}: {cls_1indexed}  (representative chosen: {rep_1indexed})")
    print()
    print(f"Minimal connected atom set (size {len(tree_1indexed)}):")
    print(f"  All atoms:       {tree_1indexed}")
    print(f"  Carbon reps:     {carbon_nodes}  ({len(carbon_nodes)})")
    print(f"  Connector atoms: {connector_nodes}  ({len(connector_nodes)})")
    for n in connector_nodes:
        symbol = molgraph.nodes[n - 1]["symbol"]
        print(f"    - atom {n}: {symbol}")

    # 0-indexed sets for drawing
    carbon_nodes_0idx = [n for n in tree_nodes if molgraph.nodes[n]["symbol"] == "C"]
    connector_nodes_0idx = [n for n in tree_nodes if molgraph.nodes[n]["symbol"] != "C"]

    image_path = molfile_path.with_suffix("").with_name(
        molfile_path.stem + "_highlighted.png"
    )
    draw_minimal_atom_set(
        mol,
        tree_nodes,
        carbon_nodes_0idx,
        connector_nodes_0idx,
        image_path,
        highlight_radius=args.highlight_radius,
    )
    print()
    print(f"Highlighted structure image saved to: {image_path}")


if __name__ == "__main__":
    main()
