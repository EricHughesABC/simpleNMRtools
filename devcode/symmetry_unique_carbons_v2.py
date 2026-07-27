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


def _render_molecule_png_bytes(mol, highlight_atoms, highlight_colors, highlight_radii, size):
    """
    Internal helper: render `mol` to PNG bytes with the given atom
    highlighting. Shared by draw_minimal_atom_set (single plot) and
    draw_combination_grid (one tile per branch combination).
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
    return drawer.GetDrawingText()


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
    highlight_atoms = list(carbon_nodes) + list(connector_nodes)
    highlight_colors = {i: carbon_color for i in carbon_nodes}
    highlight_colors.update({i: connector_color for i in connector_nodes})
    highlight_radii = {i: highlight_radius for i in highlight_atoms}

    png_bytes = _render_molecule_png_bytes(
        mol, highlight_atoms, highlight_colors, highlight_radii, size
    )
    with open(output_path, "wb") as f:
        f.write(png_bytes)
    return output_path


def is_trivial_leaf_class(molgraph, cls):
    """
    A class is 'trivial' if every member is a terminal (degree-1) atom --
    e.g. the 3 methyls of a tert-butyl, or the 3 F's of a CF3. These arise
    from free rotation of a substituent and say nothing about the molecule's
    overall branch/unit count, so they should be excluded when inferring it.
    """
    return all(molgraph.degree(atom) == 1 for atom in cls)


def detect_symmetry_unit_count(molgraph, carbon_classes, min_class_size=2, verbose=True):
    """
    Infer the number of repeating symmetric units (branches) in the molecule
    from the carbon symmetry classes, ignoring trivial terminal-leaf classes.

    The heuristic: the most common class size among non-trivial,
    multi-member classes is taken as the unit count N. This works well when
    most atoms in the molecule are "one skeletal position per branch" --
    which is the typical case -- since minority classes (extra local
    sub-symmetry, or exceptional shared/hub atoms) don't dominate the mode.

    Returns
    -------
    n_units : int or None
        The detected number of symmetric units, or None if no non-trivial
        multi-member class was found (i.e. the molecule shows no repeating
        branch symmetry worth enumerating).
    size_histogram : collections.Counter
        size -> number of non-trivial classes of that size, for transparency/
        sanity-checking the automatic detection.
    """
    from collections import Counter

    sizes = [
        len(c)
        for c in carbon_classes
        if len(c) >= min_class_size and not is_trivial_leaf_class(molgraph, c)
    ]
    size_histogram = Counter(sizes)

    if not size_histogram:
        if verbose:
            print("No non-trivial multi-member carbon classes found -- "
                  "molecule shows no repeating branch symmetry to enumerate.")
        return None, size_histogram

    n_units, n_units_count = size_histogram.most_common(1)[0]

    if verbose:
        print("Non-trivial class size histogram (size: count):")
        for size, count in sorted(size_histogram.items()):
            marker = "  <-- detected unit count" if size == n_units else ""
            print(f"  {size}: {count}{marker}")
        # flag if the detection looks ambiguous (no clear dominant size)
        total = sum(size_histogram.values())
        if n_units_count < total / 2:
            print(
                f"  Note: size {n_units} is only the most common by a small "
                "margin -- inspect the histogram above and consider passing "
                "n_units explicitly if this looks wrong."
            )

    return n_units, size_histogram


def enumerate_symmetric_representative_sets(molgraph, carbon_classes, n_units=None, anchor_atoms=None):
    """
    Discover the N distinct, self-consistent ways to pick one representative
    per carbon symmetry class, when the molecule is built from N repeating
    symmetric "branches"/units (e.g. a molecule with 6-fold symmetry).

    For each of the N branches, every class contributes whichever of its
    members is graph-closest (shortest bond path) to that branch's anchor
    atom. Atoms that are equally close to several branches (e.g. shared
    hub/core atoms, or classes smaller than N) are naturally reused across
    those branches -- this is the "singleton overlap" behavior: such an atom
    simply wins the "closest member" competition for more than one branch.

    Parameters
    ----------
    molgraph : nx.Graph
        Full heavy-atom molecular graph (weight=1 bonds).
    carbon_classes : list of lists
        Symmetry-equivalence classes (from get_carbon_symmetry_classes).
    n_units : int, optional
        Number of repeating symmetric branches/units in the molecule. If
        not provided, this is auto-detected via detect_symmetry_unit_count.
    anchor_atoms : list of int, optional
        The N atom indices to use as branch anchors -- one per branch. If
        not provided, the first class found with exactly `n_units` members
        is used automatically.

    Returns
    -------
    combinations : list of dict
        One dict per distinct branch combination, each mapping
        class_idx -> chosen representative atom idx for that branch,
        after deduplicating any branches that produced an identical result.
    """
    if n_units is None:
        n_units, _ = detect_symmetry_unit_count(molgraph, carbon_classes)
        if n_units is None or n_units < 2:
            raise ValueError(
                "Could not auto-detect a repeating unit count > 1 for this "
                "molecule -- it may not have branch symmetry, or pass "
                "n_units explicitly."
            )

    if anchor_atoms is None:
        # Prefer a non-trivial class as the anchor -- a trivial (all-terminal-
        # leaf) class of the right size, e.g. a gem-dimethyl pair, would
        # anchor the branch partitioning to a chemically meaningless axis
        # instead of the molecule's real structural symmetry.
        non_trivial_anchor_classes = [
            c
            for c in carbon_classes
            if len(c) == n_units and not is_trivial_leaf_class(molgraph, c)
        ]
        if non_trivial_anchor_classes:
            anchor_atoms = sorted(non_trivial_anchor_classes[0])
        else:
            trivial_anchor_classes = [c for c in carbon_classes if len(c) == n_units]
            if not trivial_anchor_classes:
                raise ValueError(
                    f"No carbon symmetry class with exactly {n_units} members "
                    "found to use as branch anchors -- pass anchor_atoms explicitly."
                )
            print(
                f"Warning: only trivial (terminal-leaf) classes of size {n_units} "
                "were found -- falling back to one of them as the anchor, but "
                "this may not reflect the molecule's real structural symmetry. "
                "Consider passing anchor_atoms explicitly."
            )
            anchor_atoms = sorted(trivial_anchor_classes[0])
    if len(anchor_atoms) != n_units:
        raise ValueError(f"Expected {n_units} anchor atoms, got {len(anchor_atoms)}")

    # distances[k] = {atom_idx: shortest bond-path length to anchor k}
    distances = [
        nx.single_source_shortest_path_length(molgraph, anchor) for anchor in anchor_atoms
    ]

    raw_combinations = []
    for k in range(n_units):
        chosen_reps = {}
        for ci, cls in enumerate(carbon_classes):
            # representative for this branch = class member closest to anchor k
            # (ties broken by lowest atom index for reproducibility)
            best_atom = min(cls, key=lambda a: (distances[k].get(a, float("inf")), a))
            chosen_reps[ci] = best_atom
        raw_combinations.append(chosen_reps)

    # deduplicate branches that ended up identical (can happen with partial/
    # nested symmetry, e.g. a class shared across only some branches)
    seen = set()
    combinations = []
    for combo in raw_combinations:
        key = tuple(sorted(combo.items()))
        if key not in seen:
            seen.add(key)
            combinations.append(combo)

    return combinations


def build_combination_dataframe(mol, molgraph, carbon_classes, chosen_reps):
    """
    Given one branch combination (class_idx -> representative atom idx,
    from enumerate_symmetric_representative_sets), build the connected
    minimal atom set (adding connector atoms if needed) and return it as a
    pandas DataFrame, one row per atom in the final set.
    """
    import pandas as pd

    terminals = list(chosen_reps.values())

    # check whether the chosen representatives are already connected;
    # if not, patch in connector atoms via the group Steiner tree machinery
    induced = molgraph.subgraph(terminals)
    if nx.is_connected(induced):
        tree_nodes = set(terminals)
    else:
        from networkx.algorithms.approximation.steinertree import steiner_tree

        T = steiner_tree(molgraph, terminals, weight="weight")
        tree_nodes = set(T.nodes())

    rows = []
    rep_to_class = {atom: ci for ci, atom in chosen_reps.items()}
    for atom_idx in sorted(tree_nodes):
        symbol = molgraph.nodes[atom_idx]["symbol"]
        rows.append(
            {
                "atom_idx": atom_idx,
                "atom_number_1indexed": atom_idx + 1,
                "symbol": symbol,
                "role": "carbon_representative" if atom_idx in rep_to_class else "connector",
                "class_idx": rep_to_class.get(atom_idx, None),
            }
        )
    return pd.DataFrame(rows)


def draw_combination_grid(
    mol,
    molgraph,
    carbon_classes,
    combinations,
    output_path,
    subimg_size=(350, 350),
    carbon_color=(0.4, 0.85, 0.4),
    connector_color=(1.0, 0.55, 0.15),
    highlight_radius=0.5,
    mols_per_row=None,
    label_prefix="Combination",
):
    """
    Render every branch combination (from enumerate_symmetric_representative_sets)
    as its own highlighted panel, arranged in a grid in a single image.

    Parameters
    ----------
    combinations : list of dict
        As returned by enumerate_symmetric_representative_sets.
    output_path : str or Path
        Where to save the composite grid PNG.
    subimg_size : (int, int)
        Size of each individual panel, in pixels.
    mols_per_row : int, optional
        Number of columns in the grid. Defaults to ceil(sqrt(n)) for a
        roughly square layout.
    """
    import io
    import math

    from PIL import Image, ImageDraw, ImageFont

    n = len(combinations)
    if n == 0:
        raise ValueError("No combinations to draw.")
    if mols_per_row is None:
        mols_per_row = max(1, math.ceil(math.sqrt(n)))
    n_rows = math.ceil(n / mols_per_row)

    label_height = 28
    tile_w, tile_h = subimg_size
    canvas = Image.new(
        "RGB",
        (tile_w * mols_per_row, (tile_h + label_height) * n_rows),
        "white",
    )
    draw_ctx = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.load_default()
    except Exception:
        font = None

    for idx, combo in enumerate(combinations):
        df = build_combination_dataframe(mol, molgraph, carbon_classes, combo)
        carbon_nodes = df.loc[df.role == "carbon_representative", "atom_idx"].tolist()
        connector_nodes = df.loc[df.role == "connector", "atom_idx"].tolist()

        highlight_atoms = carbon_nodes + connector_nodes
        highlight_colors = {a: carbon_color for a in carbon_nodes}
        highlight_colors.update({a: connector_color for a in connector_nodes})
        highlight_radii = {a: highlight_radius for a in highlight_atoms}

        png_bytes = _render_molecule_png_bytes(
            mol, highlight_atoms, highlight_colors, highlight_radii, subimg_size
        )
        tile_img = Image.open(io.BytesIO(png_bytes)).convert("RGB")

        row, col = divmod(idx, mols_per_row)
        x = col * tile_w
        y = row * (tile_h + label_height)
        canvas.paste(tile_img, (x, y + label_height))
        draw_ctx.text((x + 6, y + 6), f"{label_prefix} {idx}", fill="black", font=font)

    canvas.save(output_path)
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
    parser.add_argument(
        "--n-units",
        type=int,
        default=None,
        help="Number of repeating symmetric branches/units in the molecule. "
        "If omitted, this is auto-detected from the symmetry classes.",
    )
    parser.add_argument(
        "--grid-cols",
        type=int,
        default=None,
        help="Number of columns in the combinations grid image. Defaults "
        "to a roughly square layout based on the number of combinations found.",
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

    # --- Plot 1: single molecule, symmetry-highlighted minimal atom set ---
    image_path = molfile_path.with_name(molfile_path.stem + "_highlighted.png")
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

    # --- Plot 2: grid of every distinct branch combination ---
    combos = enumerate_symmetric_representative_sets(molgraph, carbon_classes, n_units=args.n_units)
    dataframes = [build_combination_dataframe(mol, molgraph, carbon_classes, c) for c in combos]
    print()
    print(f"Found {len(combos)} distinct branch combination(s).")

    grid_path = molfile_path.with_name(molfile_path.stem + "_combinations_grid.png")
    draw_combination_grid(
        mol,
        molgraph,
        carbon_classes,
        combos,
        grid_path,
        highlight_radius=args.highlight_radius,
        mols_per_row=args.grid_cols,
    )
    print(f"Combinations grid image saved to: {grid_path}")

    return dataframes


if __name__ == "__main__":
    main()
