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
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            weight=1,
            bond_order=bond.GetBondTypeAsDouble(),
        )
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


def find_graph_automorphisms(mol, molgraph, max_automorphisms=20000):
    """
    Find the graph automorphism group of the molecule (all heavy atoms,
    bond-order-aware), after pruning redundant same-parent terminal (leaf)
    siblings -- e.g. the 3 methyls of a tert-butyl, or the 3 F's of a CF3 --
    down to a single representative each. Without this pruning, free-rotating
    substituents cause a combinatorial explosion in raw automorphism count
    (thousands+) that says nothing about the molecule's real structural
    symmetry; pruning removes that noise while keeping every chemically
    meaningful symmetry element (including independent, unrelated symmetric
    sites, and any local sub-symmetry nested within a larger repeating unit).

    Returns
    -------
    automorphisms : list of dict
        Each dict maps atom_idx -> atom_idx (only for atoms retained after
        pruning; pruned leaf atoms are not keys and should be treated as
        fixed/unmapped by callers).
    """
    from networkx.algorithms.isomorphism import GraphMatcher

    groups = defaultdict(list)
    for node in molgraph.nodes():
        if molgraph.degree(node) == 1:
            parent = next(iter(molgraph.neighbors(node)))
            groups[(parent, molgraph.nodes[node]["symbol"])].append(node)

    reduced = molgraph.copy()
    for (_parent, _symbol), members in groups.items():
        if len(members) > 1:
            keep = min(members)
            for m in members:
                if m != keep:
                    reduced.remove_node(m)

    gm = GraphMatcher(
        reduced,
        reduced,
        node_match=lambda a, b: a["symbol"] == b["symbol"],
        edge_match=lambda a, b: a.get("bond_order") == b.get("bond_order"),
    )

    automorphisms = []
    for perm in gm.isomorphisms_iter():
        automorphisms.append(perm)
        if len(automorphisms) >= max_automorphisms:
            print(
                f"Warning: automorphism search capped at {max_automorphisms} -- "
                "results may be incomplete for this molecule."
            )
            break
    return automorphisms


def enumerate_all_symmetric_combinations(
    mol, molgraph, carbon_classes, baseline_tree_nodes=None, baseline_chosen_reps=None
):
    """
    Discover every distinct, structurally valid minimal atom set, by finding
    the molecule's true graph automorphism group and applying each
    automorphism to one known-good baseline solution. This correctly handles
    molecules with multiple independent symmetric sites (e.g. two unrelated
    rings each with their own local mirror) as well as nested local symmetry
    within a larger repeating branch -- both of which a single-anchor
    heuristic cannot represent.

    Parameters
    ----------
    baseline_tree_nodes, baseline_chosen_reps : optional
        A known-good minimal solution (from greedy_group_steiner_tree) to
        use as the seed. Computed automatically if not provided.

    Returns
    -------
    combinations : list of dict
        Each dict has keys:
          'chosen_reps' : class_idx -> representative atom idx
          'tree_nodes'  : full connected minimal atom set (frozenset),
                          including any non-carbon connector atoms
        Deduplicated -- automorphisms that map the baseline onto the same
        final atom set only appear once.
    """
    if baseline_tree_nodes is None or baseline_chosen_reps is None:
        baseline_tree_nodes, baseline_chosen_reps = greedy_group_steiner_tree(
            molgraph, carbon_classes
        )

    automorphisms = find_graph_automorphisms(mol, molgraph)

    atom_to_class = {}
    for ci, cls in enumerate(carbon_classes):
        for a in cls:
            atom_to_class[a] = ci

    seen = set()
    combinations = []
    for perm in automorphisms:
        mapped_tree = frozenset(perm.get(a, a) for a in baseline_tree_nodes)
        if mapped_tree in seen:
            continue
        seen.add(mapped_tree)

        chosen_reps = {}
        for a in mapped_tree:
            ci = atom_to_class.get(a)
            if ci is not None:
                chosen_reps[ci] = a

        combinations.append({"chosen_reps": chosen_reps, "tree_nodes": mapped_tree})

    return combinations


def is_trivial_leaf_class(molgraph, cls):
    """
    A class is 'trivial' if every member is a terminal (degree-1) atom --
    e.g. the 3 methyls of a tert-butyl, or the 3 F's of a CF3. These arise
    from free rotation of a substituent and are excluded from the
    side-consistency scoring below, since which specific leaf is chosen
    is chemically/visually inconsequential.
    """
    return all(molgraph.degree(atom) == 1 for atom in cls)


def score_side_consistency(mol, molgraph, carbon_classes, chosen_reps):
    """
    Score (0 to 1) how visually "one-sided" a combination looks: for every
    non-trivial, multi-member class, compute the chosen atom's 2D offset
    from that class's own local centroid (i.e. relative to just its own
    members, not the whole molecule -- this matters because a symmetric
    site positioned near the molecule's overall centroid would otherwise
    show a near-meaningless offset). Then check whether these per-class
    offset vectors, across every independent symmetric site in this
    combination, point in broadly the same 2D direction (positive pairwise
    dot product) rather than conflicting directions.

    1.0 = every site's choice points the same way (a clean, one-sided
    combination, e.g. all changes toward the bottom of the drawing).
    Lower values mean different sites point in conflicting directions (a
    "crossed over" look).

    This is a display/aesthetic heuristic, not a chemical one -- differently
    scored combinations are all equally chemically valid.
    """
    from rdkit.Chem import AllChem

    mol_for_coords = mol
    has_real_coords = mol.GetNumConformers() > 0 and not all(
        mol.GetConformer().GetAtomPosition(i).x == 0
        and mol.GetConformer().GetAtomPosition(i).y == 0
        for i in range(mol.GetNumAtoms())
    )
    if not has_real_coords:
        mol_for_coords = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol_for_coords)
    conf = mol_for_coords.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]

    vectors = []
    for ci, cls in enumerate(carbon_classes):
        if len(cls) <= 1 or is_trivial_leaf_class(molgraph, cls):
            continue
        atom = chosen_reps.get(ci)
        if atom is None:
            continue
        cx = sum(xs[a] for a in cls) / len(cls)
        cy = sum(ys[a] for a in cls) / len(cls)
        dx, dy = xs[atom] - cx, ys[atom] - cy
        if dx != 0 or dy != 0:
            vectors.append((dx, dy))

    if len(vectors) < 2:
        return 1.0  # nothing to compare against -- trivially consistent

    agreements, total = 0, 0
    for i in range(len(vectors)):
        for j in range(i + 1, len(vectors)):
            dot = vectors[i][0] * vectors[j][0] + vectors[i][1] * vectors[j][1]
            if dot != 0:
                total += 1
                if dot > 0:
                    agreements += 1

    return agreements / total if total else 1.0


def rank_combinations_by_side_consistency(mol, molgraph, carbon_classes, combinations):
    """
    Score every combination (from enumerate_all_symmetric_combinations) by
    side-consistency and sort best-first.

    Returns
    -------
    scored : list of (score, combination_dict), sorted descending by score
    """
    scored = [
        (score_side_consistency(mol, molgraph, carbon_classes, c["chosen_reps"]), c)
        for c in combinations
    ]
    scored.sort(key=lambda item: item[0], reverse=True)
    return scored


def select_preferred_combinations(mol, molgraph, carbon_classes, combinations, tol=1e-6):
    """
    From all valid combinations, return just the subset that achieves the
    best (highest) side-consistency score -- the "nicest looking" ones,
    suitable as the default display choice(s).

    Returns
    -------
    preferred : list of combination dicts, all tied for the top score
    scored : the full ranked list (score, combination_dict), for reference
    """
    scored = rank_combinations_by_side_consistency(mol, molgraph, carbon_classes, combinations)
    best_score = scored[0][0]
    preferred = [c for s, c in scored if abs(s - best_score) <= tol]
    return preferred, scored


def build_combination_dataframe(mol, molgraph, carbon_classes, chosen_reps, tree_nodes=None):
    """
    Given one combination's class_idx -> representative atom idx mapping,
    return the connected minimal atom set as a pandas DataFrame, one row
    per atom.

    Parameters
    ----------
    tree_nodes : iterable of atom idx, optional
        The full, already-known connected atom set (e.g. from
        enumerate_all_symmetric_combinations, which computes this exactly
        via graph automorphism -- pass it through to avoid recomputing).
        If not provided, connectivity is derived from `chosen_reps` alone,
        patching in connector atoms via a Steiner tree if needed.
    """
    import pandas as pd

    if tree_nodes is not None:
        tree_nodes = set(tree_nodes)
    else:
        terminals = list(chosen_reps.values())
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
    scores=None,
):
    """
    Render every combination (from enumerate_all_symmetric_combinations) as
    its own highlighted panel, arranged in a grid in a single image.

    Parameters
    ----------
    combinations : list of dict
        Each with 'chosen_reps' and 'tree_nodes' keys, as returned by
        enumerate_all_symmetric_combinations.
    output_path : str or Path
        Where to save the composite grid PNG.
    subimg_size : (int, int)
        Size of each individual panel, in pixels.
    mols_per_row : int, optional
        Number of columns in the grid. Defaults to ceil(sqrt(n)) for a
        roughly square layout.
    scores : list of float, optional
        If provided (same length/order as `combinations`), each panel's
        label includes its side-consistency score.
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
        chosen_reps = combo["chosen_reps"]
        tree_nodes = combo.get("tree_nodes")
        df = build_combination_dataframe(mol, molgraph, carbon_classes, chosen_reps, tree_nodes)
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
        label = f"{label_prefix} {idx}"
        if scores is not None:
            label += f"  (score={scores[idx]:.2f})"
        draw_ctx.text((x + 6, y + 6), label, fill="black", font=font)

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
        "--grid-cols",
        type=int,
        default=None,
        help="Number of columns in the combinations grid image. Defaults "
        "to a roughly square layout based on the number of combinations found.",
    )
    parser.add_argument(
        "--show-all-combinations",
        action="store_true",
        help="Show every distinct valid combination in the grid, instead of "
        "just the 'nicest looking' ones (those with the best side-consistency "
        "score -- i.e. every symmetric choice consistently on the same side "
        "of the molecule).",
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

    # --- Plot 2: grid of distinct combinations, preferring "same side" ones ---
    all_combos = enumerate_all_symmetric_combinations(
        mol, molgraph, carbon_classes, baseline_tree_nodes=tree_nodes, baseline_chosen_reps=chosen_reps
    )
    scored = rank_combinations_by_side_consistency(mol, molgraph, carbon_classes, all_combos)
    best_score = scored[0][0]
    preferred_count = sum(1 for s, _c in scored if abs(s - best_score) <= 1e-6)

    print()
    print(f"Found {len(all_combos)} distinct valid combination(s) total.")
    print(f"{preferred_count} of them are 'nicest looking' (best side-consistency score = {best_score:.2f}).")

    if args.show_all_combinations:
        combos_to_draw = [c for _s, c in scored]
        scores_to_draw = [s for s, _c in scored]
    else:
        combos_to_draw = [c for _s, c in scored[:preferred_count]]
        scores_to_draw = [s for s, _c in scored[:preferred_count]]

    dataframes = [
        build_combination_dataframe(mol, molgraph, carbon_classes, c["chosen_reps"], c["tree_nodes"])
        for c in combos_to_draw
    ]

    print(dataframes[-1])

    grid_path = molfile_path.with_name(molfile_path.stem + "_combinations_grid.png")
    draw_combination_grid(
        mol,
        molgraph,
        carbon_classes,
        combos_to_draw,
        grid_path,
        highlight_radius=args.highlight_radius,
        mols_per_row=args.grid_cols,
        scores=scores_to_draw,
    )
    print(f"Combinations grid image saved to: {grid_path}")

    print(type(combos_to_draw), type(combos_to_draw[0]))
    print(combos_to_draw[0].items())

    return dataframes


if __name__ == "__main__":
    dfs = main()

    for df in dfs:
        print(df)
        print()
