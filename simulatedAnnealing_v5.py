import sys
import copy

sys.path.append("../")

# import platform
import math
import random
import json
import numpy as np
import pandas as pd
import networkx as nx
from rdkit import Chem
from pathlib import Path

# import matplotlib.pyplot as plt
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol
from typing import Dict, List, Optional, Tuple

# # node color map
# color_map = {
#     0: "#FFA500",
#     1: "#98FB98",
#     2: "yellow",
#     3: "#00FFFF",
#     -1: "lightblue",  # For non-carbon atoms
#     -2: "lightgrey",
# }

from globals import SVG_DIMENSIONS as svgDimensions
from globals import NODE_COLOR_MAP as color_map


def compute_total_weight(
    graph: nx.Graph,
    mapping: Dict[int, int],
    shortest_paths: Dict[int, Dict[int, float]],
) -> float:
    """
    Compute the total weight for a given mapping.

    Args:
        graph (nx.Graph): NetworkX graph of the molecule, with nodes representing atoms and edges representing bonds.
        mapping (Dict[int, int]): Node mapping between the molecule graph and the NMR graph.
        shortest_paths (Dict[int, Dict[int, float]]): Shortest paths between all pairs of nodes in the molecule graph.

    Returns:
        float: Sum of the path lengths of COSY and HMBC edges mapped onto the molecule graph.
    """

    total_weight = 0
    for u, v, d in graph.edges(data=True):
        if d["cosy"] or d["hmbc"]:
            uu = mapping[u]
            vv = mapping[v]

            if uu == vv:
                continue

            if d.get("cosy"):
                weight = shortest_paths[uu][vv]
                total_weight += (weight - 1) ** 3
            if d.get("hmbc"):
                weight = shortest_paths[uu][vv]
                if weight < 3:
                    weight = 2
                total_weight += (weight - 2) ** 3

    return total_weight


def modify_mapping(grouped_nodes, current_mapping, nProtons_to_nodes):
    """
    Modify a mapping by randomly swapping values between two nodes.

    Parameters:
        permutations (dict): Dictionary of nodes and their possible mappings.
        current_mapping (dict): Current mapping {node_in_graph1: node_in_graph2}.

    Returns:
        dict: Updated mapping after the swap.
        bool: True if a swap was made, False otherwise.
    """
    new_mapping = current_mapping.copy()

    # Step 1: Choose a random key (node in graph1)
    random_node = random.choice(list(nProtons_to_nodes.keys()))
    nprotons_key = int(nProtons_to_nodes[random_node])

    # Step 2: Randomly select a value from its possible mappings
    possible_values = grouped_nodes[nprotons_key]
    chosen_value = random.choice(possible_values)

    # Step 3: Check if the chosen value is different from the key
    swapped = False
    if random_node != chosen_value:
        v1 = new_mapping[random_node]
        v2 = new_mapping[chosen_value]
        new_mapping[random_node] = v2
        new_mapping[chosen_value] = v1

        swapped = True

    return new_mapping, swapped


def compute_total_weight(
    graph: nx.Graph,
    mapping: Dict[int, int],
    shortest_paths: Dict[int, Dict[int, float]],
) -> float:
    """
    Compute the total weight for a given mapping.

    Args:
        graph (nx.Graph): NetworkX graph of the molecule, with nodes representing atoms and edges representing bonds.
        mapping (Dict[int, int]): Node mapping between the molecule graph and the NMR graph.
        shortest_paths (Dict[int, Dict[int, float]]): Shortest paths between all pairs of nodes in the molecule graph.

    Returns:
        float: Sum of the path lengths of COSY and HMBC edges mapped onto the molecule graph.
    """


    total_weight = 0
    for u, v, d in graph.edges(data=True):
        if d["cosy"] or d["hmbc"]:
            uu = mapping[u]
            vv = mapping[v]

            if d.get("cosy"):
                weight = shortest_paths[uu][vv]
                total_weight += (weight - 1) ** 3
            if d.get("hmbc"):
                weight = shortest_paths[uu][vv]
                if weight < 3:
                    weight = 2
                total_weight += (weight - 2) ** 3
                # total_weight += weight

    return total_weight


def simulated_annealing(
    G2: nx.Graph,
    shortest_paths: Dict[int, Dict[int, int]],
    grouped_nodes: Dict[int, List[int]],
    nProtons_to_nodes: Dict[int, int],
    current_mapping: Optional[Dict[int, int]] = None,
    randomize_mapping=True,
    initial_temp: float = 100000,
    final_temp: float = 0.001,
    cooling_rate: float = 0.99,
    max_iterations: int = 50000,
) -> Tuple[Dict[int, int], int, Dict[str, int], List[Tuple[int, int, int, float]]]:
    """
    Simulated annealing algorithm for graph mapping.

    Parameters:
    - G2: The second graph (networkx.Graph). carbon_graph
    - shortest_paths: A dictionary of shortest paths in G2.
    - grouped_nodes: A dictionary of possible permutations for each node.
    - current_mapping: The initial mapping of nodes from G1 to G2.
    - randomized_mapping: Randomize the initial mapping.
    - initial_temp: The initial temperature for the annealing process.
    - cooling_rate: The rate at which the temperature decreases.
    - max_iterations: The maximum number of iterations to perform.
    - ax: Matplotlib Axes object for plotting.
    - kwargs: Additional keyword arguments for plotting.

    Returns:
    - best_mapping: The best mapping found.
    - best_weight: The weight of the best mapping.
    - stats: Statistics about the annealing process.
    - weights: A list of tuples containing iteration, best weight, current weight, and temperature.
    """

    # Initialize mappings and weights
    if current_mapping is None:
        current_mapping = {node: node for node in G2.nodes()}
    best_mapping = current_mapping.copy()
    current_weight = best_weight = compute_total_weight(
        G2, current_mapping, shortest_paths
    )  # compute_total_weight(G2, current_mapping, shortest_paths)
    temperature = initial_temp

    stats = {
        "not_swapped": 0,
        "improvements": 0,
        "worsen_accepted": 0,
        "move_rejected": 0,
    }
    weights = []

    for iteration in range(max_iterations):
        # mapping, graph, grouped_nodes
        neighbor_mapping, swapped = modify_mapping(
            grouped_nodes, current_mapping, nProtons_to_nodes
        )  # modify_mapping(current_mapping, G2, grouped_nodes )

        if not swapped:
            stats["not_swapped"] += 1
            continue

        neighbor_weight = compute_total_weight(G2, neighbor_mapping, shortest_paths)

        # Update mapping and weights based on acceptance criteria
        if neighbor_weight < current_weight:
            stats["improvements"] += 1
            current_mapping = neighbor_mapping
            current_weight = neighbor_weight
        elif random.random() < math.exp(
            (current_weight - neighbor_weight) / temperature
        ):
            current_mapping = neighbor_mapping
            current_weight = neighbor_weight
            stats["worsen_accepted"] += 1
        else:
            stats["move_rejected"] += 1

        # Update the best solution found so far
        if current_weight < best_weight:
            best_mapping, best_weight = current_mapping, current_weight

        weights.append((iteration, best_weight, current_weight, temperature))

        # Cooling
        temperature *= cooling_rate

        # Early stopping condition
        if temperature < final_temp:
            break

    return best_mapping, best_weight, stats, weights


class NMRgraphs:
    def __init__(self, df, links):
        self.df = df
        self.graph = nx.Graph()
        self.graph.add_nodes_from(df.id)

        # add node attributes from the dataframe to the graph
        columns = df.columns
        for _, row in df.iterrows():
            id = row["id"]
            for col in columns:
                self.graph.nodes[id][col] = row[col]

        # add links to the graph
        for link in links:
            self.graph.add_edge(
                link["source"],
                link["target"],
                bond=link.get("bond", False),
                cosy=link.get("cosy", False),
                hmbc=link.get("hmbc", False),
            )


class SimulatedAnnealing2:
    def __init__(self, json_data):

        self.json_data = json_data
        self.nmr_nodes = json_data["nodes_now"]
        self.nmr_links = json_data["links"]

        # remove any links where the source equals the target
        self.nmr_links = [
            link for link in self.nmr_links if link["source"] != link["target"]
        ]

        self.nodes_offset = 0
        # adjust index offset in nmr_nodes and nmr_links
        if (
            json_data["oldjsondata"]["MNOVAcalcMethod"]["data"]["0"]
            == "NMRSHIFTDB2 Predict"
        ):
            self.nodes_offset = 0

        elif(
            json_data["oldjsondata"]["MNOVAcalcMethod"]["data"]["0"]
            == "JEOL Predict"
        ): 
            self.nodes_offset = 0
        else:
            # self.nodes_offset = 1
            self.nodes_offset = 0  # EEH 2025-sep-06


        for link in self.nmr_links:
            link["source"] = int(link["source"]) - self.nodes_offset
            link["target"] = int(link["target"]) - self.nodes_offset
            link["weight"] = 1

        for node in self.nmr_nodes:
            node["id"] = int(node["id"]) - self.nodes_offset

        mol_str = json_data["molfile"]
        self.mol = Chem.MolFromMolBlock(mol_str)

        self.mol_links = []
        # add edges from simAnneal.mol
        for idx, bond in enumerate(self.mol.GetBonds()):
            self.mol_links.append(
                {
                    "source": bond.GetBeginAtomIdx(),
                    "target": bond.GetEndAtomIdx(),
                    "bond": True,
                }
            )

        self.all_links = self.mol_links + self.nmr_links

        self.graph_df = self.create_graph_df(json_data, self.nmr_nodes)
        self.carbon_df = self.graph_df[self.graph_df["symbol"] == "C"]

        self.xy3 = self.calc_xy3_coords(self.mol)

        # update coodinates in graph_df
        for idx, (x, y) in self.xy3.items():
            self.graph_df.loc[idx, "x"] = x
            self.graph_df.loc[idx, "y"] = y

        self.nmr_graph = NMRgraphs(self.graph_df, self.all_links).graph
        self.carbon_graph = NMRgraphs(self.carbon_df, self.nmr_links).graph
        self.mol_graph = NMRgraphs(self.graph_df, self.mol_links).graph

        cosy_hmbc_edges = [
            (u, v)
            for u, v, d in self.nmr_graph.edges(data=True)
            if d.get("cosy") or d.get("hmbc")
        ]
        self.cosy_hmbc_subgraph = self.nmr_graph.edge_subgraph(cosy_hmbc_edges)

        cosy_edges = [
            (u, v) for u, v, d in self.nmr_graph.edges(data=True) if d["cosy"]
        ]
        self.cosy_subgraph = self.nmr_graph.edge_subgraph(cosy_edges)

        hmbc_edges = [
            (u, v) for u, v, d in self.nmr_graph.edges(data=True) if d["hmbc"]
        ]
        self.hmbc_subgraph = self.nmr_graph.edge_subgraph(hmbc_edges)

        bonds_edges = [
            (u, v) for u, v, d in self.nmr_graph.edges(data=True) if d["bond"]
        ]
        self.bonds_subgraph = self.nmr_graph.edge_subgraph(bonds_edges)

        self.carbon_grouped_nodes = {}
        for node in self.carbon_graph.nodes(data=True):
            # the keys are numProtons and add the node to a list of nodes with the same numProtons
            self.carbon_grouped_nodes.setdefault(node[1]["numProtons"], []).append(
                node[0]
            )

    @classmethod
    def from_json_file(cls, json_file_path: Path):
        if not json_file_path.exists():
            raise FileNotFoundError(f"File {json_file_path} not found.")
        with open(json_file_path, "r") as f:
            json_data = json.load(f)

        return cls(json_data)

    @classmethod
    def from_params(cls, nodes_now, links, molfile, oldjson_data):
        # collect all parameters internally
        # make a list of what is used from the dict

        # "nodes_now"
        # "links"
        # "molfile"
        # json_data["oldjsondata"]["allAtomsInfo"]

        json_data = {
            "nodes_now": nodes_now,
            "links": links,
            "molfile": molfile,
            "oldjsondata": oldjson_data,
        }

        return cls(json_data)

    def create_graph_df(self, json_data, nmr_nodes):

        # create a dataframe from the mol_graph nodes using the data dictionary as the columns transpose the dataframe
        graph_df = pd.DataFrame.from_dict(
            json_data["oldjsondata"]["allAtomsInfo"]["data"], orient="index"
        )
        c13predictions_df = pd.DataFrame.from_dict(
            json_data["oldjsondata"]["c13predictions"]["data"], orient="index"
        )

        # change index to integer
        graph_df.index = graph_df.index.astype(int)


        for cnode in nmr_nodes:
            for key, value in cnode.items():
                if isinstance(value, list):
                    value1 = "[" + ",".join([str(x) for x in value]) + "]"
                else:
                    value1 = value
                graph_df.at[int(cnode["id"]), key] = value1


        graph_df["iupacLabel"] = graph_df["iupacLabel"].fillna("")
        graph_df["jCouplingClass"] = graph_df["jCouplingClass"].fillna("")
        graph_df["jCouplingVals"] = graph_df["jCouplingVals"].fillna("")
        graph_df["x"] = graph_df["x"].fillna("")
        graph_df["y"] = graph_df["y"].fillna("")
        graph_df["H1_ppm"] = graph_df["H1_ppm"].fillna("")
        graph_df["ppm"] = graph_df["ppm"].fillna("")
        graph_df["ppm_calculated"] = graph_df["ppm_calculated"].fillna("")
        if "visible" in graph_df.columns:
            graph_df["visible"] = graph_df["visible"].fillna(True)
        else:
            graph_df["visible"] = True

        # add color column based on numprotons and atom type
        graph_df["color"] = graph_df.apply(
            lambda x: color_map.get(x["numProtons"], color_map.get(-1)), axis=1
        )

        # set color of non carbon atoms to lightblue
        graph_df.loc[graph_df["symbol"] != "C", "color"] = color_map.get(-1)

        return graph_df

    # def calc_xy3_coords(
    #     self,
    #     mol: Mol,
    #     molWidth: int = 1000,
    #     molHeight: int = 400,
    #     svgWidth: int = 1200,
    #     svgHeight: int = 600,
    # ) -> Tuple[str, Dict[int, Tuple[float, float]]]:
    def calc_xy3_coords(
        self,
        mol: Mol,
        molWidth: int = svgDimensions.mol_width,
        molHeight: int = svgDimensions.mol_height,
        svgWidth: int = svgDimensions.svg_width,
        svgHeight: int = svgDimensions.svg_height
    ) -> Tuple[str, Dict[int, Tuple[float, float]]]:
        """
        Create an SVG string representation of a molecule and extract atom coordinates.

        Parameters:
        mol (Mol): RDKit molecule object.
        molWidth (int): Width of the molecule drawing area in pixels. Default is 1000.
        molHeight (int): Height of the molecule drawing area in pixels. Default is 400.
        svgWidth (int): Width of the SVG canvas in pixels. Default is 1200.
        svgHeight (int): Height of the SVG canvas in pixels. Default is 600.

        Returns:
        Tuple[str, Dict[int, Tuple[float, float]]]: A tuple containing the SVG string and a dictionary of atom coordinates.
        """
        translateWidth = int((svgWidth - molWidth) / 2)
        translateHeight = int((svgHeight - molHeight) / 2)

        d2d = Draw.rdMolDraw2D.MolDraw2DSVG(molWidth, molHeight)
        d2d.drawOptions().minFontSize = 20
        d2d.DrawMolecule(mol)
        d2d.TagAtoms(mol)
        d2d.FinishDrawing()

        sss = d2d.GetDrawingText()
        sss = d2d.GetDrawingText().replace(
            f"width='{molWidth}px' height='{molHeight}px'",
            f"width={molWidth} height={molHeight}",
        )
        sss = sss.replace("fill:#FFFFFF", "fill:none").replace(
            "<svg", '<svg class="center"'
        )

        sss = sss.replace(
            f"<!-- END OF HEADER -->",
            f"<!-- END OF HEADER -->\n<g transform='translate({translateWidth}, {translateHeight})'>",
        )
        sss = sss.replace("</svg>", "</g>\n</svg>")
        sss = sss.replace(
            f"width={molWidth} height={molHeight} viewBox='0 0 {molWidth} {molHeight}'",
            f"width={svgWidth} height={svgHeight} viewBox='0 0 {svgWidth} {svgHeight}'",
        )

        idx_list = []
        xxx = []
        yyy = []

        new_xy3 = {}
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            point = d2d.GetDrawCoords(idx)
            idx_list.append(idx)
            xxx.append(point.x / molWidth)
            yyy.append(point.y / molHeight)
            new_xy3[idx] = (point.x / molWidth, point.y / molHeight)

        return new_xy3

    def initialize_mapping(
        self,
        carbon_graph: nx.Graph,
        carbon_grouped_nodes: Dict[int, List[int]],
        randomize_mapping: bool = True,
    ) -> Dict[int, int]:
        """
        Initialize a mapping of nodes in the carbon graph.

        Args:
            carbon_graph (nx.Graph): The graph representing the carbon atoms.
            carbon_grouped_nodes (Dict[int, List[int]]): A dictionary grouping nodes by the number of protons.
            randomize_mapping (bool): Whether to randomize the initial mapping.

        Returns:
            Dict[int, int]: The initialized mapping of nodes.
        """
        random_mapping = {node: node for node in carbon_graph.nodes()}
        if randomize_mapping:
            for nprotons in range(4):
                nodes1 = carbon_grouped_nodes.get(nprotons, []).copy()
                nodes2 = carbon_grouped_nodes.get(nprotons, []).copy()
                random.shuffle(nodes2)
                for i, node in enumerate(nodes1):
                    random_mapping[node] = nodes2[i]
        return random_mapping


    def map_protons_to_nodes(self, graph_df: pd.DataFrame) -> Dict[int, int]:
        """
        Map the number of protons to their corresponding nodes in the dataframe.

        Args:
            graph_df (pd.DataFrame): DataFrame containing molecular information with columns 'symbol' and 'numProtons'.

        Returns:
            Dict[int, int]: A dictionary mapping node indices to the number of protons.
        """
        nProtons_to_nodes = {}
        for idx, row in graph_df.iterrows():
            if row["symbol"] != "C":
                continue
            nProtons = row["numProtons"]
            nProtons_to_nodes[idx] = nProtons

        return nProtons_to_nodes

    def run_optimization(self, num_times: int = 1) -> None:
        """
        Run the simulated annealing optimization process multiple times.

        Args:
            num_times (int): Number of times to run the optimization process.
        """

        self.bestest_weight = float("inf")

        self.predicted_mapping = self.initialize_mapping(
            self.carbon_graph, self.carbon_grouped_nodes, randomize_mapping=False
        )

        self.predicted_weight = compute_total_weight(
            self.carbon_graph, self.predicted_mapping, self.shortest_paths
        )

        self.initial_mapping = self.initialize_mapping(
            self.carbon_graph,
            self.carbon_grouped_nodes,
            randomize_mapping=self.randomize_mapping,
        )
        self.initial_weight = compute_total_weight(
            self.carbon_graph, self.initial_mapping, self.shortest_paths
        )

        self.results = {}

        for i in range(num_times):

            self.initial_mapping = self.initialize_mapping(
                self.carbon_graph,
                self.carbon_grouped_nodes,
                randomize_mapping=self.randomize_mapping,
            )
            self.initial_weight = compute_total_weight(
                self.carbon_graph, self.initial_mapping, self.shortest_paths
            )

            best_mapping, best_weight, stats, weights = simulated_annealing(
                self.carbon_graph,
                self.shortest_paths,
                self.carbon_grouped_nodes,
                self.nProtons_to_nodes,
                self.initial_mapping,
                randomize_mapping=self.randomize_mapping,
                initial_temp=self.max_temp,
                final_temp=self.min_temp,
                cooling_rate=self.cooling_rate,
                max_iterations=self.max_iter,
            )

            if best_weight not in self.results:
                self.results[best_weight] = {}
                self.results[best_weight]["num_times"] = 0
                self.results[best_weight]["results"] = []
                self.results[best_weight]["best_mapping"] = []

            self.results[best_weight]["results"].append((best_mapping, stats, weights))
            self.results[best_weight]["num_times"] += 1
            self.results[best_weight]["best_mapping"].append(best_mapping)

            if best_weight < self.bestest_weight:
                self.bestest_weight = best_weight
                self.bestest_mapping = best_mapping
                self.bestest_stats = stats
                self.bestest_weights = weights

    def setup_run(
        self,
        max_iter: int = 100000,
        max_temp: float = 100.0,
        min_temp: float = 0.1,
        cooling_rate: float = 0.995,
        randomize_mapping: bool = False,
    ) -> None:
        """
        Set up the parameters for the simulated annealing run.

        Args:
            max_iter (int): Maximum number of iterations for the annealing process.
            max_temp (float): Maximum temperature for the annealing process.
            min_temp (float): Minimum temperature for the annealing process.
            cooling_rate (float): Cooling rate for the annealing process.
            randomize_mapping (bool): Whether to randomize the initial mapping.
        """

        self.max_iter = max_iter
        self.max_temp = max_temp
        self.min_temp = min_temp
        self.cooling_rate = cooling_rate
        self.randomize_mapping = randomize_mapping

        # initialize the mapping
        self.mapping = self.initialize_mapping(
            self.carbon_graph, self.carbon_grouped_nodes, self.randomize_mapping
        )

        self.shortest_paths = dict(nx.all_pairs_shortest_path_length(self.mol_graph))

        self.nProtons_to_nodes = self.map_protons_to_nodes(self.graph_df)

        self.predicted_mapping = self.initialize_mapping(
            self.carbon_graph, self.carbon_grouped_nodes, randomize_mapping=False
        )

        self.predicted_weight = compute_total_weight(
            self.carbon_graph, self.predicted_mapping, self.shortest_paths
        )

    def process_results(self, catoms_df, jsonGraphData):

        best_weight = self.bestest_weight

        unique_mappings = set()
        unique_mapping_dict = {}

        for mapping in self.results[best_weight]["best_mapping"]:
            v1 = list(mapping.values())
            v1_str = "".join([str(x) for x in v1])
            #  add the string to the set
            unique_mappings.add(v1_str)
            unique_mapping_dict[v1_str] = mapping

        optimized_nodes_dicts = {}
        for k, v in unique_mapping_dict.items():
            optimized_nodes = []
            for node in self.nmr_nodes:
                moved_node = copy.deepcopy(node)
                orig_id = node["id"]
                moved_id = v[orig_id]

                if orig_id != moved_id:
                    moved_node["id"] = moved_id
                    moved_node["atomNumber"] = catoms_df.loc[moved_id, "atomNumber"]
                    moved_node["ppm_calculated"] = catoms_df.loc[
                        moved_id, "ppm_calculated"
                    ]
                    moved_node["x"] = catoms_df.loc[moved_id, "x"]
                    moved_node["y"] = catoms_df.loc[moved_id, "y"]
                    moved_node["jCouplingVals"] = catoms_df.loc[
                        moved_id, "jCouplingVals"
                    ]
                    moved_node["jCouplingClass"] = catoms_df.loc[
                        moved_id, "jCouplingClass"
                    ]

                optimized_nodes.append(moved_node)

            optimized_nodes_dicts[k] = optimized_nodes

        # calc MAE and LAE for each of the unique nodes_dict
        for k, optimized_nodes in optimized_nodes_dicts.items():
            mae = 0.0
            lae_biggest = 0.0
            num_carbons = 0
            lae_atomNumber = -1
            for node in optimized_nodes:
                ppm_exptal = node["ppm"]
                ppm_calculated = node["ppm_calculated"]

                try:
                    ppm_exptal = float(ppm_exptal)
                    ppm_calculated = float(ppm_calculated)
                    lae = abs(ppm_exptal - ppm_calculated)
                    mae += lae
                    if lae > lae_biggest:
                        lae_biggest = lae
                        lae_atomNumber = node["atomNumber"]
                    num_carbons += 1
                except:
                    continue

            mae = mae / num_carbons
            optimized_nodes_dicts[k] = [
                optimized_nodes,
                mae,
                lae_biggest,
                lae_atomNumber,
            ]

        # decide which is the best solution based on lowest MAE and LAE
        best_mae = float("inf")
        best_lae = float("inf")
        best_lae_atomNumber = -1
        for k, v in optimized_nodes_dicts.items():
            mae, lae, lae_atomNumber = v[1:]
            if (mae <= best_mae) and (lae <= best_lae):
                best_mae = mae
                best_lae = lae
                best_nodes = v[0]
                best_lae_atomNumber = lae_atomNumber
                best_key = k
                best_mapping = unique_mapping_dict[k]

        optimized_nodes_eeh = best_nodes
        best_results = {
            "best_weight": best_weight,
            "best_mae": best_mae,
            "best_lae": best_lae,
            "best_lae_atomNumber": best_lae_atomNumber,
        }

        for i, row in catoms_df.iterrows():
            id = row["id"] - self.nodes_offset

        optimized_links = copy.deepcopy(self.nmr_links)

        for link in optimized_links:
            link["source"] = best_mapping[link["source"]]
            link["target"] = best_mapping[link["target"]]

        # add offset back

        for node in optimized_nodes_eeh:
            node["id"] = int(node["id"]) + self.nodes_offset

        for link in optimized_links:
            link["source"] = int(link["source"]) + self.nodes_offset
            link["target"] = int(link["target"]) + self.nodes_offset

        jsonGraphData["nodes"] = optimized_nodes_eeh
        jsonGraphData["moved_nodes"] = optimized_nodes_eeh
        jsonGraphData["links"] = optimized_links

        return jsonGraphData, best_results

    def process_results_SA_skipped(self):

        best_results = {}

        # set the best weights to the predicted weight
        best_results["best_weight"] = self.predicted_weight

        # calculated MAE and LAE from predicted values
        # decide which is the best solution based on lowest MAE and LAE
        best_mae = float("inf")
        best_lae = float("inf")

        df = self.graph_df
        num_carbons = len(self.carbon_graph.nodes())
        best_mapping = self.predicted_mapping
        ppm_rss = 0.0
        mae = 0.0
        lae = 0.0
        lae_atomNumber = -1
        for node_orig, node_moved in best_mapping.items():
            id_orig = df.loc[node_orig, "id"]
            ppm_moved = df.loc[node_moved, "ppm"]
            ppm_calculated_orig = df.loc[node_orig, "ppm_calculated"]
            atomNumber_orig = df.loc[node_orig, "atomNumber"]

            if isinstance(ppm_calculated_orig, float) and isinstance(ppm_moved, float):
                ppm_rss += (ppm_calculated_orig - ppm_moved) ** 2
                diff = abs(ppm_calculated_orig - ppm_moved)
                mae += diff
                if diff > lae:
                    lae = diff
                    lae_atomNumber = atomNumber_orig

        best_results["best_mae"] = mae / num_carbons
        best_results["best_lae"] = lae
        best_results["best_lae_atomNumber"] = lae_atomNumber

        return best_results




