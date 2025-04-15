import json
import re
import argparse
from collections import defaultdict
from augur.utils import json_to_tree
import pandas as pd
import numpy as np

def merge_branchattrs(list_of_dicts):
    """
    Merge the branch_attrs dictionaries for multiple branches
    """

    # Initialize default dict with lists
    merged = defaultdict(lambda: defaultdict(list))

    def sort_mutations(lst):
        # sort mutations numberically
        return sorted(lst, key=lambda x: int(re.search(r'\d+', x).group()))

    # Merge 'mutations' key
    for d in list_of_dicts:
        for category, mutations in d['mutations'].items():
            merged['mutations'][category].extend(mutations)

    # Sort each list of muts
    for key in merged['mutations']:
        merged['mutations'][key] = sort_mutations(merged['mutations'][key])

    # Merge 'labels' if it exists
    label_dict = defaultdict(list)
    for d in list_of_dicts:
        if 'labels' in d.keys():
            if 'aa' in d['labels']:
                for label in d['labels']['aa'].split('; '):  # Split at '; ' to handle multiple genes
                    gene, mutation = label.split(': ')
                    label_dict[gene].append(mutation)

    # Format labels correctly
    merged_labels = []
    for gene, mutations in label_dict.items():
        merged_labels.append(f"{gene}: {', '.join(mutations)}")

    merged['labels']['aa'] = '; '.join(merged_labels)

    # Convert defaultdict to dict
    merged_attrs = {k: dict(v) for k, v in merged.items()}

    return merged_attrs


def find_egg_clusters(tree_file):
    """
    Read in the tree for the virus and segment and
    find all clusters of egg-only strains
    Will return a dictionary with the node that is the base of the egg-only clade as key
    and all the egg-tips as values
    """

    with open(tree_file) as json_handle:
      tree_json = json.load(json_handle)

    tree = json_to_tree(tree_json)



    # find egg-only clusters
    egg_only_cluster_polytomy_nodes_nested = {}

    for node in tree.find_clades(terminal=False):

      descending_tips = node.get_terminals()
      # only need to edit clusters with multiple tips
      if len(descending_tips)>1:

          # find if  descending tips are all egg
          passage_type_descendants = set([x.node_attrs['passage_category']['value'] for x in descending_tips])
          if passage_type_descendants == {'egg'}:
              egg_only_cluster_polytomy_nodes_nested[node.name] = [x.name for x in descending_tips]

    # now prune the egg_only_cluster_polytomy_nodes to only have the most ancestral nodes of clusters as keys
    # currently it will have nested clades of egg-only seqs
    egg_only_cluster_polytomy_nodes = {}

    for node in tree.find_clades(terminal=False):
      # if this node is in the dict
      if node.name in egg_only_cluster_polytomy_nodes_nested.keys():
          # see whether any of its ancestors are in the dict too
          path = tree.get_path(node)
          path_names = [x.name for x in path][:-1]

          # if not, then this is the base of the egg-only clade
          if not any(item in path_names for item in list(egg_only_cluster_polytomy_nodes_nested.keys())):
              egg_only_cluster_polytomy_nodes[node.name] = egg_only_cluster_polytomy_nodes_nested[node.name]

    return egg_only_cluster_polytomy_nodes


def make_new_terminal_branches(tree_file):
    """
    For each cluster of only egg sequences,
    and aggregate the branch length and mutations to get the new terminal branches of the polytomy
    """

    # clusters of only egg seqs
    egg_only_cluster_polytomy_nodes = find_egg_clusters(tree_file)

    # read in the tree again
    with open(tree_file) as json_handle:
        tree_json = json.load(json_handle)

    tree = json_to_tree(tree_json)

    # dictionary of new terminal branches to egg seqs
    # format is node at base of clade as key (this is the node to get deleted and replaced)
    # tip info as keys (this includes name, node_attrs and branch_attrs)
    polytomy_info_by_base_node = {}


    # for each egg cluster, compute the new terminal branches by aggregating the mutations and divergence
    # that occur along the path within the clade
    for node in tree.find_clades(terminal=False):
        # looking for the base of the egg cluster
        if node.name in egg_only_cluster_polytomy_nodes.keys():
            # initiate key in dictionary to store info for tips within this cluster
            polytomy_info_by_base_node[node.name] = []

            # get all egg tips
            egg_tips = node.get_terminals()

            # then aggregate muts and div along paths within cluster
            for tip in egg_tips:
                # get path from base of cluster to egg_tip
                # get_path() excludes the start, but includes the finish
                # we want the starting node too
                cluster_path = [node]+node.get_path(tip)
                # merge the branch_attrs (mutations)
                merged_muts = merge_branchattrs([x.branch_attrs for x in cluster_path])

                # divergence is cumulative, so can just take the divergence at the tip
                # and the rest of the node_attrs at tip
                polytomy_info_by_base_node[node.name].append({"name": tip.name,
                                                              "node_attrs": tip.node_attrs,
                                                              "branch_attrs": merged_muts})

    return polytomy_info_by_base_node


def replace_clusters(tree, nodes_to_replace):
    """
    Replace egg clusters on the tree with the polytomies
    """

    # recursively go through tree to find base of egg clusters and replace them
    if "children" in tree:
        new_children = []
        for child in tree["children"]:
            if child["name"] in nodes_to_replace:
                # Replace each basal node with the provided replacement polytomies
                new_children.extend(nodes_to_replace[child["name"]])
            else:
                # Recursively process child and collect the updated subtree
                new_children.append(replace_clusters(child, nodes_to_replace))

        # Update the tree with new children
        tree["children"] = new_children

    return tree


def write_tree_with_polytomies(tree_file, output_file):
    """
    Modify and write out an auspice json tree to convert egg clusters to polytomies
    """

    # clusters to replace with egg polytomies
    polytomy_info_by_base_node = make_new_terminal_branches(tree_file)


    # read in the tree again
    with open(tree_file) as json_handle:
        auspice_json = json.load(json_handle)

    # replace all egg clusters with polytomies
    tree_json_modified = replace_clusters(auspice_json['tree'], polytomy_info_by_base_node)

    # pop this into the "tree" part of the auspice json
    auspice_json['tree'] = tree_json_modified

    # save this
    with open(output_file, "w") as file:
        json.dump(auspice_json, file)
        #indent=2


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--tree', type=str, required=True,
                        help="nextstrain tree JSON, that has egg-sequence clusters")
    parser.add_argument('--output', required=True, help="name to output the tree with egg polytomies")

    args = parser.parse_args()

    write_tree_with_polytomies(args.tree, args.output)
