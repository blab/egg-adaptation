#!/usr/bin/env python3
"""
Prune the outgroup sequence from the newick tree.
The end result will be that each HA-type-specific tree will be properly rooted (on the outgroup), 
but that the auspice tree will show only sequences from the given type
"""
import argparse
import json
from Bio import Phylo
from augur.utils import json_to_tree


def prune_outgroup(ha_type, tree_w_outgroup, auspice_w_outgroup, output_file):
    """
    Find all tips that are not of the given ha_type 
    and then prune the tree to exclude these and the branches leading to them
    """
    
    #read in the tree
    with open(auspice_w_outgroup, 'r') as f:
        tree_json = json.load(f)
                
    #put tree in Bio.phylo format
    tree = json_to_tree(tree_json)
    
    outgroup_nodes = []
    
    # find tips to prune (from outgroup)
    for node in tree.find_clades(terminal = True):
        if node.node_attrs['ha_type']['value'] != ha_type:
            outgroup_nodes.append(node.name)
            
    
    # read in the newick tree
    newick_tree = Phylo.read(tree_w_outgroup, 'newick')

    # prune out the outgroup nodes
    for x in outgroup_nodes:
        newick_tree.prune(x)


    # save pruned newick tree
    Phylo.write(newick_tree, output_file, 'newick')
        
   







if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--ha-type",
        help="HA type number.")
    parser.add_argument("--tree-w-outgroup",
        help="Path to the newick tree, including the outgroup.")
    parser.add_argument("--auspice-w-outgroup",
        help="Path to the auspice tree, including the outgroup.")
    parser.add_argument("--output-file",
        help="Output file path.")
    


    args = parser.parse_args()

    prune_outgroup(args.ha_type, args.tree_w_outgroup, args.auspice_w_outgroup, args.output_file)
