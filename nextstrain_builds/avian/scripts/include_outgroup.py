#!/usr/bin/env python3
"""
Add an outgroup sequence to the data file, so that tree will be properly rooted.
Add this sequence to the force include list.
"""
import argparse
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd


def add_outgroup(ha_type):
    # dictionary to map ha_type to what ha type should be used as an outgroup (use an HA in the other group)
    outgroup_mapper = {'H1': 'H3', 'H2': 'H3', 'H3': 'H1', 'H4': 'H1', 
                    'H5': 'H3', 'H6': 'H3', 'H7': 'H1', 'H8': 'H3', 
                    'H9': 'H3', 'H10': 'H1', 'H11': 'H3', 'H12': 'H3', 
                    'H13': 'H3', 'H14': 'H1', 'H15': 'H1', 'H16': 'H3'}
    
    # read in the fasta data for the outgroup sequences
    with open('config/outgroup_sequences.json') as f:
        outgroup_info = json.load(f)

    # initiate list to store seq records from this ha type 
    # start by adding the outgroup
    records_plus_outgroup = []

    # get outgroup sequence information for this ha type
    outgroup_for_this_ha = outgroup_mapper[ha_type]

    # there may be several outgroup sequences listed for each ha type
    outgroups = outgroup_info[outgroup_for_this_ha]

    # get the accession number or numbers for the outgroup
    outgroup_accessions = []
    for o in outgroups:
        outgroup_fasta_id = o['id']
        outgroup_fasta_seq = o['seq']
        records_plus_outgroup.append(SeqRecord(Seq(outgroup_fasta_seq), id=outgroup_fasta_id, description=''))
        outgroup_accessions.append(str(outgroup_fasta_id.split('|')[0]))


    # read in the fasta file for this HA type
    sequences_fasta = f"data/{ha_type}_ncbi_sequences.fasta"
    for record in SeqIO.parse(sequences_fasta, "fasta"):
        records_plus_outgroup.append(record)

    # write a new fasta file, including the outgroup sequence
    SeqIO.write(records_plus_outgroup, f"results/{ha_type}/sequences_w_outgroup.fasta", "fasta")


    # add the outgroup sequence to the force include list
    with open(f'config/{ha_type}/include_strains.txt', 'a') as file:
        for oa in outgroup_accessions:
            file.write(f'{oa}\n')






if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--ha-type",
        help="HA type number.")


    args = parser.parse_args()

    add_outgroup(args.ha_type)
