#!/usr/bin/env python3
"""
Write FASTA and GFF references from a Genbank file.
Assumes the genes are labeled as CDS
"""
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def make_reference_files(ha_type):
    
    fasta_record = []
    gene_locations = {}

    with open(f"config/{ha_type}/reference_{ha_type}.gb", 'r') as handle:
        for record in SeqIO.parse(handle, 'gb'):
            sequence = record.seq
            accession = record.name
            fasta_record.append(SeqRecord(Seq(sequence), id=accession, description = ''))
            for f in record.features:
                if f.type == 'source':
                    source_location = f.location
                if f.type == 'CDS':
                    location = f.location
                    gene = f.qualifiers['gene'][0]
                    gene_locations[gene] = location
            
    with open(f"config/{ha_type}/reference.fasta", 'w') as f:
        SeqIO.write(fasta_record, f, "fasta")  
            
    
    genemap_lines = [f'##sequence-region {accession} {source_location.start+1} {source_location.end}\n', 
                    f'{accession}\tfeature\tsource\t{source_location.start+1}\t{source_location.end}\t.\t+\t.\tgene=nuc\n']

    for g, l in gene_locations.items():
        genemap_lines.append(f'{accession}\tfeature\tgene\t{l.start+1}\t{l.end}\t.\t+\t.\tgene_name={g}\n')

    # write genemap.gff               
    with open(f"config/{ha_type}/genemap.gff","w") as f:
        for line in genemap_lines:
            f.write(line)
   







if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--ha-type",
        help="HA type number.")
    
    args = parser.parse_args()

    make_reference_files(args.ha_type)
