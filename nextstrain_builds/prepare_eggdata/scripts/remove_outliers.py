import argparse
from Bio import SeqIO
import pandas as pd

def remove_outliers(outliers_file, sequences_file, metadata_file, out_sequences_file, out_meta_file):

    # read in outliers
    with open(outliers_file, 'r') as f:
        outliers = set(line.strip() for line in f if line.strip())

    # filter FASTA file to remove outliers
    filtered_seqs = (record for record in SeqIO.parse(sequences_file, "fasta")
                     if record.id not in outliers)
    SeqIO.write(filtered_seqs, out_sequences_file, "fasta")

    # Filter metadata file
    metadata_df = pd.read_csv(metadata_file, sep="\t")
    filtered_metadata_df = metadata_df[~metadata_df['strain'].isin(outliers)]
    filtered_metadata_df.to_csv(out_meta_file, sep="\t", index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Removes outlier strains that should've gotten removed during augur filter but didn't",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--outliers', type=str, required=True,
                        help="text file of outliers")

    parser.add_argument('--sequences', type=str, required=True,
                        help="fasta data file")

    parser.add_argument('--metadata', type=str, required=True,
                        help="tsv metadata file")

    parser.add_argument('--out-sequences', type=str, required=True,
                        help="output sequences file")

    parser.add_argument('--out-metadata', type=str, required=True,
                        help="output tsv metadata file")


    args = parser.parse_args()

    remove_outliers(args.outliers, args.sequences, args.metadata, args.out_sequences, args.out_metadata)
