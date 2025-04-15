import pandas as pd
from Bio import SeqIO
import argparse


def prune_metadata(unpruned_metadata, filtered_seq_file, output_file_name):
    """
    Because augur parse is run before subsampling, the metadata file for the
    background sequences  will contain way more strains than will actually be used in the build.
    Prune it to just have metadata for the sequences that passed subsampling
    """

    # Fet list of strain names from the filtered sequence file
    strain_names = [record.id for record in SeqIO.parse(filtered_seq_file, "fasta")]

    # Load the full metadata file
    metadata_df = pd.read_csv(unpruned_metadata, sep='\t')

    # Filter to only include rows that passed subsampling
    pruned_df = metadata_df[metadata_df['strain'].isin(strain_names)]

    # Save the pruned metadata to a new file
    pruned_df.to_csv(output_file_name, sep='\t', index=False)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--unpruned-meta', type=str, required=True,
                        help="the metadata file produced from augur parse, before subsampling")
    parser.add_argument('--filtered-seqs', required=True, help="fasta produced after filtering")
    parser.add_argument('--output', required=True, help="pruned metadata file")

    args = parser.parse_args()

    prune_metadata(args.unpruned_meta, args.filtered_seqs, args.output)
