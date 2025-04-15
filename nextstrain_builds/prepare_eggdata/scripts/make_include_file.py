import argparse
from Bio import SeqIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Makes a list of all egg strain names so that we can get the pair of each sequence where available",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--filtered-seqs', type=str, required=True,
                        help="egg sequences after all filtering")

    parser.add_argument('--output', type=str, required=True,
                        help="output file name")

    args = parser.parse_args()

    # get all egg strains and the name of the non-egg-passaged pair
    egg_strains = []
    pair_strains = []
    with open(args.filtered_seqs) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            egg_strains.append(record.id)
            pair_strains.append(record.id.removesuffix('-egg'))

    all_strains = egg_strains + pair_strains
    # write include file
    with open(args.output, 'w') as f:
        for item in all_strains:
            f.write(item + '\n')

    
