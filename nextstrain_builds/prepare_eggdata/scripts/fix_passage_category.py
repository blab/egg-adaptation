import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-metadata', type=str, required=True, help="egg metadata file")
    parser.add_argument('--output-metadata', type=str, required=True, help="metadata file after changing all passage_category to egg")

    args = parser.parse_args()

    meta = pd.read_csv(args.input_metadata, sep='\t')

    meta['passage_category'] = 'egg'

    meta.to_csv(args.output_metadata, sep='\t', index=False)
