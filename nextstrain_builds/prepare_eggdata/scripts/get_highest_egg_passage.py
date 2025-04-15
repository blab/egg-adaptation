import re
from Bio import SeqIO
import argparse

min_length_by_segment = {
    "ha": 950,
    "na": 950,
    "pb1": 1500,
    "pb2": 1500,
    "pa": 1500,
    "np": 1000,
    "ns": 500,
    "mp": 500
}

def convert_nomenclature_to_number(egg_passage_nomenclature):
    """
    Using https://pmc.ncbi.nlm.nih.gov/articles/PMC6599686/
    as resource
    / indicates trade between institutions
    """

    # extract numbers after E
    # ex: for E6, will return 6. For E6/E2, will return 6,2, (same for E6+E2 or E6E2)
    # will ignore "passage", which is sometimes the only info
    # will also ignore cell passaging (denoted by 'C')
    # for amiguous situations like E5+1+E1, it will say 6
    E_numbers = [int(num) for num in re.findall(r'E(\d+)', egg_passage_nomenclature)]
    egg_integer= sum(E_numbers)


    return egg_integer

def find_highest_egg_passage(list_of_seqs):
  """
  """
  # convert passaging nomenclature to an integer
  list_of_seqs = [{k:convert_nomenclature_to_number(v)} for x in list_of_seqs for k,v in x.items()]

  # Find the accession greatest egg passage value
  max_passage_dict = max(list_of_seqs, key=lambda d: list(d.values())[0])

  # Extract the corresponding key (accession)
  max_passage_accession = list(max_passage_dict.keys())[0]

  return max_passage_accession

def resolve_duplicates(name_and_passage):
    """
    Find strains that have multiple different egg_passaged entries
    Find the highest egg-passage number amongst those
    Return a dict with strain name as key and accession number as the key
    where accession number maps to the sequence with highest egg passage number
    (or the *only* sequence, for the singletons)
    """

    duplicates = {}
    singletons = {}

    for k,v in name_and_passage.items():
        if len(v)>1:
            # need to find which has the highest egg-passage number
            highest_eggpassage = find_highest_egg_passage(v)
            duplicates[k] = highest_eggpassage
        else:
            # for each strain, accession is key and egg-passage is value
            singletons[k] = list(v[0].keys())[0]

    strain_to_accession = {**duplicates, **singletons}

    return strain_to_accession

def find_max_egg_passage_per_strain(virus, segment, input_seq_file, exclude):
    """
    For each strain, find the sequence that was passaged in eggs the most times
    Many strains don't have duplicates, so the max could just be the only sequence
    Since they will have the same name, select them on their accession

    Require that sequence as long as the min_length, to ensure full segment
    """
    # get list of strains to exclude
    with open(exclude, 'r') as file:
        exclude_strains = file.readlines()
        # Remove newline characters
        exclude_strains = [line.strip() for line in exclude_strains]

    # for each strain name, keep track of all
    name_and_passage = {}

    with open(input_seq_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            egg_passage = record.id.split('|')[11]
            strain_name = record.id.split('|')[0]
            egg_strain_name = strain_name+'-egg'
            accession = record.id.split('|')[4]
            # check that its not a strain we want to exclude
            if strain_name not in exclude_strains and egg_strain_name not in exclude_strains:
                # check that it meets length requirement
                if len(record.seq) >= min_length_by_segment[segment]:

                    if strain_name in name_and_passage:
                        name_and_passage[strain_name].append({accession:egg_passage})
                    else:
                        name_and_passage[strain_name]=[{accession:egg_passage}]

    # get the accession we want to use for each strain
    strain_to_accession = resolve_duplicates(name_and_passage)


    return strain_to_accession

def save_fasta(virus, segment, input_seq_file, exclude, output_seq_file):
    """
    Save a .fasta file that will include one sequence per strain (the one with highest number of egg passages)
    """

    strain_to_accession = find_max_egg_passage_per_strain(virus, segment, input_seq_file, exclude)

    accessions_to_include = list(strain_to_accession.values())

    new_records = []

    with open(input_seq_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            accession = record.id.split('|')[4]
            if accession in accessions_to_include:
                # need to edit strain name to have -egg suffix
                record_parts = record.id.split('|')
                record_parts[0] += '-egg'
                record.id = '|'.join(record_parts)
                new_records.append(record)

    with open(output_seq_file, "w") as output_handle:
        SeqIO.write(new_records, output_handle, "fasta")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--lineage', type=str, required=True,
                        help="seasonal flu lineage")
    parser.add_argument('--segment', required=True, help="flu segment")
    parser.add_argument('--input-seqs', type=str, required=True, help="file with all egg seqs")
    parser.add_argument('--exclude', type=str, required=True, help="list of outliers to exclude")
    parser.add_argument('--output-seqs', type=str, required=True, help="filename for selected seqs")

    args = parser.parse_args()

    save_fasta(args.lineage, args.segment, args.input_seqs, args.exclude, args.output_seqs)
