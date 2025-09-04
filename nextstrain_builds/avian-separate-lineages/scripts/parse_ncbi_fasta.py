#!/usr/bin/env python3
"""
Streamlines FASTA header that was downloaded from NCBI virus to have relevant info and proper headers.
Gets HA type from the genotype.
Assigns a region based on the country.
Pulls a strain name from the Genbank name.
"""
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re


def parse_ncbi_fasta(raw_fasta, geolocations_path):
    # read in the general geography rules
    general_geo_df = pd.read_csv(geolocations_path, sep='\t', header=None)[1]
    geo_locs =list(general_geo_df)
    country_to_region_mapper = {}
    for g in geo_locs:
        r = g.split('/')[0]
        c = g.split('/')[1]
        country_to_region_mapper[c] = r

    # add custom entries that don't exist
    country_to_region_mapper['Mongolia'] = 'Asia'
    country_to_region_mapper['Korea'] = 'Asia'
    country_to_region_mapper['Singapore'] = 'Asia'
    country_to_region_mapper['Bhutan'] = 'Asia'
    country_to_region_mapper['Cambodia'] = 'Asia'
    country_to_region_mapper["Cote d'Ivoire"] = 'Africa'
    country_to_region_mapper["Czechoslovakia"] = 'Europe'
    country_to_region_mapper["Kyrgyzstan"] = 'Asia'
    country_to_region_mapper["Mali"] = 'Africa'
    country_to_region_mapper["Libya"] = 'Africa'
    country_to_region_mapper["Tajikistan"] = 'Asia'
    country_to_region_mapper["Jordan"] = 'Asia'
    country_to_region_mapper["Antarctica"] = 'Antarctica'
    country_to_region_mapper["North Korea"] = 'Asia'
    country_to_region_mapper["Reunion"] = 'Africa'
    country_to_region_mapper["Greenland"] = 'Europe'

    # all influenza subtypes to care about
    ha_subtypes = ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16']


    # map influenza segment numbers to the name of the segment
    segment_mapper = {'1': 'PB1', '2': 'PB2', '3': 'PA', '4': 'HA', '5': 'NP', '6': 'NA', '7': 'MP', '8': 'NS'}

    # store new sequence records with the curated names
    # organize by ha subtype
    new_records = {haX:[] for haX in ha_subtypes}

    for record in SeqIO.parse(raw_fasta, "fasta"):
        [accession, genbank_name, organism_name, SRA_accession, submitters,
         organization, unknown_column, release_date, isolate, species, length, genotype,
         segment, publications, geo_location, country, host, isolation_source,
         collection_date] =  record.description.split('|')
        
        #remove the version number from accession
        # accession = str(accession.split('.')[0])

        # get strain name from genbank name
        # usually in the format
        # Influenza A virus (A/Northern goshawk/Nara/2912B005/2020) viral cRNA, segment 6, complete sequence
        # so find the strain name between the parentheses
        if '(' in genbank_name:
            strain_name = re.findall(r'\(.*?\)', genbank_name)[0][1:-1]
        else:
            # check whether standard flu nomenclature is included in the name
            if '/' in genbank_name:
                strain_at_front = genbank_name.split('Influenza A virus ')[1]
                split_on_slash = strain_at_front.split('/')
                year = split_on_slash[len(split_on_slash)-1][:4]
                strain_name = '/'.join(split_on_slash[:(len(split_on_slash)-1)])+f'/{year}'
            # if not, give the strain the accession number for a name
            else:
                strain_name = accession

        # get rid of the subtype designation in strain name if it is there
        if '(' in strain_name:
            strain_name = strain_name.split('(')[0]

        # get ha_type from the genotype
        ha_type = genotype.split('N')[0]

        # map segment number to segment name
        if segment in segment_mapper.keys():
            segment_name = segment_mapper[segment]
        else:
            segment_name = ''

        # manual fixes
        if country == 'Viet Nam':
            country = 'Vietnam'

        # get region from the country
        if country in country_to_region_mapper.keys():
            region = country_to_region_mapper[country]
        elif country == '':
            region = ''
        else:
            region = 'not mappable'

        # make new fasta header
        new_header = '|'.join([accession.strip(), strain_name, ha_type, segment, segment_name, collection_date,
                               genotype, species, host,
                               country, region])
        if ha_type in ha_subtypes:
            new_records[ha_type].append(SeqRecord(record.seq, id=new_header, description = ''))


    for haX in ha_subtypes:
        SeqIO.write(new_records[haX], f"data/{haX}_ncbi_sequences.fasta", "fasta")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--fasta",
        help="path to the fasta downloaded from NCBI.")
    parser.add_argument("--geolocations",
        help="path to the file listing region/country locations.")


    args = parser.parse_args()

    parse_ncbi_fasta(args.fasta, args.geolocations)
