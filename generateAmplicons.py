# Script for generating the amplicon sequences which *in theory* should
# be included within the PBEH2 PacBio library.
# Ethan Holleman
# 3/25/24

from pydna.amplify import pcr
from pydna.readers import read
import argparse
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
                    prog='Generate amplicons',
                    description='Generate PacBio sequencing amplicons'
    )
    parser.add_argument('sample_tsv', help='Path to sample tsv')
    parser.add_argument('primer_tsv', help='Path to primer list tsv file')
    parser.add_argument('plasmid_dir', help='Path to plasmid reference dir')
    parser.add_argument('output_dir', help='Output directory path')

    args = parser.parse_args()

    return args


def make_plasmid_paths(sample_row):
    """Generate a list of paths to plasmid (genbank) files. The list will have 
    at least 1 plasmid in it. Lists with many plasmids represent a mix
    of plasmids that was used in the same sample and therefore amplified with
    the same primers.

    Args:
        sample_row (pd.Series): pandas dataframe row representing the sample

    Returns:
        list: List of plasmid file paths (Path objects)
    """

    plasmid_loc = Path(sample_row['Plasmid Loc'])
    plasmid_path = PLASMID_DIR.joinpath(plasmid_loc)

    if plasmid_path.is_dir():  # if it is a directory, represents mix of plasmids
        plasmids = [
            each_path for each_path in plasmid_path.iterdir() 
            if each_path.suffix == '.gb'
            ]
        if str(plasmid_path) == 'plasmids/VRv2Mix':
            print(plasmids, '\n')
    else:
        plasmids = [plasmid_path]
    
    return plasmids



def select_primer_pair(sample_row, primer_df):

    
    def look_for_primer(primer_name):
        primer_search = primer_df.loc[primer_df['name'] == primer_name]
        try:
            assert len(primer_search) == 1
        except AssertionError:
            print('Failed primer search')
            print('Primer name:', primer_name)
            exit(1)
        return primer_search.iloc[0]
    
    
    fwd_primer = look_for_primer(sample_row.Fwd_primer)
    rev_primer = look_for_primer(sample_row.Rev_primer)


    return fwd_primer, rev_primer



def amplify_plasmids(sample_row):


    fwd_primer_seq = sample_row.primer_pair[0].seq
    rev_primer_seq = sample_row.primer_pair[-1].seq

    pcr_products = []

    def make_amplicon_name(plasmid_path):
        # Name attribute of plasmid file not always correct or very concise
        # so make some adjustments here so things are more readable and cleaner
        name = None
        if PLASMID_DIR.joinpath(Path(sample_row['Plasmid Loc'])).is_file():
            # Single plasmid sample so just use the name that is listed here
            name = Path(sample_row['Plasmid Loc']).stem
        
        else:
            name = plasmid_path.stem
        
        name = f'{name}.{sample_row.Sample_ID}'

        return name
            

    for each_plasmid in sample_row.plasmid_paths:

        plasmid_record = read(str(each_plasmid))
        amplicon = pcr(fwd_primer_seq, rev_primer_seq, plasmid_record)
        amplicon_name = make_amplicon_name(each_plasmid)
        #print(amplicon_name, len(amplicon))
        pcr_products.append((amplicon, amplicon_name))

    return pcr_products


def write_amplicons(sample_df, fa_path, gb_path):
    
    amplicon_fas = []

    for i, each_row in sample_df.iterrows():

        for each_amplicon, each_name in each_row.amplicons:

            output_fa = str(Path(fa_path).joinpath(each_name)) + '.fa'
            output_gb = str(Path(gb_path).joinpath(each_name)) + '.gb'

            each_amplicon.locus = each_amplicon.locus.replace(' ', '.')
            each_amplicon.definition = each_name

            SeqIO.write(each_amplicon, output_gb, 'gb')
            SeqIO.write(each_amplicon, output_fa, 'fasta')

            amplicon_fas.append(each_amplicon)
    
    return amplicon_fas


def main():

    args = parse_args()

    # Read in all needed data
    sample_df = pd.read_csv(args.sample_tsv, sep='\t')
    primer_df = pd.read_csv(args.primer_tsv, sep='\t')

    # Add column containing paths to files for all plasmids in the sample
    # represented by each row of the dataframe
    sample_df['plasmid_paths'] = sample_df.apply(
        lambda row: make_plasmid_paths(row), axis=1
    )

    # Add a column with the primer pair info used to amplify each sample
    sample_df['primer_pair'] = sample_df.apply(
        lambda row: select_primer_pair(row, primer_df),
        axis=1
    )

    # Add a column with list of the actual PCR products (amplicons) that
    # should in theory be included in this sample
    sample_df['amplicons'] = sample_df.apply(
        lambda row: amplify_plasmids(row),
        axis=1
    )

    amplicons_gb = Path(args.output_dir).joinpath('genbank')
    amplicons_fa = Path(args.output_dir).joinpath('fasta')

    amplicons_gb.mkdir(parents=True, exist_ok=True)
    amplicons_fa.mkdir(parents=True, exist_ok=True)

    amplicons = write_amplicons(sample_df, amplicons_fa, amplicons_gb)
    #lengths = get_sequence_lengths(amplicons)
    #plot_histogram(lengths)

    





if __name__ == '__main__':
    main()