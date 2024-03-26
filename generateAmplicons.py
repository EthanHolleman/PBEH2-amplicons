# Script for generating the amplicon sequences which *in theory* should
# be included within the PBEH2 PacBio library.
# Ethan Holleman
# 3/25/24

from pydna.amplify import pcr
import pandas as pd
from pathlib import Path


# Filepaths, these should all be included in the repo and should not change
# and are therefore hardcoded into the amplicon generation script
# ==============================================================================
SAMPLE_DF_PATH = 'PBEH2_samples.txt'  # location of sample assignment data
PRIMER_DF_PATH = 'primers.tsv'  # list of primers used (and not used) and names
PLASMID_DIR = Path('plasmids')  # location with all plasmid data files (gb)


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



def amplify_plasmids(sample_plasmids, sample_primers):


    fwd_primer_seq = sample_primers[0].seq
    rev_primer_seq = sample_primers[-1].seq

    pcr_products = []

    def make_amplicon_name(plasmid):
        return f'{plasmid}.PBEH2-{sample_id}'

        

    for each_plasmid in sample_plasmids:
        amplicon = pcr(fwd_primer_seq, rev_primer_seq, each_plasmid)
        amplicon_name = make_amplicon_name(each_plasmid)
        pcr_products.append(  # store each template and its amplicon
            (each_plasmid, amplicon, amplicon_name)
        )

    return


def make_sample_to_amplicon_dict():
    pass

    




def main():

    # Read in all needed data
    sample_df = pd.read_csv(SAMPLE_DF_PATH, sep='\t')
    primer_df = pd.read_csv(PRIMER_DF_PATH, sep='\t')

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

    





if __name__ == '__main__':
    main()