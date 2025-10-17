# PBEH2-amplicons

Script for generating the theoretical amplicons which should be included in each
sample for the PBEH2 SMRT-cells.

`generateAmplicons.py` reads in sample information from the included sample tables 
and primer sequence files (`PBEH2_samples.tsv` and `primers.tsv`) and then preforms 
in silico PCR to generate amplicons for each sample.

Some samples contain many different plasmids so the final output (which are 
both genbank and fasta files of the amplicons) are labeled with the following 
scheme

`{Plasmid name}.PBEH2-{Sample ID number}.{file type}`

These files are located in the amplicons directory.

`generateAmplicons.py` should be run from the repos root directory to avoid 
filepath issues.

## Example command 

```
python generateAmplicons.py PBEH3_samples.tsv primers.tsv plasmids PBEH3
_amplicons
```
