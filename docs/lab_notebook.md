# Project log

## TODO

### Housekeeping
- [ ] Repeat previous steps locally
- [ ] Upload and annotate all files produced so far

### Project Work
- [ ] Manually clean MSA - remove rows with gaps and remove leader regions
- [ ] Retreive pfam homologues (10k sequences)
- [ ] Use uniprot MSA as seed for pfam homologue alignment (minimum pairwise identity 60%)
- [ ] 



## Initial Setup - 27 Oct 2025
- Created GitHub repository and folders
- Installed git and basic command-line tools

## Prior to 27 Oct 2025 (Ignore, redid everything)
- Downloaded EspY3 50% identity UniRef clusters (`espy3-50.fasta`) ("https://www.uniprot.org/uniprotkb/Q8XE85/entry#similar_proteins")
- Aligned sequences with MAFFT (`espy3-50.MAFFT.fasta`)
- Ran trimAL on allignment to remove empty columns or incomplete sequences (`espy3_50.mafft.trimal_18del.fasta`)
  - Named `~.trimal_18del` because I made two files with the same name. One was the trimal output with bad parameters.
  - trimAL parameters: 0.7 0.7

## 29 Oct 2025
- Organised directory
- Re-ran MAFFT locally on 'espy3-50.fasta'
- MAFFT usage
  - "/opt/anaconda3/envs/project_2025_py313/bin/mafft"  --localpair  --maxiterate 16 --reorder "/Users/lachlanblack/Documents/GitHub/EspY3-protein-evolution/data/raw/espy3-50.fasta" > "/Users/lachlanblack/Documents/GitHub/EspY3-protein-evolution/data/cleaned/espy3-50.mafft.fasta"
  - I selected Fasta format - Sorted
  - I chose L-INS-i strategy as recommended in the MAFFT manual for alligning <200 sequences

- Used trimAL to remove spurious sequences or poorly aligned regions from 'espy3-50.mafft.fasta.
- 

- Downloaded similar sequences from Pfam (pfam-espy3.fasta)("https://www.ebi.ac.uk/interpro/protein/UniProt/Q8XE85/similar_proteins/#table")
- Aligned sequences using MAFFT
