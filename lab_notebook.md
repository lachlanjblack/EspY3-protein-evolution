# Project log

## TODO

### Housekeeping
[ ] Repeat previous steps locally
- [ ] Make environment/setup file (and update when needed)
- [ ] Upload and annotate all files produced so far

### Project Work
- [ ] Align 50% identity uniprot homologues with MAFFT
- [ ] Manually clean MSA - remove rows with gaps and remove leader regions
- [ ] Retreive pfam homologues (10k sequences)
- [ ] Use uniprot MSA as seed for pfam homologue alignment (minimum pairwise identity 60%)
- [ ] 



## Initial Setup - 27 Oct 2025
- Created GitHub repository and folders
- Installed git and basic command-line tools

## Prior to 27 Oct 2025
- Downloaded EspY3 50% identity UniRef clusters (`EspY3_50.fasta`)
- Aligned sequences with MAFFT (`EspY3_50.MAFFT.fasta`)
- Ran trimAL on allignment to remove empty columns or incomplete sequences (`EspY3_50.MAFFT.trimal_18del.fasta`)
  - Named `~.trimal_18del` because I made two files with the same name. One was the trimal output with bad parameters.
  - trimAL parameters: 0.7 0.7


