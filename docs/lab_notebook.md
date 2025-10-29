# Project log

## TODO

### Housekeeping
- [x] Repeat previous steps locally
- [x] Upload and annotate all files produced so far

### Project Work
- [x] Manually clean MSA - remove rows with gaps and remove leader regions
- [x] Retreive pfam homologues (10k sequences)
- [ ] Update tools folder and add readmes
- [ ] Figure out how to format .md files
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

- Installed trimAL v1.5.rev0 in new Tools folder
- Used trimAL to remove spurious sequences or poorly aligned regions from 'espy3-50.mafft.fasta.
- I initially used parameters -resoverlap 0.5 -seqoverlap 0.6 as my dataset has 50% sequence identity and I don't want to erase too much information. I will check the trimmer allignment in Jalview
- trimal -in /Users/lachlanblack/Documents/GitHub/EspY3-protein-evolution/data/cleaned/espy3-50.mafft.fasta -out /Users/lachlanblack/Documents/GitHub/EspY3-protein-evolution/data/cleaned/espy3-50.mafft.trimAL-0.6-0.50.fasta -keepheader -resoverlap 0.6 -seqoverlap 0.5
- It was too gappy, so I ran stricter parameters of -resoverlap 0.7 and seqoverlap 0.6
- Turns out I didn't give trimal an trimming algorithm, so I am going to run the first one again with -gappyout
- The sequences are highly conserved with truncated sequences at the ends, rather than dispersed  throughout. Therefore, I chose the parameters -resoverlap 0.7 -seqoverlap 80. This removed a majority of the truncated sequences, leaving just 3 sequences truncated by 80 or so residues which I manually removed. Similarly, there were two adjacent columns that were empty except from in a single sequence. I manually removed these too. 
- The trimal prompt:
trimal -in espy3-50.mafft.fasta -out {output.fasta} -keepheader -resoverlap 0.7 -seqoverlap 80 htmlout- {output.html}
- 

- Downloaded similar sequences from Pfam (pfam-espy3.fasta)("https://www.ebi.ac.uk/interpro/protein/UniProt/Q8XE85/similar_proteins/#table")
- Aligned sequences using MAFFT
