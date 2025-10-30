# Trimal 

Here contains data which has been processed by trimal. Outwith the archive folder lies the final processed data.
So far, this is limited to 'espy3-50.mafft.r07s08.trim.fasta'
As entailed by the name, this file was created by running trimAL on the MSA. As there were a lot of truncated sequences and they were not dispersed throughout, I set the parameters manually, `-resoverlap 0.7 and -seqoverlap 80` worked best, leaving only 3 sequences truncated by ~80 amino acids. I manually removed these, and two adjacent empty columns, in Jalview
