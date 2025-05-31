# N-terminal targeting peptides prediction

N-terminal targeting peptides were identified using the webserver of TargetP-2.0 (https://services.healthtech.dtu.dk/services/TargetP-2.0/). This will generate a gff3 file with the cleavage site (CS) reported.

To remove these peptides from alignments, 

# Molecular Evolution with HYPHY

1. label phylogeny (foreground vs background branches) using https://phylotree.hyphy.org/
See instructions https://hyphy.org/tutorials/phylotree/

2. Use interactive command line version of HYPHY
Install HYPHY ing conda
```
hyphy -i relax --alignment CI.masked.fas --tree Nmt.rooted.labeled.tre
```
