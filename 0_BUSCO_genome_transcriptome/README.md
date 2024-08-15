# Genome/Transcriptome quality assessment

1. Use the embryophyta_odb10 database from BUSCO v5.7.0 to assess the completeness of the dataset

Install BUSCO
```
singularity pull docker://ezlabgva/busco:v5.7.0_cv1
```
Then execute the command for each species under either `genome` or `transcriptome` mode depending on the species
```
singularity exec busco_v5.7.0_cv1.sif busco -i CDS/Aeginetia_trinity.Trinity.fasta -l embryophyta_odb10 -m transcriptome -c 1 -o Aeginetia
```