# Verification of gene losses

1. BUSCO evaluation of genome and transcriptome quality

Installation
```
singularity pull docker://ezlabgva/busco:v5.7.0_cv1
```
Run BUSCO
```
singularity exec busco_v5.7.0_cv1.sif busco -i CDS/Aeginetia_trinity.Trinity.fasta -l embryophyta_odb10 -m transcriptome -c 1 -o Aeginetia
singularity exec busco_v5.7.0_cv1.sif busco -i CDS/Ptr.genome.fasta -l eukaryota_odb10 -m genome -c 1 -o Ptr
```

2. Calculate pairwise DNA sequence divergence and conduct phylogenetic ANOVA test using `sequence_divergence_phyANOVA.R`. 

Input species tree `Nmt.con.rooted.tre` and Input trait values `seq_div_missing_gene.csv`.

3. Fisher's exact test `fisher_test.py` to evaluate the significance of N-mt gene losses in Orobancheae while correct for missing eukaryote BUSCOs. To calculate the corrected p-values, 10000 empirical permutation was used.