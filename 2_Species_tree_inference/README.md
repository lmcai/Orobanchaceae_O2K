# Species tree inference

1. Use the maximum inclusion algorithm from Yang and Smith (2016) to get one-to-one ortholog from the curated angiosperm Nmt dataset

```
python programs/yangya-phylogenomic_dataset_construction-489685700c2a/prune_paralogs_MI.py Nmt-angio/ .tree inf inf 5  Nmt-angio-1to1/ 
```

2. Align protein seq first and then back translate to codon, finally infer a phylogeny

```
mafft --genafpair --maxiterate 1000 in.faa > out.aln.faa
pal2nal.pl out.aln.faa cds.fas >out.codon
iqtree2 -s codon.concatenated.fas -p gene.partition -B 1000 -nt 4 --prefix Nmt.sp.tre
```