# Angiosperm Nmt survey pipeline

1. Use Arabidopsis thaliana ID in `ATH.Nmt.fas` as bait to identify the original orthogroup, output into tab-delimited format

2. Extract sequences and build a draft phylogeny

Extract seq for each species
```
python extract_gen_fam_seq.py Ortholog.tsv Ath.faa Ath.cds.fas 3 a Ath
```
After all done, align and infer a tree
```
mafft --genafpair --maxiterate 1000 in.fas > aln.fas
iqtree2 -s aln.fas -B 1000 -nt 4
```
3. Manual inspect the phylogeny to identify paralog. If paralog is present, then use the interactive phylogeny editor `TreeGraph_2.15.0-887_beta` to edit the phylogeny. Then use `output_tips_in_order_from_subtree.py` to update the gene accessions.

4. For absences, 