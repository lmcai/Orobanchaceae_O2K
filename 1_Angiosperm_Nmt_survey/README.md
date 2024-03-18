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

4. For absences, (1) HMMER search on CDS and (2) BLASTN search on genome seq was used to verify such losses.

4.1 HMMER

For a candidate gene loss in B14.5a (AT5G08060), a profile was build first based on DNA alignment
```
hmmbuild AT5G08060.hmm AT5G08060.aln.fas
```
Then search the transcriptome of the query species
```
hmmsearch AT5G08060.hmm Theobroma_cacao.cds.fas
```

4.2 BLASTN on genome, the GUI app `seq_extract_GUI.py` is very handy for extracting sequence with IDs
```
makeblastdb -in Theobroma_cacao.genome.fas -dbtype nucl -out Tca
blastn -db Tca -query AT5G08060.fas -evalue 1e-5 -outfmt 6 >AT5G08060.Tca.blast
```