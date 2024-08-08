# Species tree inference

1. Use trimAL to clean up the codon alignment and then build a draft gene tree with iqtree
```
sh trimal_iqtree.sh OG0004221
```
Then root the gene trees use the following python scripts
```
import os
files=os.listdir('.')
from ete3 import Tree
files=[i for i in files if i.endswith('.treefile')]
for f in files:
	t=Tree(f,format=1)
	tips=[node.name for node in t]
	tips=[i for i in tips if not (i.startswith('scaffol') or i.startswith('TAIR10'))]
	t.prune(tips)
	t.write(outfile=f.split('.')[0]+'.tree',format=0)
```

2. Use the maximum inclusion algorithm from Yang and Smith (2016) to get one-to-one ortholog from the curated angiosperm Nmt dataset

```
python programs/yangya-phylogenomic_dataset_construction-489685700c2a/prune_paralogs_RT.py Nmt-angio/ .tree Nmt_combined_RT 5 taxon
```

3. Concatenate one-to-one sequences and infer a phylogeny

```
import os
files=os.listdir('Nmt_combined_RT')
from ete3 import Tree
from Bio import SeqIO
files=[i for i in files if i.endswith('.ortho1.tre')]
for file in files:
	t=Tree('Nmt_combined_RT/'+file)
	valid_ids=[node.name for node in t]
	try:
		recs=SeqIO.index('ready4paml/'+file.split('.')[0]+'.trimal.fas','fasta')
		out=open(file.split('.')[0]+'.fas','w')
		for i in valid_ids:
			d=SeqIO.write(recs[i],out,'fasta')
	except ValueError:pass

```

4. Use the commands in `trimal_iqtree.sh` to infer the species tree under the GHOST mixture model.
