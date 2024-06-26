from ete3 import Tree
import os
from Bio import SeqIO

files=os.listdir('.')
files=[i for i in files if i.endswith('.treefile')]

for f in files:
	t=Tree(f)
	#reroot
	amb=[node.name for node in t if node.name.startswith('Amb')]
	if len(amb):
		root_node=t&amb[0]
		t.set_outgroup(root_node)
	else:
		cin=[node.name for node in t if node.name.startswith('Cin')]
		if len(cin):
			root_node=t&cin[0]
			t.set_outgroup(root_node)
	to_keep=[]
	shortest={}
	sp_pool=[]
	for node in t:
		sp=node.name.split('_')[0]
		if not sp in sp_pool:
			to_keep.append(node.name)
			shortest[sp]=node.get_distance(root_node)
			sp_pool.append(sp)
		else:
			if node.get_distance(root_node)<shortest[sp]:
				previous=[i for i in to_keep if i.startswith(sp)]
				to_keep.remove(previous[0])
				to_keep.append(node.name)
				shortest[sp]=node.get_distance(root_node)
			elif node.get_distance(root_node)==shortest[sp]:
				if not node.name.endswith('_2'):
					previous=[i for i in to_keep if i.startswith(sp)]
					if len(node.name) >len(previous):
						to_keep.remove(previous[0])
						to_keep.append(node.name)
						shortest[sp]=node.get_distance(root_node)
	recs=SeqIO.index(f.split('.')[0]+'.filter.fas','fasta')
	out=open(f.split('.')[0]+'.one2one.fas','w')
	for i in to_keep:
		try:
			d=SeqIO.write(recs[i],out,'fasta')
		except KeyError:
			print(i)
			pass
	out.close()