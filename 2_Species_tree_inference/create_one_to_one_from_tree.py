from ete3 import Tree
import os

files=os.listdir('.')
files=[i for i in files if i.endswith('.treefile')]

for f in files:
	t=Tree(f)
	#reroot
	amb=[node.name for node in t if node.name.startswith('Amb')]
	if len(amb):
		amb_node=t&amb[0]
		t.set_outgroup(amb_node)
	else:
		cin=[node.name for node in t if node.name.startswith('Cin')]
		if len(cin):
			cin_node=t&cin[0]
			t.set_outgroup(cin_root)
	to_keep=[]
	shortest={}
	for node in t:
		sp=node.name.split('_')[0]
		if not sp in to_keep:
			to_keep.append(node.name)
			shortest[sp]=node.get_distance(amb_node)
		else:
			if node.get_distance(amb_node)<shortest[sp]:
				previous=[i for i in to_keep if i.startswith(sp)]
				to_keep.remove(previous[0])
				to_keep.append(node.name)
				shortest[sp]=node.get_distance(amb_node)
			elif node.get_distance(amb_node)==shortest[sp]:
				if not node.name.endswith('_2'):
					previous=[i for i in to_keep if i.startswith(sp)]
					if len(node.name) >len(previous):
						to_keep.remove(previous[0])
						to_keep.append(node.name)
						shortest[sp]=node.get_distance(amb_node)
	recs=SeqIO.index(f.split('.')[0]+'.filter.fas','fasta')
	out=open(f.split('.')[0]+'.one2one.fas','w')
	for i in to_keep:
		try:
			d=SeqIO.write(recs[i],out,'fasta')
		except KeyError:
			print(i)
			pass
	out.close()