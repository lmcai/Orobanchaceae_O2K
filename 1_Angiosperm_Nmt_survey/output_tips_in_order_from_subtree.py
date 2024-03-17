from ete3 import Tree
import sys
t=Tree(sys.argv[1])
tips=[node.name for node in t]
outline=''

for i in ['Pht','Mim','Han', 'Ses', 'Ole', 'Cus', 'Ipo', 'Nic', 'Sol', 'Cof', 'Dau', 'Hel', 'Act', 'Vac', 'Mes', 'Ptr', 'Ath', 'Gos', 'Tca', 'Cit', 'Gly', 'Pru', 'Vvi', 'Aqu', 'Osa', 'Sor', 'Cin', 'Amb']:
	cur=''
	for j in tips:
		if j.startswith(i):cur=cur+', '+j[4:]
	if cur=='':
		outline=outline+'\t'
	else:
		outline=outline+'\t'+cur[2:]

print(outline+'\n')