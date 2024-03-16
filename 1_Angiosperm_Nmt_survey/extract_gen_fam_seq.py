#extract genome cds and pep sequences by entry id from poff result
#python extract_gen_fam_seq.py [orthogroup file] [genome pep] [genome cds] [index number] [family to include list] [additional taxa prefix]

from Bio import SeqIO
import sys
from Bio.SeqRecord import SeqRecord

x=open(sys.argv[1])
y=x.read().splitlines()
#ll=open(sys.argv[5]).readlines()

#ll=[l.strip() for l in ll]

aa_dict = SeqIO.index(sys.argv[2], "fasta")
na_dict = SeqIO.index(sys.argv[3], "fasta")
seq_count = 0
for line in y[1:]:
        z=line.split('\t')
        if True:
	#if z[0] in ll:
		id=z[int(sys.argv[4])].split(", ")
		aa_filename= z[0]+'.aa.fas'
		na_filename= z[0]+'.na.fas'
		if id[0]!='':
			try:
				for j in id:
					SeqIO.write(SeqRecord(na_dict[j].seq,id=sys.argv[6]+'_'+j), open(na_filename,'a'), "fasta")
					SeqIO.write(SeqRecord(aa_dict[j].seq,id=sys.argv[6]+'_'+j), open(aa_filename,'a'), "fasta")
			except KeyError as e:
				print e
			
x.close()

