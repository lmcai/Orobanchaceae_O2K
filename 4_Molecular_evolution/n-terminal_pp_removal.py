from Bio import SeqIO
files=open('/Users/limingcai/Downloads/s.txt').readlines()
#files=open('n-pep_position.txt').readlines()

def find_nth_non_gap_index(sequence, n):
    count = 0
    for i, char in enumerate(sequence):
        if char != '-':
            count += 1
            if count == n:
                return i
    return len(sequence)  # If fewer than n non-gap characters


def process_alignment(fasta_file, cutoff_site):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    ath_seq = next((r for r in sequences if r.id.startswith('ArabidopsisThaliana')), None)
    if cutoff_site is None:
        return sequences  # Return unmodified sequences
    cutoff_index = find_nth_non_gap_index(ath_seq, cutoff_site)
    for rec in sequences:
        rec.seq = rec.seq[cutoff_index + 1:]
    return sequences	

for l in files:
	og_id, value = l.strip().split('\t')
	cutoff_site = None if value == 'NA' else int(value) * 3
	trimmed = process_alignment(og_id+'.one2one.fas', cutoff_site)
	output_file = f"{og_id}.trimmed.fas"
	SeqIO.write(trimmed, output_file, "fasta")