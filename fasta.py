from Bio import SeqIO
seq_dict = {rec.id : rec.seq for rec in SeqIO.parse("data.txt", "fasta")}

sqs = []

for i in seq_dict:
    sqs.append(str(seq_dict[i]))
