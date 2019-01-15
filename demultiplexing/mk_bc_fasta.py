# script to convert reference fasta and barcode list to barcode fasta for alignment
# Calla Martyn 01-14-18

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

# read in barcode list
gen_barcodes = pd.read_csv('barcode_list.csv', header = None)

# read in GBC reference
GBC = SeqIO.read('GBC_ref.fa', 'fasta')
# locate region of reference that has Ns (to replace with barcodes)
bc_start = GBC.seq.find('N')
bc_end = bc_start + 18

# names of barcodes are in first column of dataframe
bc_names = gen_barcodes[1]


records = []
for i in range(len(bc_names)):
    # creating a sequence with each barcodes reverse complement substituted for the Ns
    sequence = GBC.seq[0:bc_start] + Seq(gen_barcodes[0].iloc[i], generic_dna).reverse_complement() + GBC.seq[bc_end:]
    # saving sequence and name as a SeqRecord object
    records.append(SeqRecord(sequence, str(bc_names[i])))

# writing output to a fasta file
SeqIO.write(records, 'gen_barcodes.fasta', 'fasta')
