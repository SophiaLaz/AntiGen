from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import pandas as pd

input_path = 'files/'
output_path = 'clustering/'
input_file = 'all_sequences.csv'
output_file = 'clustering/DbForClustering.fasta'
df = pd.read_csv(input_path + input_file)

all_seq = []
for ind in tqdm(range(df.shape[0])):
    seq_id = str(ind) + '_' + df['v_call_heavy'][ind]
    seq_desc = df['Isotype_heavy'][ind]
    seq = Seq(df['sequence_aa_heavy'][ind])
    seq_r = SeqRecord(seq, id=seq_id, description=seq_desc)
    all_seq.append(seq_r)

SeqIO.write(all_seq, output_path + output_file, 'fasta')
