from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import pandas as pd


file_path = 'files/'
cluster_path = 'clustering/98_cluster/'
dataset_file = 'all_sequences.csv'
new_dataset_file = 'clustering/DbForClustering.fasta'
cluster_file = 'resDB_seq_cluster.csv'
count_cluster = 'resDB_clu_count.csv'
small_dataset_file = 'DbForClusteringSmall.fasta'


def create_db():
    df = pd.read_csv(file_path + dataset_file)
    all_seq = []
    for ind in tqdm(range(df.shape[0])):
        seq_id = str(ind) + '_' + df['v_call_heavy'][ind]
        seq_desc = df['Isotype_heavy'][ind]
        seq = Seq(df['sequence_aa_heavy'][ind])
        seq_r = SeqRecord(seq, id=seq_id, description=seq_desc)
        all_seq.append(seq_r)

    SeqIO.write(all_seq, cluster_path + new_dataset_file, 'fasta')


def create_cluster_fasta_files():
    df = pd.read_csv(file_path + dataset_file)

    cluster_df = pd.read_csv(cluster_path + cluster_file)
    cluster_names = cluster_df['cluster'].tolist()
    members = cluster_df['members'].tolist()

    for i in tqdm(range(len(cluster_names))):
        all_seqs = []
        seqs = members[i].split('|')
        if len(seqs) == 1:
            continue
        if len(seqs) <= 10:
            cur_dir = '10/'
        elif len(seqs) <= 20:
            cur_dir = '20/'
        elif len(seqs) <= 50:
            cur_dir = '50/'
        else:
            cur_dir = '100/'

        indexes = [int(seq.split('_')[0]) for seq in seqs]

        for j, ind in enumerate(indexes):
            seq_id = seqs[j]
            seq_desc = df['Isotype_heavy'][ind]
            seq = Seq(df['sequence_aa_heavy'][ind])
            seq_r = SeqRecord(seq, id=seq_id, description=seq_desc)
            all_seqs.append(seq_r)
        SeqIO.write(all_seqs, cluster_path + cur_dir + cluster_names[i] + '.fasta', 'fasta')


if __name__ == '__main__':
    # create_db()
    create_cluster_fasta_files()
