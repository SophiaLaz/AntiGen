from typing import List, Tuple
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO


cluster_path = 'clustering/98_cluster/'
dataset_path = 'files/'

filename = 'resDB_clu.tsv'
count_filename = 'resDB_clu_count.csv'
new_filename = 'resDB_seq_cluster.csv'
dataset_filename = 'all_sequences.csv'
lightchain_filename = 'count_vgene.csv'

ref_heavy_file = 'resDB_rep.fasta'
final_dataset_name = 'dataset.csv'


def load_clustering_result() -> dict:
    df = pd.read_csv(cluster_path + filename, sep='\t', header=None)
    res_name, res_count = {}, {}

    for i in range(df.shape[0]):
        cluster, name = df[0][i], df[1][i]
        if cluster in res_name:
            res_name[cluster].append(name)
            res_count[cluster] += 1
        else:
            res_name[cluster] = [name]
            res_count[cluster] = 1

    # сохраним в отдельные файлы
    dump_file(sorted(res_name.items(), key=lambda x: len(x[1]), reverse=True),
              sorted(res_count.items(), key=lambda x: x[1], reverse=True))
    return res_name


def dump_file(res_name: List[Tuple[str, List[str]]], res_count: List[Tuple[str, int]]) -> None:
    df = pd.DataFrame(data={'cluster': [x[0] for x in res_count], 'count': [x[1] for x in res_count]})
    df.to_csv(cluster_path + count_filename, columns=df.columns, index=False)
    df2 = pd.DataFrame(data={'cluster': [x[0] for x in res_name], 'members': ["|".join(x[1]) for x in res_name]})
    df2.to_csv(cluster_path + new_filename)


def analyze_clusters() -> None:
    try:
        df = pd.read_csv(cluster_path + count_filename)
    except Exception as err:
        print(f"Can't load the {count_filename} file! {err}")
        return
    size = list(df['count'])
    print(f'Максимальные размеры кластеров: {size[:5]}')
    print('Количество кластеров разного размера:')
    print(f'\t>100:\t\t\t{len([i for i in size if i > 100])}')
    print(f'\t50...100:\t{len([i for i in size if 50 < i <= 100])}')
    print(f'\t20...50:\t{len([i for i in size if 20 < i <= 50])}')
    print(f'\t10...20:\t{len([i for i in size if 10 < i <= 20])}')
    print(f'\t2..10:\t\t{len([i for i in size if 2 <= i <= 10])}')
    print(f'\t1:\t\t\t\t{len([i for i in size if i == 1])}')
    analyze_current_cluster(df['cluster'][0])


def analyze_current_cluster(curr_cluster: str):
    try:
        df = pd.read_csv(cluster_path + new_filename)
    except Exception as err:
        print(f"Can't load the {new_filename} file! {err}")
        return
    members = list(df['members'][df['cluster'] == curr_cluster])[0]
    members = members.split("|")
    ind = [int(member.split('_')[0]) for member in members]
    seqs = load_current_seq(ind)
    length = np.array([len(seq) for seq in seqs])
    uniq_len = np.unique(length, return_counts=True)

    print(f'\n{"-" * 50}')
    print(f'Проанализируем кластер {curr_cluster} размера {len(members)}:')
    for len_i, count_i in zip(uniq_len[0], uniq_len[1]):
        print(f"\tпоследовательностей длины {len_i} - {count_i} шт.")


def analyze_light_chain_in_cluster():
    try:
        df_members = pd.read_csv(cluster_path + new_filename)
    except Exception as err:
        print(f"Can't load the {new_filename} file! {err}")
        return
    try:
        df_count = pd.read_csv(cluster_path + count_filename)
    except Exception as err:
        print(f"Can't load the {count_filename} file! {err}")
        return
    try:
        df = pd.read_csv(dataset_path + dataset_filename)
    except Exception as err:
        print(f"Can't load the {dataset_filename} file! {err}")
        return

    lc_dict = {}
    res = []
    for i in tqdm(range(len(df_members))):
        if df_count['count'][i] == 1:
            continue
        if df_members['cluster'][i] != df_count['cluster'][i]:
            print("Incorrect dataset files...")
        members = df_members['members'][i].split("|")
        ind = [int(member.split('_')[0]) for member in members]
        v_call_light = np.array([df['v_call_light'][i] for i in ind])
        # v_call_light = list(np.unique([v.split("*")[0] for v in v_call_light]))
        v_call_light = list(np.unique(v_call_light))

        for vi in v_call_light:
            if vi not in lc_dict:
                lc_dict[vi] = 0
            lc_dict[vi] += 1

        v_call_light_count = len(v_call_light)
        size = df_count['count'][i]
        res.append((df_members['cluster'][i], size, v_call_light_count, "|".join(v_call_light)))

    df2 = pd.DataFrame(data={
        'cluster_name': [x[0] for x in res],
        'cluster_size': [x[1] for x in res],
        'v_call_light_count': [x[2] for x in res],
        'v_call_light': [x[3] for x in res]})
    df2.to_csv(cluster_path + lightchain_filename)

    lc_dict = lc_dict.items()
    df3 = pd.DataFrame(data={
        'v_call_light': [x[0] for x in lc_dict],
        'count': [x[1] for x in lc_dict]
    })
    df3.to_csv(cluster_path + 'frequency_' + lightchain_filename)


def load_current_seq(ind: List[int]) -> List[str]:
    try:
        df = pd.read_csv(dataset_path + dataset_filename)
    except Exception as err:
        print(f"Can't load the {dataset_filename} file! {err}")
        return []

    return [df['sequence_aa_heavy'][i] for i in ind]


def create_final_dataset():
    try:
        df_members = pd.read_csv(cluster_path + new_filename)
    except Exception as err:
        print(f"Can't load the {new_filename} file! {err}")
        return
    try:
        df_count = pd.read_csv(cluster_path + count_filename)
    except Exception as err:
        print(f"Can't load the {count_filename} file! {err}")
        return
    try:
        df = pd.read_csv(dataset_path + dataset_filename)
    except Exception as err:
        print(f"Can't load the {dataset_filename} file! {err}")
        return

    def numbering_light_chain(light_seqs_ind: List[int]) -> List[str]:
        regions = []
        for light_seq_ind in light_seqs_ind:
            count = 0
            regions_i = ["0"]
            for name_col in ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]:
                count += df[name_col + "_length_light"][light_seq_ind]
                regions_i.append(str(count))
            regions.append("|".join(regions_i))
        return regions

    # Загрузка репрезентативных последовательностей из fasta-файла:
    ref_heavy = {}
    with open(cluster_path + ref_heavy_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ref_heavy[record.id] = record.seq

    res = []
    for i in tqdm(range(len(df_members))):
        if df_count['count'][i] == 1:
            continue
        if df_members['cluster'][i] != df_count['cluster'][i]:
            print("Incorrect dataset file...")
        members = df_members['members'][i].split("|")
        ind = [int(member.split('_')[0]) for member in members]
        v_call_light = np.array([df['v_call_light'][i] for i in ind])
        # v_call_light = list([v.split("*")[0] for v in v_call_light])
        light_seqs = np.array([df['sequence_aa_light'][i] for i in ind])
        regions = numbering_light_chain(ind)

        v_call_unique = list(np.unique(v_call_light))
        v_call_light_count = len(v_call_unique)
        size = df_count['count'][i]

        if v_call_light_count > 1:
            for j in range(size):
                for k in range(size):
                    if j != k:
                        res.append((df_members['cluster'][i], size, ref_heavy[df_members['cluster'][i]],
                                    light_seqs[j], regions[j],
                                    v_call_light[k], light_seqs[k], regions[k]))

    df2 = pd.DataFrame(data={
        'cluster_name': [x[0] for x in res],
        'cluster_size': [x[1] for x in res],
        'seq_ref_heavy': [x[2] for x in res],
        'seq_input_light': [x[3] for x in res],
        'seq_input_numb': [x[4] for x in res],
        'v_call_light': [x[5] for x in res],
        'seq_output_light': [x[6] for x in res],
        'seq_output_numb': [x[7] for x in res]})
    df2.to_csv(dataset_path + final_dataset_name)


if __name__ == '__main__':
    # cluster_dict = load_clustering_result()
    # analyze_clusters()
    # analyze_light_chain_in_cluster()
    create_final_dataset()
