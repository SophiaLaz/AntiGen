from typing import List, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


cluster_path = 'clustering/'
dataset_path = 'files/'
filename = 'clusterRes_cluster.tsv'
count_filename = 'clusterRes_count.csv'
new_filename = 'clusterRes_seq_cluster.csv'
dataset_filename = 'all_sequences.csv'


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


def load_current_seq(ind: List[int]) -> List[str]:
    try:
        df = pd.read_csv(dataset_path + dataset_filename)
    except Exception as err:
        print(f"Can't load the {dataset_filename} file! {err}")
        return []

    return [df['sequence_aa_heavy'][i] for i in ind]


if __name__ == '__main__':
    # cluster_dict = load_clustering_result()
    analyze_clusters()
