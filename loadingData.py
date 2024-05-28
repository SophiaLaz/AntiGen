import pandas as pd
import numpy as np
from tqdm import tqdm
import os


def length_region(dfj: pd.DataFrame, chain: str, name: str):
    dfj[name + '_length_' + chain] = dfj[name + '_aa_' + chain].apply(lambda x: len(x))
    return dfj


def concatenate_regions(dfj: pd.DataFrame, chain: str):
    dfj['sequence_aa_' + chain] = dfj['fwr1_aa_' + chain] + dfj['cdr1_aa_' + chain] + \
                                dfj['fwr2_aa_' + chain] + dfj['cdr2_aa_' + chain] + \
                                dfj['fwr3_aa_' + chain] + dfj['cdr3_aa_' + chain] + dfj['fwr4_aa_' + chain]
    return dfj


parth = 'dataset/'
filenames = os.listdir(parth)
print(f'Количество файлов в каталоге: {len(filenames)}')
df = pd.DataFrame()
countChains = 0
countShortedFw, countShortedCdr = np.array([[0, 0], [0, 0], [0, 0], [0, 0]]), np.array([[0, 0], [0, 0], [0, 0]])
columns = ['v_call_heavy', 'Isotype_heavy',
           'fwr1_aa_heavy', 'cdr1_aa_heavy', 'fwr2_aa_heavy', 'cdr2_aa_heavy',
           'fwr3_aa_heavy', 'cdr3_aa_heavy', 'fwr4_aa_heavy',
           'v_identity_heavy', 'ANARCI_status_heavy',
           'v_call_light', 'Isotype_light',
           'fwr1_aa_light', 'cdr1_aa_light', 'fwr2_aa_light', 'cdr2_aa_light',
           'fwr3_aa_light', 'cdr3_aa_light', 'fwr4_aa_light',
           'v_identity_light', 'ANARCI_status_light']

for file_i in tqdm(filenames[:1]):
    try:
        dfi = pd.read_csv(parth+file_i, skiprows=1, usecols=columns).drop_duplicates()
        countChains += dfi.shape[0]
        for chain_i in ['heavy', 'light']:
            dfi = concatenate_regions(dfi, chain_i)
            for name_i in ['fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4']:
                dfi = length_region(dfi, chain_i, name_i)
        dfi['ANARCI_status_heavy'] = dfi['ANARCI_status_heavy'].apply(lambda x: x[1:-1].split(sep='|'))
        dfi['ANARCI_status_light'] = dfi['ANARCI_status_light'].apply(lambda x: x[1:-1].split(sep='|'))
        maskCh = dfi['ANARCI_status_heavy'].apply(lambda x: x[2] == '')
        maskf1h = dfi['ANARCI_status_heavy'].apply(lambda x: 'fw1' in x[3])
        maskf2h = dfi['ANARCI_status_heavy'].apply(lambda x: 'fw2' in x[3])
        maskf3h = dfi['ANARCI_status_heavy'].apply(lambda x: 'fw3' in x[3])
        maskf4h = dfi['ANARCI_status_heavy'].apply(lambda x: 'fw4' in x[3])
        maskc1h = dfi['ANARCI_status_heavy'].apply(lambda x: 'cdr1' in x[3])
        maskc2h = dfi['ANARCI_status_heavy'].apply(lambda x: 'cdr2' in x[3])
        maskc3h = dfi['ANARCI_status_heavy'].apply(lambda x: 'cdr3' in x[3])
        maskCl = dfi['ANARCI_status_light'].apply(lambda x: x[2] == '')
        maskf1l = dfi['ANARCI_status_light'].apply(lambda x: 'fw1' in x[3])
        maskf2l = dfi['ANARCI_status_light'].apply(lambda x: 'fw2' in x[3])
        maskf3l = dfi['ANARCI_status_light'].apply(lambda x: 'fw3' in x[3])
        maskf4l = dfi['ANARCI_status_light'].apply(lambda x: 'fw4' in x[3])
        maskc1l = dfi['ANARCI_status_light'].apply(lambda x: 'cdr1' in x[3])
        maskc2l = dfi['ANARCI_status_light'].apply(lambda x: 'cdr2' in x[3])
        maskc3l = dfi['ANARCI_status_light'].apply(lambda x: 'cdr3' in x[3])
        dfi = dfi[maskCh & maskCl]
        df = pd.concat([df, dfi], ignore_index=True)
        countShortedFw += [[sum(maskf1h), sum(maskf1l)], [sum(maskf2h), sum(maskf2l)], [sum(maskf3h), sum(maskf3l)], [sum(maskf4h), sum(maskf4l)]]
        countShortedCdr += [[sum(maskc1h), sum(maskc1l)], [sum(maskc2h), sum(maskc2l)], [sum(maskc3h), sum(maskc3l)]]

        print(file_i)
        print()
    except Exception as except_i:
        print(f'File {file_i} is not correct: {except_i}')

print(f'\ncountShortedFw = {countShortedFw}, \ncountShortedCdr = {countShortedCdr}')
print(f'Количество пар цепей изначально {countChains} в датафрейме {df.shape[0]}, размер датафрейма: {df.shape}')
df.to_csv('files/all_sequences.csv', columns=df.columns, index=False)
