import pandas as pd
import matplotlib.pyplot as plt


parth = 'files/'
file_name = 'all_sequences.csv'
df = pd.read_csv(parth + file_name)

columns = []
for chain_i in ['heavy', 'light']:
    for name_i in ['fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3']:
        columns.append(name_i + '_length_' + chain_i)

for i in range(2):
    df[columns[i * 6: (i + 1) * 6]].hist(bins=25)
plt.show()
