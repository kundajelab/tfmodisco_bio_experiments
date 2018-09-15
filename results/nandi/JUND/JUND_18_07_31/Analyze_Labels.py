
# read in header (list of experiments)
import gzip

with gzip.open('labels.txt.gz', 'r') as f:
    head = f.readline()
    head = head.split()[1:]

print(len(head))

# read in validation set labels to df (one column per experiment)
import pandas as pd

df = pd.read_csv("splits/valid.txt.gz", sep='\t', index_col=0, header=None, compression='gzip')
df.columns = head
print(df.shape)

# generate label counts for amb, neg, and pos [-1,0,1] for each column
counts_list = []
for c in range(df.shape[1]):
    counts = df.iloc[:, c].value_counts().sort_index() # [-1,0,1]
    counts_list.append(counts)
counts_df = pd.concat(counts_list, axis=1)
print(counts_df.index)
print(counts_df)

# get percentage of amb, neg, pos
sums = counts_df.sum()
print(sums)
percent_df = counts_df / sums
print(percent_df)

all = pd.concat([counts_df, percent_df])
all.transpose().to_csv("counts.tsv", sep='\t')
all.columns = ["amb","neg","pos","%amb","%neg","%pos"]
