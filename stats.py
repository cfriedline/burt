# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
fastq_files = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe -name "*.fastq"

# %%
len(fastq_files)

# %%
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from ipyparallel import Client

# %%
rc = Client(profile="sge")

# %%
dv =rc[:]
lv = rc.load_balanced_view()
len(dv)

# %%
with dv.sync_imports():
    from Bio.SeqIO.QualityIO import FastqGeneralIterator


# %%
def count_records(f):
    recs = 0
    for name, seq, qual in FastqGeneralIterator(open(f)):
        recs += 1
    return f, recs


# %%
dv['count_records'] = count_records

# %%
jobs = lv.map_async(count_records, fastq_files)

# %%
jobs.progress

# %%
for j in jobs:
    print(j)

# %%
import pandas as pd
import numpy as np
import os

# %%
counts = pd.DataFrame([x for x in jobs], columns = ['file', 'reads'])

# %%
counts["sample_name"] = counts.file.apply(lambda x: os.path.basename(x).replace(".R1.fastq", ""))
counts['population'] = counts.sample_name.apply(lambda x: "-".join(x.split("-")[0:-1]))
counts['species'] = counts.sample_name.apply(lambda x: (x.split("-")[0]))

# %%
fmt = lambda x: '%.2f' % x
fmt2 = lambda x: int(np.round(x))
counts[['reads', 'species']].groupby("species").describe().applymap(fmt2)

# %%
counts.to_csv("counts.txt", sep="\t", header=True, index=False)

# %%
counts = pd.read_csv("counts.txt", sep="\t")

# %%
counts.sort_values("reads")

# %%
counts.groupby("species")['population'].describe()

# %%
counts

# %%
