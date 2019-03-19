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
import os, sys
import rpy2
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import rpy2.robjects as ro
import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib inline
import seaborn as sns
import numpy as np
import dill
import random
import vcf
import statsmodels.api as sm
import statsmodels.formula.api as smf
import operator
import traceback
# %load_ext rpy2.ipython
from rpy2.robjects import pandas2ri as p2r
p2r.activate()
r = ro.r
import shutil
from utils import read_df, save_df
from pathlib import Path, PurePath
from ipyparallel import Client
from collections import Counter, defaultdict, namedtuple, OrderedDict
from scipy.stats import mannwhitneyu, ks_2samp, f_oneway
import tables
import ujson
import pickle
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from IPython.display import display

# %%
analysis_dir = "/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3"

# %%
spp = ["E", "G", "P", "T"]

# %%
z12_df = {}
for s in spp:
    d = os.path.join(analysis_dir, s)
    print(d)
    if not "beagle" in d:
        z12_df[d] = read_df(d, "z12_swapped")

# %%
test_key = '/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/E'

# %%
z12_df[test_key].head()

# %%
for k in z12_df:
    v = z12_df[k]
    v['population'] = v.apply(lambda x: "-".join(x.name.split("-")[0:-1]), axis=1)
    v = v.replace(-1, 9)
    z12_df[k] = v

# %%
pd.read_csv("pheno.txt", sep='\t', index_col=0).join(z12_df[test_key], how="inner").head(30)

# %%

# %%
