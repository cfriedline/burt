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
z12_df['/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/E'].head()

# %%
for k in z12_df:
    v = z12_df[k]
    v['population'] = v.apply(lambda x: "-".join(x.name.split("-")[0:-1]), axis=1)
    v = v.replace(-1, 9)
    z12_df[k] = v

# %%

# %%
snpmat = {}
for k in z12_df:
    snpmat[k] = z12_df[k].ix[:,:-1]

# %%
locus_names = {}
for k in z12_df:
    locus_names[k] = z12_df[k].columns[:-1]

# %%
pop_names = {}
for k in z12_df:
    pop_names[k] = z12_df[k]['population']

# %%
for k in snpmat:
    print(k)
    s = snpmat[k]
    l = pd.Series(locus_names[k])
    p = pd.Series(pop_names[k])
    s.to_csv(os.path.join(k, "snpmat"), sep="\t", header=False, index=False)
    l.to_csv(os.path.join(k, "locus_names"), sep="\t", header=False, index=False)
    p.to_csv(os.path.join(k, "pop_names"), sep="\t", header=False, index=False)
    

# %% {"language": "R"}
# library(OutFLANK)
# library(data.table)

# %%
len(pop_names[k].unique())

# %%
for k in snpmat:
    scr = """library(OutFLANK)
        library(data.table)
        pop_names = unlist(fread("pop_names", data.table=F, header=F), use.names=F)
        locus_names = unlist(fread("locus_names", data.table=F, header=F), use.names=F)
        snpmat = fread("snpmat", header=F, data.table=F)
        fstdata = MakeDiploidFSTMat(snpmat,locus_names,pop_names)
        outflank_res = OutFLANK(fstdata, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples, qthreshold=0.05)
        saveRDS(outflank_res, "outflank_res.rds")"""
        
    with open(os.path.join(k, "outflank.R"), "w") as o:
        o.write("setwd('{}')\n".format(k))
        o.write("NumberOfSamples={}\n".format(len(pop_names[k].unique())))
        o.write(scr)
    

# %%
outflank_res_files = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3 -name "outflank*.rds"
for o in outflank_res_files:
    var = "outflank_res_{}".format(os.path.basename(os.path.dirname(o)))
    r("{} = readRDS('{}')".format(var, o))

# %%
outflank_res = {}
for s in spp:
    var = "outflank_res_{}$results".format(s)
    outflank_res[s] = r(var)

# %%

# %%
outliers = {}
for k in outflank_res:
    d = outflank_res[k]
    o = pd.DataFrame(d[d.OutlierFlag == 1])
    if len(o) > 0:
        outliers[k] = o

# %%
sets = []
for k in outliers:
    sets.append(set(outliers[k].LocusName.tolist()))

# %%
set.intersection(*sets)

# %%
sets

# %% {"language": "R"}
# library(OutFLANK)
# library(data.table)

# %% {"language": "R"}
# setwd("/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3")
# files=c("./E/outflank_res.rds",
# "./P/outflank_res.rds",
# "./G/outflank_res.rds",
# "./T/outflank_res.rds")
# names = c("Slash", "Longleaf", "TMP", "Loblolly")
# for (i in 1:length(files)) {
#     outflank_res = readRDS(files[i])
#     print(paste(files[i], outflank_res$numberHighFstOutliers))
#     OutFLANKResultsPlotter(outflank_res,
#                        withOutliers=TRUE, 
#                        NoCorr=TRUE, 
#                        Hmin=0.1, 
#                        binwidth=0.005, 
#                        Zoom=FALSE,
#                        RightZoomFraction=0.05, 
#                        titletext=names[i])
# }
#

# %%
