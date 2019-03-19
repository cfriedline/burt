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
import sys

sys.path.append("../include_utils/")

from ipyparallel import Client
import ipyparallel as ipp
import os, time
import include_utils as u
import pandas as pd
import numpy as np
import scipy as sp
import numbers
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import vcf
from sklearn import preprocessing
from subprocess import Popen, PIPE
import seaborn as sns
from IPython.display import FileLink
#import urllib2
import urllib.request as urllib2
import urllib
import dill
import traceback
from pandas import Series, DataFrame
import gzip
import warnings
warnings.filterwarnings('ignore',category=pd.io.pytables.PerformanceWarning)
# %config InlineBackend.figure_format = 'retina'
from Bio import SeqIO
import pysam
from collections import OrderedDict, namedtuple
import operator
import multiprocessing as mp
import pickle
from IPython.display import FileLink, FileLinks, display

samtools = "/home/cfriedline/gpfs/src/samtools-1.3/samtools"
bcftools = "/home/cfriedline/gpfs/src/bcftools-1.3/bcftools"
picard = "/home/cfriedline/gpfs/src/broadinstitute-picard-03a1d72/dist/picard.jar"
java = "/home/cfriedline/g/src/jdk1.8.0_60/bin/java"
perl = "/home/cfriedline/gpfs/opt/ActivePerl-5.18/bin/perl"

vcfutils = "perl /home/cfriedline/g/src/bcftools-1.3/vcfutils.pl"
vcftools = "/home/cfriedline/bin/vcftools"
bcftools = "/home/cfriedline/gpfs/src/bcftools-1.3/bcftools"
tabix = "/home/cfriedline/gpfs/src/htslib-1.3/tabix"
bgzip = "/home/cfriedline/gpfs/src/htslib-1.3/bgzip"


def setup_r():
    os.environ['R_HOME'] = '/home/cfriedline/g/R3/lib64/R'
    os.environ['LD_LIBRARY_PATH'] = "%s/lib:%s:%s" % (os.environ['R_HOME'], 
                                                   os.environ['LD_LIBRARY_PATH'],
                                                     "/home/cfriedline/lib64")

setup_r()
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
r = robjects.r

# %reload_ext autoreload
# %autoreload 2
# %matplotlib inline
# %reload_ext rpy2.ipython

# %%
species = ["T", "E", "P", "G"]
ni_dir = "/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3"
imp_dir = "/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40" 

# %%
vcfs = []
for spp in species:
    for d in [ni_dir, imp_dir]:   
        d = os.path.join(d, spp)
        vcfs.append(os.path.join(d, "isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz"))

# %%
vcfs

# %%
for v in vcfs:
    !$vcftools --gzvcf $v --012 --out $v

# %%
z12s = ["%s.012" % x for x in vcfs]

# %%
z12s

# %%
rc = Client(profile="sge")


# %%
def get_z12_df(z12_file):
    import numpy as np
    import pandas as pd
    print(z12_file)
    indv_file = "%s.indv" % z12_file
    pos_file = "%s.pos" % z12_file
    z12_data = []
    for i, line in enumerate(open(z12_file)):
        line = line.strip()
        line = [int(x) for x in line.split("\t")]
        z12_data.append(np.array(line))
    z12_data = np.array(z12_data)
    p = pd.read_csv(pos_file, sep="\t", names=['contig', 'pos'])
    i = pd.read_csv(indv_file, names=['sample_name'])
    df = pd.DataFrame(z12_data)
    df = df.drop(0, axis=1)
    df.columns = p.apply(lambda x: "%s_%s" % (x.contig, x.pos), axis=1)
    df.index = [x for x in i.sample_name]
    return z12_file, df
#z12_dfs = [get_z12_df(x) for x in z12s]
#z12_dfs = [x[keep_snps.index] for x in z12_dfs]


# %%
dv, lv = u.get_views(rc)
len(dv)

# %%
dv['get_z12_df'] = get_z12_df

# %%
z12_jobs = lv.map_async(get_z12_df, z12s)

# %%
z12_jobs.completed

# %%
z12_dfs = {}
for j in z12_jobs:
    z12_dfs[j[0]] = j[1]

# %%
test_key = '/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/E/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz.012'

# %%
z12_dfs[test_key]


# %%
def get_pheno():
    pheno = pd.read_csv("/home/cfriedline/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/pheno/evolution2016_sample_info.txt", sep="\t")
    pheno.index = pheno['name']
    pheno = pheno.drop("name", axis=1)
    return pheno
pheno = get_pheno()
pheno['state'] = pheno.population.apply(lambda x: x.split("-")[1])

# %%
[(k, v.shape) for k, v in z12_dfs.items()]


# %%
def get_correction(n):
    #for finite sample size
    return (2*n)/(2*n-1)

def get_allele_freqs(locus, debug):
    c = locus[locus != -1].value_counts()
    total_alleles = 2.0*sum(c)
    num_individuals = sum(c)
    P = 0
    Q = 0
    PQ = 0
    if 0 in c:
        P = 2*c[0]
    if 2 in c:
        Q = 2*c[2]
    if 1 in c:
        PQ = c[1]
    P += PQ
    Q += PQ
    if total_alleles == 0:
        return None
    p = P/total_alleles
    q = Q/total_alleles
    assert p + q == 1.0
    He = 2 * p * q * get_correction(num_individuals)
    Ho = PQ*1.0/num_individuals
    Fis = 1 - (Ho/He)
    #print p, q, He, Ho, Fis
    
        
    ret = pd.Series({"p":p, 
                      "q":q,
                      "P":P,
                      "Q":Q,
                      "He":He,
                      "Ho":Ho, 
                      "Fis":Fis})
    if debug:
        print(ret)
    return ret

# %%
ni_stats = pickle.load(open(os.path.join(ni_dir, "combined_stats.pkl"), "rb"))
imp_stats = pickle.load(open(os.path.join(imp_dir, "combined_stats.pkl"), "rb"))


# %%
def get_pos(key):
    d = pd.read_csv("{}.pos".format(key), header=None, sep="\t", names=["ctg", "position"])
    d['snp_name'] = d.apply(lambda x: "{}-{}".format(x.ctg, x.position), axis=1)
    return d

def get_stat(key, col):
    stat = ni_stats
    spp = os.path.basename(os.path.dirname(key))
    if "beagle40" in key:
        stat = imp_stats
    snp_pos = get_pos(key)
    return stat[spp][0][col]


# %%
ni_keys = [x for x in z12_dfs if not 'beagle40' in x]
for key in ni_keys:
    key2 = key.replace("samtools1.3", "samtools1.3/beagle40")
    mafs0 = get_stat(key, "MAF")
    mafs1 = get_stat(key2, "MAF")
    j = pd.concat([mafs0, mafs1], join="inner", axis=1)
    j.columns = ["maf_ni", "maf_imp"]
    plt.scatter(j['maf_ni'], j['maf_imp'])
    plt.title("{} MAF".format(os.path.basename(os.path.dirname(key))))
    plt.xlabel("not imputed")
    plt.ylabel("imputed")
    plt.show()


# %%
def convert_locus(name):
    return "-".join(name.rsplit("_", 1))

def swap_alleles(locus, af):
    locus_id = convert_locus(locus.name)
    freqs = af.ix[locus_id]
    if freqs['MAF'] == freqs["A1_freq"]:
        return locus.replace({0:2,2:0})
    return locus


# %%
z12_swapped = {}
for k, v in z12_dfs.items():
    print(k)
    allele_freqs = get_stat(k, ['A1_freq', "MAF"])
    z12_swapped[k] = v.apply(swap_alleles, args=(allele_freqs,))

# %%
z12_swapped[test_key].ix[:,0:5].head()

# %%
z12_dfs[test_key].ix[:,0:5].head()

# %%
pop_id = {}
i = 1
for p in sorted(pheno.population.unique()):
    pop_id[p] = i
    i+=1

# %%
state_id = {}
i = 1
for p in sorted(pheno.state.unique()):
    state_id[p] = i
    i+=1

# %%
state_id


# %%
def assign_popid(series):
    p = series.name.rsplit("-", 1)[0]
    series['popid'] = pop_id[p]
    return series


# %%
z12_swapped_pop = {key: value.apply(assign_popid, axis=1) for (key, value) in z12_swapped.items()}


# %%
def center_and_standardize_value(val, u, var):
    if val == -1:
        return 0.0
    return (val-u)/np.sqrt(var)

def center_and_standardize(locus, af):
    if "_" in locus.name:
        locus_id = convert_locus(locus.name)
        maf = af.ix[locus_id]
        var = maf*(1-maf)
        u = np.mean([x for x in locus if x != -1])
        return locus.apply(center_and_standardize_value, args=(u, var))
    return locus


# %%
pca_std = {}
for k, df in z12_swapped_pop.items():
    print(k)
    allele_freqs = get_stat(k, "MAF")
    pca_std[k] = df.apply(center_and_standardize, args=(allele_freqs,))

# %%
len(pca_std)

# %%
pca_std[test_key].head()

# %%
pca_std_data = {key: value.ix[:,:-1] for (key, value) in pca_std.items()}

# %%
for k in pca_std_data:
    outdir = os.path.dirname(k)
    fname = os.path.join(outdir, "pca_std_data.txt")
    print(fname)
    pca_std_data[k].to_csv(fname, header=True, index=True, sep="\t")

# %%
pca_data_files = sorted(list(pca_std_data.keys()))
pca_data_files = [os.path.join(os.path.dirname(x), "pca_std_data.txt") for x in pca_data_files]

# %%
pca_data_files

# %%
# %R -i pca_data_files

# %% [markdown]
# ## Run PCA

# %% {"language": "R"}
# library(data.table)
# run_pca = function(data_file) {
#     print(data_file)
#     d = fread(data_file, sep="\t", data.table=F)
#     rownames(d) = d$V1
#     drops = c("V1")
#     d = d[,!(names(d) %in% drops)]
#     res = prcomp(d, scale=F, center=F)
#     rownames(res$x) = rownames(d)
#     fname = 'pca_res.rds'
#     out = file.path(dirname(data_file), fname)
#     print(out)
#     saveRDS(res, out)
#     return(out)
# }
# results = lapply(pca_data_files, run_pca)

# %%
pca_result_files = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3 -name "pca_res.rds"

# %%
pca_results = {}
for i, f in enumerate(pca_result_files):
    print(f)
    var = "res{}".format(i)
    r("res{}=readRDS('{}')".format(i, f))
    pca_results[f] = var

# %%
pca_results


# %%
def get_pca_x(res):
    x = pd.DataFrame(pandas2ri.ri2py(res.rx2("x")))
    x.index = res.rx2("x").names[0]
    x.columns = res.rx2("x").names[1]
    return x


# %%
summary = r('summary')

# %%
pca_x = {}
for key, var in pca_results.items():
    pca_x[key] = get_pca_x(r[var])


# %%
def imputed_name(key):
    if 'beagle' in key:
        return "imputed"
    return "not_imputed"

sns.set_style("white")
def plot_pca(key, pca_std, pca_std_data, pca_x, prcomp_res, pheno, color_by):
    pop_dict = pop_id
    pop_key = color_by
    if color_by == "state":
        pop_dict = state_id
    
    joined = pd.concat([pca_std_data, pca_x, pheno], join="inner", axis=1)
    legend_dict = {}
    legend_idx = 0
    for elem in sorted(joined[color_by].unique()):
        legend_dict[elem] = legend_idx
        legend_idx+=1
        
    norm = mcolors.Normalize(0, len(legend_dict))

    legend = {}
    
    for row in joined.iterrows():
        pop =row[1][pop_key]
        n = norm(legend_dict[pop])
        color=cm.rainbow(n)
        legend[pop] = color
        plt.scatter(row[1].PC1, 
                    row[1].PC2, 
                    s=50, 
                    c=color)
    fig = plt.gcf()
    ax = plt.gca()
    cmap = plt.get_cmap()
    fig.set_size_inches(10,8)
    
    imp = summary(prcomp_res).rx("importance")[0]
    plt.xlabel("PC1 (%g)" % imp.rx(2,1)[0])
    plt.ylabel("PC2 (%g)" % imp.rx(2,2)[0])

    handles = []
    for pop in sorted(legend):
        handles.append(mpatches.Patch(color=legend[pop], label=pop))
    plt.legend(handles=handles)
    
    out_file = "{}-{}-{}-{}.pdf".format(os.path.basename(os.path.dirname(key)),
                                        imputed_name(key),
                                        os.path.basename(key),
                                        color_by)
    
    plt.title("PCA of n={} {} samples on {} loci ({})".format(len(joined), 
                                                              joined.species[0].lower(), 
                                                              len(pca_std_data.columns), 
                                                              imputed_name(key)))
    sns.despine()
    plt.savefig(out_file)
    plt.show()
    return out_file


# %%
def get_pca_std(pca_results_key):
    for k in pca_std:
        if os.path.dirname(k) == os.path.dirname(pca_results_key):
            return pca_std[k]
    return None

def get_pca_std_data(pca_results_key):
    for k in pca_std_data:
        if os.path.dirname(k) == os.path.dirname(pca_results_key):
            return pca_std_data[k]
    return None

for key, var in pca_results.items():
    for color_by in ['state', 'population']:
        f = plot_pca(key, get_pca_std(key), get_pca_std_data(key), pca_x[key], r[var], pheno, color_by)
        display(FileLink(f))


# %%
def save_df(dirname, fname, df):
    f = os.path.join(dirname, "%s.txt" % fname) 
    df.to_csv(f, 
              header=True,
              index=True,
              sep="\t")
    print("saved %s" % f)


# %%
pca_std_data_files = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3 -name "pca_std_data.txt"

# %%
pca_std_data_files

# %%
# %R -i pca_std_data_files

# %%
pwd

# %% {"language": "R"}
# library(data.table)
# source("tw_calc.R")
# test=read.table("twtable", header=F)
#
# run_tw = function(pca_data_file) {
#     d = fread(pca_data_file, sep="\t", data.table=F)
#     rownames(d) = d$V1
#     drops = c("V1")
#     d = d[,!(names(d) %in% drops)]
#     return(TWcalc(as.matrix(d),20))
# }
# tw_results = lapply(pca_std_data_files, run_tw)

# %%
tws = {}
for i, f in enumerate(pca_data_files):
    ri = i + 1
    tws[f] = r('tw_results[[{}]][[2]]'.format(ri))


# %% {"language": "R"}
# ls()

# %%
#tws = [r("tw_ni[[2]]"), r("tw_imp[[2]]")]

# %%
def get_sig_tracywidom(tw_p):
    ps = []
    for i, p in enumerate(tw_p):
        if p > 0.05:
            #print(i, p)
            break
        else:
            ps.append(p)
    return len(ps), ps
    


# %%
tw_nums = {}
for key in sorted(tws):
    pvals = tws[key]
    t = get_sig_tracywidom(pvals)[0]
    print(key, t)
    tw_nums[key] = t

# %% [markdown]
# ### Tracy-Widom
#
# | Sample| TW |
# |-------|---------|
# | /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/P/pca_std_data.txt | 1 |
# | /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/P/pca_std_data.txt | 2 |
# | /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/G/pca_std_data.txt | 11 |
# | /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/T/pca_std_data.txt | 9 |
# | /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/E/pca_std_data.txt | 2 |
# | /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/T/pca_std_data.txt | 9 |
# | /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/E/pca_std_data.txt | 3 |
# | /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/G/pca_std_data.txt | 11 |
#

# %%
for key, df in z12_swapped.items():
    d = os.path.dirname(key)
    save_df(d, "z12_swapped", df)

# %%
pca_cov = {}
for k, v in tw_nums.items():
    pca_x_key = os.path.join(os.path.dirname(k), "pca_res.rds")
    x = pca_x[pca_x_key]
    pca_cov[k] = x.ix[:,0:v]

# %%
for key in sorted(pca_cov):
    out = os.path.join(os.path.dirname(key), "{}-{}-pca_cov.txt".format(os.path.basename(os.path.dirname(key)), 
                                                                        imputed_name(key)))
    print(out)
    pca_cov[key].to_csv(out, sep="\t", header=True, index=True)
