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
import urllib.request as urllib2
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
import dill
from scipy import stats
from IPython.display import display
import geopy

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
# ni_dir = "/home/cfriedline/eckertlab/gypsy_indiv/raw_demult/analysis/samtools1.3_masurca3/ni"
# imp_dir = "/home/cfriedline/eckertlab/gypsy_indiv/raw_demult/analysis/samtools1.3_masurca3/beagle40"
# notimputed_vcf_gz = os.path.join(ni_dir, "isect_snps.recode.vcf.gz_sorted.vcf.gz")
# imputed_vcf_gz = os.path.join(imp_dir, "isect_snps.recode.vcf.gz_sorted.vcf.gz")
# vcf_files = [notimputed_vcf_gz, imputed_vcf_gz]

# %%
hierf_trans = {0:11, 1:12, 2:22, -1:'NA'}
def apply_hierf_trans(series):
    return [hierf_trans[x] if x in hierf_trans else x for x in series]


# %%
def read_df(dirname, fname):
    f = os.path.join(dirname, "%s.txt" % fname)
    return pd.read_csv(f, sep="\t", index_col=0)


# %%
z12_swapped = {}
for spp in species:
    for indir in [ni_dir, imp_dir]:
        d = os.path.join(indir, spp)
        z12_swapped[d] = read_df(d, "z12_swapped")

# %%
z12_swapped['/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/E'].head()

# %% {"run_control": {"marked": false}}
hierf_df = {}
for k, df in z12_swapped.items():
    print(k)
    hierf_df[k] = df.apply(apply_hierf_trans)


# %%
def get_pheno():
    pheno = pd.read_csv("/home/cfriedline/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/pheno/evolution2016_sample_info.txt", sep="\t")
    pheno.index = pheno['name']
    pheno = pheno.drop("name", axis=1)
    return pheno
pheno = get_pheno()
pheno['state'] = pheno.population.apply(lambda x: x.split("-")[1])

# %%
pop_id = {}
i = 1
for p in sorted(pheno.population.unique()):
    pop_id[p] = i
    i+=1


# %%
def assign_popid(series):
    p = series.name.rsplit("-", 1)[0]
    series['popid'] = pop_id[p]
    return series


# %%
z12_swapped_pop = {key: value.apply(assign_popid, axis=1) for (key, value) in z12_swapped.items()}

# %%
hierf_df = {}
for k, v in z12_swapped_pop.items():
    hierf_df[k] = v.ix[:,-1:].join(v.ix[:,:-1])

# %%
for outdir, df in hierf_df.items():
    outfile = os.path.join(outdir, "isect_hierfstat.txt")
    print(outfile)
    df.to_csv(outfile, header=True, index=False, sep="\t")

# %%
scr = """library(hierfstat)
library(data.table)
data = fread("isect_hierfstat.txt", header=T, sep="\\t", data.table=F)
levels = data.frame(data$popid)
loci = data[,2:ncol(data)]
bs = basic.stats(data)
saveRDS(bs, "isect_hierfstat_basic_stats.rds")
res = varcomp.glob(levels=levels, loci=loci, diploid=T)
saveRDS(res, "isect_hierfstat_varcomp.rds")
"""
for outdir, df in hierf_df.items():
    with open(os.path.join(outdir, "hierf.R"), "w") as o:
        o.write("{}\n".format(scr))
        
    with open(os.path.join(outdir, "hierf.q"), "w") as o:
        o.write("""#$ -S /bin/bash
#$ -N hierf
#$ -cwd
#$ -V
#$ -j y
cd {}
/home/cfriedline/g/R3/bin/R --vanilla < hierf.R
""".format(outdir))

# %% [markdown]
# # Put into R (because it can be slow)
#
# ```
# find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3 -name "hierf.q" -exec qsub {} \;
# ```
#
# ### Also compute pairwise Fst (qrsh)
# ```
# rm(list=ls())
# library(hierfstat)
# library(data.table)
# library(snow)
# data = fread("isect_hierfstat.txt", header=T, sep="\t", data.table=F)
# levels = data.frame(data$popid)
# loci = data[,2:ncol(data)]
#
# run_varcomp = function(idx) {
#     i = args[[idx]]$i
#     j = args[[idx]]$j
#     key = paste(i, j, sep="-")
#     d = copy(data)
#     d = subset(d, d$popid == i | d$popid == j)
#     levels = data.frame(d$popid)
#     loci = d[,2:ncol(d)]
#     return(varcomp.glob(levels=levels, loci=loci, diploid=T))
# }
#
# args = list()
# for (i in 1:6) {
#     for (j in 1:i) {
#         if (i != j) {
#             args[[length(args)+1]] = list(i=i, j=j)
#         }
#     }
# }
#
# hosts = rep("localhost", length(args))
# cl = makeSOCKcluster(hosts)
# clusterExport(cl, "data")
# clusterEvalQ(cl, library(hierfstat))
# clusterEvalQ(cl, library(data.table))
# clusterExport(cl, "args")
# clusterExport(cl, "run_varcomp")
# pairwise_res = parLapply(cl, 1:length(args), "run_varcomp")
# saveRDS(pairwise_res, "isect_hierfstat_pairwise.rds")
# saveRDS(args, "isect_hierfstat_pairwise_args.rds")
# stopCluster(cl)
#
# ```

# %%
hierf_rds = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3 -name "isect*hierf*.rds"

# %%
assert len(hierf_rds) == 16

# %%
hierf_results = {}
for h in hierf_rds:
    key = os.path.dirname(h)
    if not key in hierf_results:
        hierf_results[key] = {}
    if "basic_stats" in h:
        hierf_results[key]['bs'] = h
    else:
        hierf_results[key]['varcomp'] = h

# %%
h = 0
for key, vals in hierf_results.items():
    vals['index'] = h
    r("bs{} = readRDS('{}')".format(h, vals['bs']))
    r("vc{} = readRDS('{}')".format(h, vals['varcomp'])) 
    h+=1


# %% {"language": "R"}
#
# ls()

# %%
def get_r_series(key):
    s = pd.Series(get_r(key))
    s.index = get_r("names(%s)" % key)
    return s

def get_r_df(key):
    df = pd.DataFrame(get_r(key))
    try:
        rname = get_r("rownames(%s)" % key)
        df.index = rname
    except:
        pass
    
    try:
        cname = get_r("colnames(%s)" % key)
        df.columns = cname
    except:
        pass
    
    return df

def get_r(key):
    return r(key)


# %% {"language": "R"}
# print(vc1)

# %%
hierf_keys = sorted(hierf_results.keys())
for key in hierf_keys:
    vals = hierf_results[key]
    print(key)
    print(get_r_df('vc{}$F'.format(vals['index'])))
    print()

# %%
perloc_not = get_r_df("bs_not$perloc")
Ho_not = get_r_df("bs_not$Ho")
Hs_not = get_r_df("bs_not$Hs")
Fis_not = get_r_df("bs_not$Fis")
overall_not = get_r_series("bs_not$overall")
n_ind_samp_not = get_r_df("bs_not$n.ind.samp")

# %%
perloc_imp = get_r_df("bs_imp$perloc")
Ho_imp = get_r_df("bs_imp$Ho")
Hs_imp = get_r_df("bs_imp$Hs")
Fis_imp = get_r_df("bs_imp$Fis")
overall_imp = get_r_series("bs_imp$overall")
n_ind_samp_imp = get_r_df("bs_imp$n.ind.samp")


# %%
def save_df(dirname, fname, df):
    f = os.path.join(dirname, "%s.txt" % fname) 
    df.to_csv(f, 
              header=True,
              index=True,
              sep="\t")
    print("saved %s" % f)


# %%
loc_df_not = get_r_df('varcomp_not$loc')
F_df_not = get_r_df('varcomp_not$F')
overall_df_not = get_r_df('varcomp_not$overall')

# %%
loc_df_imp = get_r_df('varcomp_imp$loc')
F_df_imp = get_r_df('varcomp_imp$F')
overall_df_imp = get_r_df('varcomp_imp$overall')

# %%
loc_df_imp.head()

# %%
F_df_not

# %%
F_df_imp


# %%
def compute_fst(series):
    Hs = series[0]
    Ht = sum(series)
    return Hs/Ht


# %%
loci_fst_not = loc_df_not.apply(compute_fst, axis=1)
loci_fst_imp = loc_df_imp.apply(compute_fst, axis=1)

# %%
loci_fst_not.describe()

# %%
loci_fst_imp.describe()

# %%
for i, d in enumerate([ni_dir, imp_dir]):
    if i == 0:
        save_df(d, 'perloc', perloc_not)
        save_df(d, 'Ho', Ho_not)
        save_df(d, "Hs", Hs_not)
        save_df(d, 'Fis', Fis_not)
        save_df(d, 'overall', overall_not)
        save_df(d, 'n_ind_samp', n_ind_samp_not)
        save_df(d, 'loc_df', loc_df_not)
        save_df(d, 'F_df', F_df_not)
        save_df(d, 'varcomp_overall', overall_df_not)
        save_df(d, 'loci_fst', loci_fst_not)
    else:
        save_df(d, 'perloc', perloc_imp)
        save_df(d, 'Ho', Ho_imp)
        save_df(d, "Hs", Hs_imp)
        save_df(d, 'Fis', Fis_imp)
        save_df(d, 'overall', overall_imp)
        save_df(d, 'n_ind_samp', n_ind_samp_imp)
        save_df(d, 'loc_df', loc_df_imp)
        save_df(d, 'F_df', F_df_imp)
        save_df(d, 'varcomp_overall', overall_df_imp)
        save_df(d, 'loci_fst', loci_fst_imp)

# %%
plt.hist(loci_fst_not, bins=50)
plt.title("not imputed n=%d mean=%.4f +/- %.4f [%.4f, %.4f]" % (len(loci_fst_not), 
                                                    np.mean(loci_fst_not), 
                                                    np.std(loci_fst_not),
                                                    np.min(loci_fst_not), 
                                                    np.max(loci_fst_not)))
plt.xlabel(r"$F_{ST}$")
plt.show()

# %%

sns.set(font_scale=1.5)

sns.set_style("white")

loci_fst_imp = read_df(imp_dir, "loci_fst")
plt.hist(loci_fst_imp["0"], bins=50)
plt.title("imputed n=%d mean=%.4f +/- %.4f [%.4f, %.4f]" % (len(loci_fst_imp), 
                                                    np.mean(loci_fst_imp), 
                                                    np.std(loci_fst_imp),
                                                    np.min(loci_fst_imp), 
                                                    np.max(loci_fst_imp)))
plt.xlabel(r"$F_{ST}$")
fig = plt.gcf()
fig.set_size_inches(10,8)
sns.despine()
out = "fst_hist.pdf"
plt.savefig(out)
plt.show()
display(FileLink(out))

# %%
pwd

# %%
popid_map = {}
for population, data in z12_swapped[1].groupby("population"):
    popid = data['popid'].unique()[0]
    print(population, popid)
    popid_map[popid] = population

# %%
perloc = read_df(imp_dir, "perloc")
loc_df = read_df(imp_dir, "loc_df")
Ho = read_df(imp_dir, "Ho")
Ho.columns = [popid_map[int(x)] for x in Ho.columns]
overall = read_df(imp_dir, "overall")
n_ind_samp = read_df(imp_dir, "n_ind_samp")
n_ind_samp.columns = [[popid_map[int(x)] for x in n_ind_samp.columns]]
fis = read_df(imp_dir, "Fis")
fis.columns = [popid_map[int(x)] for x in fis.columns]

# %%
Ho.to_csv("Ho_labelled.txt", 
          index=True, 
          header=True, 
          sep="\t")

# %%
from IPython.display import FileLink
FileLink("Ho_labelled.txt")

# %%
Ho.head()

# %%
Ho.columns

# %%
stats.f_oneway(Ho['NC'], Ho['NY'], Ho['QC32'], Ho['QC93'], Ho['VA1'], Ho['VA2'])

# %%
sns.set_context("talk")

# %%
sns.boxplot(data=Ho);

# %% {"language": "R"}
# bs_imp$pop.freq$ctg7180005039298_50

# %%
pop_freq = r('bs_imp$pop.freq')

pop_freq_names = pandas2ri.ri2py(pop_freq.names)


# %%
def compute_He(elem):
    He = 2
    for x in elem:
        He *= x
    return He

He_dict = {}

for i, name in enumerate(pop_freq_names):
    af = pandas2ri.ri2py_dataframe(pop_freq.rx2(name))
    af.columns = [x+1 for x in af.columns]
    He_dict[name] = af.apply(compute_He).to_dict()
    
    if i % 10000 == 0:
        print("at %d" % i)

# %%
He = pd.DataFrame(He_dict).T

# %%
He.columns = [popid_map[x] for x in He.columns]

# %%
Ho_He = Ho.join(He, lsuffix = "_Ho", rsuffix = "_He")

# %%
means = pd.DataFrame(Ho_He.apply(np.mean)).T

# %%
means

# %%
diffs = Ho-He

# %%
sns.boxplot(data=diffs)
plt.ylabel(r"$Ho - He$")

# %%
sns.boxplot(data=np.abs(diffs))
plt.ylabel(r"$|Ho-He|$")

# %% {"language": "R"}
# pairwise = readRDS("/gpfs_fs/home/eckertlab/gypsy_indiv/raw_demult/analysis/samtools1.3_masurca3/beagle40/isect_hierfstat_pairwise.rds")
# pairwise_args = readRDS("/gpfs_fs/home/eckertlab/gypsy_indiv/raw_demult/analysis/samtools1.3_masurca3/beagle40/isect_hierfstat_pairwise_args.rds")
#

# %% {"language": "R"}
# pairwise[1]

# %%
pairwise = r("pairwise")
pairwise_args = r("pairwise_args")

# %%
pairwise_fst = np.zeros((7,7))
for arg in range(len(pairwise_args)):
    arg = arg+1
    i = pairwise_args.rx2(arg).rx2("i")[0]
    j = pairwise_args.rx2(arg).rx2("j")[0]
    print(i, j)
    F = pandas2ri.ri2py_dataframe(pairwise.rx2(arg).rx2("F"))
    pairwise_fst[i, j] = F.ix[0,0]

# %%
pd.DataFrame(pairwise_fst)

# %%
pop_dict = {}
for idx, data in z12_swapped[1][['population', 'popid']].iterrows():
    pop_dict[data.popid] = data.population

# %%
np.max(pairwise_df.max())

# %%
pairwise_df = pd.DataFrame(pairwise_fst)
pairwise_df = pairwise_df.drop(0, axis=1).drop(6, axis=1).drop(0)
pairwise_df.columns = [pop_dict[x] for x in pairwise_df]
pairwise_df.index = [pop_dict[x] for x in pairwise_df.index]
pairwise_df.replace(0, False)

# %%
sns.set_context("talk")
cmap = sns.cubehelix_palette(light=1, as_cmap=True)
sns.heatmap(pairwise_df, vmin=0, vmax=pairwise_df.max().max(), cmap=cmap, annot=True)
plt.title("Pairwise multilocus %s" % r'$F_{ST}$')
out = "mfst.pdf"
plt.gcf().set_size_inches(10, 10)
plt.savefig()
plt.show()

# %%
latlon = read_df(imp_dir, "bioclim_df")[["lat", "lon"]]

# %%
latlon

# %%
from geopy.distance import vincenty
dist = {}
for i in range(len(latlon.index)):
    ipop = latlon.index[i]
    icoord = (latlon.ix[ipop].lat, latlon.ix[ipop].lon)
    if not ipop in dist:
        dist[ipop] = {}
    for j in range(i):
        jpop = latlon.index[j]
        jcoord = (latlon.ix[jpop].lat, latlon.ix[jpop].lon)
        d = vincenty(icoord, jcoord).km
        print(ipop, jpop, d)
        dist[ipop][jpop] = d
        if not jpop in dist:
            dist[jpop] = {}
        dist[jpop][ipop] = d

# %%
dist_df = pd.DataFrame(dist)

# %%
dist_df = dist_df.fillna(0)

# %%
full_pairwise = pd.DataFrame(index=pairwise_df.index, columns=pairwise_df.index)

# %%
pairwise_df

# %%
for col in pairwise_df:
    for idx in pairwise_df.index:
        if pairwise_df.ix[idx, col] == 0 and col is not idx:
            pairwise_df.ix[idx, col] = pairwise_df.ix[col, idx]

# %%
pairwise_df = pairwise_df.join(pd.DataFrame(pairwise_df.ix["VA2"])).fillna(0)

# %%
for i, idx in enumerate(pairwise_df.index):
    for j, jdx in enumerate(pairwise_df.index):
        if pairwise_df.ix[idx, jdx] == 0:
            print(idx, jdx)

# %%
fst_gps = {}
for p in pairwise_df.index:
    for q in pairwise_df.index:
        if not p == q:
            print(p, q)
            fst = -1
            try:
                fst = pairwise_df.ix[p, q]
            except:
                fst = pairwise_df.ix[q, p]
            key = "-".join(sorted([p, q]))
            fst_gps[key] = {"vincenty": dist_df.ix[p,q],
                                         "fst": fst}

# %%
fst_gps_df = pd.DataFrame(fst_gps).T

# %%
fst_gps_df

# %%
g = sns.lmplot("fst", "vincenty", data=fst_gps_df, scatter_kws={"s": 100})
g.fig.set_size_inches(10, 10)
for i in fst_gps_df.index:
    plt.annotate(i, 
                 xy=(fst_gps_df.ix[i, "fst"], fst_gps_df.ix[i, "vincenty"]), 
                xytext = (-10, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom', 
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
plt.ylabel("Vincenty (km)")
plt.xlabel(r'Pairwise $F_{ST}$')
plt.show()
#bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
#arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0')

# %%
from statsmodels.formula.api import ols

# %%
ols("fst~vincenty", data=fst_gps_df).fit().summary()

# %%
dist_df

# %%
# %R -i pairwise_df -i fst_gps_df -i dist_df

# %% {"language": "R"}
# pairwise_df

# %% {"language": "R"}
# dist_df

# %% {"language": "R"}
# library(vegan)
# mantel(pairwise_df, dist_df, permutations=10000)

# %%
