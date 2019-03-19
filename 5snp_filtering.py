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

sys.path.append("../include_utils")

#from IPython.parallel import Client
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
from ipyparallel import Client
import pysam


# %%
def setup_r():
    os.environ['R_HOME'] = '/home/cfriedline/g/R3/lib64/R'
    os.environ['LD_LIBRARY_PATH'] = "%s/lib:%s:%s" % (os.environ['R_HOME'], 
                                                   os.environ['LD_LIBRARY_PATH'],
                                                     "/home/cfriedline/lib64")


# %%
setup_r() #skip on mac

# %%
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
r = robjects.r


# %%
# %reload_ext autoreload
# %autoreload 2
# %matplotlib inline
# %reload_ext rpy2.ipython

# %%
def convert_GQ_to_p(q):
    return pow(10,(q/-10.0))


# %%
vcfutils = "perl /home/cfriedline/g/src/bcftools-1.3/vcfutils.pl"
vcftools = "/home/cfriedline/bin/vcftools"
bcftools = "/home/cfriedline/gpfs/src/bcftools-1.3/bcftools"
tabix = "/home/cfriedline/gpfs/src/htslib-1.3/tabix"
bgzip = "/home/cfriedline/gpfs/src/htslib-1.3/bgzip"
java  = "/home/cfriedline/g/src/jdk1.8.0_60/bin/java"
plink = "/home/cfriedline/g/src/plink-1.07-x86_64/plink --noweb"
plink2 = "/home/cfriedline/g/src/plink_beta_3.29/plink"

# For Mac
# vcfutils = "perl /Users/chris/src/bcftools-1.3/vcfutils.pl"
# vcftools = "/Users/chris/bin/vcftools"
# bcftools = "/Users/chris/src/bcftools-1.3/bcftools"
# tabix = "/Users/chris/src/htslib-1.3/tabix"
# bgzip = "/Users/chris/src/htslib-1.3/bgzip"

# %%
analysis_dir = '/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3'
vcf_file = os.path.join(analysis_dir, "concat2.vcf.gz")
assert os.path.exists(vcf_file)
vcf_file

# %%
!$vcftools --remove-indels \
--max-missing 0.5 \
--min-alleles 2 \
--max-alleles 2 \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--gzvcf \
$vcf_file \
--out $vcf_file

# %% [markdown]
# ```
# VCFtools - 0.1.14
# (C) Adam Auton and Anthony Marcketta 2009
#
# Parameters as interpreted:
# 	--gzvcf /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/concat2.vcf.gz
# 	--recode-INFO-all
# 	--max-alleles 2
# 	--min-alleles 2
# 	--max-missing 0.5
# 	--out /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/concat2.vcf.gz
# 	--recode
# 	--remove-filtered-all
# 	--remove-indels
#
# Using zlib version: 1.2.8
# After filtering, kept 768 out of 768 Individuals
# Outputting VCF file...
# After filtering, kept 643619 out of a possible 36629450 Sites
# Run Time = 13095.00 seconds
# ```

# %%
vcf_filtered = "%s.recode.vcf" % vcf_file
vcf_filtered_gz = "%s.gz" % vcf_filtered

# %%
!$bgzip -c $vcf_filtered > {vcf_filtered_gz}
!$tabix {vcf_filtered_gz}


# %%
def get_vcf_stats(args):
    vcftools, vcf_gz, stat = args
    res = !$vcftools --gzvcf $vcf_gz --out $vcf_gz {"--%s" % stat} 
    return stat


# %%
samples = []
for x in pysam.VariantFile(vcf_filtered_gz):
    samples = list(x.samples)
    break

with open("evolution2016_samples.txt", "w") as o:
    for s in samples:
        o.write("{}\n".format(s))

# %%
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/T
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/E
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/P
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/G

# %%
spp_dict = {}
for s in samples:
    spp = s[0]
    if not spp in spp_dict:
        spp_dict[spp] = []
    spp_dict[spp].append(s)

# %%
filedir = "/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3"
cmds = []
for spp in spp_dict:
    spp_dir = os.path.join(filedir, spp)
    with open(os.path.join(spp_dir, "spp.txt"), "w") as o:
        for elem in spp_dict[spp]:
            o.write("{}\n".format(elem))
    out_vcf = os.path.join(spp_dir, "{}-snps.vcf".format(spp))
    cmd = "{3} --keep {0} --remove-filtered-all --recode --recode-INFO-all --gzvcf {1} --out {2}".format(o.name, vcf_filtered_gz, out_vcf, vcftools)
    cmds.append(cmd)
        

# %%
cmds[1:]

# %%
rc = Client(profile="sge")

# %%
dv, lv = u.get_views(rc)
len(dv)


# %%
def run_cmd(cmd):
    res  = !$cmd
    return res


# %%
dv['run_cmd'] = run_cmd

# %%
by_spp = lv.map_async(run_cmd, cmds[1:])

# %%
by_spp.progress

# %%
stats = ['depth',
            'site-depth',
            'site-mean-depth',
            'site-quality',
            'missing-indv',
            'missing-site',
            'freq',
            'counts',
            'hardy',
            'het']

# %%
dv['get_vcf_stats'] = get_vcf_stats

# %%
stat_jobs = []
for s in stats:
    stat_jobs.append(lv.apply_async(get_vcf_stats, (vcftools, vcf_filtered_gz, s)))

# %%
split_vcfs = !find $filedir -name "*-snps.vcf.recode.vcf"

# %%
split_vcfs

# %%
split_stat_args = []
for v in split_vcfs:
    for s in stats:
        split_stat_args.append((vcftools, v, s))

# %%
split_stat_jobs = lv.map_async(get_vcf_stats, split_stat_args)

# %%
split_stat_jobs.progress, len(split_stat_jobs)

# %%
pd.set_option('display.max_columns', 100)

def get_MAF(row):
    import numpy as np
    try:
        return np.min([row.A1_freq, row.A2_freq])
    except:
        traceback.print_exc()
        
def get_correction(n):
    #for finite sample size
    return (2*n)/(2*n-1)

def calculate_Fis(vals):
    import numpy as np
    try:
        data = [float(x) for x in vals.split("/")]
        assert len(data) == 3
        num_individuals = np.sum(data)
        total_alleles = 2*num_individuals
        a1_count = 2*data[0]
        a2_count = 2*data[2]
        het_count = data[1]
        a1_count += het_count
        a2_count += het_count
        a1_freq = a1_count/total_alleles
        a2_freq = a2_count/total_alleles
        assert a1_freq + a2_freq == 1.0
        He = 2 * a1_freq * a2_freq * get_correction(num_individuals)
        Ho = het_count/num_individuals
        Fis = 1 - (Ho/He)
        return Fis
    except:
        return -9

def combine_vcf_stats(args):
    filedir, prefix = args
    print(filedir, prefix)
    import pandas as pd
    import numpy as np
    hardy_files = !ls {filedir}/{prefix}*.hwe
    print(hardy_files)
    hardy = pd.read_csv(hardy_files[0], sep="\t")

    hardy.columns = ['CHROM', 'POS', 'OBS(HOM1/HET/HOM2)', 'E(HOM1/HET/HOM2)', 'ChiSq_HWE',
       'P_HWE', 'P_HET_DEFICIT', 'P_HET_EXCESS']
    hardy.index = hardy.apply(lambda x: "%s-%d" % (x.CHROM, x.POS), axis=1)
    
    loci_files = !ls {filedir}/{prefix}*.l* | grep -v log
    print(loci_files)
    loci_df = pd.concat([pd.read_csv(x, sep="\t", skiprows=0) for x in loci_files], axis=1)
    chrom_pos = loci_df.ix[:,0:2]
    
    frq_files = !ls {filedir}/{prefix}*.frq* | grep -v count
    print(frq_files)
    frq_data = []
    h = open(frq_files[0])
    header = h.readline().strip().split()
    for line in h:
        frq_data.append(line.strip().split('\t'))

    header = ['CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'A1_FREQ', "A2_FREQ"]
    frq_df = pd.DataFrame(frq_data)
    print(frq_df.columns)
    #frq_df = frq_df.drop([6,7],axis=1)
    frq_df.columns = header
    frq_df.index = frq_df.apply(lambda x: "%s-%s" % (x.CHROM, x.POS), axis=1)
    
    loci_df = loci_df.drop(['CHROM','CHR','POS'], axis=1)
    loci_df = pd.concat([chrom_pos, loci_df], axis=1)
    loci_df.index = loci_df.apply(lambda x: "%s-%d" % (x.CHROM, x.POS), axis=1)
    
    loci_df = pd.concat([loci_df, frq_df, hardy], axis=1)
    loci_df["A1_allele"] = loci_df.apply(lambda row: row.A1_FREQ.split(":")[0], axis=1)
    loci_df["A2_allele"] = loci_df.apply(lambda row: row.A2_FREQ.split(":")[0], axis=1)
    
    loci_df["A1_freq"] = loci_df.apply(lambda row: float(row.A1_FREQ.split(":")[1]), axis=1)
    loci_df["A2_freq"] = loci_df.apply(lambda row: float(row.A2_FREQ.split(":")[1]), axis=1)
    
    loci_df['MAF'] = loci_df.apply(get_MAF, axis=1)
    loci_df = loci_df.drop(['CHROM', 'POS'], axis=1)
    
    loci_df['Fis'] = loci_df['OBS(HOM1/HET/HOM2)'].apply(calculate_Fis)
    
    return loci_df, frq_df, hardy


# %%
dv['combine_vcf_stats'] = combine_vcf_stats
dv['calculate_Fis'] = calculate_Fis
dv['get_MAF'] = get_MAF
dv['get_correction'] = get_correction

# %%
spp_dirs = ["T", "E", "P", "G"]

# %%
combined_stats_jobs = {}
for spp in spp_dirs:
    print(spp)
    d = os.path.join(analysis_dir, spp)
    combined_stats_jobs[spp] = lv.apply_async(combine_vcf_stats, (d, spp))

# %%
for k, v in combined_stats_jobs.items():
    print(v.ready())

# %%
combined_stats = {}
for spp in combined_stats_jobs:
    combined_stats[spp] = combined_stats_jobs[spp].r

# %%
combined_stats["T"][0].head()

# %%
import pickle

# %%
pickle.dump(combined_stats, open(os.path.join(analysis_dir, "combined_stats.pkl"), "wb"), pickle.HIGHEST_PROTOCOL)

# %%
#loci_df, frq_df, hardy = combine_vcf_stats(analysis_dir, "samtools")

# %% [markdown]
# ## Impute genotypes with beagle
#
# ```bash
# cd /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3
# mkdir beagle40
# cd beagle40
# ln -s ../concat2.vcf.gz.recode.vcf.gz
#
# > cat beagle.q
# #$ -S /bin/bash
# #$ -cwd
# #$ -V
# #$ -N beagle
# #$ -pe smp 64
# #$ -o beagle.out
# #$ -e beagle.err
# #$ -q godel199@godel97
#
# ~/g/src/jdk1.8.0_92/bin/java -jar ~/g/src/BEAGLE4/beagle.r1399.jar \
# gl=concat2.vcf.gz.recode.vcf.gz \
# out=/tmp/cfriedline/beagle40 \
# nthreads=64 \
# phase-its=10 \
# burnin-its=10 \
# impute-its=10
#
# cp /tmp/cfriedline/beagle40* /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40
# ```

# %%
beagle_dir = os.path.join(analysis_dir, "beagle40")

# %%
beagle_vcf_gz = os.path.join(beagle_dir, "beagle40.vcf.gz")

# %%
assert os.path.exists(beagle_vcf_gz)

# %%
# !vcftools --gzvcf {beagle_vcf_gz}

# %%
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/T
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/E
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/P
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/G

# %%
imp_cmds = []
for spp in spp_dict:
    spp_dir = os.path.join("/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40", spp)
    with open(os.path.join(spp_dir, "spp.txt"), "w") as o:
        for elem in spp_dict[spp]:
            o.write("{}\n".format(elem))
    out_vcf = os.path.join(spp_dir, "{}-snps.vcf".format(spp))
    cmd = "{3} --keep {0} --remove-filtered-all --recode --recode-INFO-all --gzvcf {1} --out {2}".format(o.name, beagle_vcf_gz, out_vcf, vcftools)
    imp_cmds.append(cmd)

# %%
imp_cmds

# %%
imp_by_spp = lv.map_async(run_cmd, imp_cmds)

# %%
imp_split_vcfs = !find $filedir/beagle40 -name "*-snps.vcf.recode.vcf"

# %%
imp_split_stat_args = []
for v in imp_split_vcfs:
    for s in stats:
        imp_split_stat_args.append((vcftools, v, s))

# %%
imp_split_stat_jobs = lv.map_async(get_vcf_stats, imp_split_stat_args)

# %%
imp_split_stat_jobs.progress, len(imp_split_stat_jobs)

# %%
imp_combined_stats_jobs = {}
for spp in spp_dirs:
    print(spp)
    d = os.path.join("{}/beagle40".format(analysis_dir), spp)
    imp_combined_stats_jobs[spp] = lv.apply_async(combine_vcf_stats, (d, spp))

# %%
for k, v in imp_combined_stats_jobs.items():
    print(v.ready())

# %%
imp_combined_stats = {}
for spp in imp_combined_stats_jobs:
    imp_combined_stats[spp] = imp_combined_stats_jobs[spp].r

# %% [markdown]
# ### Combined stats
#
# For each species, combined_stats is:
#
# 1. loci_df
# 2. freq_df
# 3. hardy

# %%
pickle.dump(imp_combined_stats, open(os.path.join("{}/beagle40".format(analysis_dir), "combined_stats.pkl"), "wb"), 
            pickle.HIGHEST_PROTOCOL)

# %%
for spp in combined_stats:
    loci_df = combined_stats[spp][0]
    chroms = sorted(set([x.split("-")[0] for x in loci_df.index]))
    out_dir = os.path.join(analysis_dir, spp)
    with open(os.path.join(out_dir, "chrom_map.txt"), "w") as o:
        for i, c in enumerate(chroms):
            o.write("%s\t%d\n" % (c, i))

# %%
for spp in imp_combined_stats:
    loci_df = imp_combined_stats[spp][0]
    chroms = sorted(set([x.split("-")[0] for x in loci_df.index]))
    out_dir = os.path.join("{}/beagle40".format(analysis_dir), spp)
    with open(os.path.join(out_dir, "chrom_map.txt"), "w") as o:
        for i, c in enumerate(chroms):
            o.write("%s\t%d\n" % (c, i))


# %%
def write_plink_files(args):
    vcftools, vcf_gz, chrom_map = args
    cmd = "{0} --gzvcf {1} --out {1} --plink --chrom-map {2}".format(vcftools, vcf_gz, chrom_map)
    return cmd


# %%
plink_args = []
for spp in combined_stats:
    v = os.path.join("{}/{}".format(analysis_dir, spp), "{}-snps.vcf.recode.vcf.gz".format(spp))
    c = os.path.join("{}/{}".format(analysis_dir, spp), "chrom_map.txt")
    plink_args.append(write_plink_files((vcftools, v, c)))

# %%
imp_plink_args = []
for spp in imp_combined_stats:
    v = os.path.join("{}/{}".format("{}/beagle40".format(analysis_dir), spp), "{}-snps.vcf.recode.vcf.gz".format(spp))
    c = os.path.join("{}/{}".format("{}/beagle40".format(analysis_dir), spp), "chrom_map.txt")
    imp_plink_args.append(write_plink_files((vcftools, v, c)))

# %%
plink_args

# %%
imp_plink_args

# %%
plink_jobs = lv.map_async(run_cmd, plink_args)

# %%
imp_plink_jobs = lv.map_async(run_cmd, imp_plink_args)

# %%
plink_jobs.progress

# %%
imp_plink_jobs.progress


# %%
def write_plink_recode(args):
    plink, vcf_gz = args
    cmd = "{0} --recodeA --tab --file {1} --out {1}_recodeA".format(plink, vcf_gz)
    return cmd


# %%
peds = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3  -maxdepth 2 -name "*.ped"

# %%
imp_peds = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40  -maxdepth 2 -name "*.ped"

# %%
imp_peds

# %%
recode_args = []
for p in peds:
    v = p.replace(".ped", "")
    recode_args.append(write_plink_recode((plink, v)))

# %%
imp_recode_args = []
for p in imp_peds:
    v = p.replace(".ped", "")
    imp_recode_args.append(write_plink_recode((plink, v)))

# %%
imp_recode_args

# %%
recode_jobs = lv.map_async(run_cmd, recode_args)

# %%
recode_jobs.progress

# %%
imp_recode_jobs = lv.map_async(run_cmd, imp_recode_args)

# %%
imp_recode_jobs.progress

# %%
for spp in combined_stats:
    loci_df = combined_stats[spp][0]
    print(spp, loci_df.SUM_DEPTH.describe())
    loci_df.to_csv(os.path.join(analysis_dir, "{}_loci_stats.txt".format(spp)),
              sep="\t",
              index=False)

# %%
for spp in combined_stats:
    loci_df = combined_stats[spp][0]
    print(spp)
    print(len(loci_df[loci_df.Fis == -9]))
    print(len(loci_df[loci_df.QUAL >= 10]) - len(loci_df[loci_df.QUAL >= 20]))
    print(len(loci_df[loci_df.QUAL < 20]), len(loci_df[loci_df.QUAL < 10]))
    print(len(loci_df[loci_df.Fis >= 0.5]), len(loci_df[loci_df.Fis <= -0.5]), len(loci_df[loci_df.MAF < 0.01]))


# %%
def filter_snps(df, imputed=False):
    if imputed:
        return df[(df.MAF >= 0.01) & 
                  (df.Fis < 0.5) & 
                  (df.Fis > -0.5)]
    else:
        return df[(df.SUM_DEPTH >= 50) & 
                  (df.SUM_DEPTH < 1500) & 
                  (df.QUAL >= 20) & 
                  (df.MAF >= 0.01) & 
                  (df.Fis < 0.5) & 
                  (df.Fis > -0.5)]


# %%
loci_stage1 = {}
for spp in spp_dirs:
    loci_stage1[spp] = filter_snps(combined_stats[spp][0])
    print(spp, loci_stage1[spp].shape)

# %%
beagle_stage1 = {}
for spp in spp_dirs:
    beagle_stage1[spp] = filter_snps(imp_combined_stats[spp][0], imputed=True)
    print(spp, beagle_stage1[spp].shape)

# %%
for spp in spp_dirs:
    with open(os.path.join("{}/{}".format(analysis_dir, spp), "stage1_positions.txt"), "w") as o:
        for elem in loci_stage1[spp].index:
            o.write("%s\n" % "\t".join(elem.split("-")))

    with open(os.path.join("{}/{}".format(beagle_dir, spp), "stage1_positions.txt"), "w") as o:
        for elem in beagle_stage1[spp].index:
            o.write("%s\n" % "\t".join(elem.split("-")))
    

# %% {"run_control": {"marked": false}}
good_args = []
for spp in spp_dirs:
    for d in [analysis_dir, beagle_dir]:
        d = os.path.join(d, spp)
        v = os.path.join(d, "{}-snps.vcf.recode.vcf.gz".format(spp))
        p = os.path.join(d, "stage1_positions.txt")
        out = os.path.join(d, "good_snps.vcf")
        cmd = "{} --gzvcf {} --remove-filtered-all --recode --recode-INFO-all --positions {} --out {}".format(vcftools, v, p, out)
        good_args.append(cmd)

# %%
good_jobs = lv.map_async(run_cmd, good_args)

# %%
good_jobs.progress, len(good_jobs)


# %% [markdown]
# ### Zip and index good snps
#
# ```
# find . -name "good*.recode.vcf.gz" | parallel tabix {}
# find . -name "good*.recode.vcf" | parallel bgzip -c {} \> {}.gz
# ```

# %%
def get_intersection(imp, ni):
    return set.intersection(set(ni.index), set(imp.index))


# %%
isect = {}
for spp in spp_dirs:
    isect[spp] = get_intersection(beagle_stage1[spp], loci_stage1[spp])
    isect[spp] = sorted(isect[spp])

# %%
for spp in isect:
    print(spp, len(loci_stage1[spp].index), len(beagle_stage1[spp].index), len(isect[spp]))

# %%
for spp in spp_dirs:
    for d in [analysis_dir, beagle_dir]:
        d = os.path.join(d, spp)
        with open(os.path.join(d, "isect_positions.txt"), "w") as o:
            for elem in isect[spp]:
                o.write("%s\n" % "\t".join(elem.split("-")))


# %%
isect_args = []
for spp in spp_dirs:
    for d, vcf_gz in zip([analysis_dir, beagle_dir], [vcf_filtered_gz, beagle_vcf_gz]):
        d = os.path.join(d, spp)
        v = os.path.join(d, "good_snps.vcf.recode.vcf.gz")
        p = os.path.join(d, "isect_positions.txt")
        o = os.path.join(d, "isect_snps")
        cmd = "{} --gzvcf {} --remove-filtered-all --recode --recode-INFO-all --positions {} --out {}".format(vcftools, v, p, o)
        isect_args.append(cmd)

# %%
isect_args

# %%
isect_jobs = lv.map_async(run_cmd, isect_args)

# %%
isect_jobs.progress, len(isect_jobs)

# %% [markdown]
# ### zip and index isect snps
#
# ```
# find . -name "isect_snps.recode.vcf" | parallel --bar bgzip -c {} \> {}.gz
# find . -name "isect_snps.recode.vcf.gz" | parallel --bar tabix
# ```
#
# ### sort, zip, index
#
# ```
# find . -name "isect_snps.recode.vcf.gz" | parallel --bar vcf-sort {} \> {}_sorted.vcf
# find . -name "isect_snps.recode.vcf.gz_sorted.vcf" | parallel --bar bgzip -c {} \> {}.gz
# find . -name "isect_snps.recode.vcf.gz_sorted.vcf.gz" | parallel --bar tabix
# ```
#
# ### thin, zip, index
#
# ```
# find . -name "isect_snps.recode.vcf.gz_sorted.vcf.gz" | parallel --bar vcftools --gzvcf {} --remove-filtered-all --recode --recode-INFO-all --thin 50 --out {}_thin
# find . -name "*thin.recode.vcf" | parallel --bar bgzip -c {} \> {}.gz
# find . -name "*thin.recode.vcf.gz" | parallel --bar tabix
# ```
#

# %%
plink_file_cmds = []
plink_recode_cmds = []
for spp in spp_dirs:
    for d in [analysis_dir, beagle_dir]:
        d = os.path.join(analysis_dir, spp)
        v = os.path.join(d, "isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz")
        assert os.path.exists(f)
        c = os.path.join("{}/{}".format(analysis_dir, spp), "chrom_map.txt")
        assert os.path.exists(c)
        plink_file_cmds.append(write_plink_files((vcftools, v, c)))
        plink_recode_cmds.append(write_plink_recode((plink, v)))
#     write_plink_files(f)
#     write_plink_recode(f)

# %%
plink_file_cmds

# %%
plink_recode_cmds

# %%
file_jobs = lv.map_async(run_cmd, plink_file_cmds)

# %%
file_jobs.progress, len(file_jobs)

# %%
recode_jobs = lv.map_async(run_cmd, plink_recode_cmds)

# %%
recode_jobs.progress, len(recode_jobs)


# %% [markdown]
# ### output 012 files
#
# ```
# find . -name "isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz" | parallel --bar vcftools --gzvcf {} --012 --out {}
# ```

# %%
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/merged
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/merged

# %%
def write_merge_command(args):
    bcftools, outdir, vcfs = args
    cmd = "{} merge {} -Oz -o {} --threads 16".format(bcftools, " ".join(vcfs), os.path.join(outdir, "merged.vcf.gz"))
    return cmd


# %%
# merge non imputed files

ni_spp_files = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3 -maxdepth 2 -name "isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz" 
imp_spp_files = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40 -name "isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz" 
print(write_merge_command((bcftools, "/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/merged", ni_spp_files)))
print()
print(write_merge_command((bcftools, "/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle/merged", imp_spp_files)))


# %% [markdown]
# ### merge commands
#
# ```
# /home/cfriedline/gpfs/src/bcftools-1.3/bcftools merge \
# /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/E/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz \
# /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/P/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz \
# /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/G/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz \
# /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/T/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz \ 
# -Oz -o /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/merged/merged.vcf.gz \
# --threads 16
# ```
#
# ```
# /home/cfriedline/gpfs/src/bcftools-1.3/bcftools merge \
# /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/T/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz \
# /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/E/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz \
# /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/P/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz \
# /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/G/isect_snps.recode.vcf.gz_sorted.vcf.gz_thin.recode.vcf.gz \
# -Oz -o /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/beagle40/merged/merged.vcf.gz \
# --threads 16
# ```
#
# ### create 012 files
#
# ```
# find . -name "merged.vcf.gz" | parallel --bar tabix {}
# find . -name "merged.vcf.gz" | parallel --bar vcftools --gzvcf {} --thin 50 --recode --recode-INFO-all --out {}_thin
# find . -name "merged.vcf.gz_thin*" | parallel bgzip -c {} \> {}.gz
# find . -name "merged.vcf.gz_thin*.gz" | parallel tabix {}
# find . -name "merged.vcf.gz_thin*.gz" | parallel --bar vcftools --gzvcf {} --012 --out {}
# ```

# %%
merged_files = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3 -name "merged.vcf.gz_thin.recode.vcf.gz"

# %%
merged_stat_args = []
for m in merged_files:
    for s in stats:
        merged_stat_args.append((vcftools, m, s))

# %%
merged_stat_jobs = lv.map_async(get_vcf_stats, merged_stat_args)

# %%
merged_stat_jobs.progress, len(merged_stat_jobs)

# %%
merged_combined_stats_args = []
for m in merged_files:
    merged_combined_stats_args.append((os.path.dirname(m), "merged"))

# %%
merged_combined_stats_args

# %%
