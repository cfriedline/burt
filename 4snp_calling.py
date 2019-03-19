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
from subprocess import Popen, PIPE, STDOUT, check_output, call
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
import glob
import tempfile
from ipyparallel import Client

# %%
cd "~/eckertlab/projects/burt/seq/dedupe/"

# %%
bam_dir = "."
analysis_dir = os.path.join(bam_dir, "samtools1.3")
if not os.path.exists(analysis_dir):
    os.makedirs(analysis_dir)
assert os.path.exists(analysis_dir)

# %%
bam_files = !ls *mapped.bam
bam_files = [os.path.abspath(x) for x in bam_files]

# %%
len(bam_files) == 768

# %%
samtools = "/home/cfriedline/bin/samtools"
bcftools = "/home/cfriedline/bin/bcftools"
picard = "/home/cfriedline/gpfs/src/broadinstitute-picard-03a1d72/dist/picard.jar"
java = "/home/cfriedline/g/src/jdk1.8.0_92/bin/java"
perl = "/home/cfriedline/gpfs/opt/ActivePerl-5.18/bin/perl"

# %%
assembly = "/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/conitgs.fa_mapped.fasta"

# %%
!$samtools faidx {assembly}


# %%
def create_split_beds(nodes, bed):
    lines = 0
    for line in open(bed):
        lines += 1
    print(lines, lines//nodes)
    per_bed = lines//nodes
    cmd = "split -a 3 -d -l %d %s contig.bed." % (per_bed, bed)
    call(cmd.split())
create_split_beds(150, "contigs.bed")

# %%
beds = !ls contig.bed.*
beds = [os.path.abspath(x) for x in beds]

# %%
rc = u.get_client(profile="sge")

# %%
dv, lv = u.get_views(rc)
len(dv)

# %%
with dv.sync_imports():
    from subprocess import Popen, PIPE, STDOUT, check_output, call
    import os, sys, socket, glob, tempfile, shutil


# %%
def create_parallel_bams(args):
    samtools, bam_file, bed_file = args
    num = bed_file.split(".")[-1]
    out = "%s.%s" % (bam_file, num)
    t = tempfile.NamedTemporaryFile(delete=False, dir="/tmp")
    cmd = "%s view -L %s -b %s -o %s" % (samtools, bed_file, bam_file, t.name)
    call(cmd.split())
    shutil.copy(t.name, out)
    os.remove(t.name)
    return out


# %%
dv['create_parallel_bams'] = create_parallel_bams

# %%
jobs = []
args = []
for bam in bam_files:
    for bed in beds:
        a = [samtools, bam, bed]
        args.append(a)

# %%

# %%
len(bam_files), len(bam_files)*len(beds)

# %%
jobs = dv.map_async(create_parallel_bams, args)

# %%
jobs.progress

# %% [markdown]
# ### Move `*.bam.*` and `contig.bed.*` files to `./work`

# %%
pwd

# %%
par_bams = !ls *.bam.*

# %%
len(par_bams)

# %%
job_map = {}
for b in par_bams:
    num = b.split(".")[-1]
    if not num in job_map:
        job_map[num] = []
    job_map[num].append(b)

# %%
for num in job_map:
    job_map[num] = sorted(job_map[num])
    assert len(job_map[num]) == 768

# %%
snp_args = []
for num in sorted(job_map):
    bam_files = job_map[num]
    bed_file = "contig.bed.%s" % num
    a = (samtools,
         bed_file,
         assembly, 
         bam_files, 
         bcftools, 
         "samtools_1.3.vcf.gz.%s" % num, 
         'samtools1.3')
    snp_args.append(a)

# %%
snp_args[0]


# %%
def call_snps(args):
    import socket, os, stopwatch
    #print(socket.gethostname())
    timer = stopwatch.Timer()
    samtools, bed, reference, bam_sorted, bcftools, raw_vcf, out_dir = args 
    if not out_dir:
        out_dir = os.environ['TMPDIR']
    raw_vcf = os.path.join(out_dir, os.path.basename(raw_vcf))
    pileup = "%s mpileup -t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -Iugf %s %s | %s call -f GP,GQ -vmO z -o %s" % (samtools, 
                                                                     reference, 
                                                                     ' '.join(bam_sorted), 
                                                                     bcftools,                                                                
                                                                     raw_vcf) 
    
    #print(pileup)
    #!$pileup
    timer.stop()
    return pileup, timer.elapsed


# %%
cmds = []
for a in snp_args:
    cmds.append(call_snps(a)[0])

# %%
with open("jobs.sh", "w") as j:
    j.write("#!/bin/bash\n")
    for i, a in enumerate(snp_args):
        scr = "run%d.sh" % i
        j.write("qsub %s\n" % scr)
        with open(scr, "w") as o:
            header = """#!/bin/bash

#$ -S /bin/bash
#$ -N psnp%d
#$ -cwd
#$ -V
#$ -o /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/psnp%d.out
#$ -e /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3/psnp%d.err
#$ -q all.q""" % (i, i, i)
            
            o.write("%s\n" % header)
            o.write("%s\n" % call_snps(a)[0])

# %% [markdown]
# ## Run on SGE
# ```bash
# cd /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work
# chmod +x *.sh
# ./jobs.sh
# ```

# %%
cd work/samtools1.3/

# %%
vcfs = !ls samtools_1.3.vcf.gz.* | grep -v tbi | grep -v sorted
vcfs = [os.path.abspath(x) for x in vcfs]

# %%
perl = "/home/cfriedline/gpfs/opt/ActivePerl-5.18/bin/perl"
vcf_concat = "{} /home/cfriedline/g/src/vcftools-0.1.14/src/perl/vcf-concat".format(perl)
vcf_sort = "{} /home/cfriedline/g/src/vcftools-0.1.14/src/perl/vcf-sort".format(perl)
bgzip = "/home/cfriedline/g/src/htslib-1.3/bgzip"
tabix = "/home/cfriedline/g/src/htslib-1.3/tabix"

# %% [markdown]
# ### index parallel vcf files
#
# ```
# $ tabix --version                                                                                                       
# tabix (htslib) 1.3.1                                                                                                     
# Copyright (C) 2016 Genome Research Ltd.
# ```
#
# ```
# cd /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3
# ls samtools_1.3.vcf.* | parallel --bar tabix -f {}
# ```

# %%
pwd

# %%
with open("concat2.sh", "w") as o:
    o.write("{} concat --threads 50 -Oz -o concat2.vcf.gz {}\n".format(bcftools, " ".join(vcfs)))

# %% [markdown]
# ### concatenate vcf files
#
# ```
# cd /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/work/samtools1.3
# chmod +x concat2.sh
# ./concat2.sh
# ```
#
# ### sort concat file
# ```
# vcf-sort concat2.vcf.gz > concat2_sorted.vcf.gz
# ```

# %%
