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
sys.path.append("/home/cfriedline/ipynb/include_utils")
import include_utils as u
import gzip
import shutil
import tempfile
from ipyparallel import Client
import scandir
import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pickle
import Levenshtein as lv
from Bio import SeqIO
from subprocess import Popen, PIPE, call, check_output
import numpy as np

# %%
rc = Client(profile="sge")

# %%
dview = rc[:]
lview = rc.load_balanced_view()
len(dview)

# %%
rootdir = "/home/cfriedline/eckertlab/projects/burt/seq"

# %%
# organize files for demultiplexing with GBSX
for seqdir in ['160525', '160520']:
    os.chdir('{}/{}'.format(rootdir, seqdir))
    for i in range(8):
        d = "Burt{}".format(i)
        if not os.path.exists(d):
            os.mkdir(d)
        f = glob.glob("{}*.gz".format(d))
        if f:
            assert len(f) == 1
            shutil.move(f[0], d)

# %%
cd $rootdir

# %%
src_files = []
for seqdir in ['160525', '160520']:
    os.chdir('{}/{}'.format(rootdir, seqdir))
    files = !find . -name 'Burt*.gz'
    files = [os.path.abspath(x) for x in files]
    for x in files:
        src_files.append(x)
src_files = sorted(src_files)

# %%
src_files


# %%
#repair beginning Ns before mapping
def check_bc(seq, bc_len):
    bc = seq[0:bc_len]
    print(bc_len, bc)
    min_dist = 100
    min_bc = None
    for b in bc_lens[bc_len]:
        dist = lv.hamming(b, bc)
        if dist < min_dist:
            min_dist = dist
            min_bc = b
    return bc_len, min_dist, min_bc

def convert_ascii(qual):
    return [(ord(x)-33) for x in qual]

for s in src_files:
    print(s)
    for name, seq, qual in FastqGeneralIterator(gzip.open(s, "rt")):
        if seq.startswith("N"):
            print(seq)
            print(qual)
            
            res = []
            for i in range(8, 11):
                res.append(check_bc(seq, i))
            res = sorted(res, key=lambda x: x[1])
            print(res)
            seq2 = res[0][2] + seq[res[0][0]:]
            qual2 = "I"+qual[1:] #I = ASCII 73
            print(seq2)
            print(qual2)
            print(convert_ascii(qual2))
            break
    break

# %%
for s in src_files:
    print(s)
    for rec in SeqIO.parse(gzip.open(s, "rt"), "fastq"):
        print(rec)
        print(rec.letter_annotations)
        break
    break


# %%
def build_gbsx_cmd(fastq, bc, enz):
    cmd = "/home/cfriedline/g/src/jdk1.8.0_92/bin/java -jar /home/cfriedline/g/src/GBSX/GBSX_v1.2.jar --Demultiplexer"
    return "{} -f1 {} -i {} -gzip true -rad true -mb 2 -me 1 -ea {}".format(cmd, fastq, bc, enz), "gbsx"

def write_qsub(workdir, cmd, label, run, cmd_label):
    with open(os.path.join(workdir, "run_{}.sh".format(cmd_label)), "w") as o:
        o.write("""#!/bin/bash
#$ -N {4}{0}
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -e {4}_{3}_burt_{0}.err
#$ -o {4}_{3}_burt_{0}.out
cd {1}
{2}
""".format(label, workdir, cmd, run, cmd_label))

for s in src_files:
    run = os.path.basename(os.path.dirname(os.path.dirname(s)))
    label = os.path.basename(s).split("Burt")[1].split("_")[0]
    bc_file = os.path.join(rootdir, "barcode_{}_gbsx.txt".format(label))
    enz_file = os.path.join(rootdir, "ecori.txt")
    workdir = os.path.dirname(s)
    gbsx_cmd, gbsx_label = build_gbsx_cmd(s, bc_file, enz_file)
    write_qsub(workdir, gbsx_cmd, label, run, gbsx_label)

# %%
cd $rootdir

# %% [markdown]
# ### Submit jobs to SGE
# ```
# cd /gpfs_fs/home/eckertlab/projects/burt/seq
# find . -name "run_gbsx.sh" | xargs chmod +x
# find . -name "run_gbsx.sh" -exec qsub {} \;
# ```

# %%
cd $rootdir

# %%
fastq_files = !find . -name "*.fastq.gz" | grep -v undet
fastq_files = [os.path.abspath(x) for x in fastq_files]
for s in src_files:
    fastq_files.remove(s)
len(fastq_files)

# %%
pwd

# %%
# !mkdir /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe

# %%
cd dedupe

# %%
fastq_dict = {}
for f in fastq_files:
    name = os.path.basename(f)
    if not name in fastq_dict:
        fastq_dict[name] = []
    fastq_dict[name].append(f)

# %%
assert len(fastq_dict) == 768 # (8 lanes * 96/plate)


# %%
def combine_fastq(args):
    name, fastq_list = args
    out_dir = "/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe"
    out_file = os.path.join(out_dir, name)
    cmd = "zcat {} | /home/cfriedline/bin/bgzip -c > {}".format(" ".join(fastq_list), out_file)
    return cmd


# %%
def run_cmd(cmd):
    res = !$cmd
    return res


# %%
dview['run_cmd'] = run_cmd
dview['combine_fastq'] = combine_fastq

# %%
with dview.sync_imports():
    import os

# %%
jobs = []
for k, v in fastq_dict.items():
    cmd = combine_fastq((k, v))
    jobs.append(lview.apply_async(run_cmd, cmd))

# %%
np.sum([x.ready() for x in jobs])

# %%
pwd

# %%
fastq_dedupe = !find '/gpfs_fs/home/eckertlab/projects/burt/seq/dedupe' -name "*.fastq.gz"

# %%
len(fastq_dedupe) == 768

# %%
with open("fastq_files.txt", "w") as o:
    for f in fastq_dedupe:
        o.write("{}\n".format(f))

# %%
