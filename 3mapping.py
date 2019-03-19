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
from ipyparallel import Client
import matplotlib.pyplot as plt
# %matplotlib inline
from subprocess import Popen, PIPE
from Bio import SeqIO
import pandas as pd
import pickle
import scandir
import numpy as np
import tempfile
import socket
from collections import Counter

# %%
sys.path.append("/home/cfriedline/ipynb/include_utils")

# %%
import include_utils as u

# %%
root = "/home/cfriedline/eckertlab/projects/burt/seq/dedupe/"

# %%
cd $root

# %%
pwd

# %%
fastq_files = !ls *.fastq
fastq_files = sorted([os.path.abspath(x) for x in fastq_files])

# %%
len(fastq_files) == 768

# %%
assembly = "/gpfs_fs/home/eckertlab/loblolly2/conitgs.fa"

#indexed as: ~/g/src/bowtie2-2.2.9/bowtie2-build --threads 20 conitgs.fa conitgs.fa -offrate 9 --large-index --threads 40 (b/c it's so big)

# %%
# --very-fast-local
# Same as: -D 5 -R 1 -N 0 -L 25 -i S,1,2.00

# --fast-local
# Same as: -D 10 -R 2 -N 0 -L 22 -i S,1,1.75

# --sensitive-local
# Same as: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default in --local mode)

# --very-sensitive-local
# Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

#@lview.remote()
def run_bowtie2(args):
    import os, stopwatch, multiprocessing, socket
    timer = stopwatch.Timer()
    cpus = 64
    assembly, reads, outdir = args
    sam = os.path.join(outdir, "{}.sam".format(os.path.basename(reads)))
    sam = tempfile.NamedTemporaryFile(delete=False, dir="/tmp")
    cmd = "/home/cfriedline/g/src/bowtie2-2.2.9//bowtie2 --local -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p %d -x %s -U %s -S %s" % (8,
                                                               assembly,
                                                               reads,
                                                               sam.name)
    res = None
    res = cmd
#     if not os.path.exists(sam):
#         res = !$cmd
    timer.stop()
    return assembly, sam.name, cmd, timer.elapsed, res


# %%
# ?tempfile.NamedTemporaryFile

# %%
sam_outdir = "/gpfs_fs/home/eckertlab/projects/burt/bowtie2"

# %%
hosts = []
cpus = {}
qhost = !qhost | grep godel
for q in qhost:
    q = q.split()
    mem = float(q[4][:-1])
    host = q[0]
    if mem > 60:
        hosts.append(host)
        cpus[host] = int(q[2])


# %%
def create_bt_jobs(files):
    bt_jobs = {}
    host_id = 0
    for f in files:
        if host_id == len(hosts):
            host_id = 0
        host = hosts[host_id]
        res = run_bowtie2((assembly, f, sam_outdir))
        if not host in bt_jobs:
            bt_jobs[host] = {'cmds':[], 'outs':[]}
        bt_jobs[host]['cmds'].append(res[2])
        bt_jobs[host]['outs'].append((res[1], f+".sam"))
        host_id += 1
    return bt_jobs


# %%
bt_jobs = create_bt_jobs(fastq_files)    

# %%
len(bt_jobs) == len(hosts)


# %%
def write_copy_cmds(key, jobs):
    s = ""
    for elem in jobs[key]['outs']:
        s += "cp -f {0} {1}\nrm {0}\n".format(elem[0], elem[1])
    return s

def write_jobs(key, jobs):
    with open("{}_jobs".format(key), "w") as o:
        for elem in jobs[key]['cmds']:
            o.write("{}\n".format(elem))
    
    script = """#!/bin/bash
    
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o {0}.out
#$ -e {0}.err
#$ -N bt
#$ -q *@{0}
cat {0}_jobs | parallel -j3 --progress
{1}
""".format(key, write_copy_cmds(key, jobs))
    
    with open("{}_job.sh".format(key), "w") as o:
        o.write("{}\n".format(script))
    return o.name


# %%
with open("runbowtie.sh", "w") as o:
for k in bt_jobs:
    o.write("qsub {}\n".format(write_jobs(k, bt_jobs)))

# %%
fqs = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/ -name "*.fastq"

# %%
sams = !find /gpfs_fs/home/eckertlab/projects/burt/seq/dedupe/ -name "*.sam"

# %%
len(fqs), len(sams)

# %%
missed = []
for f in fqs:
    if not os.path.exists(f + ".sam"):
        missed.append(f)

# %%
missed_jobs = create_bt_jobs(missed)

# %%
len(missed_jobs) == len(missed)

# %%
with open("runbowtie_missed.sh", "w") as o:
    for k in missed_jobs:
        o.write("qsub {}\n".format(write_jobs(k, missed_jobs)))

# %%
# for h in missed_jobs:
#     print(h)
#     !ssh $h pkill -9 bash
#     !ssh $h pkill -9 perl
#     !ssh $h pkill -9 bowtie2-align-l 

# %%
rc = Client(profile="sge")

# %%
dv, lv = u.get_views(rc)
len(dv)


# %%
@lv.remote()
def convert_sam_to_bam(sam):
    import stopwatch, multiprocessing, os
    timer = stopwatch.Timer()
    cpus = multiprocessing.cpu_count()
    bam = sam.replace(".sam", ".bam")
    bam_sorted = "%s_sorted.bam" % bam.replace(".bam", "")
    if not os.path.exists(bam):
        !/home/cfriedline/gpfs/src/samtools-1.3/samtools view -b $sam -o $bam
        !/home/cfriedline/gpfs/src/samtools-1.3/samtools sort -@ $cpus $bam -o $bam_sorted
        !/home/cfriedline/gpfs/src/samtools-1.3/samtools index $bam_sorted
    timer.stop()
    return bam, bam_sorted, timer.elapsed


# %%
len(fastq_files)

# %%
demult_dir = "/gpfs_fs/home/eckertlab/gypsy_indiv/raw_demult_gbsx"

# %%
pwd

# %%
sam_files = !find . -type f -name "*.sam"
sam_files = [os.path.abspath(x) for x in sam_files]
assert len(sam_files) == len(fastq_files)

# %%
sam_files

# %%
sam_bam_jobs = []
for f in sam_files:
    sam_bam_jobs.append(convert_sam_to_bam(os.path.abspath(f)))

# %%
u.get_async_progress(sam_bam_jobs)

# %%
sorted_bams = !find . -type f -name '*sorted.bam'
sorted_bams = [os.path.abspath(x) for x in sorted_bams if 'bam' in x]
assert len(sorted_bams) == len(fastq_files)


# %%
@lv.remote()
def get_lane_info(bam):
    res = !/home/cfriedline/g/src/samtools-1.3/samtools view $bam | tail -n1
    return bam, res[0].split("\t")[0]


# %%
rg_info = []
for f in sorted_bams:
    rg_info.append(get_lane_info(f))

# %%
u.get_async_progress(rg_info)

# %%
lane_info = [x.r for x in rg_info]

# %%
rg_dict = {}
for bam, header in lane_info:
    sample = os.path.basename(bam).split(".")[0]
    instr, run, flowcell, lane, tile, x, y = header.split(":")
    rg_dict[bam] = {"id": "{}.{}.{}".format(flowcell, lane, sample),
                   "pl": "ILLUMINA",
                   "lb": "{}.{}".format(flowcell, lane),
                   "sm": sample}

# %%
dv['rg_dict'] = rg_dict

# %%
hv = u.get_single_host_lview(rc, "all")

# %%
len(hv)


# %%
@hv.remote()
def add_rg_info_to_bam(bam):
    import os
    cmd = "java -jar /home/cfriedline/gpfs/src/picard-tools-1.112/AddOrReplaceReadGroups.jar"
    bam_rg = bam.replace(".bam", "_rg.bam")
    info = rg_dict[bam]
    rg_string = "RGID={0} RGLB={1} RGPL=ILLUMINA RGPU={1} RGSM={2}".format(info['id'],
                                                                           info['lb'],
                                                                           info['sm'])
    cmd = "{} INPUT={} OUTPUT={} {} CREATE_INDEX=true".format(cmd,
                                                              bam,
                                                              bam_rg,
                                                              rg_string)
#     if not os.path.exists(bam_rg):
    !$cmd
    return bam_rg, rg_string, cmd


# %%
add_rg = []
for f in sorted_bams:
    add_rg.append(add_rg_info_to_bam(f))

# %%
u.get_async_progress(add_rg)

# %%
rg_bams = !find . -name "*rg.bam"
rg_bams = sorted([os.path.abspath(x) for x in rg_bams if 'rg.bam' in x])
assert len(rg_bams) == len(fastq_files)

# %%
len(rg_bams)


# %%
def get_mapped(bam):
    import os
    out = "%s_mapped.bam" % bam.split(".")[0]
    if not os.path.exists(out):                                       
        cmd = "%s view -b -F 4 %s > %s" % (samtools, bam, out)
        res = !$cmd
    index_bam(out)
    return bam


# %%
def index_bam(bam):
    cmd = "%s index %s" % (samtools, bam)
    !$cmd
    return bam


# %%
samtools = "/home/cfriedline/g/src/samtools-1.3/samtools"
dv['samtools'] = samtools
dv['index_bam'] = index_bam
dv['get_mapped'] = get_mapped

# %%
jobs = []
for b in rg_bams:
    jobs.append(lv.apply_async(get_mapped, b))

# %%
u.get_async_progress(jobs)


# %%
def get_contigs(bam):
    contigs = set()
    cmd = "/home/cfriedline/g/src/samtools-1.3/samtools view %s" % bam
    sys.stderr.write("%s: %s\n" % (socket.gethostname(), cmd))
    p = Popen(cmd, stdout=PIPE, shell=True)
    for line in p.stdout:
        d = line.decode().split("\t")
        contigs.add(d[2])
    return contigs


# %%
dv['get_contigs'] = get_contigs

# %%
mapped = !ls *mapped.bam
mapped = [os.path.abspath(x) for x in mapped]

# %%
len(mapped)

# %%
with dv.sync_imports():
    import os, sys, socket
    from subprocess import Popen, PIPE

# %%
mapped[0]

# %%
jobs = []
for b in mapped:
    jobs.append(lv.apply_async(get_contigs, b))

# %%
u.get_async_progress(jobs)

# %%
contig_counts = Counter()
for j in jobs:
    for contig in j.r:
        contig_counts[contig] += 1

# %%
contig_counts_sav = "contig_counts.pkl"
pickle.dump(contig_counts, open(contig_counts_sav, "wb"), pickle.HIGHEST_PROTOCOL)

# %%
assembly = "/home/cfriedline/eckertlab/loblolly2/conitgs.fa"

# %%
count = !grep -c ">" $assembly

# %%
count = 2855700

# %%
len(contig_counts), len(contig_counts)/count

# %%
cd $root

# %%
mapped_fasta = os.path.join(root, "%s_mapped.fasta" % os.path.basename(assembly))

# %%
mapped_fasta

# %%
faSomeRecords = "/home/cfriedline/g/src/kentUtils/bin/faSomeRecords"

# %%
with open("all_contigs.txt", "w") as o:
    for c in contig_counts:
        o.write("{}\n".format(c))

# %%
## Run this in a terminal
"{} {} {} {}".format(faSomeRecords, assembly, "all_contigs.txt", mapped_fasta)

# %%
with open("contigs.bed", "w") as o:
    for rec in SeqIO.parse(mapped_fasta, "fasta"):
        o.write("%s\t%d\t%d\n" % (rec.name, 0, len(rec)))
