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
import os
import time

# %%
qhost = !qhost | grep godel

# %%
procs_per_host = 4

# %%
hosts = []
cpus = {}
for q in qhost:
    q = q.split()
    mem = float(q[4][:-1])
    host = q[0]
    if mem > 60:
        hosts.append(host)
        cpus[host] = int(q[2])


# %%
def write_numfile(f, n):
    with open(f, "w") as o:
        for i in range(n):
            o.write("{}\n".format(i))
            
            
def write_qsub_file(host, cpu):
    cpu= int(cpu*.7)
    if cpu > 20:
        cpu = 20
        
    cpu = 10
    numfile = "/home/cfriedline/ipython/num{}".format(cpu)
    if not os.path.exists(numfile):
        write_numfile(numfile, cpu)
        
    cmd = """#!/bin/bash
#$ -cwd
#$ -V
#$ -N eng
#$ -q *@{0}
#$ -pe smp 1
#$ -e {0}.err
#$ -o {0}.out
#$ -S /bin/bash
source /home/cfriedline/anaconda/bin/activate py34
/home/cfriedline/bin/parallel -j {1} --delay 2 ipengine --profile=sge :::: /home/cfriedline/ipython/num{1}
""".format(host, cpu)
    
    with open("qsub_{}.sh".format(host), "w") as o:
        o.write("{}\n".format(cmd))
    !chmod +x {o.name}
    return o.name

scripts = []
# !rm qsub*.sh
# !rm godel*.err
# !rm godel*.out
for h in hosts:
    if h != 'godel96'
    scripts.append(write_qsub_file(h, cpus[h]))

# %%
len(scripts)

# %%
for s in scripts:
    !qsub $s
    time.sleep(0.5)

# %%
