# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.1.5
#   kernelspec:
#     display_name: burt
#     language: python
#     name: burt
# ---

# %%
import glob
import os 
import multiprocessing as mp
from pathlib import Path
import shutil
import pandas as pd 
import hashlib

# %%
plate_map_dir = "/home/cfriedline/eckertlab/BURT/plate_maps/"

# %%
plate_files = glob.glob(f"{plate_map_dir}/*.xls")


# %%
def read_plate(f):
    df = pd.read_excel(f)
    df.index = df.iloc[:,0]
    df = df.drop(df.columns[0], axis=1)
    df.index.name = None
    df = df.reset_index().rename(columns={"index":"row"})
    df["source"] = Path(f).name
    return df.melt(id_vars=["source", "row"], var_name="col", value_name="sample_name")


# %%
all_plates = pd.concat([read_plate(x) for x in plate_files])

# %%
all_plates = all_plates.reset_index(drop=True)

# %%
all_plates.head()


# %%
def split_name(name):
    try:
        data = name.split("-")
        return pd.Series(dict(species=data[0], 
                             state=data[1], 
                             popn=data[2],
                             ind="-".join(data[3:]),
                             full_pop="-".join(data[0:3])))
    except:
        return pd.Series()
all_plates = all_plates.join(all_plates.sample_name.apply(split_name))

# %%
all_plates = all_plates[~all_plates.sample_name.isna()].copy()


# %%
def derive_library_name(name):
    if "layout" in name:
        return f'Burt{name.split("_")[1]}'
    return f'Burt{name.split(".")[0].replace("plate", "")}'
all_plates["library"] = all_plates.source.apply(derive_library_name)

# %%
sorted(all_plates.library.unique())

# %%
all_plates.head()

# %%
all_plates.to_csv(Path(plate_map_dir, "all_plates.txt"), sep="\t", index=False)

# %%
all_plates[all_plates.sample_name=="G-VA-1-15"]

# %%
all_plates.groupby("species").count()


# %%
def make_df(gz_list, seq_center):
    df = pd.DataFrame(gz_list, columns=["fastq"])
    df["sample_name"] = df.fastq.apply(lambda path: path.name.split(".R1")[0])
    sample_name_split = df.sample_name.str.split("-")
    df["species"] = [x[0] for x in sample_name_split]
    df["state"] = [x[1] for x in sample_name_split]
    df["ind"] = [x[3] for x in sample_name_split]
    df["popn"] = ["-".join(x[0:3]) for x in sample_name_split]
    df["spp_popn"] = ["-".join(x[1:3]) for x in sample_name_split]
    df["run"] = [x.parent.parent.name for x in df.fastq]
    if seq_center == "Novogene":
        df["run"] = "Novogene1"
    df["lib"] = [x.parent.name for x in df.fastq]
    df["seq_center"] = seq_center
    return df


# %%
round1_root = Path("/home/cfriedline/eckertlab/projects/burt/seq")

# %%
round1_gz = list(round1_root.glob("**/*-*.R1.*.gz"))

# %%
round1_gz = [x for x in round1_gz if x.parent.parent.name in ["160520", "160525"]]

# %%
df1 = make_df(round1_gz, "NARF")

# %%
df1.head()

# %%
round2_root = Path("/home/cfriedline/eckertlab/Novogene/burt/demult")

# %%
round2_gz = list(round2_root.glob("*BURT*/*-*.R1.*.gz"))

# %%
df2 = make_df(round2_gz, "Novogene")

# %%
df_all = pd.concat([df1, df2])


# %%
def md5(fname):
    res = !md5sum {fname}
    return res[0].split()[0]


# %%
def progress(l):
    ready = sum([x.ready() for x in l])
    return f"{ready}/{len(l)} ({ready*100/len(l)}%)"


# %%
jobs = []
pool = mp.Pool()
for f in df_all.fastq:
    jobs.append(pool.apply_async(md5, (f,)))
pool.close()
# pool.join()

# %%
progress(jobs)

# %%
df_all["md5"] = [x.get() for x in jobs]

# %%
pool.join()

# %%
df_all["RGPL"] = "ILLUMINA"
df_all["RGSM"] = df_all.sample_name
df_all["RGLB"] = df_all.lib
df_all["RGID"] = df_all.sample_name + "." + df_all.run + "." + df_all.lib

# %%
df_all.head()

# %%
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip


# %%
def get_PU(fastq_file):
    with gzip.open(fastq_file, "rt") as f:
        for name, seq, qual in FastqGeneralIterator(f):
            data = name.split(":")
            instrument = data[0]
            run_number = data[1]
            flowcell_id = data[2]
            lane = data[3]
            return [flowcell_id, lane]
flowcell_lane = df_all.fastq.apply(get_PU)

# %%
df_all["flowcell_lane"] = [".".join(x for x in y) for y in flowcell_lane]

# %%
df_all["RGPU"] = df_all.flowcell_lane + "." + df_all.sample_name


# %%
def get_read_num(fastq_file):
    count = 0
    try:
        with gzip.open(fastq_file, "rt") as f:
            for name, seq, qual in FastqGeneralIterator(f):
                count+=1
    except:
        return -1
    return count


# %%
pool = mp.Pool()
jobs = []
for f in df_all.fastq:
    jobs.append(pool.apply_async(get_read_num, (f,)))
pool.close()

# %%
progress(jobs)

# %%
df_all["num_reads"] = [x.get() for x in jobs]

# %%
df_all[df_all.num_reads==-1].fastq.values

# %% [markdown]
# ## Failed fastq files
# These files failed. Will uncompress and recompres as follows.
#
# ```
# array([PosixPath('/home/cfriedline/eckertlab/Novogene/burt/demult/BURT_10/G-VA-3-5.R1.fastq.gz'),
#        PosixPath('/home/cfriedline/eckertlab/Novogene/burt/demult/BURT_11/G-VA-3-14.R1.fastq.gz'),
#        PosixPath('/home/cfriedline/eckertlab/Novogene/burt/demult/BURT_12/G-VA-3-6.R1.fastq.gz'),
#        PosixPath('/home/cfriedline/eckertlab/Novogene/burt/demult/BURT_12/G-VA-3-10.R1.fastq.gz')],
#       dtype=object)
# ```
#
# ```
# zcat /home/cfriedline/eckertlab/Novogene/burt/demult/BURT_10/G-VA-3-5.R1.fastq.gz | gzip > /home/cfriedline/eckertlab/Novogene/burt/demult/BURT_10/G-VA-3-5.R1.fastq.gz_new.gz
#
# mv /home/cfriedline/eckertlab/Novogene/burt/demult/BURT_10/G-VA-3-5.R1.fastq.gz_new.gz /home/cfriedline/eckertlab/Novogene/burt/demult/BURT_10/G-VA-3-5.R1.fastq.gz
# ```
#       

# %%
for f in df_all[df_all.num_reads==-1].fastq.values:
    print(f)
    old = f
    bak = f"{f}.bak"
    new = f"{f}.temp"
    !cp {old} {bak}
    !zcat {old} | gzip > {new}
    !mv -f {new} {old}

# %%
fixed_counts = df_all[df_all.num_reads==-1].fastq.apply(get_read_num)

# %%
for f in fixed_counts.index:
    df_all.loc[f,"num_reads"] = fixed_counts[f]

# %%
df_all.to_csv("/gpfs_fs/home/eckertlab/BURT/seq/burt_fastq_data.txt", sep="\t")

# %%
df_all

# %%
# !cp /gpfs_fs/home/eckertlab/BURT/seq/burt_fastq_data.txt .

# %%
