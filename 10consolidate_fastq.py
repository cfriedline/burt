# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.0.3
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
dest_dir = "/gpfs_fs/home/eckertlab/BURT/"

# %%
first_runs = "/home/cfriedline/eckertlab/projects/burt/seq"

# %%
fastq_files1 = []
fastq_map1 = {}
for root, dirs, files in os.walk(first_runs):
    for f in files:
        p = Path(root, f)
        if "fastq.gz" in p.name and ".R1." in p.name:
            sample_name = p.name.split(".")[0]
            library = p.parent.name
            fastq_files1.append(dict(sample_name=sample_name,
                                    fastq_path=str(p),
                                    library=library))
            
            if sample_name not in fastq_map1:
                fastq_map1[sample_name] = []
            fastq_map1[sample_name].append(str(p))

# %%
fastq_data1 = []
for k, v in fastq_map1.items():
    fastq_data1.append(dict(sample_name=k, fastq_files=v, 
                            library=list(set([Path(x).parent.name for x in v]))[0],
                           processed_fastq=f"{Path(first_runs, 'dedupe', Path(v[0]).name)}"))

fastq_df1 = pd.DataFrame(fastq_data1)


# %%
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
fastq_df1.md5 = fastq_df1.processed_fastq.apply(md5)


# %%
def copy_file(s, d):
    shutil.copy(s, d)

pool = mp.Pool(20)
jobs = []
for f in fastq_files1:
    s = Path(f)
    d = Path(dest_dir, s.name)
    jobs.append(pool.apply_async(copy_file, (s, d)))
pool.close()

# %%
sum([x.ready() for x in jobs])

# %%
pool.join()

# %%
second_library = "/home/cfriedline/eckertlab/Novogene/burt/merged"

# these files were the result of a merge of the failed NARF libraries and the good novogene libraries using 
# a script that trevor wrote ~/eckertlab/Novogene/burt/merge_fastq.sh

# %%
fastq_files2 = glob.glob(f"{second_library}/*.fastq.gz")

# %%
len(fastq_files2)

# %%
map2 = {}
for root, dirs,files in os.walk(Path(second_library).parent):
    for f in files:
        if "undetermined" not in f:
            if f.endswith("fastq.gz"):
                p = Path(root, f)
                if p.name not in map2:
                    map2[p.name] = []
                if "BURT" in p.parent.name:
                    map2[p.name].append(p)

# %%
map2

# %%
