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
import pandas as pd
from scipy.spatial.distance import hamming
import Levenshtein as lv
import pickle

# %%
cd /home/cfriedline/eckertlab/projects/burt/seq

# %%
barcodes = pd.read_excel("barcodes.xls")

# %%
plates = !ls plate*.xls


# %%
def convert_plate(plate_dict):
    d = {}
    for row, inner in plate_dict.items():
        for col, sample in inner.items():
            key = "{}{}".format(row, col)
            d[key] = sample
    return pd.DataFrame(d, index=["sample"]).T
test = convert_plate(p)


# %%
def write_barcode(df, out):
     df[['sample', 'well', 'barcode1', 'barcode2']].to_csv(out, sep="\t", index=False)


# %%
def write_gbsx_barcode(df, out):
    d = df.copy()
    d['none'] = ""
    d['gbsx_bc'] = d.barcode2.apply(lambda x: x.replace("CTCTTTCCCTACACGACGCTCTTCCGATCT", "").upper())
    d['gbsx_enz'] = "EcoRI"
    d[['sample', 'gbsx_bc', 'gbsx_enz']].to_csv(out, sep="\t", index=False, header=False)


# %%
for plate in plates:
    p = pd.read_excel(plate, index_col=0).T.to_dict()
    p = convert_plate(p)
    d = p.merge(barcodes, left_index=True, right_on="well")
    lib =plate.split("_")[1] 
    gbsx = "barcode_{}_gbsx.txt".format(lib)
    d = d.sort_values("well")
    write_gbsx_barcode(d, gbsx)

# %%
bcs = barcodes.barcode2.apply(lambda x: x.replace("CTCTTTCCCTACACGACGCTCTTCCGATCT", "")[:-1].upper())

# %%
bc_lens = {}
for b in bcs:
    l = len(b)
    if not l in bc_lens:
        bc_lens[l] = []
    bc_lens[l].append(b)

# %%
pickle.dump(bc_lens, open("bc_lens.pkl", "wb"))

# %%
