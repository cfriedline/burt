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
import os

# %%
counts = pd.read_csv("counts.txt", sep="\t")

# %%
counts.sort_values("reads").head(54-15).to_csv("reseq.txt", sep="\t", header=True, index=False)

# %%
counts[counts['sample_name'] == 'T-AL-1-24']

# %%
counts.sort_values("reads").head(60).to_csv("reseq_60.txt", sep="\t", header=True, index=False)

# %%
