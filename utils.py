import os
import pandas as pd


def save_df(dirname, fname, df):
    f = os.path.join(dirname, "%s.txt" % fname)
    df.to_csv(f,
              header=True,
              index=True,
              sep="\t")
    print("saved %s" % f)


def read_df(dirname, fname):
    f = os.path.join(dirname, "%s.txt" % fname)
    return pd.read_csv(f, sep="\t", index_col=0)
