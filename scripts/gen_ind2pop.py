import pandas as pd
import numpy as np
import click

@click.command()
@click.option('--metapath', help='path to meta file')
@click.option('--popspath', help='path to pops file')
@click.option('--out', help='path to out')
def gen_ind2pop(metapath, popspath, out):
    """
    """
    pops = np.loadtxt(popspath, dtype=str).tolist()
    meta_df = pd.read_table(metapath, sep=",")
    idx = meta_df[meta_df["clst"].isin(pops)].index
    meta_df.iloc[idx][["clst"]].to_csv(out, sep="\t", index=False, header=False)

if __name__ == '__main__':
    gen_ind2pop()
