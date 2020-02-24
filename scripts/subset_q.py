import pandas as pd
import numpy as np
import click

@click.command()
@click.option('--qpath', help='path to admixture Q file')
@click.option('--metapath', help='path to meta file')
@click.option('--popspath', help='path to pops file')
@click.option('--qout', help='path to output admixture Q file')
def subset_q(qpath, metapath, popspath, qout):
    """
    """
    pops = np.loadtxt(popspath, dtype=str).tolist()
    meta_df = pd.read_table(metapath, sep=",")   
    idx = meta_df[meta_df["clst"].isin(pops)].index
    q_df = pd.read_table(qpath, sep=" ", header=None)
    q_sub_df = q_df.iloc[idx]
    q_sub_df.to_csv(qout, sep=" ", index=False, header=False)

if __name__ == '__main__':
    subset_q()
