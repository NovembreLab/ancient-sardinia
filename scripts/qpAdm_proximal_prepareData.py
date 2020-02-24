# Proximal Modeling 
# Created 2019-06-14 by Chi-Chun

import os
import subprocess
import shutil
import numpy as np
import pandas as pd
import itertools as it


# change directory to repo
os.chdir('/project/jnovembre/jhmarcus/ancient-sardinia')


### Define outgroups
#     - a15: 15 outgroup population
#     - a15plus: Nuragic Sardinians are added as an additional outgroup

a15 = ["Mota", "Ust_Ishim", "Kostenki14", "GoyetQ116-1", 
       "Vestonice16", "MA1", "ElMiron", "Villabruna", 
       "EHG", "CHG", "Iran-N", "Natufian", "Levant_N", 
       "WHG", "Anatolia_N"]

a15plus = a15 + ["Sar-Nur"]


### Some Functions to prepare the files for the Analysis
# make a random seed here

input_folder = "output/eigenstrat/full/"

post_Nur_inds = ["MSR002", "MSR003",
                 "VIL004", 'VIL006', 'VIL007', 'VIL009', 'VIL010', "VIL011",
                 "AMC001", "AMC005", "AMC014",
                 "COR001", "COR002",
                 "SNN001", "SNN002", "SNN003", "SNN004"]
mod_Sar = ["Cag", "Car", "Cam", "Ori", "Ogl", "Nuo", "Sas", "Olb"]

def replace_inds(df, indlist=[], poplabel="default"):
    """Replace Population Labels in df for indlist with poplabel"""
    idcs =  df[0].isin(indlist)
    if not idcs.sum():
        print(f"Found {idcs.sum()} / {len(indlist)} individuals for {poplabel}")
    df.loc[idcs, 2] = poplabel
    return df
    
def ind_as_pop(ind_path=input_folder, ind="qpall.ind", seed= 1234):
    ''' set individual label as population label (See Harald's original notebook files for comparison)
    '''
    df = pd.read_csv((ind_path+ind), index_col=False, delim_whitespace=True, header=None)
    
    # set population labels as individual IDs 
    # Post Nuragic ancients
    for i in post_Nur_inds:
        df = replace_inds(df, [i], i)
    
    # Modern provinces
    np.random.seed(seed)
    for m in mod_Sar:
        ids = np.where(df[2] == m)[0]
        np.random.shuffle(ids)
        ids = ids[:5]
        new_ids = [m + str(i) for i in range(1,6)] # Create the new Labels
        assert(len(ids) == len(new_ids))
        df.loc[ids, 2] = new_ids
    
    # Save to a ind file
    fname = ind.split(".")[0]
    savepath = ind_path + fname + "_proximal.ind"
    df.to_csv(savepath, index=False, header=None, sep=" ")
    print("Saved files to: %s" % savepath)
    
    df.loc[:, 0] = pd.Series(range(df.shape[0]))
    fname = ind.split(".")[0]
    savepath = ind_path + fname + "_proximal_uniqID.ind"
    df.to_csv(savepath, index=False, header=None, sep=" ")
    print("Saved files to: %s" % savepath)
    
    
ind_as_pop()


def make_dir_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def print_time():
    print('\n')
    now = datetime.datetime.now()
    print(now.strftime("%Y-%m-%d %H:%M:%S"))


# In[6]:

make_dir_if_not_exists('output/qpAdm_rev/NuragicSource')
make_dir_if_not_exists('output/qpAdm_rev/NuragicOutgroup')


# ### Additional preparation for submitting jobs
#     - creating a dictionary for source population pools
#     - creating a proper subset

sources_dict = {'Sard.Nuragic': ['Sar-Nur'],
           'N.African': ['Morocco_LN', 'Morocco_EN', 'Guanche'],
           'E.Mediterranean': ['Levant-BA', 'Iran-CA', 'Canaanite',
                               'Anatolia-BA', 'Minoan-BA', 'Myc-BA'],
           'Beaker': ['Italy_North-Bk', 'CE-EBA', 'CE-LBA',
                      'Iberia-BA', 'Balkans-BA'],
           'Sicily.BronzeAge': ['Sicily-Bk']}
pools = sources_dict.keys()


# In[28]:
def get_relevant_population(mod_Sar = mod_Sar, post_Nur_inds = post_Nur_inds, sources_dict = sources_dict, outgroups = a15):
    ''' relevant populations: targets: anc+mod Sar, sources: sources_dict, outgroups: a15
    '''
    pops = []
    for p in pools:
        pops.extend(sources_dict[p])
    mod_Sar_pops = [p + str(i) for p in mod_Sar for i in range(1,6)]
    pops.extend(mod_Sar_pops)
    pops.extend(post_Nur_inds)
    pops.extend(outgroups)
    return list(pops)

relevant_pops = get_relevant_population()
with open('output/eigenstrat/subset/proximal_populations.txt', 'w') as file:
    for l in relevant_pops:
        file.write(l + '\n')

# Run the below manually in shell preferablly.
#os.system('sbatch /project/jnovembre/jhmarcus/ancient-sardinia/scripts/eigenstrat_subset_proximal.sbatch')

