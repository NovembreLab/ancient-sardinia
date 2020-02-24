'''
Created on October 19th, 2018
Python File for runs on the cluster
@author: Harald Ringbauer
'''
import allel
import h5py
import numpy as np
import pandas as pd
#import pysam
import os
import itertools as it
import pickle as pkl
import warnings


###################################
meta_path = "output/meta/meta_rev_final.csv"
snp_path = "output/indices_rev/idx_dict.pkl"
data_path = "output/h5_rev/mod_reich_sardinia_ancients_rev_mrg_dedup_3trm_anno.h5"

# Load some parameters
param_path = "data/meta/parameters/ancsard_params.csv"
df_prs = pd.read_csv(param_path)

#ancsard_ind = df_prs["n_ancsard"][0] # The highest ancsard index
#anc_i = df_prs["n_anc"][0]    # The highest ancient index
nr_snps = df_prs["snp_cutoff"][0]
#cutoff1, cutoff2  = df_prs.iloc[0, 4:].values # Cutoffs of the periods
#periods = list(df_prs)[3:]  # Name of the periods NEED TO WORK THAT IN


####################################
####################################

def saveSNP(data, snp_list, path="output/eigenstrat/original/all.snp"):
    """Saves the SNP file"""
    chrom = data["variants/CHROM"][:][snp_list].astype("int")
    pos = data["variants/POS"][:][snp_list].astype("int")
    snp_id = data["variants/ID"][:][snp_list]
    ref = data["variants/REF"][:][snp_list]
    alt = data["variants/ALT"][:][snp_list]
    p = chrom.shape[0]
    
    snp_df = pd.DataFrame({"snp_id": snp_id, "chrom": chrom, "cm": [0 for _ in range(p)], "pos": pos, "ref": ref, "alt": alt})
    snp_df[["snp_id", "chrom", "cm", "pos", "ref", "alt"]].to_csv(path, header=None, index=False, sep="\t")

def saveInd(iids, clsts, path="output/eigenstrat/original/all.ind"):
    """Save the Ind file
    iids: List of Individual IDs
    clsts: List of Inidivual Cluster Labels
    """
    assert(len(iids) == len(clsts))
    iid_df = pd.DataFrame({"iid": iids, "sex": ["F" for _ in range(len(iids))], "clst": clsts})
    iid_df[["iid", "sex", "clst"]].to_csv(path, header=None, index=False, sep="\t")
    
def saveGeno(G, path="output/eigenstrat/original/all.geno"):
    """Save the Geno file"""
    np.savetxt(path, G, fmt="%.0f", delimiter="")

def produce_eigenstrat(meta_df, data, snp_list, pops=[]):
    """Produce Eigenstrat File from meta_df, hdf5 and SNP-List
    meta_df: Meta Data Dataframe
    data: hdf5 file
    snp_list: list of SNPs to save from hdf5
    pops: List of populations in meta_df to save!"""
    
    
    assert(len(set(pops)) == len(pops)) # Sanity Check to confirm that every population occurs only once.


    ### Extract Populations:
    target_pops_idx = meta_df["clst"].isin(pops)
    clsts = meta_df[target_pops_idx]["clst"].tolist()
    iids = meta_df[target_pops_idx]["iid"].tolist()
    
    # Some Output
    print("Nr of Populations: %i" % len(pops))
    print("Nr of SNPs to extract: %i" % len(snp_list))
    print("Nr of extracted Inds: %i" % len(iids))
    
    print("Extracting Genotypes...")
    H = data["calldata/GT"][:, target_pops_idx, :][snp_list, :, :]
    D = np.sum(H, axis=-1) # Calculate the Genotype Matrix
    
    print("Correcting the Genotypes...")
    ### Rewrite the Entries in Eigenstrat Format:
    G = np.ones(np.shape(D)) * 9 # Empty Matrix with default missing
    G[D==2] = 0
    G[D==1] = 1
    G[D==0] = 2 
    
    print("Saving Files...")
    # Save the files  
    saveGeno(G)
    saveInd(iids, clsts)
    saveSNP(data, snp_list)
    
    print("Completed!")
    

####################################
####################################

def extract_eigenstrat():
	"""Simple Wrapper to call"""
	anc_pops = ['Anatolia-BA', 'Anatolia-CA', 'Anatolia-N', 'Balkans-BA','Balkans-EN', 'Balkans-IA',
         'Balkans-MNCA', 'CE-EBA', 'CE-EN', 'CE-LBA', 'CE-MNCA', 'EHG-HG', 'France-N', 'Greece-EN', 'Iberia-BA',
         'Iberia-ECA', 'Iberia-EN', 'Iberia-LCA',  'Iran-CA', 'Iran-LN', 'Iran-N', 'Iron_Gates-HG', 'Italy_North-Bk',
	 'Levant-BA', 'Levant-N', 'Minoan-BA', 'Myc-BA', 'Natufian-HG', 'Netherlands-BA','Poland-BA', 'SHG-HG',
         'Sicily-Bk', 'Steppe-EMBA','Steppe-EN', 'Steppe-IA', 'Steppe-MLBA', 'GB-EBA','GB-EN', 'GB-LBA',
         'GB-LN', 'WHG-HG', 'Ukraine-N', "Morocco_EN", 'Morocco_Iberomaurusian', 'Morocco_LN', 'Guanche', 'Canaanite']

	ho_pops = ["Ami", "Biaka", "Bougainville", "Chukchi", "Eskimo", "Han", 
		"Ju_hoan_North", "Karitiana", "Kharia", "Mbuti", "Onge", 
		"Papuan", "She", "Ulchi", "Yoruba", "Bergamo", "Libyan_Jew", 
                "Tunisian", "Saharawi", "Algerian", 'Egyptian', 'Mozabite', "Cypriot", "Turkish", "Lebanese", "Greek", "ItyS", 
                "Italian_South", "Maltese",  "Tuscan", "Sicilian", "Spanish", "Basque",
                "French", "Norwegian", "Czech", 'Tunisian_Jew', 'Turkish_Jew', "Moroccan_Jew", "BedouinA", 
                "Canary_Islanders", 'Druze', 'Palestinian', 'Jordanian']

	sard_pops = ["Olb", "Sas", "Nuo", "Ogl", "Ori", "Cam", "Car", "Cag"]

	### Load the Metafile
	meta_df = pd.read_csv(meta_path)  # Load The Meta file
	anc_sard_df = meta_df[meta_df["study"] == "Marcus et al. 2018"]
	anc_sard_df = anc_sard_df[(anc_sard_df["n_cov_snp"] > nr_snps) & (anc_sard_df["include_alt"] > 0)]
	ancsard_pops = list(anc_sard_df.clst.unique())

	pops = ancsard_pops + anc_pops + ho_pops + sard_pops

	### Load the SNPs
	idx_dict = pkl.load(open(snp_path, "rb"))
	#ho_snp_fil = np.where(idx_dict["ho"]["snp"])[0]
	ho_anc_fil = np.where(idx_dict["hq"]["snp"])[0]

	### Load the hdf5
	data = h5py.File(data_path, mode='r')

	###### Run the actual File
	print("Starting the run...")
	produce_eigenstrat(meta_df, data, snp_list=ho_anc_fil, pops=pops)

####################################
####################################
# Only run if called from this file:
if __name__ == "__main__":
	extract_eigenstrat()

