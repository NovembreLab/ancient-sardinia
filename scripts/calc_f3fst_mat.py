'''
Created on November 10th, 2018
Calculate F_ST and F_3 matrices.
Saves them to .csvs.
@author: Harald Ringbauer
'''

import allel
import h5py
import numpy as np
import pandas as pd
import os
import itertools as it
import _pickle as pkl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from collections import Counter # Do do count of arrays

#h5path = "output/h5/mod_reich_sardinia_ancients_mrg_dedup_3trm_anno.h5" 
h5path = "output/h5_rev/mod_reich_sardinia_ancients_rev_mrg_dedup_3trm_anno.h5"
#meta_path = "output/meta/meta_rev_final.csv"    # The original Meta Path
meta_path = "output/meta/meta_rev_final_pn_sites.csv"

meta_df = pd.read_csv(meta_path)
print(f"Using Meta from {meta_path}")
#idx_dict = pkl.load(open("output/indices/idx_dict.pkl", "rb"))
idx_dict = pkl.load(open("output/indices_rev/idx_dict.pkl", "rb"))

foldercsv = "output/data-csvs/"  # Folder to save csvs
figf3_path = "output/figures/Supplemental/fstf3/pw_outgroup-f3.pdf"
figfst_path = "output/figures/Supplemental/fstf3/pw_fst.pdf"

### populations that are extracted:
comps = ['Iberia-EN', "France-N", "GB-EN", "CE-EN", "Balkans-EN", "Anatolia-N", "Steppe-EMBA", "Steppe-MLBA", "EHG-HG", 
         "Iron_Gates-HG", "WHG-HG", "CE-MNCA", 'Iberia-ECA', 'Iberia-LCA', 'Iberia-BA', 'GB-LN', 'GB-EBA', 'CE-EBA', "Iran-CA", "Minoan-BA", 
         "Myc-BA", "Canaanite", "Levant-BA", "Morocco_EN", "Morocco_LN", "French", "Basque", "Spanish",  "Tuscan", "Maltese", "Sicilian", "Cypriot", "Greek", "Libyan_Jew", "Tunisian_Jew", "Turkish_Jew", "Moroccan_Jew",
	'Druze', 'Palestinian', 'Jordanian', "Tunisian", "Saharawi", "BedouinA", "Turkish", "Lebanese", "Ogl"]

comps_sard =  ["Cag", "Car", "Cam", "Ori", "Nuo", "Sas", "Olb"]  # Provence labels that get grouped
ancsard_comps = ["Sar-MN", "Sar-ECA", "Sar-EMBA", "Sar-Nur", "Sar-VIL", "Sar-MSR", "Sar-AMC", "Sar-SNN"]   # Labels of ancient Sardinia in Plots. Also for searching their.

for c in comps:
	l=np.sum(meta_df["clst"]==c)
	print(f"{c} : {l} individuals")

### Load parameters
param_path = "data/meta/parameters/ancsard_params.csv"
df_prs = pd.read_csv(param_path)

ancsard_ind = df_prs["n_ancsard"][0] # The highest ancsard index
anc_i = df_prs["n_anc"][0]    # The highest ancient index
nr_snps = df_prs["snp_cutoff"][0]


### Do the cutoffs now by naming the periods
#cutoff1, cutoff2, cutoff3  = df_prs.iloc[0, 4:].values # Cutoffs of the periods
#periods = list(df_prs)[3:]  # Name of the periods NEED TO WORK THAT IN


#############################################
#############################################
### Prepare the objects:
## Load HDF5
f = h5py.File(h5path, "r")

snp_fil = np.where(idx_dict["ho"]["snp"] * idx_dict["sard"]["snp"] * idx_dict["hq"]["snp"])[0]
snp_fil_nog = np.where(idx_dict["ho"]["snp"] * idx_dict["sard"]["snp"] * idx_dict["hq"]["snp"] * idx_dict["nog"]["snp"])[0]

###### Output:
print("\nNr all SNPs: %i" % np.shape(f["calldata/GT"])[0])
print("Filtered SNPs: %i" % len(snp_fil))
print("Nr of filtered no g: %i" % len(snp_fil_nog))

assert(np.shape(f["calldata/GT"])[1] == len(meta_df))  # Sanity Check

######
# Prepare the Ancestral Allele File
# Do the pre-processing 
aa=f["variants"]["AA"][:]
aa_c = np.array([s.strip('|').upper() for s in aa]) # Remove all unneccessary elements, and transform to uppercase
print("\nAncestral Alleles:")
print(Counter(aa_c))  # Ref is 0; Alt is 1.
assert(np.shape(f["calldata"]["GT"])[0]==len(aa_c))  # Sanity Check
p_ref = np.zeros((len(aa_c), 2), dtype="int8") # Create empty allele Frequency vector
p_ref[:,0] = 1 # Set everything to Ref.

#set(alt=f["variants"]["ALT"][:]) 
alts = (f["variants"]["ALT"][:]==aa_c)  # Extract the indices where ALT is ancestral
print("Nr of Alleles where ALT is Ancestral: %i" % np.sum(alts))
p_ref[alts,0] = 0
p_ref[alts,1] = 1 


#############################################
#############################################
### Some helper functions:
def calc_all_freqs(gts, ls, snps_okay=[]):
    """Calcualtes the mean of sublist. NOT USED. LEGACY"""
    if len(snps_okay)==0:
        snps_okay=range(np.shape(gts)[0]) # Slice out all SNPs
    pt=gts[:,ls,:][snps_okay,:,:] # Slice out the right stuff
    pt=np.ma.masked_less(pt, 0, copy=False) # Mask out invalid snps with values <0!
    pt=np.mean(pt, axis=1)
    pt=np.mean(pt, axis=1)
    return pt

def calc_all_counts(gts, ls, snps_okay=[]):
    """Calcualtes the mean of sublist.
    anc_i: Up until which index to divide by two for pseudo haplotypes"""
    if len(snps_okay)==0:
        snps_okay=range(np.shape(gts)[0]) # Slice out all SNPs
    pt=gts[:,ls,:][snps_okay,:,:] # Slice out the right stuff
    
    # Apply the correction factor for ancestral individuals (Pseudo Haplotypes)
    anc_ind = np.array(ls) < anc_i
    pt[:, anc_ind, 1] = -1 # Set the extra allele as missing.
    
    
    # Calculate the Allele Counts.
    # I.e. nx2 Vector of the counts
    c_ref=np.sum(pt==0, axis=1) # Sum the Counts over all Individuals
    c_ref=np.sum(c_ref, axis=1) # Sum over both haplotypes
    
    c_alt=np.sum(pt==1, axis=1) # Sum the Counts over all Individuals
    c_alt=np.sum(c_alt, axis=1) # Sum over both haplotypes
    
    # Double 0,0 no problem, goes to NaN and is then caught by allel
    return np.column_stack((c_ref, c_alt)) # Return the nx2 Allele COunts
    
def f3_l(lt, ls1, ls2, snps_okay=None, gts=None):
    """Calculate f3 for Individuals in Target and List1 and List2
    lt, ls1, ls2: Set or Lists of Individuals indices
    snps_okay: Which SNPs to actually use. if none all
    gts: The Full Genotype Matrix. If none, take default
    """
    # Load Genotypes:
    if gts==None:
        gts=f["calldata/GT"]
        
    # Calculate the Mean:
    pt=calc_all_counts(gts, lt, snps_okay)
    p1=calc_all_counts(gts, ls1, snps_okay)
    p2=calc_all_counts(gts, ls2, snps_okay)
    
    #print("Finished Mean Allele Frequency Calculation!")
    #f3 = np.mean((pt-p1)*(pt-p2))
    f3 = allel.average_patterson_f3(pt, p1, p2, blen=1000, normed=False)
    return [f3[0],f3[1],f3[2]]  # f4, se, z

def f3_freqs(pt, p1, p2, snps_okay=None, gts=None, blen=1000):
    """Calculate f3 for Allele Counts (lx2 arrays)
    snps_okay: Which SNPs to actually use. If none use all
    gts: The Full Genotype Matrix. If none, take default
    blen: Block Nr for Bootstrap
    """
    # Load Genotypes:
    if gts==None:
        gts=f["calldata/GT"]
    #f3 = np.mean((pt-p1)*(pt-p2))
    f3 = allel.average_patterson_f3(pt, p1, p2, blen=blen, normed=False)
    return [f3[0], f3[1], f3[2]]  # f4, se, z

def f4_freqs(p1, p2, p3, p4, snps_okay=None, gts=None, blen=1000):
    """Calculate f4 for Allele Counts (lx2 arrays)
    snps_okay: Which SNPs to actually use. if none use all
    gts: The Full Genotype Matrix. If none, take default
    blen: Block Nr for Bootstrap
    """
    # Load Genotypes:
    if gts==None:
        gts=f["calldata/GT"]
    f4 = allel.average_patterson_d(p1, p2, p3, p4, blen=blen)
    return [f4[0], f4[1], f4[2]]  # f4, se, z

def fst_counts(p1, p2, blen=1000):
    """Calculate f3 for Allele Counts (lx2 arrays)
    blen: Block Nr for Bootstrap
    A sim wrapper, so later on different methods can be implemented.
    Return fst, se, z value (based on jackkniving)
    """
    res = allel.average_patterson_fst(p1, p2, blen=blen)
    f4, se = res[0], res[1]
    z = f4 / se # Calculate the z-Value
    return [res[0], res[1], z]  # f4, se, z

def fst_counts_hudson(p1, p2, blen=1000):
    """Calculate f3 for Allele Counts (lx2 arrays)
    blen: Block Nr for Bootstrap
    A sim wrapper, so later on different methods can be implemented.
    Return fst, se, z value (based on jackkniving)
    """
    res = allel.average_hudson_fst(p1, p2, blen=blen)
    f4, se = res[0], res[1]
    z = f4 / se # Calculate the z-Value
    return [res[0], res[1], z]  # f4, se, z

def plot_2d(pops, f, ste, title="Test", vrange=[], fsl=4, mpl=100, full=False, show=False, savename="f4.pdf", reverse=False): #fsl8
    """Create 2D colored plot of all combinations of pops.
    pops: Populations, f, ste Values and errors
    title: What to take as title.
    v_range: What range for the color values. If none, choose automatically
    fsl: Fontsize
    mpl: Multiplicator of f value for better visualition
    full: Whether to show the full Matrix
    show: Whether to plot the picture
    reverse: Whether to reverse color bar"""

    if len(vrange)>0:
        vmin, vmax = vrange[0], vrange[1]
    else:
        vmin, vmax = np.min(f), np.max(f)
    k = len(pops)
    print("Nr. of populations: %i" % k)
    
    #il1 = np.tril_indices(k, k=-1)  # Indices of lower triangular matrix
    il1=np.array(list(it.combinations(np.arange(k),2)))
    il1 = (il1[:,0], il1[:,1])   # Split for indexing
    assert(len(f) == (len(pops) * (len(pops) - 1) / 2) and (len(f)==len(ste)))
    
    f_mat, ste_mat = np.zeros((k,k)), np.zeros((k,k))     # Transform f, ste into matrices
    
    # Set the values
    f_mat[il1] = f
    f_mat[(il1[1], il1[0])] = f # Set the other triangular matrix 
    
    if full==False:
        mask = np.logical_not(np.tri(k, k=-1)).T  # Mask out triangular matrix
        f_mat = np.ma.array(f_mat, mask=mask)  # mask out the lower triangle
        
    elif full==True:
        mask = np.diag(np.ones(k))
        f_mat = np.ma.array(f_mat, mask=mask)   # Mask out the Diagonal
        f = f_mat.flatten().compressed() # So that list of entries for plotting values
        
    ste_mat[il1] = ste
    
    fs=24
    plt.figure(figsize=(16,14))
    plt.axis("equal")
    cmap = "viridis"
    if reverse == True:
       cmap = "viridis_r"

    c = plt.pcolor(f_mat, vmin=vmin, vmax=vmax, alpha=0.9, cmap=cmap)  # cmap RdBu
    ax = plt.gca()  # Get the axis object
    plt.xticks(np.arange(k) + 0.5, pops, rotation='vertical', fontsize=12)
    plt.yticks(np.arange(k) + 0.5, pops, rotation='horizontal', fontsize=12)
    
    def show_values(pc, fmt=" %.2f"):    # Comment out %.3f \n" + r"$\pm$" + r"%.4f"
        pc.update_scalarmappable()
        #ax = pc.get_axes()
        for p, color, f0 in zip(pc.get_paths(), pc.get_facecolors(), f):  # pc.get_array()
            x, y = p.vertices[:-2, :].mean(0)
            if np.all(color[1] > 0.5):
                color = (0.0, 0.0, 0.0)
            else:
                color = (1.0, 1.0, 1.0)
            ax.text(x, y, fmt % (f0 * mpl), ha="center", va="center", color=color, fontsize=fsl)

    show_values(c)
    plt.colorbar()
    plt.title(title, fontsize=fs)
    
    if full==True: 
        low=0.0    
    else: low=1.0
        
    plt.xlim([low, k])
    plt.ylim([0.0,k])
    plt.savefig(savename, bbox_inches = 'tight', pad_inches = 0)
    if show:
        plt.show()


###################################################################
###################################################################
def f3_mat():
	"""Calculate and save f3 Matrix"""

	gts=f["calldata/GT"] # Where to find the raw data
	snps_okay = snp_fil # With G Sites. Without G Sites: snp_fil_nog!!
	print("Nr of SNPS: %i" % len(snps_okay))

	# First calculate Allele Frequency of Outgroup: 
	t_label="Ancestral"  # With wich population to compare  Yoruba
	print("Outgroup: %s" % t_label)

	pt = int(1e6) * p_ref[snps_okay,:] # 100 Ancestral Individuals

	#### Get list of indices for all of the single populations 
	inds = [np.where((meta_df.clst == clst) & (meta_df["include_alt"] == 1))[0] for clst in comps]

	# Get Ancestral Individuals for first position
	anc_sard_meta = meta_df.iloc[:ancsard_ind]
	anc_sard_meta = anc_sard_meta[(anc_sard_meta["include_alt"] == 1) & (anc_sard_meta["n_cov_snp"] >= nr_snps)]

	#as_med_samples = anc_sard_meta[(anc_sard_meta["age"] < cutoff1)].index.tolist()
	#as_late_ba_samples =  anc_sard_meta[(anc_sard_meta["age"] >= cutoff1) & (anc_sard_meta["age"] < cutoff2)].index.tolist()
	#as_middle_ba_samples = anc_sard_meta[(anc_sard_meta["age"] >= cutoff2) & (anc_sard_meta["age"] < cutoff3)].index.tolist()
	#as_neo_samples = anc_sard_meta[(anc_sard_meta["age"] >= cutoff3)].index.tolist()

	anc_inds = [np.where((meta_df.clst == clst) & (meta_df["include_alt"] == 1))[0] for clst in ancsard_comps]

	# Get Modern Sardinian Inidividuals for last position
	mod_inds = [np.where(meta_df.clst.isin(comps_sard))[0]]

	inds = anc_inds + inds + mod_inds
	fcomps = ancsard_comps + comps + ["non_Ogl"]  # Append the name as well

	#### Calculate all allele frequencies
	#print(list(map(len, inds)))
	#print(fcomps)
	assert(np.min(list(map(len, inds)))>0)
	assert(len(inds)==len(fcomps)) # Sanity Check

	ps = []
	for i in range(len(fcomps)):
	    print("Calculating AF for: %s. Sample Size: %i" % (fcomps[i], len(inds[i])))
	    p = calc_all_counts(gts, inds[i], snps_okay) # Calculate all allele Frequencies
	    ps.append(p) 
	    
	print("Finished pre-processing all freqs. Starting f_st calculation!")

	########## Calculate all f3s
	res = []  # Vector for the Results
	pops1, pops2 = [], [] # Vector for the populations pairs

	k = len(fcomps)
	for i1, i2 in it.combinations(np.arange(k),2):
	    # Save for the population list
	    pops1.append(fcomps[i1])
	    pops2.append(fcomps[i2])    
	    
	    # Do the f3 Calculation
	    res_new = f3_freqs(pt, ps[i1], ps[i2], snps_okay=snps_okay)
	    res.append(res_new)
	    
	print("Finished calculation!\n")


	### The plotting part:
	res = np.array(res)
	t_labels = [t_label for _ in res] # Create the Outgroup Label

	df4 = pd.DataFrame({"t": t_labels, "s1": pops1, "s2": pops2,
		            "f4": res[:,0], "se": res[:,1], 
		            "z": res[:,2]})

	### Sort to have things on top. For printing if needed
	df4_top = df4.sort_values(by=['f4'], ascending=False)

	########### Do the Plot
	plot_2d(fcomps, df4.f4, df4.se, title="F3(Ancestral-Pop x; Ancestral-Pop y)", vrange=[0.217, 0.225], fsl=3.6, mpl=100, full=True, savename=figf3_path)

	########### Do the saving
	### Save the f3-Matrix as well as the Countries:
	df4.to_csv(foldercsv + 'pwf3.csv', index=False) # Save the Values
	np.savetxt(foldercsv + 'pops.csv', fcomps, fmt="%s")  # Save the Population

	#### save relevent matricies ####
	k = len(fcomps)
	    
	il1 = np.array(list(it.combinations(np.arange(k),2)))
	il1 = (il1[:,0], il1[:,1])   # Split for indexing    
	f_mat, ste_mat, z_mat = np.zeros((k,k)), np.zeros((k,k)),  np.zeros((k,k)) # Transform f, ste, z-scores into matrices
	    
	# Set the values
	f_mat[il1] = df4.f4
	f_mat[(il1[1],il1[0])] = df4.f4 # Set the other triangular matrix 
	np.savetxt(foldercsv + "pwf3_re.mat", f_mat)

	ste_mat[il1] = df4.se
	ste_mat[(il1[1],il1[0])] = df4.se # Set the other triangular matrix 
	np.savetxt(foldercsv + "pwf3_se.mat", ste_mat)

	z_mat[il1] = df4.z
	z_mat[(il1[1], il1[0])] = df4.z # Set the other triangular matrix 
	np.savetxt(foldercsv + "pwf3_z.mat", z_mat)

	print("Saving to %s complete!\n" % foldercsv)


###################################################################
###################################################################


def fst_mat():
	"""Calculate and save fst Matrix"""
	gts=f["calldata/GT"] # Where to find the raw data
	snps_okay = snp_fil # With G Sites. Without G Sites: snp_fil_nog!!
	print("Nr of SNPS: %i" % len(snps_okay))


	#### Get list of indices for all of the single populations 
	inds = [np.where((meta_df.clst == clst) & (meta_df["include_alt"] == 1))[0] for clst in comps]

	# Get Ancestral Individuals for first position
	anc_sard_meta = meta_df.iloc[:ancsard_ind]
	anc_sard_meta = anc_sard_meta[(anc_sard_meta["include_alt"] == 1) & (anc_sard_meta["n_cov_snp"] >= nr_snps)]

	#as_med_samples = anc_sard_meta[(anc_sard_meta["age"] < cutoff1)].index.tolist()
	#as_late_ba_samples =  anc_sard_meta[(anc_sard_meta["age"] >= cutoff1) & (anc_sard_meta["age"] < cutoff2)].index.tolist()
	#as_middle_ba_samples = anc_sard_meta[(anc_sard_meta["age"] >= cutoff2) & (anc_sard_meta["age"] < cutoff3)].index.tolist()
	#as_neo_samples = anc_sard_meta[(anc_sard_meta["age"] >= cutoff3)].index.tolist()

	anc_inds = [np.where((meta_df.clst == clst) & (meta_df["include_alt"] == 1))[0] for clst in ancsard_comps]

	# Get Modern Sardinian Inidividuals for last position
	mod_inds = [np.where(meta_df.clst.isin(comps_sard))[0]]

	inds = anc_inds + inds + mod_inds
	fcomps = ancsard_comps + comps + ["non_Ogl"]  # Append the name as well

	##############################################################
	assert(np.min(list(map(len,inds)))>0) # Sanity Check if all Populations have Individuals
	assert(len(inds)==len(fcomps)) # Sanity Check

	ps = []  # Vector for allele Frequencies
	for i in range(len(fcomps)):
	    print("Calculating AF for: %s. Sample Size: %i" % (fcomps[i], len(inds[i])))
	    p = calc_all_counts(gts, inds[i], snps_okay) # Calculate all allele Frequencies
	    ps.append(p) 
	    
	print("Finished pre-processing all freqs. Starting f3 calculation!")

	# Calculate all f3s
	res, res_h = [], []  # Vector for the Results
	pops1, pops2 = [], [] # Vector for the populations pairs

	k = len(fcomps)
	for i1, i2 in it.combinations(np.arange(k),2):
	    # Save for the population list
	    pops1.append(fcomps[i1])
	    pops2.append(fcomps[i2])    
	    ### Do the fst Calculation
	    res_new = fst_counts(ps[i1], ps[i2])
	    res.append(res_new)
	    #res_new_h = fst_counts_hudson(ps[i1], ps[i2]) # Do the Hudson Estimator 
	    #res_h.append(res_new_h)
	    
	print("Finished calculation!")
	
	##############################################################

	res= np.array(res)
	#res = np.array(res_h) # Do the Hudson

	df4 = pd.DataFrame({"s1": pops1, "s2": pops2,
		       "f4": res[:,0], "se": res[:,1], "z": res[:,2]})

	### Sort to have things on top. For printing if needed
	df4_top = df4.sort_values(by=['f4'], ascending=False)

	### Do the Plot
	plot_2d(fcomps, df4.f4, df4.se, title="Patterson F_ST x100 (Pop x; Pop y)", vrange=[0.001, 0.075], fsl=4.6, mpl=100, full=True, savename=figfst_path, reverse=True)

	##############################################################
	# Save the Figure

	df4.to_csv(foldercsv + 'pw_fst.csv', index=False) # Save the Values
	np.savetxt(foldercsv + 'pops.csv', fcomps, fmt="%s")  # Save the Population

	#### save relevent matricies ####

	k = len(fcomps)
	    
	il1 = np.array(list(it.combinations(np.arange(k),2)))
	il1 = (il1[:,0], il1[:,1])   # Split for indexing    
	f_mat, ste_mat, z_mat = np.zeros((k,k)), np.zeros((k,k)),  np.zeros((k,k)) # Transform f, ste, z-scores into matrices
	    
	# Set the values
	f_mat[il1] = df4.f4
	f_mat[(il1[1],il1[0])] = df4.f4 # Set the other triangular matrix 
	np.savetxt(foldercsv + "pw_fst_re.mat", f_mat)

	ste_mat[il1] = df4.se
	ste_mat[(il1[1],il1[0])] = df4.se # Set the other triangular matrix 
	np.savetxt(foldercsv + "pw_fst_se.mat", ste_mat)

	z_mat[il1] = df4.z
	z_mat[(il1[1], il1[0])] = df4.z # Set the other triangular matrix 
	np.savetxt(foldercsv + "pw_fst_z.mat", z_mat)

	print("Saving to %s complete!\n" % foldercsv)


####################################
####################################
# Only run if called from this file:
if __name__ == "__main__":
	f3_mat()    # Run and save f3 Matrix
	fst_mat()   # Run and save fst Matrix

