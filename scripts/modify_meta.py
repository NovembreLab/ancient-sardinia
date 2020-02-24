'''
Created on November 10th, 2018
Python file to modify the clsts of the meta_file
@author: Harald Ringbauer
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  # For Plotting
import json 

# For identifying countries
import requests
from shapely.geometry import mapping, shape
from shapely.prepared import prep
from shapely.geometry import Point


#meta_path = "output/meta/meta.csv"
#meta_path_new =  "output/meta/meta_final.csv"

meta_path = "output/meta/meta_rev.csv"
meta_path_new =  "output/meta/meta_rev_final.csv"


json_path = "data/geojson/countries.geojson"
cont_path = "data/meta/parameters/exclude_samples.csv"
rel_path = "data/meta/parameters/samples_relatives.csv"
param_path = "data/meta/parameters/ancsard_params.csv"

################################### Load the Parameters
df_prs = pd.read_csv(param_path)

anc_sard_ind = df_prs["n_ancsard"][0] # The highest ancsard index
anc_ind = df_prs["n_anc"][0]    # The highest ancient index
snp_cutoff = df_prs["snp_cutoff"][0]
period_cutoffs = df_prs.iloc[0,3:].values # Cutoffs of the periods
periods = list(df_prs)[3:]  # Name of the periods

###
contaminated = np.loadtxt(cont_path, dtype=np.str) # Load list of contaminated Individuals
related = np.loadtxt(rel_path, dtype=np.str) # Load list of related individuals (one of each pair)

####################################
####################################


def modifiy_meta():
	"""Reads out meta_path, modify it.
	Save it to meta_path_new"""
	
	full_df = pd.read_csv(meta_path) # The original one
	meta_df = pd.read_csv(meta_path)
	meta_df = meta_df[:anc_ind] # Only deal with ancient samples
	clst_new = meta_df.clst.values	

	### Define Countriee
	with open(json_path) as f:
	    data = json.load(f)
	#data = requests.get("../../data/geojson/countries.geojson").json()
	    
	cts = {}
	for feature in data["features"]:
	    geom = feature["geometry"]
	    country = feature["properties"]["ADMIN"]
	    cts[country] = prep(shape(geom))

	print("Countries loaded: %i" % len(cts))

	def get_country(lon, lat):
	    """Return country for given lon and lat"""
	    point = Point(float(lon), float(lat))
	    for country, geom in cts.items():
	        if geom.contains(point):
	            return country

	    return "unknown"

	print(get_country(10.0, 47.0)) # Sanity Test. 


	### Assign each sample to it's location
	lat, lon = meta_df.lat.values, meta_df.lon.values

	countries = [get_country(lon[i], lat[i]) for i in range(anc_ind)]
	meta_df["Country"]= countries

	# Set the island samples to their Nations:
	for s in ["Scotland-N", "Scotland-EBA", "Wales-MBA", "England-EBA", "Scotland-MBA", "Britain-Bk",
		"UK-LN", "UK-LBA", "UK-EBA", "UK-EN"]:
	    meta_df.loc[meta_df.clst==s, "Country"] = "United Kingdom"

	meta_df.loc[meta_df.iid=="I5071", "Country"]="Croatia"
	meta_df.loc[meta_df.iid=="I5072", "Country"]="Croatia"


	### Turkey. 30 samples. Problem is 1 chalcolithic sample (5775 y old). 26 Neolithic. 3 BA
	# Leave

	def set_labels(df, country, label_new, clst_new, age=[]):
	    """Set Cluster in clst_new (array) within age (if given)
	    to label_new (string). Age: Array of length 2 """
	    assert(len(df)==len(clst_new)) # Sanity Check
	    
	    if isinstance(country, list):
	        inds = (df["age"] > age[0]) & (df["age"]< age[1]) & (df["Country"].isin(country))
	    else:
	        inds = (df["age"] > age[0]) & (df["age"]< age[1]) & (df["Country"]==country)
		
	    print("Reset to %s: %i" % (label_new, np.sum(inds)))
	    clst_new[inds] = label_new
	    return clst_new
	    
	    
	### Austria 
	#clst_new[meta_df.clst=="LBK_Austria-N"] = "Germany-EN" # Rename the Austria samples

	### Balkans-Countries:
	balkan_cts = ["Republic of Serbia", "Bulgaria", "Romania", "Macedonia", "Hungary", "Croatia"]
	set_labels(meta_df, balkan_cts, "Balkans-IA", clst_new, age=[1500,3500]) # One sample only!!
	set_labels(meta_df, balkan_cts, "Balkans-BA", clst_new, age=[4000,5500])
	set_labels(meta_df, balkan_cts, "Balkans-MNCA", clst_new, age=[5500,7600])
	#set_labels(meta_df, balkan_cts, "Balkans-LN", clst_new, age=[7000,7600])
	set_labels(meta_df, balkan_cts, "Balkans-EN", clst_new, age=[7599,8000]) # Balkan HG are reset below
	clst_new[meta_df.clst=="Trypillia-CA"] = "Balkans-MNCA" # north Rumanian/south Ukrainian culture

	### France
	set_labels(meta_df, "France", "France-N", clst_new, age=[5000,7000])
	set_labels(meta_df, "France", "CE-EBA", clst_new, age=[4100,5001]) # Only 1 low cov. sample <6000
	set_labels(meta_df, "France", "CE-LBA", clst_new, age=[3001,4101])        
	#clst_new[meta_df.clst=="Southern_France-Bk"] = "S_France-BA" # Rename the South French samples

	### Germany
	german_cts = ["Germany", "Czech Republic", "Austria",  "Slovakia", "Switzerland", "Poland"]
	set_labels(meta_df, german_cts, "CE-EN", clst_new, age=[6700,7400])
	set_labels(meta_df, german_cts, "CE-MNCA", clst_new, age=[4800,6700]) # Only one sample >5800
	#set_labels(meta_df, german_cts, "CE-CA", clst_new, age=[4800,5300]) # Only one sample >5800
	set_labels(meta_df, german_cts, "CE-EBA", clst_new, age=[4100,4800])
	set_labels(meta_df, german_cts, "CE-LBA", clst_new, age=[3000,4100])

	### Poland 
	set_labels(meta_df, "Poland", "Ignore--9", clst_new, age=[4500,10000]) # All faulty samples

	### Greece Rename 5 Pelopennes-N samples
	clst_new[meta_df.clst=="Peloponnese-N"] = "Greece-N"

	### Italy 
	clst_new[meta_df.clst=="Northern_Italy-Bk"] = "Italy_North-Bk"

	### Netherlands
	set_labels(meta_df, "Netherlands", "Netherlands-BA", clst_new, age=[3700,4400]) # Four samples

	### Spain/Iberia
	set_labels(meta_df, "Spain", "Iberia-EN", clst_new, age=[6000,7500])
	set_labels(meta_df, "Spain", "Iberia-ECA", clst_new, age=[4800,6000])
	set_labels(meta_df, "Spain", "Iberia-LCA", clst_new, age=[4500,4800])
	set_labels(meta_df, "Spain", "Iberia-BA", clst_new, age=[3600,4501])

	set_labels(meta_df, "Portugal", "Iberia-LCA", clst_new, age=[4500,4800])
	set_labels(meta_df, "Portugal", "Iberia-BA", clst_new, age=[3600,4501])

	### Serbia
	#set_labels(meta_df, "Portugal", "Portugal-BA", clst_new, age=[3600,4800])

	### Switzerland 
	#set_labels(meta_df, "Switzerland", "France-BA", clst_new, age=[4000,4800]) # Four samples

	### UK
	set_labels(meta_df, "United Kingdom", "GB-EN", clst_new, age=[5500,7000])
	set_labels(meta_df, "United Kingdom", "GB-LN", clst_new, age=[4800,5500])
	set_labels(meta_df, "United Kingdom", "GB-EBA", clst_new, age=[3600,4800])
	set_labels(meta_df, "United Kingdom", "GB-LBA", clst_new, age=[2000,3600])
	#clst_new[meta_df.clst=="WHG-HG"] = "WHG-HG"  # To get the Orkey island samples as well

	### Ukraine
	set_labels(meta_df, "Ukraine", "Ukraine-HG", clst_new, age=[8100,12000])
	set_labels(meta_df, "Ukraine", "Ukraine-N", clst_new, age=[6300, 8100])
	set_labels(meta_df, "Ukraine", "Ukraine-CA", clst_new, age=[4700,6300])

	####################################
	### Reset some labels for consistency
	clst_new[meta_df.clst=="Mycenaean-BA"] = "Myc-BA"        # To save space
	clst_new[meta_df.clst=="Anatolia-ChL"] = "Anatolia-CA"   # CA
	clst_new[meta_df.clst=="Greece-EN"] = "Greece-N"         # Time consistency
	clst_new[meta_df.clst=="Iran-ChL"] = "Iran-CA"           # CA
	clst_new[meta_df.clst=="Armenia-ChL"] = "Armenia-CA"     # CA


	# Reset the Hunter-Gatherers:
	clst_new[meta_df.clst=="WHG-HG"] = "WHG-HG"
	clst_new[meta_df.clst=="Iron_Gates-HG"] = "Iron_Gates-HG"
	clst_new[meta_df.clst=="Romania-HG"] = "Iron_Gates-HG" # One sample only
	clst_new[meta_df.clst=="Ignore--9"] = "Ignore--9"

	print("Reduced from %i to %i unique Labels." % (len(set(meta_df.clst)), len(set(clst_new))))

	####################################

	ancsard_df = full_df.loc[:anc_sard_ind].copy()

	### Set the three age clusters
	ancsard_df = full_df[:anc_sard_ind]

	### Set the age clsts based on parameter filed
	print("Reseting periods to:")
	print(periods)
	print(period_cutoffs)	
	
	for i in range(len(periods)):
		lower = period_cutoffs[i]
		if i<len(periods)-1:
			upper = period_cutoffs[i+1]
		else:
			upper=10000 # Default upper value

		ancsard_df.loc[(ancsard_df.age < upper) & (ancsard_df.age >= lower), 'clst'] = periods[i]
	
	#print(ancsard_df[:20])
	### Do the Samples from Serra Crabiles that are unclear. LEGACY after we got C14 ages.
	# ancsard_df.loc[meta_df["iid"].str[:3]=="SEC", 'clst'] = "Sar-Unknown"

	### Reset the inclusion:
	ancsard_df.loc[:, 'include_alt'] = 1   # Default: Include everyone
	ancsard_df.loc[ancsard_df["n_cov_snp"] < snp_cutoff, 'include_alt'] = 0

	####################################
	### Save the relabeling
	save_df = full_df.copy()
	save_df.loc[:anc_ind-1, "clst"] = clst_new
	save_df.loc[:anc_sard_ind-1,:] = ancsard_df.loc[:anc_sard_ind,:] # Overwrite the ancient Sardinians


	# Rename the Batch 3 South Spanish Data:
	save_df.loc[save_df.clst=="Iberia_EN0", "clst"] = "IberiaS_EN"
	save_df.loc[save_df.clst=="Canary_Islands_Guanche_mummy.SG", "clst"] = "Guanche"
	save_df.loc[save_df.clst=="Lebanon_Canaanite", "clst"] = "Canaanite"

	
	# Last check for bad samples
	bad_inds = save_df["period_alt"]=="-9"
	
	print("Setting %i bad individuals to Ignore--9" % np.sum(bad_inds))
	save_df.loc[bad_inds, "include_alt"] = 0    # Exclude bad samples
	save_df.loc[save_df["include_alt"].isnull(), "include_alt"] = 1 # Set missing include_alt values to 1 (all moderns)
	save_df.loc[bad_inds, "clst"] = "Ignore--9"

	### Run the explicit Lists:
	print("Setting %i related Samples to include_alt = 2" % len(related))
	print(related)
	save_df.loc[save_df.iid.isin(related), 'include_alt'] = 2   # Set related individuals

	print("Setting %i contaminated Samples to include_alt = 0" % len(contaminated))
	print(contaminated)	
	save_df.loc[save_df.iid.isin(contaminated), 'include_alt'] = 0   # Contaminated ind
	save_df.loc[save_df["include_alt"]==0, "clst"] = "Sar--9"
	
	print(save_df.loc[:anc_sard_ind-1,:].groupby(["clst", "include_alt"]).size())  # Give out overview of Nr used Sardinian Samples

	
	assert(len(save_df)==len(full_df)) # Sanity Check
	save_df.to_csv(meta_path_new, index=False)
	print("Saved to: %s" % meta_path_new)

####################################
####################################
# Only run if called from this file:
if __name__ == "__main__":
	modifiy_meta()

