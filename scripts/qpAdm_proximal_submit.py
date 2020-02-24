# Proximal Modeling 
# Created 2019-06-14 by Chi-Chun

import os
import subprocess
import time
import datetime
import shutil
import tempfile
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


# make a random seed here
input_folder = "output/eigenstrat/subset/"

post_Nur_inds = ["MSR002", "MSR003",
                 "VIL004", 'VIL006', 'VIL007', 'VIL009', 'VIL010', "VIL011",
                 "AMC001", "AMC005", "AMC014",
                 "COR001", "COR002",
                 "SNN001", "SNN002", "SNN003", "SNN004"]
mod_Sar = ["Cag", "Car", "Cam", "Ori", "Ogl", "Nuo", "Sas", "Olb"]



# ### Driver script for qpAdm submission
def qpAdm_run(leftpops, rightpops, output_file, input_folder = "output/eigenstrat/subset", 
              input_file="qpall_proximal_subset",
              par_file_folder = "scripts/parfiles_AT/", input_ind_suff="",
              output_folder ="output/qpAdm_rev/", 
              submit_qpAdm = "sbatch scripts/run_qpAdm_lowmem.sbatch",
              all_snps=False, output=True):
    ''' qpAdm submission on midway (See Harald's original notebook files for comparison)
    '''
    tmp_suffix = next(tempfile._get_candidate_names())
    parfile_path = par_file_folder + "parfile_qpAdm" + leftpops[0] + tmp_suffix
    left_path, right_path = par_file_folder + "left" + tmp_suffix, par_file_folder + "right" + tmp_suffix
    
    while any([os.path.isfile(f) for f in [parfile_path, left_path, right_path]]):
        tmp_suffix = next(tempfile._get_candidate_names())
        parfile_path = par_file_folder + "parfile_qpAdm" + leftpops[0] + tmp_suffix
        left_path, right_path = par_file_folder + "left" + tmp_suffix, par_file_folder + "right" + tmp_suffix

    ### Create the parfile:
    with open(parfile_path, 'w') as f:
        f.write("%s\n" % ("DIR: " + input_folder))
        f.write("%s\n" % ("S1: " + input_file))
        indline = "indivname: DIR/S1" + input_ind_suff + ".ind"
        f.write("%s\n" % indline)
        f.write("%s\n" % "snpname: DIR/S1.snp")
        f.write("%s\n" % "genotypename: DIR/S1.geno")
        f.write("%s\n" % ("popleft: " + left_path))
        f.write("%s\n" % ("popright: " + right_path))
        f.write("%s\n" % "details: YES")   
        if all_snps:
            f.write("%s\n" % "allsnps: YES")
    
    ### Write leftpops rightpops:       
    with open(left_path, 'w') as f:
        f.write("\n".join(leftpops))
        
    with open(right_path, 'w') as f:
        f.write("\n".join(rightpops))
      
    ### Run qpAdm
    output_path = output_folder + output_file
    os.system(' '.join([submit_qpAdm, parfile_path, output_path]))
    return 0

def make_dir_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def print_time():
    print('\n')
    now = datetime.datetime.now()
    print(now.strftime("%Y-%m-%d %H:%M:%S"))


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



### n-way model

mod_Sar_inds = [p + str(i) for p in mod_Sar for i in range(1,6)]

def count_midway_job():
    cmd = "squeue --partition=jnovembre | wc -l"
    ps = subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    n_jobs = int(ps.communicate()[0].decode('ascii').strip()) - 1
    return n_jobs

def produce_pool_combination(n, pools = pools):
    return  list(it.combinations(pools, n))

def run_n_way_model(n, target_list = post_Nur_inds + mod_Sar_inds, sources_dict = sources_dict, sleep = 10, max_job = 150):
    '''
    n: n-way
    '''
    start = time.time()
    pool_combination = produce_pool_combination(n)
    
    #count =  count_midway_job()
    count = 1    
    for c in pool_combination:
        if 'Sard.Nuragic' in c:
            out_dir = "NuragicSource/"
            rightpops = a15
            sources_list = [sources_dict[c[i]] for i in range(n)]
            sources_list = list(it.product(*sources_list))
            for t in target_list:
                for s in sources_list:
                    leftpops = [t] + list(s)
                    output_file = out_dir + '.'.join(leftpops) + ".log"
                    if not os.path.exists(output_file):
                        print(output_file)
                        if not count % 100:
                            print('sleeping...')
                            time.sleep(180)
                        else:
                            time.sleep(0.5)
                        qpAdm_run(leftpops, rightpops, output_file, all_snps=True)
                    count += 1
        else:
            out_dir = "NuragicOutgroup/"
            rightpops = a15plus
            sources_list = [sources_dict[c[i]] for i in range(n)]
            sources_list = list(it.product(*sources_list))
            for t in target_list:    
                for s in sources_list:
                    leftpops = [t] + list(s)
                    output_file = out_dir + '.'.join(leftpops) + ".log"
                    if not os.path.exists(output_file):
                        print(output_file)
                        if not count % 100:
                            print('sleeping...')
                            time.sleep(180)
                        else:
                            time.sleep(0.5)
                        qpAdm_run(leftpops, rightpops, output_file, all_snps=True)
                    count += 1
    end = time.time()
    print('Doing {}-way model'.format(n))
    print('elapsed: {} mins'.format((end - start)/60))
    print('submitted {} jobs'.format(count))

# test run
#run_n_way_model(5)
#run_n_way_model(4)
#run_n_way_model(3)
#run_n_way_model(2)
run_n_way_model(1)
