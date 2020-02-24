#!/bin/bash

#source activate anc-sard
module load python/3.6.1+intel-16.0 python/3.6.1+intel-16.0
module load gsl # For Admixtools
module load openblas # For Admixtools

python3 -m snakemake \
    -kp \
    --ri \
    -j 65 \
    --max-jobs-per-second 5 \
    --cluster-config cluster.json \
    -c "sbatch \
        --time={cluster.time} \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=jnovembre \
        --job-name={cluster.name} \
        --mail-user={cluster.email} \
        --mail-type={cluster.emailtype} \
	--output={cluster.logfile}" \
    $*
