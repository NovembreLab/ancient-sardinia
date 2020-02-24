#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --partition=jnovembre
#SBATCH --account=pi-jnovembre
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30G
#SBATCH --job-name=nb-tunnel
#SBATCH --output=nb-log-%J.out
#SBATCH --error=nb-log-%J.err

# activate the conda env 
source activate anc-sard

ip=$(/sbin/ip route get 8.8.8.8 | awk '{print $NF;exit}')
port=$((10000+ $RANDOM % 20000))

echo "http://"$ip":"$port"/"
jupyter lab --no-browser --ip=$ip --port=$port --log-level='ERROR'
