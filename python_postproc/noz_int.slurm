#!/bin/bash
#SBATCH -J noz_int_post
#SBATCH -p wholenode
#SBATCH -o noz_int_post.%j.out
#SBATCH -e noz_int_post.%j.err
#SBATCH -N 1 #56 #50 #28 #48 #40 
#SBATCH -n 128  #7168 #6400 #3584 #6144 #5120
#SBATCH -t 3:00:00 #RUN time (hh:mm:ss)
#SBATCH --mail-user=sd0126@stanford.edu
#SBATCH --mail-type=ALL
##SBATCH --exclude=a485,a[421-430],a513,a958,a952 
date
#configure module
source ~/.bashrc
setupCharles
conda activate post_proc_env

echo $pwd
mybalance

python3 /home/x-sdai/ONR_SCRIPTS/python_postproc/noz_int_post.py

echo "post proc done"
mybalance
