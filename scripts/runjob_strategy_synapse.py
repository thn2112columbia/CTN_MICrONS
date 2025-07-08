import os
import argparse
import time
from sys import platform

from importlib import reload


def runjobs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", "-t", type=int, default=0)
    parser.add_argument("--cluster_", help=" String", default='burg')
    parser.add_argument("--ver", "-v", help="materialization version",type=int, default=1412)
    parser.add_argument("--unpf_samp", "-n", help="number of unproofread cells to sample",type=int, default=1000)
    parser.add_argument("--seed", "-s", help="RNG seed",type=int, default=0)
    parser.add_argument('--gb', '-g', help='number of gbs per cpu',type=int, default=2)
    args = vars(parser.parse_args())
    cluster = str(args["cluster_"])
    ver = int(args["ver"])
    unpf_samp = int(args["unpf_samp"])
    seed = int(args["seed"])
    gb = int(args['gb'])

    if (args["test"]):
        print ("testing commands")
    
    currwd = os.getcwd()

    #--------------------------------------------------------------------------
    # Ofiles folder

    user = os.environ["USER"]

    if cluster=='burg':
        path_2_package="/burg/theory/users/"+user+"/CTN_MICrONS/scripts"
        ofilesdir = path_2_package + "/Ofiles/"
        resultsdir = path_2_package + "/results/"
        
    elif cluster=='axon':
        path_2_package="/home/"+user+"/CTN_MICrONS/scripts"
        ofilesdir = path_2_package + "/Ofiles/"
        resultsdir = path_2_package + "/results/"

    if not os.path.exists(ofilesdir):
        os.makedirs(ofilesdir)

    if not os.path.exists(resultsdir):
        os.makedirs(resultsdir)

    time.sleep(0.2)
    
    #--------------------------------------------------------------------------
    # Make SBTACH
    inpath = currwd + "/strategy_synapse.py"
    jobprms = f"{inpath}  -v {ver} -n {unpf_samp} -s {seed}"
    jobname = f"pf_strat_v={ver}_n={unpf_samp}_s={seed}"

    if not args["test"]:
        jobnameDir=os.path.join(ofilesdir, jobname)
        text_file=open(jobnameDir, "w");
        os. system("chmod u+x "+ jobnameDir)
        text_file.write("#!/bin/sh \n")
        if cluster=='burg':
            text_file.write("#SBATCH --account=theory \n")
        text_file.write("#SBATCH --job-name="+jobname+ "\n")
        text_file.write("#SBATCH -t 0-11:59  \n")
        text_file.write("#SBATCH --mem-per-cpu={:d}gb \n".format(gb))
        # text_file.write("#SBATCH --gres=gpu\n")
        text_file.write("#SBATCH -c 1 \n")
        text_file.write("#SBATCH -o " + ofilesdir + "/%x.%j.o # STDOUT \n")
        text_file.write("#SBATCH -e " + ofilesdir + "/%x.%j.e # STDERR \n")
        text_file.write("python  -W ignore " + jobprms +" \n")
        text_file.write("echo $PATH  \n")
        text_file.write("exit 0  \n")
        text_file.close()

        if cluster=='axon':
            os.system(f"sbatch -p burst " + jobnameDir);
        else:
            os.system(f"sbatch " + jobnameDir);
    else:
        print (jobprms)

if __name__ == "__main__":
    runjobs()
