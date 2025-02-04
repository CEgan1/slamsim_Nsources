#!/usr/bin/env python3

##libraries
import os
from random import random
import argparse
import subprocess
import numpy as np

## Arguments
parser = argparse.ArgumentParser(description='Launches SLiM simulations from parameters.')
parser.add_argument('--p',
                    default = str(os.getcwd()) + '/SLAMSim.parfile',
                    dest = 'parfile',
                    help = 'Path to SLAMSim parfile. Default is SLAMSim.parfile in current directory. See readme for more details on parfile format.')
parser.add_argument('--compute_sumstats',
                    default = str(os.getcwd()) + '/compute_sumstats_n_sources.py',
                    dest='compute_sumstats',
                    help='Path to programm to compute summary statistics. Default is named compute_sumstats_from_vcf.py and is in current directory.')
parser.add_argument('--transform_vcf',
                    default = str(os.getcwd()) + '/vcf_change.py',
                    dest='transform_vcf',
                    help='Path to programm to ad invariant site to vcf in slim output. Default is named add_invariant_sites_to_vcf.py and is in current directory.')
parser.add_argument('--finalize_sumstats',
                    default = str(os.getcwd()) + '/finalize_sumstats.R',
                    dest='finalize_sumstats',
                    help='Path to programm to finalize computation summary statistics. Default is named finalize_sumstats.R and is in current directory.')
values = parser.parse_args()

#current directory
cwd = str(os.getcwd()) + '/'
print(cwd)

#Read parameters from parfile:
parfile = open(values.parfile, 'r')
slim_sim = parfile.readline()[:-1].split(': ')[1]
print('slim simulation is ' + slim_sim)
simu = parfile.readline()[:-1].split(': ')[1]
print('simulation parameters are in ' + simu)
parallel = parfile.readline()[:-1].split(': ')[1]
print(parallel + ' jobs are run in parallel')
N_sources = int(parfile.readline().split(': ')[1])
print(str(N_sources) + ' source populations')
sample = int(parfile.readline().split(': ')[1])
print(str(sample) + ' individuals sampled to compute sumstats')
bp = int(parfile.readline().split(': ')[1])
print('size of autosomal genome is ' + str(bp))
vcf = parfile.readline()[:-1].split(': ')[1]
print('autosomes vcf is ' + vcf)
fAM = parfile.readline()[:-1].split(': ')[1]
print('fAM range is ' + fAM)
sigma = parfile.readline()[:-1].split(': ')[1]
print('sigma range is ' + sigma)
rec_map = parfile.readline()[:-1].split(': ')[1]
print('Recombination map is ' + rec_map)
asd_path = parfile.readline()[:-1].split(': ')[1]
print('Path to ASD software is ' + asd_path)
activateX = parfile.readline()[:-1].split(': ')[1] == 'YES'
print('Simulate chromosome X? ' + str(activateX))
activateY = parfile.readline().split(': ')[1] == 'YES'
print('Simulate chromosome Y? ' + str(activateY))
if activateX:
    vcfX = parfile.readline()[:-1].split(': ')[1]
    print('X vcf is ' + vcfX)
    bpX = int(parfile.readline().split(': ')[1])
    print('size of X chromosome is ' + str(bpX))
else :
    bpX = 0
    vcfX = None
if activateY:
    vcfY = parfile.readline()[:-1].split(': ')[1]
    print('Y vcf is ' + vcfY)
    bpY = int(parfile.readline().split(': ')[1])
    print('size of Y chromosome is ' + str(bpY))
else :
    bpY = 0
    vcfY = None
parfile.close()


#Set ranges of variations for fAM and sigma
min_fAM, max_fAM = float(fAM.split('-')[0]), float(fAM.split('-')[1])
min_sigma, max_sigma = float(sigma.split('-')[0]), float(sigma.split('-')[1])

#Arrange variables depending on sexual chromosomes activation
fullVCF = vcf + '_all'
vcf_all = open(fullVCF, 'w')
VCF = open(vcf, 'r')
vcf_all.writelines(VCF.readlines())
VCF.close()
if activateX:
    VCF = open(vcfX, 'r')
    vcf_all.writelines([line for line in VCF.readlines() if line[0] != '#'])
    VCF.close()
if activateY:
    VCF = open(vcfY, 'r')
    vcf_all.writelines([line for line in VCF.readlines() if line[0] != '#'])
    VCF.close()
vcf_all.close()


##Write commands in a command file to be later run in parallel.
# Get list of all simulation folders where parfiles are
simu_files = [dirname for dirname in os.listdir(simu)]
# Create command file
command_file = open(simu + "command_file", "w")
fAM_file = open(simu + "fAM_file", "w")
sigma_file = open(simu + "sigma_file", "w")

#Extend parfiles
print('Running extend_parfiles.py to make compatible parfiles')
subprocess.run('./extend_parfiles.py --path ' + simu, shell = True)

# One directory for each simulation
for dirname in simu_files:

    #ignore eventual command files and fAM_files left from other runs
    if (dirname == "command_file" or dirname == "fAM_file" or dirname == "sigma_file"):
        continue
    #get  path to parfile
    par_file = dirname + ".extpar"
    #output directory
    path_to_directory = simu + dirname + "/"

    par_lines = open(path_to_directory + par_file, "r")
    parameters = np.array([line.split() for line in par_lines])
    par_lines.close()
    
    #keep memory of ids in each population
    for s in range(N_sources):
        s_ind = open(path_to_directory + "s"+ str(s+1) +"_ind", "w")
        s_ind.writelines(["i" + str(i) + "\n" for i in range(s*sample, (s+1)*sample)])
        s_ind.close()
    adm_ind = open(path_to_directory + "adm_ind", "w")
    adm_ind.writelines(["i" + str(i) + "\n" for i in range(N_sources*sample, (N_sources+1)*sample)])
    adm_ind.close()

    #output vcf
    vcf_out = dirname + ".vcf_out"
    sumstat = dirname + ".sumstat"
    
    #Choose fAM randomly (uniform) and register it.
    fAM = (max_fAM - min_fAM)*random() + min_fAM
    fAM_file.write(dirname + '\t' + str(fAM) + '\n')
    
    #Choose sigma randomly (uniform) and register it.
    sigma = (max_sigma - min_sigma)*random() + min_sigma
    sigma_file.write(dirname + '\t' + str(sigma) + '\n')

    #This command will allow to run a SLIM simulation, to compute summary statistics on the VCF output, then to remove the VCF and other intermediary files.

    #go to sim directory
    command1 = 'cd ' + path_to_directory

    # run slim sim for given parameters
    if activateX and activateY:
        command2 = ' && slim -d "r=' + "'" + cwd + rec_map + "'" + '" -d bp=' + str(bp) + ' -d bpX=' + str(bpX) + ' -d bpY=' + str(bpY) + ' -d "parfile=' + "'" + par_file + "'" + '" -d "vcf_out=' + "'" + vcf_out + "'" + '" -d "vcf=' + "'" + cwd + vcf + "'" + '" -d "vcfX=' + "'" + cwd + vcfX + "'" + '" -d "vcfY=' + "'" + cwd + vcfY + "'" + '" -d fAM=' + str(fAM) + ' -d sigma=' + str(sigma) + ' -d sample_size=' + str(sample) + ' ' + cwd + slim_sim
    elif activateX:
        command2 = ' && slim -d "r=' + "'" + cwd + rec_map + "'" + '" -d bp=' + str(bp) + ' -d bpX=' + str(bpX) + ' -d bpY=0 -d "parfile=' + "'" + par_file + "'" + '" -d "vcf_out=' + "'" + vcf_out + "'" + '" -d "vcf=' + "'" + cwd + vcf + "'" + '" -d "vcfX=' + "'" + cwd + vcfX + "'" + '" -d fAM=' + str(fAM) + ' -d sigma=' + str(sigma) + ' -d sample_size=' + str(sample) + ' ' + cwd + slim_sim
    elif activateY:
        command2 = ' && slim -d "r=' + "'" + cwd + rec_map + "'" + '" -d bp=' + str(bp) + ' -d bpY=' + str(bpY) + ' -d bpX=0 -d "parfile=' + "'" + par_file + "'" + '" -d "vcf_out=' + "'" + vcf_out + "'" + '" -d "vcf=' + "'" + cwd + vcf + "'" + '" -d "vcfY=' + "'" + cwd + vcfY + "'" + '" -d fAM=' + str(fAM) + ' -d sigma=' + str(sigma) + ' -d sample_size=' + str(sample) + ' ' + cwd + slim_sim
    else:
        command2 = ' && slim -d "r=' + "'" + cwd + rec_map + "'" + '" -d bp=' + str(bp) + ' -d bpY=0 -d bpX=0 -d "parfile=' + "'" + par_file + "'" + '" -d "vcf_out=' + "'" + vcf_out + "'" + '" -d "vcf=' + "'" + cwd + vcf + "'" + '" -d fAM=' + str(fAM) + ' -d sigma=' + str(sigma) + ' -d sample_size=' + str(sample) + ' ' + cwd + slim_sim

    # transform output vcf to add invariant sites
    command3 = ' && ' + values.transform_vcf + ' --vcf_slim_output ' + vcf_out + ' --vcf_original ' +  cwd + fullVCF + ' --output ' + vcf_out + '_transformed --sample_size ' + str(sample)

    # compute summary statistics on output vcf (with invariant sites)
    command4 = ' && ' + values.compute_sumstats + ' --N_source ' + str(N_sources) + ' --vcf ' + vcf_out + '_transformed_A --Rscript-path  ' + values.finalize_sumstats + ' --asd-path  ' + asd_path + ' --out_final ' + sumstat
    # remove intermediary files
    command5 = ' && rm ' + vcf_out + '*transformed* && rm ' + vcf_out + ' && rm *ind\n'

    #write command in command file
    command_file.write(command1 + command2 + command3 + command4 + command5)
command_file.close()
fAM_file.close()
sigma_file.close()

## Run Eidos script with different sets of parameters
print('Running : parallel -j ' + parallel + ' --bar < ' + simu + 'command_file')
subprocess.run('parallel -j ' + parallel + ' --bar < ' + simu + 'command_file', shell = True, executable='/bin/bash')
