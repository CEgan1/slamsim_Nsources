#!/usr/bin/env python3

##libraries
import os
from random import random
import argparse
import subprocess
import numpy as np

## Arguments
parser = argparse.ArgumentParser(description='Launches SLiM simulations from parameters.')
parser.add_argument('--slim_sim',
                    type=str,
                    required=True,
                    help='Path to slim simulations script',
                    dest='slim_path')
parser.add_argument('--directory',
                    type=str,
                    required=True,
                    help='Path to directory in where the par files are',
                    dest='directory')
parser.add_argument('--N_source',
                    type=int,
                    default=2,
                    help='Number of source populations',
                    dest='N_s')
parser.add_argument('--vcf',
                    type=str,
                    required=True,
                    help='Path to VCF for source population. /!\ It MUST contain ids from source pop 1 and source pop 2 in that order.',
                    dest='VCF')
parser.add_argument('--sample',
                    type=int,
                    default = 100,
                    help='size of samples of s1, s2 and adm on which to compute summary statistics.',
                    dest='sample')
parser.add_argument('--bp',
                    type=int,
                    required=True,
                    help='Length of simulated genomes in bp. Must not contradict VCF files.',
                    dest='L')
parser.add_argument('--fAM',
                    type=str,
                    required=True,
                    help='Range of variation of fAM, format is min-max. If set value x, write x-x.',
                    dest='fAM')
parser.add_argument('--sigma',
                    type=str,
                    required=True,
                    help='Range of variation of sigma, format is min-max. If set value x, write x-x.',
                    dest='sigma')
parser.add_argument('--R',
                    default='default_rec_map',
                    help='Recombination map. Default is 10^(-8) for all genome, see default_rec_map.',
                    dest='rec_rate')
parser.add_argument('--j',
                    default = 1,
                    dest='parallel_jobs',
                    help='Number of simulations to run simultaneously (parallel). Default is 1.')
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
parser.add_argument('--asd_path',
                    required = True,
                    dest='asd_path',
                    help='Path to allelic sharing distance computation software')

parser.add_argument('--activate_X', \
                    help='Set this flag to simulate X chromosomes', \
                    action='store_true', \
                    default=False)
parser.add_argument('--activate_Y', \
                    help='Set this flag to simulate Y chromosomes', \
                    action='store_true', \
                    default=False)
parser.add_argument('--bpX',
                    type=int,
                    help='Length of simulated genomes in bp. Must not contradict VCF files.',
                    dest='L_x')
parser.add_argument('--bpY',
                    type=int,
                    help='Length of simulated genomes in bp. Must not contradict VCF files.',
                    dest='L_y')
parser.add_argument('--vcfX',
                    type=str,
                    help='Path to VCF for source population. /!\ It MUST contain ids from source pop 1 and source pop 2 in that order.',
                    dest='VCF_x')
parser.add_argument('--vcfY',
                    type=str,
                    help='Path to VCF for source population. /!\ It MUST contain ids from source pop 1 and source pop 2 in that order.',
                    dest='VCF_y')

values = parser.parse_args()

#current directory
cwd = str(os.getcwd()) + '/'
print(cwd)

#Set ranges of variations for fAM and sigma
min_fAM, max_fAM = float(values.fAM.split('-')[0]), float(values.fAM.split('-')[1])
min_sigma, max_sigma = float(values.sigma.split('-')[0]), float(values.sigma.split('-')[1])

#Arrange variables depending on sexual chromosomes activation
fullVCF = values.VCF + '_all'
vcf_all = open(fullVCF, 'w')
vcf = open(values.VCF, 'r')
vcf_all.writelines(vcf.readlines())
vcf.close()
if values.activate_X:
    vcf = open(values.VCF_x, 'r')
    vcf_all.writelines([line for line in vcf.readlines() if line[0] != '#'])
    vcf.close()
if values.activate_Y:
    vcf = open(values.VCF_y, 'r')
    vcf_all.writelines([line for line in vcf.readlines() if line[0] != '#'])
    vcf.close()
vcf_all.close()


##Write commands in a command file to be later run in parallel.
# Get list of all simulation folders where parfiles are
simu_files = [dirname for dirname in os.listdir(values.directory)]
# Create command file
command_file = open(values.directory + "command_file", "w")
fAM_file = open(values.directory + "fAM_file", "w")
sigma_file = open(values.directory + "sigma_file", "w")

#Extend parfiles
print('Running extend_parfiles.py to make compatible parfiles')
subprocess.run('./extend_parfiles.py --path ' + values.directory, shell = True)

# One directory for each simulation
for dirname in simu_files:

    #ignore eventual command files and fAM_files left from other runs
    if (dirname == "command_file" or dirname == "fAM_file" or dirname == "sigma_file"):
        continue
    #get  path to parfile
    par_file = dirname + ".extpar"
    #output directory
    path_to_directory = values.directory + dirname + "/"

    par_lines = open(path_to_directory + par_file, "r")
    parameters = np.array([line.split() for line in par_lines])
    par_lines.close()
    
    #keep memory of ids in each population
    for s in range(values.N_s):
        s_ind = open(path_to_directory + "s"+ str(s+1) +"_ind", "w")
        s_ind.writelines(["i" + str(i) + "\n" for i in range(s*values.sample, (s+1)*values.sample)])
        s_ind.close()
    adm_ind = open(path_to_directory + "adm_ind", "w")
    adm_ind.writelines(["i" + str(i) + "\n" for i in range(values.N_s*values.sample, (values.N_s+1)*values.sample)])
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
    if values.activate_X and values.activate_Y:
        command2 = ' && slim -d "r=' + "'" + cwd + values.rec_rate + "'" + '" -d bp=' + str(values.L) + ' -d bpX=' + str(values.L_x) + ' -d bpY=' + str(values.L_y) + ' -d "parfile=' + "'" + par_file + "'" + '" -d "vcf_out=' + "'" + vcf_out + "'" + '" -d "vcf=' + "'" + cwd + values.VCF + "'" + '" -d "vcfX=' + "'" + cwd + values.VCF_x + "'" + '" -d "vcfY=' + "'" + cwd + values.VCF_y + "'" + '" -d fAM=' + str(fAM) + ' -d sigma=' + str(sigma) + ' -d sample_size=' + str(values.sample) + ' ' + cwd + values.slim_path
    elif values.activate_X:
        command2 = ' && slim -d "r=' + "'" + cwd + values.rec_rate + "'" + '" -d bp=' + str(values.L) + ' -d bpX=' + str(values.L_x) + ' -d bpY=0 -d "parfile=' + "'" + par_file + "'" + '" -d "vcf_out=' + "'" + vcf_out + "'" + '" -d "vcf=' + "'" + cwd + values.VCF + "'" + '" -d "vcfX=' + "'" + cwd + values.VCF_x + "'" + '" -d fAM=' + str(fAM) + ' -d sigma=' + str(sigma) + ' -d sample_size=' + str(values.sample) + ' ' + cwd + values.slim_path
    elif values.activate_Y:
        command2 = ' && slim -d "r=' + "'" + cwd + values.rec_rate + "'" + '" -d bp=' + str(values.L) + ' -d bpY=' + str(values.L_y) + ' -d bpX=0 -d "parfile=' + "'" + par_file + "'" + '" -d "vcf_out=' + "'" + vcf_out + "'" + '" -d "vcf=' + "'" + cwd + values.VCF + "'" + '" -d "vcfY=' + "'" + cwd + values.VCF_y + "'" + '" -d fAM=' + str(fAM) + ' -d sigma=' + str(sigma) + ' -d sample_size=' + str(values.sample) + ' ' + cwd + values.slim_path
    else:
        command2 = ' && slim -d "r=' + "'" + cwd + values.rec_rate + "'" + '" -d bp=' + str(values.L) + ' -d bpY=0 -d bpX=0 -d "parfile=' + "'" + par_file + "'" + '" -d "vcf_out=' + "'" + vcf_out + "'" + '" -d "vcf=' + "'" + cwd + values.VCF + "'" + '" -d fAM=' + str(fAM) + ' -d sigma=' + str(sigma) + ' -d sample_size=' + str(values.sample) + ' ' + cwd + values.slim_path

    # transform output vcf to add invariant sites
    command3 = ' && ' + values.transform_vcf + ' --vcf_slim_output ' + vcf_out + ' --vcf_original ' +  cwd + fullVCF + ' --output ' + vcf_out + '_transformed --sample_size ' + str(values.sample)

    # compute summary statistics on output vcf (with invariant sites)
    command4 = ' && ' + values.compute_sumstats + ' --N_source ' + str(values.N_s) + ' --vcf ' + vcf_out + '_transformed_A --Rscript-path  ' + values.finalize_sumstats + ' --asd-path  ' + values.asd_path + ' --out_final ' + sumstat
    # remove intermediary files
    command5 = ' && rm ' + vcf_out + '*transformed* && rm ' + vcf_out + ' && rm *ind\n'

    #write command in command file
    command_file.write(command1 + command2 + command3 + command4 + command5)
command_file.close()
fAM_file.close()
sigma_file.close()

## Run Eidos script with different sets of parameters
print('Running : parallel -j ' + values.parallel_jobs + ' --bar < ' + values.directory + 'command_file')
subprocess.run('parallel -j ' + values.parallel_jobs + ' --bar < ' + values.directory + 'command_file', shell = True, executable='/bin/bash')
