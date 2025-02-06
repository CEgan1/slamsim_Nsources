#!/usr/bin/env python3

import sys
import argparse
import distutils.spawn
import subprocess
import shlex
import os
import glob

def parse_args() :
    vcftools_path = distutils.spawn.find_executable('vcftools')
    R_path = distutils.spawn.find_executable('Rscript')
    asd_path = distutils.spawn.find_executable('asd')
    parser = argparse.ArgumentParser(description='Computes summary statistics from a VCF file', \
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', \
                        metavar='file.vcf[.gz]', \
                        help='The VCF file to process (possibly gziped). This file should only contain samples from the admixed and source populations', \
                        required=True,
                        dest='vcf')
    parser.add_argument('--N_source',
                        type=int,
                        default=2,
                        help='Number of source populations',
                        dest='N_s')
    parser.add_argument('--threads', \
                        metavar='N', \
                        help='Number of threads used to compute ASD matrix', \
                        default=1, \
                        type=int)
    parser.add_argument('--Rscript-path', \
                        metavar='PATH/TO/finalize_sumstats.R', \
                        help='Path to R script "finalize_sumstats.R"', \
                        default='./R_codes/finalize_sumstats.R')
    parser.add_argument('--vcftools-path', \
                        metavar='PATH', \
                        default=vcftools_path, \
                        help='Path to vcftools')
    parser.add_argument('--R-path', \
                        metavar='PATH', \
                        default=R_path, \
                        help='Path to Rscript')
    parser.add_argument('--asd-path', \
                        metavar='PATH', \
                        default=asd_path, \
                        help='Path to asd')
    parser.add_argument('--keep-files', \
                        help='Set this flag to keep all intermediate files produced by vcftools, etc...', \
                        action='store_true', \
                        default=False)
    parser.add_argument('--out_final', \
                        help='Path to ouput (summary stats)', \
                        required=True)
    tmp = parser.parse_args()
    vcf, N_s, vcftools_path, R_path, asd_path, threads, keep_files, out_final = tmp.vcf, tmp.N_s, tmp.vcftools_path, tmp.R_path, tmp.asd_path, tmp.threads, tmp.keep_files, tmp.out_final
    return vars(tmp)

# tries to determine if a file is gzip (or bgzip) compressed
# uses gzip Magic Bytes, so not garantied to work 100% of the times
def is_gzip_file(path) :
    f = open(path, 'rb')
    magic = f.read(3)
    f.close()
    if magic == b'\x1f\x8b\x08' :
        return True
    return False

# runs vcftools to compute allele frequencies on each population
def compute_allele_frequencies(args) :
    print('Computing allele frequencies...', file=sys.stderr)
    # ADM POP
    command = f"{args['vcftools_path']} {args['vcftools_flag']} {args['vcf']} --keep adm_ind --freq2 --out adm"
    command = shlex.split(command)
    p = subprocess.run(command, stderr=subprocess.PIPE)
    if p.returncode != 0 :
        print('Something went wrong while computing allele frequencies on adm population', file=sys.stderr)
    # Source POP
    for pop in range(args['N_s']):
        command = f"{args['vcftools_path']} {args['vcftools_flag']} {args['vcf']} --keep s" + str(pop + 1) + "_ind --freq2 --out s" + str(pop + 1)
        command = shlex.split(command)
        p = subprocess.run(command, stderr=subprocess.PIPE)
        if p.returncode != 0 :
            print('Something went wrong while computing allele frequencies on s' + str(pop + 1) + ' population', file=sys.stderr)


# runs vcftools to compute inbreeding coefficient on each population
def compute_inbreeding(args) :
    print('Computing inbreeding coefficients...', file=sys.stderr)
    # ADM POP
    command = f"{args['vcftools_path']} {args['vcftools_flag']} {args['vcf']} --keep adm_ind --het --out adm"
    command = shlex.split(command)
    p = subprocess.run(command, stderr=subprocess.PIPE)
    if p.returncode != 0 :
        print('Something went wrong while computing inbreeding coefficients on adm population', file=sys.stderr)
    # Source POP
    for pop in range(args['N_s']):
        command = f"{args['vcftools_path']} {args['vcftools_flag']} {args['vcf']} --keep s" + str(pop + 1) + "_ind --het --out s" + str(pop + 1)
        command = shlex.split(command)
        p = subprocess.run(command, stderr=subprocess.PIPE)
        if p.returncode != 0 :
            print('Something went wrong while computing inbreeding coefficients on s' + str(pop + 1) + ' population', file=sys.stderr)
    

# changes data format and computes ASD matrix
def compute_asd_matrix(args) :
    # put data in plink TPED format
    print('Reformatting data to plink format...', file=sys.stderr)
    command = f"{args['vcftools_path']} {args['vcftools_flag']} {args['vcf']} --plink-tped --out data"
    command = shlex.split(command)
    p = subprocess.run(command, stderr=subprocess.PIPE)
    if p.returncode != 0 :
        print('Something went wrong while reformatting data', file=sys.stderr)
    # compute asd matrix
    print('Computing ASD matrix...', file=sys.stderr)
    command = f"{args['asd_path']} --tped data.tped --tfam data.tfam --threads {args['threads']} --out data"
    print(command)
    command = shlex.split(command)
    p = subprocess.run(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    if p.returncode != 0 :
        print('Something went wrong while computing ASD matrix', file=sys.stderr)

# computes all pairwise Fst
def compute_fst(args) :
    print('Computing all pairwise Fst...', file=sys.stderr)
    # Fst for all pairs of populations
    Fst = []
    for pop1 in range(args['N_s'] - 1):
        command = f"{args['vcftools_path']} {args['vcftools_flag']} {args['vcf']} --weir-fst-pop s" + str(pop1 + 1) + "_ind --weir-fst-pop adm_ind -c"
        command = shlex.split(command)
        p = subprocess.run(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='ascii')
        if p.returncode != 0 :
            print('Something went wrong while computing Fst between s' + str(pop1 + 1) + ' and adm.', file=sys.stderr)
        tmp = p.stderr
        print(tmp)
        tmp = [i for i in tmp.split('\n') if 'weighted' in i][0]
        Fst.append(tmp.split()[-1])
        for pop2 in range(pop1 + 1, args['N_s']):
            command = f"{args['vcftools_path']} {args['vcftools_flag']} {args['vcf']} --weir-fst-pop s" + str(pop1 + 1) + "_ind --weir-fst-pop s" + str(pop2 + 1) + "_ind -c"
            command = shlex.split(command)
            p = subprocess.run(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='ascii')
            if p.returncode != 0 :
                print('Something went wrong while computing Fst between s' + str(pop1 + 1) + ' and ' + str(pop2 + 1), file=sys.stderr)
            tmp = p.stderr
            tmp = [i for i in tmp.split('\n') if 'weighted' in i][0]
            Fst.append(tmp.split()[-1])
    command = f"{args['vcftools_path']} {args['vcftools_flag']} {args['vcf']} --weir-fst-pop s" + str(args['N_s']) + "_ind --weir-fst-pop adm_ind -c"
    command = shlex.split(command)
    p = subprocess.run(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='ascii')
    if p.returncode != 0 :
        print('Something went wrong while computing Fst between s' + str(args['N_s']) + ' and adm.', file=sys.stderr)
    tmp = p.stderr
    print(tmp)
    tmp = [i for i in tmp.split('\n') if 'weighted' in i][0]
    Fst.append(tmp.split()[-1])
    f = open('result.fst', 'w')
    for fst in Fst:
        print(fst, file=f)
    f.close()

# runs the R script to finalize sumstats computations
def run_r_script(args) :
    print('Finalizing sumstats computation...', file=sys.stderr)
    command = f"{args['R_path']} {args['Rscript_path']} {args['out_final']} {args['N_s']}"
    print(command)
    command = shlex.split(command)
    p = subprocess.run(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='ascii')
    if p.returncode != 0 :
        print('Something went wrong while finalizing sumstats computation using R code', file=sys.stderr)

# removes files created during sumstats computations
def remove_temporary_files() :
    try :
        os.remove('data.tped')
        os.remove('data.tfam')
        os.remove('data.error')
        os.remove('data.asd.dist')
        os.remove('data.log')
        os.remove('result.fst')
        for f in glob.glob("*het"):
            os.remove(f)
        for f in glob.glob("*frq"):
            os.remove(f)
    except :
        pass
        
def main() :
    args = parse_args()
    if is_gzip_file(args['vcf']) :
        args['vcftools_flag'] = '--gzvcf'
    else :
        args['vcftools_flag'] = '--vcf'
    compute_allele_frequencies(args)
    compute_inbreeding(args)
    compute_asd_matrix(args)
    compute_fst(args)
    run_r_script(args)
    if not args['keep_files'] :
        remove_temporary_files()

    print('Sumstats computed', file=sys.stderr)

if __name__ == '__main__' :
    main()