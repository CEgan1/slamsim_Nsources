#!/usr/bin/env python3

##Libraries
import argparse

##Arguments

parser = argparse.ArgumentParser(description='Complete slim output vcf with invariant sites.')
parser.add_argument('--vcf_slim_output',
                    type=str,
                    required=True,
                    help='Path to slim slim output vcf',
                    dest='vcf_slim')
parser.add_argument('--vcf_original',
                    type=str,
                    required=True,
                    help='Path to slim input vcf',
                    dest='vcf_original')
parser.add_argument('--output',
                    type=str,
                    required=True,
                    help='Path to output transformed VCF.',
                    dest='output_vcf')
parser.add_argument('--sample_size',
                    type=int,
                    required=True,
                    help='size of sample from each pop.',
                    dest='N')

values = parser.parse_args()

##Transformation
#Collect all the sites that are referenced in the original vcf
all_sites = []
vcf_original = open(values.vcf_original, "r")
for line in vcf_original.readlines():
    if line[0] != "#":
        site = line.split('\t')[1]
        all_sites.append(site)
vcf_original.close()

#Collect variant sites that are referenced in the vcf in output of the slim simulation
#Also save the header
header = []
sites_in_slim = []
vcf_slim_output = open(values.vcf_slim, "r")
for line in vcf_slim_output.readlines():
    #header lines start with #
    if line[0] != "#":
        site = line.split('\t')[1]
        sites_in_slim.append(site)
    else:
        header.append(line)
vcf_slim_output.close()


#Add invariant sites in a new output vcf
output_vcf = open(values.output_vcf + '_intermediary', "w")
#Write header
output_vcf.writelines(header)
#read the vcf in output of the slim simulation ans skip header
vcf_slim_output = open(values.vcf_slim, "r")
for line in range(len(header)):
    vcf_slim_output.readline()
#One line per site
for site in all_sites:
    #if the site is in the site is variant (in the slim output vcf), write it directly
    if site in sites_in_slim:
        output_vcf.write(vcf_slim_output.readline())
    #if the site is invariant in the sample, write a line with only 0|0
    else:
        new_line = ['\n'] + ["1\t", site, "\t.\tA\tT\t1000\tPASS\tMID=3;S=0;DOM=0.5;PO=-1;TO=1;MT=1;AC=4;DP=1000\tGT"] + ["\t0|0" for n in range(3*values.N)]
        output_vcf.write(new_line)
vcf_slim_output.close()
output_vcf.close()

#Matching chromosomes between input and output
input = open(values.vcf_original,'r')
output = open(values.output_vcf + '_intermediary', 'r')
new_output = open(values.output_vcf , 'w')
lines_input = input.readlines()
input_header = 0
lines_output = output.readlines()
output_header = 0
for i in range(len(lines_output)):
    if lines_output[i][0] == '#':
        new_output.write(lines_output[i])
        output_header += 1
    if i < len(lines_input):
        if lines_input[i][0] == '#':
            input_header += 1
    if lines_output[i][0] != '#':
        new_output.write('\t'.join((lines_input[i + input_header - output_header]).split('\t')[0:9] + (lines_output[i]).split('\t')[9:]))
input.close()
output.close()
new_output.close()

#Separating autosomes, X and Y and formatting hemizygosity
input = open(values.output_vcf, 'r')
autosome_output = open(values.output_vcf + '_A', 'w')
X_output = open(values.output_vcf + '_X', 'w')
Y_output = open(values.output_vcf + '_Y', 'w')
for line in input.readlines():
    if line[0] == '#':
        X_output.write(line)
        Y_output.write(line)
        autosome_output.write(line)
    if line[0] == 'X':
        X_output.write(line)
    if line[0] == 'Y':
        Y_output.write(line)
    else:
        autosome_output.write(line)
input.close()
X_output.close()
Y_output.close()
