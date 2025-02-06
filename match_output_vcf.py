#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Match chromosome names on output vcf to input vcf.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf_input',
                    metavar='file.vcf[.gz]',
                    dest='vcf_input',
                    help='input vcf',
                    required=True)
parser.add_argument('--vcf_output',
                    metavar='file.vcf[.gz]',
                    dest='vcf_output',
                    help='output vcf',
                    required=True)
values = parser.parse_args()

input = open(values.vcf_input,'r')
output = open(values.vcf_output, 'r')
new_output = open(values.vcf_output + '_modified', 'w')
for i in range(len(output.readlines())):
    new_output.write('\t'.join((input.readlines()[i]).split('\t')[0:9] + (output.readlines()[i]).split('\t')[9:]))
input.close()
output.close()
new_output.close()