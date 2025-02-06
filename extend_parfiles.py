#!/usr/bin/env python3

import sys
import argparse
import random
import os.path
import os

parser = argparse.ArgumentParser(description='Extends parameters for \'slamsim\' simulations.')
parser.add_argument('--path',
                    required=True,
                    help='path to parfiles',
                    dest='path')
                     
values = parser.parse_args()

for i in range(len(os.listdir(values.path)) - 3):
    f_read = open(values.path + 'simu_' + str(i + 1) + '/simu_' + str(i + 1) + '.par', 'r')
    f_write = open(values.path + 'simu_' + str(i + 1) + '/simu_' + str(i + 1) + '.extpar', 'w')
    line = f_read.readline()
    line = line[:-1].split('\t')
    while len(line) < 10:
        n = int((len(line)-2)/2)
        line.insert(n+1, 'N' + str(n+1))
        line.append('c' + str(n+1))
    f_write.write('\t'.join(line) + '\n')
    for line in f_read.readlines():
        line = line[:-1].split('\t')
        while len(line) < 10:
            n = int((len(line)-2)/2)
            line.insert(n+1, '0')
            line.append('0')
        f_write.write('\t'.join(line) + '\n')
    f_read.close()
    f_write.close()


