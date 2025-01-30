#!/usr/bin/env python3

import sys
import argparse
import random
import os.path
import os

global max_nb_try
max_nb_try = 10000

# Add new functions following the existing pattern for s1 and s2
def exit_on_contrib_s3_parsing_fail(contrib_str) :
    print('Failed parsing contrib_s3 parameter: {0}'.format(contrib_str), file = sys.stderr)
    print('Exiting', file = sys.stderr)
    sys.exit(1)

def exit_on_contrib_s4_parsing_fail(contrib_str) :
    print('Failed parsing contrib_s4 parameter: {0}'.format(contrib_str), file = sys.stderr)
    print('Exiting', file = sys.stderr)
    sys.exit(1)

def exit_on_contrib_s5_parsing_fail(contrib_str) :
    print('Failed parsing contrib_s5 parameter: {0}'.format(contrib_str), file = sys.stderr)
    print('Exiting', file = sys.stderr)
    sys.exit(1)

# creates parser for command line arguments and parses them
def parse_args() :
    parser = argparse.ArgumentParser(description='Generates parameters for \'slamsim\' simulations.')
    parser.add_argument('-S', '--nb_simulation',
                        type=int,
                        default=1,
                        help='Number of parameters sets to create',
                        dest='nb_simulation')
    parser.add_argument('-N', '--nb_generation',
                        type=int,
                        required=True,
                        help='Number of generations to simulate',
                        dest='nb_generation')
    parser.add_argument('-P', '--prefix',
                        required=True,
                        help='Prefix for parameters files names',
                        dest='prefix')
    parser.add_argument('--N1',
                        required=True,
                        dest='N1',
                        help='Description of admixed population N1 evolution over time. See documentation for help about syntax')
    parser.add_argument('--N2',
                        required=True,
                        dest='N2',
                        help='Description of admixed population N2 evolution over time. See documentation for help about syntax')
    parser.add_argument('--N3',
                        required=True,
                        dest='N3',
                        help='Description of admixed population N3 evolution over time. See documentation for help about syntax')
    parser.add_argument('--N4',
                        required=True,
                        dest='N4',
                        help='Description of admixed population N4 evolution over time. See documentation for help about syntax')
    parser.add_argument('--N5',
                        required=True,
                        dest='N5',
                        help='Description of admixed population N5 evolution over time. See documentation for help about syntax')
    parser.add_argument('--Nadm',
                        required=True,
                        dest='Nadm',
                        help='Description of admixed population Nadm evolution over time. See documentation for help about syntax')
    parser.add_argument('--contrib_s1',
                        required=True,
                        dest='contrib_s1',
                        help='Description of s1 population contribution to admixed population over time. See documentation for help about syntax')
    parser.add_argument('--contrib_s2',
                        required=True,
                        dest='contrib_s2',
                        help='Description of s2 population contribution to admixed population over time. See documentation for help about syntax')
    # New arguments for s3, s4, s5
    parser.add_argument('--contrib_s3',
                        required=True,
                        dest='contrib_s3',
                        help='Description of s3 population contribution to admixed population over time. See documentation for help about syntax')
    parser.add_argument('--contrib_s4',
                        required=True,
                        dest='contrib_s4',
                        help='Description of s4 population contribution to admixed population over time. See documentation for help about syntax')
    parser.add_argument('--contrib_s5',
                        required=True,
                        dest='contrib_s5',
                        help='Description of s5 population contribution to admixed population over time. See documentation for help about syntax')
    parser.add_argument('--force-rewrite',
                        default=False,
                        action='store_true',
                        dest='force_rewrite',
                        help='Should we rewrite a file if it already exists? Default=False')
    
                        
    values = parser.parse_args()
    return (values.nb_simulation, values.nb_generation, values.prefix, 
            values.N1, values.N2, values.N3, values.N4, values.N5, values.Nadm, 
            values.contrib_s1, values.contrib_s2, 
            values.contrib_s3, values.contrib_s4, values.contrib_s5, 
            values.force_rewrite)

# Add functions for print_contrib_interpretation for s3, s4, s5
def print_contrib_s3_interpretation(contrib, contrib_id=3) :
    print('Contribution {0} interpreted as:'.format(contrib_id), file = sys.stderr)
    print('\tFounding contribution: {0}'.format(contrib[0]), file = sys.stderr)
    print('\tContribution scheme: {0}'.format(contrib[1]), file = sys.stderr)
    print ('\tInitial contribution range: {0}'.format(contrib[2]), file = sys.stderr)
    print('\tFinal contribution range: {0}'.format(contrib[3]), file = sys.stderr)

def print_contrib_s4_interpretation(contrib, contrib_id=4) :
    print('Contribution {0} interpreted as:'.format(contrib_id), file = sys.stderr)
    print('\tFounding contribution: {0}'.format(contrib[0]), file = sys.stderr)
    print('\tContribution scheme: {0}'.format(contrib[1]), file = sys.stderr)
    print ('\tInitial contribution range: {0}'.format(contrib[2]), file = sys.stderr)
    print('\tFinal contribution range: {0}'.format(contrib[3]), file = sys.stderr)

def print_contrib_s5_interpretation(contrib, contrib_id=5) :
    print('Contribution {0} interpreted as:'.format(contrib_id), file = sys.stderr)
    print('\tFounding contribution: {0}'.format(contrib[0]), file = sys.stderr)
    print('\tContribution scheme: {0}'.format(contrib[1]), file = sys.stderr)
    print ('\tInitial contribution range: {0}'.format(contrib[2]), file = sys.stderr)
    print('\tFinal contribution range: {0}'.format(contrib[3]), file = sys.stderr)

# Modify parse_contrib_str to use new print_contrib interpretation functions
def parse_contrib_str(contrib_str, contrib_id) :
    # Content remains the same as original, just add more print interpretation functions
    if contrib_id == 3:
        return parse_contrib_str_with_check(contrib_str, contrib_id, print_contrib_s3_interpretation)
    if contrib_id == 4:
        return parse_contrib_str_with_check(contrib_str, contrib_id, print_contrib_s4_interpretation)
    if contrib_id == 5:
        return parse_contrib_str_with_check(contrib_str, contrib_id, print_contrib_s5_interpretation)
    # Existing code for s1 and s2...

def parse_contrib_str_with_check(contrib_str, contrib_id, print_interpretation_func):
    # Special case for the 'Pulse' scheme
    if 'Pulse' in contrib_str :
        contrib = parse_pulse_contrib(contrib_str, contrib_id)
        return contrib

    # Rest of existing parse_contrib_str code...
    contrib_splitted = contrib_str.split('/')
    if 'Con' in contrib_str :
        if len(contrib_splitted) != 3 :
            exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    else :
        if len(contrib_splitted) != 4 :
            exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    if contrib_splitted[0] == 'default' :
        contrib_0 = 'default'
    else :
        try :
            contrib_0 = float(contrib_splitted[0])
        except :
            exit_on_contrib_parsing_fail(contrib_str, contrib_id)

    if contrib_splitted[1] not in ('Con','Inc','Dec','All', 'Pulse') :
        exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    if contrib_splitted[2] == 'default' :
        contrib_range_init = (0.0, 1.0)
    else :
        try :
            tmp = contrib_splitted[2].split('-')
            if len(tmp) != 2 :
                raise ValueError
            contrib_range_init = float(tmp[0]), float(tmp[1])
            contrib_range_init = tuple(sorted(contrib_range_init))
        except :
            exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    if contrib_splitted[1] != 'Con' :
        if contrib_splitted[3] == 'default' :
            contrib_range_final = (0.0, 1.0)
        else :
            try :
                tmp = contrib_splitted[3].split('-')
                if len(tmp) != 2 :
                    raise ValueError
                contrib_range_final = float(tmp[0]), float(tmp[1])
                contrib_range_final = tuple(sorted(contrib_range_final))
            except :
                exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    else :
        contrib_range_final = contrib_range_init        
    contrib = (contrib_0, contrib_splitted[1], contrib_range_init, contrib_range_final)
    print_interpretation_func(contrib)
    return contrib

# Modify check_contributions_consistency to handle s3, s4, s5
def check_contributions_consistency(p1, p2, p3, p4, p5) :
    # Existing checks
    if p1[0] != 'default' and p2[0] != 'default' :
        if p1[0] + p2[0] > 1 :
            exit_on_compatibility_failure()
    
    # Add checks for additional sources
    contribution_sources = [p1, p2, p3, p4, p5]
    
    for i in range(len(contribution_sources)):
        for j in range(i+1, len(contribution_sources)):
            pi = contribution_sources[i]
            pj = contribution_sources[j]
            
            if pi[1] == 'Pulse' or pj[1] == 'Pulse':
                continue
            
            if pi[0] != 'default' and pj[0] != 'default' :
                if pi[0] + pj[0] > 1 :
                    exit_on_compatibility_failure()
            
            if pi[2][0] + pj[2][0] > 1 :
                exit_on_compatibility_failure()
            
            if pi[3][0] + pj[3][0] > 1 :
                exit_on_compatibility_failure()

# Modify generate_c0s to handle additional sources
def generate_c0s(contribs):
    # Find how many sources have 'default' contribution
    default_indices = [i for i, contrib in enumerate(contribs) if isinstance(contrib[0], str)]
    
    if len(default_indices) == 0:
        return contribs
    
    total_known_contrib = sum(contrib[0] for contrib in contribs if not isinstance(contrib[0], str))
    
    # Distribute remaining contribution equally among default sources
    remaining_contrib = 1 - total_known_contrib
    default_contrib = remaining_contrib / len(default_indices)
    
    for idx in default_indices:
        contribs[idx] = list(contribs[idx])
        contribs[idx][0] = default_contrib
        contribs[idx] = tuple(contribs[idx])
    
    return contribs

# Modify print_real_parameters to include s3, s4, s5
def print_real_parameters(prefix, idx_simu, u_values, contribs, uN1, N1s) :
    header = ['N1.0', 'N1.1', 'N1.N', 'N1.u']
    params = [str(N1s[0]), str(N1s[1]), str(N1s[-1]), str(uN1)]
    
    for i in range(1, 6):  # s1 to s5
        header.append(f's{i}.0')
        params.append(str(contribs[i-1][0]))
        
        u_value = u_values[i-1]
        if isinstance(u_value, tuple):
            times, intensities = u_value
            nb_pulse = len(times)
            for j in range(nb_pulse):
                header.append(f't{i}.{j+1}')
                params.append(str(times[j]))
                header.append(f'i{i}.{j+1}')
                params.append(str(intensities[j]))
        else:
            header.extend([f's{i}.1', f's{i}.N', f's{i}.u'])
            params.extend([str(contribs[i-1][1]), str(contribs[i-1][-1]), str(u_value)])
    
    f = open('{0}/simu_{1}/simu_{1}.txt'.format(prefix, idx_simu), 'w')
    print('\t'.join(header), file = f)
    print('\t'.join(params), file = f)
    f.close()

# Modify main to handle s3, s4, s5
def main() :
    global max_nb_try
    # parsing command line
    nb_simulation, nb_generation, prefix, N1_str, N2_str, N3_str, N4_str, N5_str, Nadm_str, contrib_s1_str, contrib_s2_str, contrib_s3_str, contrib_s4_str, contrib_s5_str, force_rewrite = parse_args()

    # parsing patterns
    N1_pattern = parse_N1_str(N1_str)
    N2_pattern = parse_N2_str(N2_str)
    N3_pattern = parse_N3_str(N3_str)
    N4_pattern = parse_N4_str(N4_str)
    N5_pattern = parse_N5_str(N5_str)
    Nadm_pattern = parse_Nadm_str(Nadm_str)
    s1_pattern = parse_contrib_str(contrib_s1_str, 1)
    s2_pattern = parse_contrib_str(contrib_s2_str, 2)
    s3_pattern = parse_contrib_str(contrib_s3_str, 3)
    s4_pattern = parse_contrib_str(contrib_s4_str, 4)
    s5_pattern = parse_contrib_str(contrib_s5_str, 5)
    
    # sanity checks + adjusting range
    N1_pattern = make_N1_sanity_check(N1_pattern)
    N2_pattern = make_N2_sanity_check(N2_pattern)
    N3_pattern = make_N3_sanity_check(N3_pattern)
    N4_pattern = make_N4_sanity_check(N4_pattern)
    N5_pattern = make_N5_sanity_check(N5_pattern)
    Nadm_pattern = make_Nadm_sanity_check(Nadm_pattern)
    s1_pattern = make_contrib_sanity_check(s1_pattern, 1)
    s2_pattern = make_contrib_sanity_check(s2_pattern, 2)
    s3_pattern = make_contrib_sanity_check(s3_pattern, 3)
    s4_pattern = make_contrib_sanity_check(s4_pattern, 4)
    s5_pattern = make_contrib_sanity_check(s5_pattern, 5)
                                      
    for i in range(1, nb_simulation+1) :
        uN1, new_N1s = generate_N1s(N1_pattern, nb_generation)
        uN2, new_N2s = generate_N2s(N2_pattern, nb_generation)
        uN3, new_N3s = generate_N3s(N3_pattern, nb_generation)
        uN4, new_N4s = generate_N4s(N4_pattern, nb_generation)
        uN5, new_N5s = generate_N5s(N5_pattern, nb_generation)
        uNadm, new_Nadms = generate_Nadms(Nadm_pattern, nb_generation)
        contribs_compatible = False
        nb_try = 0
        while not contribs_compatible :
            if nb_try > max_nb_try :
                exit_on_nb_try(max_nb_try)
            nb_try += 1
            u1, contribs_s1 = generate_contribution(s1_pattern, nb_generation)
            u2, contribs_s2 = generate_contribution(s2_pattern, nb_generation)
            u3, contribs_s3 = generate_contribution(s3_pattern, nb_generation)
            u4, contribs_s4 = generate_contribution(s4_pattern, nb_generation)
            u5, contribs_s5 = generate_contribution(s5_pattern, nb_generation)
            contribs_s1, contribs_s2, contribs_s3, contribs_s4, contribs_s5 = generate_c0s(contribs_s1, contribs_s2, contribs_s3, contribs_s4, contribs_s5)

        # printing result to file
        filename = '{0}/simu_{1}/simu_{1}.par'.format(prefix, i)
        if os.path.isfile(filename) and force_rewrite == False :
            print('File {0} exists. Doing nothing. Check the force-rewrite argument'.format(filename), file = sys.stderr)
            continue
        try :
            os.mkdir('{0}/simu_{1}'.format(prefix, i))
        except :
            pass
        f = open(filename, 'w')
        print('g\tN1\tN2\tN3\tN4\tN5\tNadm\tc1\tc2\tc3\tc4\tc5', file = f)
        for j in range(nb_generation+1) :
            print('\t'.join([str(j), str(new_N1s[j]), str(new_N2s[j]), str(new_N3s[j]), str(new_N4s[j]), str(new_N5s[j]), str(new_Nadms[j]), str(contribs_s1[j]), str(contribs_s2[j]), str(contribs_s3[j]), str(contribs_s4[j]), str(contribs_s5[j])]), file = f)
        f.close()
        # printing real parameters we want to estimate latter
        print_real_parameters(prefix, i, u1, contribs_s1, u2, contribs_s2, uN1, new_N1s)

if __name__ == '__main__' :
    main()