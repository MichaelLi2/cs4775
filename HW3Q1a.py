#!/usr/bin/env python

'''Script for computing GC-rich and GC-poor intervals in a given sequence.
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states
    -out: file to output intervals to (1 interval per line)

Outputs:
    File with list of intervals (a_i, b_i) such that bases a_i to b_i are
    classified as GC-rich.
Example Usage:
    python HW3Q1a.py -f hmm-sequence.fa -mu 0.01 -out viterbi-intervals.txt
'''

import argparse
import numpy as np


'''Reads the fasta file and outputs the sequence to analyze.
Arguments:
	filename: name of the fasta file
Returns:
	s: string with relevant sequence
'''
def read_fasta(filename):
    with open(filename, "rb") as f:
        s = ""
        for l in f.readlines()[1:]:
            s += l.strip()
    return s


''' Outputs the Viterbi decoding of a given observation.
Arguments:
	obs: observed sequence of emitted states (list of emissions)
	trans_probs: transition log-probabilities (dictionary of dictionaries)
	emiss_probs: emission log-probabilities, formatted in the reverse order
                     of the posterior probability (dictionary of dictionaries)
	init_probs: initial log-probabilities for each hidden state (dictionary)
Returns:
	l: list of most likely hidden states at each position (list of hidden
           states)
	p: log-probability of the returned hidden state sequence
'''
def viterbi(obs, trans_probs, emiss_probs, init_probs):
    vit = [0.0]*(len(obs) * 2)
    traceback = ['n']*(len(obs) * 2)
    states = ['n']*len(obs)
    prob = 0.0
    vit[0] = init_probs['h']
    vit[1] = init_probs['l']

    for i in range (1, len(obs)):
        prev_max_h = 0.0
        prev_max_l = 0.0
        prev_h2h = vit[i*2 - 2] + trans_probs['h']['h']
        prev_l2h = vit[i*2 - 1] + trans_probs['l']['h']
        prev_h2l = vit[i*2 - 2] + trans_probs['h']['l']
        prev_l2l = vit[i*2 - 1] + trans_probs['l']['l']
        if prev_h2h > prev_l2h:
            prev_max_h = prev_h2h
            traceback[i*2] = 'h'
        else:
            prev_max_h = prev_l2h
            traceback[i*2] = 'l'

        if prev_h2l > prev_l2l:
            prev_max_l = prev_h2l
            traceback[i*2 + 1] = 'h'
        else:
            prev_max_l = prev_l2l
            traceback[i*2 + 1] = 'l'

        v_h = emiss_probs['h'][obs[i: i+1]] + prev_max_h
        v_l = emiss_probs['l'][obs[i: i+1]] + prev_max_l

        vit[i*2] = v_h
        vit[i*2 + 1] = v_l

    trace_pointer = -1
    if v_h > v_l:
        states[len(obs)-1] = 'h'
        prob = v_h
        trace_pointer = 0
    else:
        states[len(obs)-1] = 'l'
        prob = v_l
        trace_pointer = 1

    for i in range (len(obs)-1, 0, -1):
        states[i-1] = traceback[i*2 + trace_pointer]
        if traceback[i*2 + trace_pointer] == 'h':
            trace_pointer = 0
        else:
            trace_pointer = 1
    #print states
    #print traceback
    #print prob
    return states, prob


''' Returns a list of non-overlapping intervals describing the GC rich regions.
Arguments:
	sequence: list of hidden states
Returns:
	intervals: list of tuples (i, j), 1 <= i <= j <= len(sequence), that
                   describe GC rich regions in the input list of hidden states.
'''
def find_intervals(sequence):
    intervals = []
    start = -1
    end = -1
    for i in range (1, len(sequence)):
        if sequence[i] == 'h' and start == -1:
            start = i+1
        if sequence[i] == 'l' and sequence[i-1] == 'h':
            end = i
            intervals.append((start, end))
            start = -1
            end = -1
    if sequence[len(sequence)-1] == 'h':
        intervals.append((start, len(sequence)))
    #print intervals
    return intervals



def main():
    parser = argparse.ArgumentParser(
        description='Parse a sequence into GC-rich and GC-poor regions using Viterbi.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-mu', action="store", dest="mu", type=float, required=True)
    parser.add_argument('-out', action="store", dest="out", type=str, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    mu = args.mu
    intervals_file = args.out

    obs_sequence = read_fasta(fasta_file)
    transition_probabilities = {
        'h': {'h': np.log(1 - mu), 'l': np.log(mu)},
        'l': {'h': np.log(mu), 'l': np.log(1 - mu)}
    }
    emission_probabilities = {
        'h': {'A': np.log(0.13), 'C': np.log(0.37), 'G': np.log(0.37), 'T': np.log(0.13)},
        'l': {'A': np.log(0.32), 'C': np.log(0.18), 'G': np.log(0.18), 'T': np.log(0.32)}
    }
    initial_probabilities = {'h': np.log(0.5), 'l': np.log(0.5)}
    sequence, p = viterbi(obs_sequence, transition_probabilities,
                        emission_probabilities, initial_probabilities)
    intervals = find_intervals(sequence)
    with open(intervals_file, "w") as f:
        f.write("\n".join([("(%d,%d)" % (start, end)) for (start, end) in intervals]))
        f.write("\n")


if __name__ == "__main__":
    main()
