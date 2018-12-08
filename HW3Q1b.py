#!/usr/bin/env python

'''Script for computing posterior probabilities of hidden states at each
   position of a given sequence.
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states

Outputs:
    posteriors.csv - a KxN matrix outputted as a CSV with the posterior
                     probability of each state at each position

Example Usage:
    python HW3Q1b.py -f hmm-sequence.fa -mu 0.01
'''

import argparse
import numpy as np
from numpy import log1p
from math import exp


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


''' Sum of log probabilities
Arguments:
    a: log of first probability to be summed
    b: log of second probability to be summed
Returns:
    s: log of sum of probabilities
'''
def sumLogProb(a, b):
    if a > b:
        return a + log1p(exp(b-a))
    else:
        return b + log1p(exp(a-b))

''' Outputs the forward and backward probabilities of a given observation.
Arguments:
	obs: observed sequence of emitted states (list of emissions)
	trans_probs: transition log-probabilities (dictionary of dictionaries)
	emiss_probs: emission log-probabilities, formatted in the reverse order
                     of the posterior probability (dictionary of dictionaries)
	init_probs: initial log-probabilities for each hidden state (dictionary)
Returns:
	F: matrix of forward probabilities
        likelihood_f: P(obs) calculated using the forward algorithm
	B: matrix of backward probabilities
        likelihood_b: P(obs) calculated using the backward algorithm
	R: matrix of posterior probabilities
'''
def forward_backward(obs, trans_probs, emiss_probs, init_probs):
    forw = [0.0]*(len(obs) * 2)
    back = [0.0]*(len(obs) * 2)
    prob_b = 0.0
    post = []

    forw[0] = init_probs['h'] + emiss_probs['h'][obs[0: 1]]
    forw[1] = init_probs['l'] + emiss_probs['l'][obs[0: 1]]

    for i in range (1, len(obs)):
        prev_h = sumLogProb(forw[i*2 - 2] + trans_probs['h']['h'], forw[i*2 - 1] + trans_probs['l']['h'])
        prev_l = sumLogProb(forw[i*2 - 2] + trans_probs['h']['l'], forw[i*2 - 1] + trans_probs['l']['l'])
        f_h = emiss_probs['h'][obs[i: i+1]] + prev_h
        f_l = emiss_probs['l'][obs[i: i+1]] + prev_l
        forw[i*2] = f_h
        forw[i*2 + 1] = f_l
    prob_f = sumLogProb(forw[len(forw)-1], forw[len(forw)-2])
    print prob_f

    back[len(forw)-2] = np.log(1.0)
    back[len(forw)-1] = np.log(1.0)

    for j in range (len(obs)-2, -1, -1):
        b_h = sumLogProb(back[j*2 + 2] + trans_probs['h']['h'] + emiss_probs['h'][obs[j+1: j+2]],
                         back[j*2 + 3] + trans_probs['h']['l'] + emiss_probs['l'][obs[j+1: j+2]])
        b_l = sumLogProb(back[j*2 + 2] + trans_probs['l']['h'] + emiss_probs['h'][obs[j+1: j+2]],
                         back[j*2 + 3] + trans_probs['l']['l'] + emiss_probs['l'][obs[j+1: j+2]])
        back[j*2] = b_h
        back[j*2 + 1] = b_l
    prob_b = sumLogProb(back[0] + emiss_probs['h'][obs[0: 1]] + init_probs['h'],
                        back[1] + emiss_probs['l'][obs[0: 1]] + init_probs['l'])
    print prob_b

    post_row1 = [0.0]*(len(obs))
    post_row2= [0.0]*(len(obs))
    for k in range (0, len(obs)):
        post_h = forw[k*2] + back[k*2]
        post_l = forw[k*2 + 1] + back[k*2 + 1]
        p = sumLogProb(post_h, post_l)
        post_row1[k] = exp(post_h - p)
        post_row2[k] = exp(post_l - p)
    post.append(post_row1)
    post.append(post_row2)

    return forw, prob_f, back, prob_b, post


def main():
    parser = argparse.ArgumentParser(
        description='Compute posterior probabilities at each position of a given sequence.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-mu', action="store", dest="mu", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    mu = args.mu

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
    F, like_f, B, like_b, R = forward_backward(obs_sequence,
                                              transition_probabilities,
                                              emission_probabilities,
                                              initial_probabilities)
    np.savetxt("posteriors.csv", R, delimiter=",", fmt='%.4e')


if __name__ == "__main__":
    main()
