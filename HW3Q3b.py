#!/usr/bin/env python

'''Script for Question 3b.
Arguments:
    -f: file containing the sequence (fasta file)
    -mu: the probability of switching states

Outputs:
    prints out a list of 50 state paths generated by random sampling and a list
    of the corresponding likelihoods

Example Usage:
    python HW3Q3b.py -f hmm-sequence.fa -mu 0.01 > sequences3b.txt
'''

import argparse
import numpy as np
from numpy import log1p
from math import exp
from random import random


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
    #print prob_f

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
    #print prob_b

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

''' Samples 50 state paths while calculating the log likelihood of the sampled path
Arguments:
    ford: list of probabilities generated by the forward algorithm
    trans_probs: transition log-probabilities (dictionary of dictionaries)
Returns:
    paths: an array of 50 state paths create by random sampling
    likelis: an array of 50 likelihoods corresponding to the state paths
'''
def sample(ford, trans_probs):
    likelis = [0.0]*(50)
    paths = [""]*(50)
    leng = len(ford)/2

    for i in range (0, 50):
        st = 'n'
        prob_h = exp(ford[(leng-1)*2] - sumLogProb(ford[(leng-1)*2], ford[leng*2 - 1]))
        if random() < prob_h:
            st = 'h'
            likeli = ford[(leng-1)*2]
        else:
            st = 'l'
            likeli = ford[leng*2 - 1]
        path = st

        for k in range (leng-2, -1, -1):
            numer = ford[k*2] + trans_probs['h'][st]
            denom = sumLogProb(ford[k*2] + trans_probs['h'][st], ford[k*2 + 1] + trans_probs['l'][st])
            prob_h = exp(numer - denom)
            if random() < prob_h:
                st = 'h'
                likeli = likeli + np.log(prob_h)
            else:
                st = 'l'
                likeli = likeli + np.log(1-prob_h)
            path = st + path
        paths[i] = path
        likelis[i] = likeli

    return paths, likelis


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
    pths, likli = sample(F, transition_probabilities)
    print "State paths:"
    print pths
    print "Likelihoods:"
    print likli

    print "observations at each position:"
    for i in range (0, 1000):
        num_h = 0
        num_l = 0
        for j in range (0,50):
            if pths[j][i: i+1] == 'h':
                num_h = num_h+1
            else:
                num_l = num_l+1
        print str(i) + "," + str(num_h) + ", " + str(num_l)


if __name__ == "__main__":
    main()
