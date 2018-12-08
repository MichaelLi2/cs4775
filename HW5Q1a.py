#!/usr/bin/env python

''' Estimate mu, theta_l, and theta_h via EM.

Arguments:
    -f: sequence to read in
    -mu, -theta_h, -theta_l: initializations for parameter values
Outputs:
    EM estimates for theta_h, theta_l, and mu
    em_<mu>_<theta_h>_<theta_l>.png - file containing plot of log likelihoods
                                      over the EM iterations
                                      (see ```saveplot```)

Example Usage:
    python HW5Q1a.py -f hmm-sequence.fa -mu 0.05 -theta_h 0.6 -theta_l 0.4
'''

import argparse
import numpy as np
#import matplotlib.pyplot as plt
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


''' Generates the transition and emission probabilities table given the
    parameters.
Arguments:
    mu, theta_h, theta_l: parameters as described in the problem
Returns:
    transition_probabilities, emission_probabilities
        (both are dictionaries of dictionaries)
'''
def get_probabilities(mu, theta_h, theta_l):
    transition_probabilities = {
        'h': {'h': np.log(1 - mu), 'l': np.log(mu)},
        'l': {'h': np.log(mu), 'l': np.log(1 - mu)}
    }
    emission_probabilities = {
        'h': {'A': np.log(.5 * (1 - theta_h)), 'C': np.log(.5 * theta_h),
              'G': np.log(.5 * theta_h), 'T': np.log(.5 * (1 - theta_h))},
        'l': {'A': np.log(.5 * (1 - theta_l)), 'C': np.log(.5 * theta_l),
              'G': np.log(.5 * theta_l), 'T': np.log(.5 * (1 - theta_l))},
    }
    return transition_probabilities, emission_probabilities


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
	trans_probs: transition log probabilities (dictionary of dictionaries)
	emiss_probs: emission log probabilities (dictionary of dictionaries)
	init_probs: initial log probabilities for each hidden state (dictionary)
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


''' Performs 1 EM step.
Arguments:
    fb_values: relevant output variables from forward-backward
    obs: the sequence under analysis
    tp: transition probabilities in log space
    ep: emission probabilities in log space
    ip: initalization probabilities in log space
Returns:
    tp: updated transition probabilities, in log space
    ep: updated emission probabilities, in log space
'''
def em(fb_output, obs, tp, ep, ip):
    #e step
    ekb = [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
    akl = [[0.0, 0.0], [0.0, 0.0]]

    for i in range (0, len(obs)):
        if obs[i:i+1] == 'A':
            ekb[0][0] += fb_output[4][0][i]
            ekb[1][0] += fb_output[4][1][i]
        elif obs[i:i+1] == 'C':
            ekb[0][1] += fb_output[4][0][i]
            ekb[1][1] += fb_output[4][1][i]
        elif obs[i:i+1] == 'G':
            ekb[0][2] += fb_output[4][0][i]
            ekb[1][2] += fb_output[4][1][i]
        else:
            ekb[0][3] += fb_output[4][0][i]
            ekb[1][3] += fb_output[4][1][i]

        if i == len(obs)-1:
            continue

        akl[0][0] += exp(fb_output[0][i*2] + tp['h']['h'] + ep['h'][obs[i+1:i+2]] + fb_output[2][i*2+2]
                    -sumLogProb(fb_output[0][len(obs)*2-2], fb_output[0][len(obs)*2-1]))
        akl[0][1] += exp(fb_output[0][i*2] + tp['h']['l'] + ep['l'][obs[i+1:i+2]] + fb_output[2][i*2+3]
                    -sumLogProb(fb_output[0][len(obs)*2-2], fb_output[0][len(obs)*2-1]))
        akl[1][0] += exp(fb_output[0][i*2+1] + tp['l']['h'] + ep['h'][obs[i+1:i+2]] + fb_output[2][i*2+2]
                    -sumLogProb(fb_output[0][len(obs)*2-2], fb_output[0][len(obs)*2-1]))
        akl[1][1] += exp(fb_output[0][i*2+1] + tp['l']['l'] + ep['l'][obs[i+1:i+2]] + fb_output[2][i*2+3]
                    -sumLogProb(fb_output[0][len(obs)*2-2], fb_output[0][len(obs)*2-1]))

    #m step
    ekb_htotal = 0.0
    ekb_ltotal = 0.0
    akl_htotal = 0.0
    akl_ltotal = 0.0

    for a in range (0,4):
        ekb_htotal += ekb[0][a]
        ekb_ltotal += ekb[1][a]
    for b in range (0,2):
        akl_htotal += akl[0][b]
        akl_ltotal += akl[1][b]

    for c in range (0,4):
        ekb[0][c] = ekb[0][c] / ekb_htotal
        ekb[1][c] = ekb[1][c] / ekb_ltotal
    for d in range (0,2):
        akl[0][d] = akl[0][d] / akl_htotal
        akl[1][d] = akl[1][d] / akl_ltotal

    new_tp = {
        'h': {'h': np.log(akl[0][0]), 'l': np.log(akl[0][1])},
        'l': {'h': np.log(akl[1][0]), 'l': np.log(akl[1][1])}
    }
    new_ep = {
        'h': {'A': np.log(ekb[0][0]), 'C': np.log(ekb[0][1]),
              'G': np.log(ekb[0][2]), 'T': np.log(ekb[0][3])},
        'l': {'A': np.log(ekb[1][0]), 'C': np.log(ekb[1][1]),
              'G': np.log(ekb[1][2]), 'T': np.log(ekb[1][3])},
    }
    return new_tp, new_ep


''' Helper function to save plot of log likelihoods over iterations to file for
    visualization.
Arguments:
    log_likelihoods: list of log likelihoods over iterations
    init_mu, init_theta_h, init_theta_l: the initial values of parameters used
        (for naming the file containing the plot)
Outputs:
    plot of log likelihoods to file
'''
def saveplot(log_likelihoods, mu, theta_h, theta_l):
    plt.title("EM log likelihoods with initialization %.2f, %.2f, %.2f" % (mu, theta_h, theta_l))
    plt.xlabel("Iteration")
    plt.ylabel("Log likelihood")
    plt.plot(range(len(log_likelihoods)), log_likelihoods, 'r-')
    plt.savefig("em_%.2f_%.2f_%.2f.png" % (mu, theta_h, theta_l))


''' Uses EM to infer the parameters ```mu, theta_h, theta_l```, iterating until
    a valid stopping condition is reached.
Arguments:
    sequence: sequence data to train on
    mu: the value of mu to use for initializing the transition probabilities
    theta_h, theta_l: parameters of the emission probability distributions
Returns:
    mu: parameter of trained transition probability distribution
    theta_h, theta_l: parameters of trained emission probability distribution
'''
def train(sequence, mu, theta_h, theta_l, stop_diff = 0.0001):
    init_mu, init_theta_h, init_theta_l = mu, theta_h, theta_l
    trans_probs, emiss_probs = get_probabilities(mu, theta_h, theta_l)
    init_probs = {'h': np.log(0.5), 'l': np.log(0.5)}
    log_likelihoods = [] # list of log likelihoods from each iteration

    while len(log_likelihoods) < 3 or log_likelihoods[len(log_likelihoods)-1] - log_likelihoods[len(log_likelihoods)-2] > stop_diff:
    #for i in range (0,100):
        outpfb = forward_backward(sequence, trans_probs, emiss_probs, init_probs)
        trans_probs, emiss_probs = em(outpfb, sequence, trans_probs, emiss_probs, init_probs)
        log_likelihoods.append(sumLogProb(outpfb[0][len(sequence)*2-2], outpfb[0][len(sequence)*2-1]))
        #print trans_probs, emiss_probs
        #print exp(trans_probs['h']['l']), 2*exp(emiss_probs['h']['C']), 2*exp(emiss_probs['l']['C'])
        print sumLogProb(outpfb[0][len(sequence)*2-2], outpfb[0][len(sequence)*2-1])

    #saveplot(log_likelihoods, init_mu, init_theta_h, init_theta_l)
    return ((exp(trans_probs['h']['l'])+exp(trans_probs['l']['h']))/2, exp(emiss_probs['h']['C'])+exp(emiss_probs['h']['G']),
            exp(emiss_probs['l']['C'])+exp(emiss_probs['l']['G']))


def main():
    parser = argparse.ArgumentParser(description='Compute mu, theta_l, and theta_h via EM.')
    parser.add_argument('-f', action="store", dest="f", type=str, default='hmm-sequence.fa')
    parser.add_argument('-mu', action="store", dest="mu", type=float, default=0.05)
    parser.add_argument('-theta_h', action="store", dest="theta_h", type=float, default=0.6)
    parser.add_argument('-theta_l', action="store", dest="theta_l", type=float, default=0.4)

    args = parser.parse_args()
    sequence = read_fasta(args.f)
    mu = args.mu
    theta_h = args.theta_h
    theta_l = args.theta_l

    mu, theta_h, theta_l = train(sequence, mu, theta_h, theta_l)
    print "theta_h: %.5f\ntheta_l: %.5f\nmu: %.5f" % (theta_h, theta_l, mu)

if __name__ == "__main__":
    main()
