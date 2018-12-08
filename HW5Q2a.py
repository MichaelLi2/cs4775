#!/usr/bin/env python

''' Performs Gibbs sampling on the specified set of sequences to identify motifs

Arguments:
    -f: set of sequences to identify a motif in
    -k: the length of the desired motif to find
Outputs:
    - a list of start positions for each sequence, along with the corresponding
      motif identified
    - the consensus motif (based on majority base at each position)
    - the position weight matrix

Example Usage:
    python HW5Q2a.py -f motif1.fa -k 10
'''

import argparse
import numpy as np
from random import randint
import random
from numpy import log1p
from math import exp


''' Reads in the sequences from the motif files.

Arguments:
    filename: which filename to read sequences from
Returns:
    output: list of sequences in motif file
'''
def read_fasta(filename):
    with open(filename, "rb") as f:
        output = []
        s = ""
        for l in f.readlines():
            if l.strip()[0] == ">":
                # skip the line that begins with ">"
                if s == "": continue
                output.append(s)
                s = ""
            # keep appending to the current sequence when the line doesn't begin
            # with ">"
            else: s += l.strip()
        output.append(s)
        return output


''' Returns the majority vote (consensus) motif sequence.

Arguments:
    motif_sequences: list of all motif instances (list of strings)
Returns:
    majority sequence: consensus motif sequence (string)
'''
def majority(motif_sequences):
    base_ordering = {c: i for i, c in enumerate("ACGT")}
    freqs = np.zeros((len(motif_sequences[0]), 4))
    for m_s in motif_sequences:
        for i, b in enumerate(m_s): freqs[i, base_ordering[b]] += 1
    return ''.join(['ACGT'[np.argmax(row)] for row in freqs])


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

''' Performs Gibbs sampling to estimate the start positions of the motif in
    each sequence and the motif model.

Arguments:
    sequences: list of sequences to find motifs in
    k: length of motif
    epsilon: pseudocounts for each base
Returns:
    starts: list of start positions for each sequence
    pi: corresponding motif model (incorporating pseudocounts).
'''
def gibbs_sampling(sequences, k = 10, epsilon = 0.5):
    #random start sites
    starts = []
    for b in range(0, len(sequences)):
        starts.append(randint(0, len(sequences[b])-k+1))

    starts[0]=421
    starts[1]=858
    starts[2]=1030
    starts[4]=358
    starts[5]=615
    starts[8]=770

    #pi matrix
    probs = []
    for a in range(0,k):
        probs.append([0.0,0.0,0.0,0.0])

    #likelihood matrix
    likes = []
    for c in range(0, len(sequences)):
        likes.append([0.0]*(len(sequences[c])-k+1))

    #new pi matrix
    new_probs = []
    for m in range(0,k):
        new_probs.append([0.0,0.0,0.0,0.0])
    #new start matrix
    new_starts = [0]*len(sequences)

    old_like = -9999999999999

    #iterations
    for z in range(0, 10):
        #if z%100 ==0
        #    print z
        for j in range (0, len(sequences)):
            #pi matrix calculation
            for i in range(0,k):
                num_A = epsilon
                num_T = epsilon
                num_C = epsilon
                num_G = epsilon
                for y in range (0, len(sequences)):
                    if j != y:
                        if sequences[y][starts[y]+i:starts[y]+i+1] == 'A':
                            num_A += 1
                        elif sequences[y][starts[y]+i:starts[y]+i+1] == 'T':
                            num_T += 1
                        elif sequences[y][starts[y]+i:starts[y]+i+1] == 'C':
                            num_C += 1
                        else:
                            num_G += 1
                probs[i][0] = num_A / (len(sequences) - 1 + 4*epsilon)
                probs[i][1] = num_T / (len(sequences) - 1 + 4*epsilon)
                probs[i][2] = num_C / (len(sequences) - 1 + 4*epsilon)
                probs[i][3] = num_G / (len(sequences) - 1 + 4*epsilon)

            #likelihood at each start position
            for r in range(0, len(sequences[j])-k+1):
                for s in range(r, r+k):
                    if sequences[j][s:s+1] == 'A':
                        likes[j][r] += np.log(probs[s-r][0]) - np.log(0.25)
                    elif sequences[j][s:s+1] == 'T':
                        likes[j][r] += np.log(probs[s-r][1]) - np.log(0.25)
                    elif sequences[j][s:s+1] == 'C':
                        likes[j][r] += np.log(probs[s-r][2]) - np.log(0.25)
                    else:
                        likes[j][r] += np.log(probs[s-r][3]) - np.log(0.25)
            total_like = likes[j][0]
            for w in range(1, len(sequences[j])-k+1):
                total_like = sumLogProb(total_like, likes[j][w])
            for t in range(0, len(sequences[j])-k+1):
                likes[j][t] = exp(likes[j][t] - total_like)

            #random sampling
            ran = random.random()
            count = 0
            while ran > 0:
                #print ran
                ran = ran-likes[j][count]
                count += 1
            new_starts[j] = count
        #new pi matrix calculation
        for f in range(0,k):
            new_num_A = epsilon
            new_num_T = epsilon
            new_num_C = epsilon
            new_num_G = epsilon
            for n in range (0, len(sequences)):
                if sequences[n][new_starts[n]+f:new_starts[n]+f+1] == 'A':
                    new_num_A += 1
                elif sequences[n][new_starts[n]+f:new_starts[n]+f+1] == 'T':
                    new_num_T += 1
                elif sequences[n][new_starts[n]+f:new_starts[n]+f+1] == 'C':
                    new_num_C += 1
                else:
                    new_num_G += 1
            new_probs[f][0] = new_num_A / (len(sequences) + 4*epsilon)
            new_probs[f][1] = new_num_T / (len(sequences) + 4*epsilon)
            new_probs[f][2] = new_num_C / (len(sequences) + 4*epsilon)
            new_probs[f][3] = new_num_G / (len(sequences) + 4*epsilon)

        #new likelihood calculation
        new_like = 0.0
        for p in range(0, len(sequences)):
            for q in range(0, k):
                if sequences[p][new_starts[p]+q:new_starts[p]+q+1] == 'A':
                    new_like += np.log(new_probs[q][0]) - np.log(0.25)
                elif sequences[p][new_starts[p]+q:new_starts[p]+q+1] == 'T':
                    new_like += np.log(new_probs[q][1]) - np.log(0.25)
                elif sequences[p][new_starts[p]+q:new_starts[p]+q+1] == 'C':
                    new_like += np.log(new_probs[q][2]) - np.log(0.25)
                else:
                    new_like += np.log(new_probs[q][3]) - np.log(0.25)
        #update start position if new likelihood is better
        if new_like > old_like:
            old_like = new_like
            for g in range(0, len(starts)):
                starts[g] = new_starts[g]

    #print starts
    return starts, probs


def main():
    parser = argparse.ArgumentParser(description='Estimate the start positions and motif model via Gibbs sampling.')
    parser.add_argument('-f', action="store", dest="f", type=str, default='motif1.fa')
    parser.add_argument('-k', action="store", dest="k", type=int, default=10)

    args = parser.parse_args()
    sequences = read_fasta(args.f)
    k = args.k
    ''' Must define an epsilon value. '''
    epsilon = 0.5

    starts, pi = gibbs_sampling(sequences, k=k, epsilon=epsilon)
    motif_sequences = [s[starts[i]:starts[i]+k] for i, s in enumerate(sequences)]
    print '\n'.join([("Sequence %d: " % i) + m for i, m in enumerate(motif_sequences)])
    print "Consensus motif: %s" % majority(motif_sequences)
    print pi



if __name__ == '__main__':
    main()
