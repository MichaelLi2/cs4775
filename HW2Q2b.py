#!/usr/bin/env python

'''Script for computing sequence alignments using Needleman-Wunsch with
   linear gap penalties.
Arguments:
    f - FASTA file with sequences in FASTA format.
    s - JSON with the score matrix for alignment.
    d - The gap penalty for the alignment.

Outputs:
    Prints alignment to console.

Example Usage:
    python 2b.py -f sequences.fasta -s score_matrix.json -d 100
'''

import argparse
import json

'''Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    t: the traceback matrix
Returns:
    a_x: the string for the alignment of x's sequence
    a_y: the string for the alignment of y's sequence
'''
def traceback(x, y, t):
    x = '-' + x
    y = '-' + y
    lenx = len(x)
    leny = len(y)
    i = lenx
    j = leny
    a_x = ""
    a_y = ""

    while i>1 or j>1:
        if t[(j-1)*lenx + i - 1] == 'g':
            a_x = x[i-1: i] + a_x
            a_y = y[j-1: j] + a_y
            i = i-1
            j = j-1
            #print "g"
        elif t[(j-1)*lenx + i - 1] == 'u':
            a_x = '-' + a_x
            a_y = y[j-1: j] + a_y
            j = j-1
            #print "u"
        elif t[(j-1)*lenx + i - 1] == 'l':
            a_y = '-' + a_y
            a_x = x[i-1: i] + a_x
            i = i-1
            #print "l"

    if a_x[:1] == '-' and a_y[:1] == '-':
        a_x = a_x[1:]
        a_y = a_y[1:]

    return a_x, a_y


'''Computes the score and alignment of two strings.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening/extension penalty
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''
def sequence_alignment(x, y, s, d):
    lenx = len(x)+1
    leny = len(y)+1
    f = [0]*(lenx*leny) #Recurrence matrix
    t = ['n']*(lenx*leny) #Traceback matrix

    for i in range (1, lenx):
        f[i] = f[i-1] - d
        t[i] = 'l'
    for j in range (1, leny):
        f[j*lenx] = f[(j-1)*lenx] - d
        t[j*lenx] = 'u'

    for j in range (1, leny):
        for i in range (1, lenx):
            score = s[x[i-1: i]][y[j-1: j]]

            v_g = f[(j-1)*lenx + i-1] + score
            v_l = f[j*lenx + i-1] - d
            v_u = f[(j-1)*lenx + i] - d

            if v_g >= v_l and v_g >= v_u:
                f[j*lenx + i] = v_g
                t[j*lenx + i] = 'g'
            elif v_u >= v_g and v_u >= v_l:
                f[j*lenx + i] = v_u
                t[j*lenx + i] = 'u'
            else:
                f[j*lenx + i] = v_l
                t[j*lenx + i] = 'l'


    a_x, a_y = traceback(x, y, t)
    return f[len(f)-1], (a_x, a_y)


'''Prints two aligned sequences formatted for convenient inspection.
Arguments:
    a_x: the first sequence aligned
    a_y: the second sequence aligned
Outputs:
    Prints aligned sequences (80 characters per line) to console
'''
def print_alignment(a_x, a_y):
    assert len(a_x) == len(a_y), "Sequence alignment lengths must be the same."
    for i in range(1 + (len(a_x) / 80)):
        start = i * 80
        end = (i + 1) * 80
        print a_x[start:end]
        print a_y[start:end]
        print


def main():
    parser = argparse.ArgumentParser(
        description='Calculate sequence alignments for two sequences with a linear gap penalty.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d

    with open(fasta_file) as f:
        _, x, _, y = [line.strip() for line in f.readlines()]
    with open(score_matrix_file) as f:
        s = json.loads(f.readlines()[0])

    score, (a_x, a_y) = sequence_alignment(x, y, s, d)
    print "Alignment:"
    print_alignment(a_x, a_y)
    print "Score: " + str(score)

if __name__ == "__main__":
    main()
