#!/usr/bin/env python

'''Script for computing sequence alignments using Needleman-Wunsch with
   affine gap penalties.
Arguments:
    f - FASTA file with sequences in FASTA format.
    s - JSON with the score matrix for alignment.
    d - The gap opening penalty for the alignment.
    e - The gap extension penalty for the alignment.

Outputs:
    Prints alignment to console.

Example Usage:
    python 2c.py -f sequences.fasta -s score_matrix.json -d 430 -e 30
'''

import argparse
import json


'''Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    t: the traceback matrix, which stores values that point to which
       prior matrix was used to reach a given location in each of the
       3 matrices.
    start: value indicating the starting matrix (that had the optimal value)
Returns:
    a_x: the string for the alignment of x's sequence
    a_y: the string for the alignment of y's sequence
'''
def traceback(x, y, t, start):
    x = '-' + x
    y = '-' + y
    lenx = len(x)
    leny = len(y)
    i = lenx
    j = leny
    a_x = ""
    a_y = ""
    matrix = start
    print start


    while i>1 or j>1:
        if matrix == 0:
            a_x = x[i-1: i] + a_x
            a_y = y[j-1: j] + a_y
        elif matrix == 1:
            a_y = '-' + a_y
            a_x = x[i-1: i] + a_x
        elif matrix == 2:
            a_x = '-' + a_x
            a_y = y[j-1: j] + a_y

        if t[(j-1)*lenx + i - 1][matrix] == 'm':
            #print 'm'
            if matrix == 0:
                i = i-1
                j = j-1
            elif matrix == 1:
                i = i-1
            else:
                j = j-1
            matrix = 0

        elif t[(j-1)*lenx + i - 1][matrix] == 'x':
            #print 'x'
            if matrix == 1:
                i = i-1
            elif matrix == 0:
                i = i-1
                j = j-1
            else:
                j = j-1
            matrix = 1

        elif t[(j-1)*lenx + i - 1][matrix] == 'y':
            #print 'y'
            if matrix == 2:
                j = j-1
            elif matrix == 0:
                i = i-1
                j = j-1
            else:
                i = i-1
            matrix = 2

    if a_x[:1] == '-' and a_y[:1] == '-':
        a_x = a_x[1:]
        a_y = a_y[1:]

    return a_x, a_y


'''Computes the score and alignment of two strings using an affine gap penalty.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening penalty
    e: the gap extension penalty
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''
def affine_sequence_alignment(x, y, s, d, e):
    lenx = len(x)+1
    leny = len(y)+1
    m = [0]*(lenx*leny) #Recurrence matrix M
    i_x = [0]*(lenx*leny) #Recurrence matrix Ix
    i_y = [0]*(lenx*leny) #Recurrence matrix Iy
    t = [('n','n','n','n','n','n')]*(lenx*leny) #Traceback matrix

    i_x[1] = -d
    i_y[1] = -float('inf')
    m[1] = -float('inf')
    t[1] = ('x','x','x')

    i_y[lenx] = -d
    i_x[lenx] = -float('inf')
    m[lenx] = -float('inf')
    t[lenx] = ('y','y','y')

    for i in range (2, lenx):
        i_x[i] = i_x[i-1] - e
        i_y[i] = -float('inf')
        m[i] = -float('inf')
        t[i] = ('x','x','x')
    for j in range (2, leny):
        i_y[j*lenx] = i_y[(j-1)*lenx] - e
        i_x[j*lenx] = -float('inf')
        m[j*lenx] = -float('inf')
        t[j*lenx] = ('y','y','y')


    for j in range (1, leny):
        for i in range (1, lenx):
            score = s[x[i-1: i]][y[j-1: j]]

            m_m = m[(j-1)*lenx + i-1] + score
            m_x = i_x[(j-1)*lenx + i-1] + score
            m_y = i_y[(j-1)*lenx + i-1] + score
            i_x_m = m[j*lenx + i-1] - d
            i_x_x = i_x[j*lenx + i-1] - e
            i_y_m = m[(j-1)*lenx + i] - d
            i_y_y = i_y[(j-1)*lenx + i] - e

            trace_m = ''
            trace_x = ''
            trace_y = ''

            #score for m
            if m_m >= m_x and m_m >= m_y:
                m[j*lenx + i] = m_m
                trace_m = 'm'
            elif m_x > m_m and m_x >= m_y:
                m[j*lenx + i] = m_x
                trace_m = 'x'
            else:
                m[j*lenx + i] = m_y
                trace_m = 'y'
            #score for x
            if i_x_m >= i_x_x:
                i_x[j*lenx + i] = i_x_m
                trace_x = 'm'
            else:
                i_x[j*lenx + i] = i_x_x
                trace_x = 'x'
            #score for y
            if i_y_m >= i_y_y:
                i_y[j*lenx + i] = i_y_m
                trace_y = 'm'
            else:
                i_y[j*lenx + i] = i_y_y
                trace_y = 'y'

            #traceback
            t[j*lenx + i] = (trace_m, trace_x, trace_y)

    start = -1
    final_score = 0.0
    score_m = m[len(m)-1]
    score_x = i_x[len(i_x)-1]
    score_y = i_y[len(i_y)-1]
    if score_m >= score_x and score_m >= score_y:
        start = 0
        final_score = score_m
    elif score_x >= score_m and score_x >= score_y:
        start = 1
        final_score = score_x
    else:
        start = 2
        final_score = score_y

    a_x, a_y = traceback(x, y, t, start)
    return final_score, (a_x, a_y)


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
        description='Calculate sequence alignments for two sequences with an affine gap penalty.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)
    parser.add_argument('-e', action="store", dest="e", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d
    e = args.e

    with open(fasta_file) as f:
        _, x, _, y = [line.strip() for line in f.readlines()]
    with open(score_matrix_file) as f:
        s = json.loads(f.readlines()[0])

    score, (a_x, a_y) = affine_sequence_alignment(x, y, s, d, e)
    print "Alignment:"
    print_alignment(a_x, a_y)
    print "Score: " + str(score)


if __name__ == "__main__":
    main()
