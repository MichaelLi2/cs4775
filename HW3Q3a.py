'''Script for Question 3a.
Arguments:

Outputs:
    Generates 50 sequencess according to the HMM and prints the list out, along
    with the corresponding state paths and likelihoods.  Then computes the
    mean and variance of the number of transitions

Example Usage:
    python HW3Q3a.py > sequences3a.txt
'''

import numpy as np
from random import random

''' Computes the mean and variance of an array of numbers.
Arguments:
    ar - the array of numbers
Returns:
    expect - the expectation
    variance - the variance
'''
def compute_stats(ar):
    sum_ex = 0.0
    sum_exsq = 0.0
    for i in ar:
        sum_ex = sum_ex + i
        sum_exsq = sum_exsq + i**2

    expect = sum_ex/len(ar)
    variance = sum_exsq/len(ar) - expect**2

    return expect, variance



def main():
    seqs = [""]*(50)
    hidden_states = [""]*(50)
    transitions = [0]*(50)
    likelihoods = [0.0]*(50)

    for j in range (1, 51):
        state = 'n'
        if random() < 0.5:
            state = 'h'
        else:
            state = 'l'

        seq = ""
        states = state
        transi = 0
        likeli = np.log(0.5)
        for i in range (1, 1001):
            #emission
            ran = random()
            if state == 'h':
                if ran < 0.13:
                    seq = seq + "A"
                    likeli = likeli + np.log(0.13)
                elif ran < 0.5:
                    seq = seq + "C"
                    likeli = likeli + np.log(0.37)
                elif ran < 0.87:
                    seq = seq + "G"
                    likeli = likeli + np.log(0.13)
                else:
                    seq = seq + "T"
                    likeli = likeli + np.log(0.13)
            else:
                if ran < 0.32:
                    seq = seq + "A"
                    likeli = likeli + np.log(0.32)
                elif ran < 0.5:
                    seq = seq + "C"
                    likeli = likeli + np.log(0.18)
                elif ran < 0.68:
                    seq = seq + "G"
                    likeli = likeli + np.log(0.18)
                else:
                    seq = seq + "T"
                    likeli = likeli + np.log(0.32)
            #transition
            if random() < 0.01 and i != 1000:
                transi = transi+1
                likeli = likeli + np.log(0.01)
                if state == 'h':
                    state = 'l'
                else:
                    state = 'h'
            else:
                likeli = likeli + np.log(0.99)
            states = states + state
        #print seq
        #print states
        seqs[j-1] = seq
        hidden_states[j-1] = states
        transitions[j-1] = transi
        likelihoods[j-1] = likeli

    print "Sequences:"
    print seqs
    print "States:"
    print states
    print "Likelihoods:"
    print likelihoods
    print "Mean, Variance:"
    print compute_stats(transitions)
    #end


if __name__ == "__main__":
    main()
