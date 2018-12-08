#!/usr/bin/env python

'''Script for computing h_{Y_n} given n and the probabilities of observing 1-6.
Arguments:
    n - integer representing the number of observations
    p - sequence of 6 float values summing to 1 representing pi_i, i in [1,6]

Outputs:
    h_<n>_dp.csv - a file the probabilities h_{Y_n}(y), y in [n,6n] in CSV
                   format, for use in later parts of question 1.
Example Usage:
    python 1c.py -n 50 -p .1 .2 .2 .2 .1 .2

    This would output a file called h_dp_50.csv.
'''

import argparse
import numpy as np


'''Computes h_{Y_n}(y) for all y in {n, ..., 6n} for a given n and pi.
Arguments:
    n - the number of random variables being summed
    pi - the probabilities of 1-6 for a single observation

Returns:
    a vector of probabilities for y in {n, ..., 6n}.
        index 0 should correspond to y = n, index 1 to n+1, etc.
'''
def h_Y(n, pi):
    probs = pi
    for i in range (2,n+1):
        new_prob = []
        for j in range (i, 6*i+1):
            p = 0
            for k in range (1, 7):
                if j-k > i-2 and j-k < len(probs)+i-1:
                    p += probs[j-k-i+1]*pi[k-1]
            new_prob.append(p)
        probs = new_prob
    return probs


''' Returns the minimum ten probabilities of a given array.
Arguments:
    n - the number of random variables being summed
    probs - the probabilities of [n, 6n]
Returns:
    a vector of the values of y in [n, 6n] that correspond to the minimum 10
        probabilities in probs priority is given to higher indices in case of
        ties.
'''
def min10(n, probs):
    arr = np.array(probs)
    return arr.argsort()[:10] + n


''' Returns the maximum ten probabilities of a given array.
Arguments:
    n - the number of random variables being summed
    probs - the probabilities of [n, 6n]
Returns:
    a vector of the values of y in [n, 6n] that correspond to the minimum 10
        probabilities in probs priority is given to higher indices in case of
        ties.
'''
def max10(n, probs):
    arr = np.array(probs)
    return np.fliplr([arr.argsort()[-10:] + n])[0]


''' Computes the mean and variance of an array of probabilities for [n, 6n].
Arguments:
    n - the value of n defining the range [n, 6n]
    probs - the array of probabilities
Returns:
    expect - the expectation
    variance - the variance
'''
def compute_stats(n, probs):
    # redo with correct stuf
    sum_ex = 0.0
    sum_exsq = 0.0
    num = n
    for i in probs:
        sum_ex = sum_ex + i*num
        sum_exsq = sum_exsq + i*num*num
        num = num+1
    expect = sum_ex

    variance = sum_exsq - sum_ex**2

    return expect, variance


def main():
    parser = argparse.ArgumentParser(
        description='Calculate discrete convolutions of sum of three RVs,'
        ' each ranging from 1 to 6.')
    parser.add_argument('-n', action="store", dest="n", type=int, required=True)
    parser.add_argument("-pi", action="store", dest = "pi", nargs=6,
        metavar=('p1', 'p2', 'p3', 'p4', 'p5', 'p6'),
        help="The probabilities of observing 1 through 6",
        type=float, required=True)

    args = parser.parse_args()
    n = args.n
    pi = args.pi
    assert(sum(pi) == 1.0)

    h = h_Y(n, pi)
    #for y in h:
    #    print (y)
    err = 1.0 - sum(h)
    assert(err * err <= 10 ** -10)
    min10_values = min10(n, h)
    max10_values = max10(n, h)

    print "Min probabilities:"
    for y in min10_values:
        print (str(y) + " %.4g" % h[y - n])
    print "\nMax probabilities:"
    for y in max10_values:
        print (str(y) + " %.4g" % h[y - n])
    with open("h_%d_dp.csv" % n, "wb") as fil:
        for i in range(n, n + len(h)):
            s = str(i) + ",%.4e\n" % h[i-n]
            fil.write(s)
    mean, var = compute_stats(n, h)
    print "\nMean is %.4g" % mean
    print "Variance is %.4g" % var


if __name__ == "__main__":
    main()
