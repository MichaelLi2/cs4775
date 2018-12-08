''' Performs STRUCUTRE algorithm on the specified set of sequences to separate population into groups

Arguments:
    -f: set of sequences to process
    -k: the number of populations desired
Outputs:
##################################################################
    DON'T FORGET TO DO THIS
##################################################################

Example Usage:
    python my_STRUCTURE.py -f strucutre_input.txt -k 2
'''

import numpy as np
import argparse
import random
from numpy import log1p
from math import exp

''' Reads in the allele sequence from text files.

Arguments:
    filename: which filename to read sequences from
Returns:
    output: list of sequences in allele file
    ids: list of identification for each sequence corresponding to output list
'''
def read_allele(filename):
    with open(filename, "rb") as f:
        output = []
        ids = []
        s = ""
        first = 0
        for l in f.readlines():
            #skip first line
            if first == 0:
                first = 1
                continue
            # input each line
            else:
                #isolate IDs
                name = l.split(' ', 1)[0]
                ids.append(name)
                #print name
                #add sequence to output list
                output.append(''.join(l.split())[len(name):])#l.replace("" "", "")[len(name):])
                #print ''.join(l.split())[len(name):]
        return output, ids

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

''' Performs STRUCUTRE to predict the populations structure for the
 sequences provided.

Arguments:
    seqs: list of sequences
    name: list of names associated with sequences
    k: number of population structures to separate into
Returns:
    pstructure: list of list of sequences in each population structure (list of k
    populations)
'''
def struct(seqs, name, k):
    pstructure = []
    for a in seqs:
        pstructure.append([])
    #ancestral population matrix for each allele z[sample][locus]
    z1 = []
    z2 = []
    #matrix of allele frequency at each locus in each population p[pop][locus][allele]
    p = []
    #counts of alleles at each locus in each population n[pop][locus][allele]
    n = []
    for e in range(k):
        freqs = []
        nums = []
        for f in range(len(seqs[0])):
            freqs.append([0.0,0.0])
            nums.append([0.0,0.0])
        p.append(freqs)
        n.append(nums)
    #list of admixture proportions for each sample q[samp][pop]
    q = []
    #counts of alleles in each sample origionated in each population m[samp][pop]
    m = []
    for g in seqs:
        props = []
        nums = []
        for h in range(k):
            props.append(0.0)
            nums.append(0.0)
        q.append(props)
        m.append(nums)

    #random initialization of Z
    for b in seqs:
        rand_init1 = []
        rand_init2 = []
        for c in range(len(b)):
            ran1 = random.random()
            for d in range(k):
                ran1 = ran1 - 1.0/k
                if ran1 < 0:
                    rand_init1.append(d)
                    break
        z1.append(rand_init1)
        for c2 in range(len(b)):
            ran2 = random.random()
            for d2 in range(k):
                ran2 = ran2 - 1.0/k
                if ran2 < 0:
                    rand_init2.append(d2)
                    break
        z2.append(rand_init2)

    #probability matrix for z, pis[sample][locus][populations]
    pis1 = []
    pis2 = []
    for b1 in seqs:
        samps1 = []
        samps2 = []
        for b2 in range(len(b1)):
            pops1 = []
            pops2 = []
            for b5 in range(k):
                pops1.append(0.0)
                pops2.append(0.0)
            samps1.append(pops1)
            samps2.append(pops2)
        pis1.append(samps1)
        pis2.append(samps2)
    #temporary z values
    z1temp = []
    z2temp = []
    for b3 in seqs:
        temp1 = []
        temp2 = []
        for b4 in range(len(b)):
            temp1.append(0)
            temp2.append(0)
        z1temp.append(temp1)
        z2temp.append(temp2)

    likelihood = -999999999999999
    #iterations
    for i in range(300):
        if i%10 == 0: print i
        #zero out n and m
        for i1 in range(k):
            for i2 in range(len(seqs[0])):
                n[i1][i2] = [0.5,0.5] #lambda
            for i3 in range(len(seqs)):
                m[i3][i1] = 1 #alpha

        #figure out counts for n and m
        for j in range(len(seqs)):
            for l in range(len(seqs[j])):
                #sample j has 1 more allele from populations z1[j][l] and z2[j][l]
                m[j][z1[j][l]] = m[j][z1[j][l]] + 1
                m[j][z2[j][l]] = m[j][z2[j][l]] + 1
                #if both 1st allele, count of first allele at locus l and in populations
                #z1[j][l] and z2[j][l] is increased by 1
                if seqs[j][l:l+1] == "0":
                    n[z1[j][l]][l][0] = n[z1[j][l]][l][0] + 1
                    n[z2[j][l]][l][0] = n[z2[j][l]][l][0] + 1
                #if one 1st allele, count of first allele at locus l and in populations
                #z1[j][l] and z2[j][l] is increased by 1, assuming first?
                elif seqs[j][l:l+1] == "1":
                    n[z1[j][l]][l][0] = n[z1[j][l]][l][0] + 1
                    #n[z2[j][l]][l][0] = n[z2[j][l]][l][0] + 0.5
                    #n[z1[j][l]][l][1] = n[z1[j][l]][l][1] + 0.5
                    n[z2[j][l]][l][1] = n[z2[j][l]][l][1] + 1
                #else count of first allele at locus l and in populations
                #z1[j][l] and z2[j][l] is increased by 1
                elif seqs[j][l:l+1] == "2":
                    n[z1[j][l]][l][1] = n[z1[j][l]][l][1] + 1
                    n[z2[j][l]][l][1] = n[z2[j][l]][l][1] + 1
                #otherwise data might be 9 but that just means missing info

        #sample p from dirchlet distribution
        for r in range(k):
            for s in range(len(seqs[0])):
                p[r][s] = np.random.dirichlet(n[r][s])
        #sample q from dirchlet distribution
        for t in range(len(seqs)):
            q[t] = np.random.dirichlet(m[t])
        #print q
        #print p

        #create probability matix for alleles to sample z from
        for u in range(len(seqs)):
            for v in range(len(seqs[u])):
                total_pi1 = 0.0
                total_pi2 = 0.0
                for w in range(k):
                    #if allele in sample u, locus v is the first allele
                    if seqs[u][v] == "0":
                        pis1[u][v][w] = q[u][w] * p[w][v][0]
                        total_pi1 = total_pi1 + q[u][w] * p[w][v][0]
                        pis2[u][v][w] = q[u][w] * p[w][v][0]
                        total_pi2 = total_pi2 + q[u][w] * p[w][v][0]
                    #if allele in sample u, locus v is the second allele
                    elif seqs[u][v] == "2":
                        pis1[u][v][w] = q[u][w] * p[w][v][1]
                        total_pi1 = total_pi1 + q[u][w] * p[w][v][1]
                        pis2[u][v][w] = q[u][w] * p[w][v][1]
                        total_pi2 = total_pi2 + q[u][w] * p[w][v][1]
                    #if there is one of each at the locus, allele could be either, assume first?
                    elif seqs[u][v] == "1":
                        pis1[u][v][w] = q[u][w] * p[w][v][0]
                        total_pi1 = total_pi1 + q[u][w] * p[w][v][0]
                        pis2[u][v][w] = q[u][w] * p[w][v][1]
                        total_pi2 = total_pi2 + q[u][w] * p[w][v][1]
                    #if unkown, average?
                    else:
                        pis1[u][v][w] = q[u][w] * (p[w][v][1] + p[w][v][0])/2
                        total_pi1 = total_pi1 + q[u][w] * (p[w][v][1] + p[w][v][0])/2
                        pis2[u][v][w] = q[u][w] * (p[w][v][1] + p[w][v][0])/2
                        total_pi2 = total_pi2 + q[u][w] * (p[w][v][1] + p[w][v][0])/2

                for x in range(k):
                    pis1[u][v][x] = pis1[u][v][x] / total_pi1
                    pis2[u][v][x] = pis2[u][v][x] / total_pi2

        #sample z from pi and calculate likelihood
        new_likelihood = 0.0
        for y in range(len(seqs)):
            for aa in range(len(seqs[y])):
                ran1 = random.random()
                for bb in range(k):
                    ran1 = ran1 - pis1[y][aa][bb]
                    if ran1 < 0:
                        #z1[y][aa] = bb
                        z1temp[y][aa] = bb
                        new_likelihood = new_likelihood + np.log(pis1[y][aa][bb])
                        break
                ran2 = random.random()
                for cc in range(k):
                    ran2 = ran2 - pis2[y][aa][cc]
                    if ran2 < 0:
                        #z2[y][aa] = cc
                        z2temp[y][aa] = cc
                        new_likelihood = new_likelihood + np.log(pis2[y][aa][cc])
                        break

        #if likelihood is iproved, update z
        if new_likelihood > likelihood:
            print new_likelihood
            likelihood = new_likelihood
            for dd in range(len(seqs)):
                for ee in range(len(seqs[y])):
                    z1[dd][ee] = z1temp[dd][ee]
                    z2[dd][ee] = z2temp[dd][ee]


    #print pis2
    #print p[0]
    #print "----------------------"
    #print p[1]
    return z1,z2, q

def main():
    parser = argparse.ArgumentParser(description='something')
    parser.add_argument('-f', action="store", dest="f", type=str, default='structure_input.txt')
    parser.add_argument('-k', action="store", dest="k", type=int, default=10)

    args = parser.parse_args()
    sequences, names = read_allele(args.f)
    k = args.k

    structure,structure2, proportions = struct(sequences, names, k)
    total1 = 0
    total2 = 0
    for i in range(len(structure)):
        if i == 48:
            print "dwarf from pop1: " + str(total1)
            print "dwarf from pop2: " + str(total2)
            print "--------------------------"
            total1 = 0
            total2 = 0
        print proportions[i]
        if proportions[i][0] < 0.5:
            total1 = total1 + 1
        else:
            total2 = total2 + 1
    print "normal from pop1: " + str(total1)
    print "normal from pop2: " + str(total2)


if __name__ == '__main__':
    main()
