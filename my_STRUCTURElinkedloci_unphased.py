''' Performs STRUCUTRE algorithm on the specified set of sequences to separate population into groups

Arguments:
    -f: set of sequences to process
    -k: the number of populations desired
Outputs:
##################################################################
    DON'T FORGET TO DO THIS
##################################################################

Example Usage:
    python my_STRUCTURElinkedloci.py -f strucutre_input.txt -k 2
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
    qtemp = []
    #counts of alleles in each sample origionated in each population m[samp][pop]
    m = []
    for g in seqs:
        props = []
        propst = []
        nums = []
        for h in range(k):
            props.append(0.0)
            propst.append(0.0)
            nums.append(0.0)
        q.append(props)
        qtemp.append(propst)
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
    pis = []
    for b1 in seqs:
        samps1 = []
        for b2 in range(len(b1)):
            pops1 = []
            for b5 in range(k):
                pops1.append([0.0,0.0])
            samps1.append(pops1)
        pis.append(samps1)
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

    #rate for the linked loci model
    r = 0.001
    #likelihood of population assignments, initially set to min number
    likelihood = -999999999999999
    q_likelihood = -999999999999999
    #iterations
    for i in range(60):
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
        for r1 in range(k):
            for s in range(len(seqs[0])):
                p[r1][s] = np.random.dirichlet(n[r1][s])
        #sample q from dirchlet distribution
        new_q_likeli = 0
        for t in range(len(seqs)):
            q[t] = np.random.dirichlet(m[t])#qtemp?
            for t2 in range(k):
                new_q_likeli = new_q_likeli + np.log(q[t][t2])#qtemp?
        #metropolis hastings update for q
        #if new_q_likeli >= q_likelihood:
        #    for t3 in range(len(seqs)):
        #        for t4 in range(k):
        #            q[t3][t4] = qtemp[t3][t4]
        #else:
        #    rand = random.random()
        #    if rand < new_q_likeli/q_likelihood:
        #        for t3 in range(len(seqs)):
        #            for t4 in range(k):
        #                q[t3][t4] = qtemp[t3][t4]


        #create probability matix for alleles to sample z from
        for u in range(len(seqs)):
            #value of betas, beta[locus][pop][allele]
            beta = []
            for lcs in range(len(seqs[0])):
                betak = []
                for po in range(k):
                    betak.append([0.0,0.0])
                beta.append(betak)
            #forward step, compute all the betas
            for v in range(len(seqs[u])):
                #first loci betas calculated differently
                if v == 0:
                    for w in range(k):
                        for y in range(k):
                            #if allele in sample u, locus v is the first allele
                            if seqs[u][v:v+1] == "0":
                                beta[v][w][y] = q[u][w] * p[w][v][0] * q[u][y] * p[y][v][0]
                            #if allele in sample u, locus v is the second allele
                            elif seqs[u][v:v+1] == "2":
                                beta[v][w][y] = q[u][w] * p[w][v][1] * q[u][y] * p[y][v][1]
                            #if there is one of each at the locus, allele could be either, assume first/second
                            elif seqs[u][v:v+1] == "1":
                                beta[v][w][y] = q[u][w] * p[w][v][0] * q[u][y] * p[y][v][1]
                            #if unkown, average?
                            else:
                                beta[v][w][y] = q[u][w] * (p[w][v][1] + p[w][v][0]) * q[u][y] * (p[y][v][1] + p[y][v][0])/4
                else:
                    for w2 in range(k):
                        for y2 in range(k):
                            #if allele in sample u, locus v is the first allele
                            for ks in range(k):
                                for ks2 in range(k):
                                    if w2 == ks:
                                        p11 = exp(r*(-1)) + (1-exp(r*(-1))) * q[u][ks]
                                    else:
                                        p11 = (1-exp(r*(-1))) * q[u][ks]
                                    if w2 == ks2:
                                        p12 = exp(r*(-1)) + (1-exp(r*(-1))) * q[u][ks2]
                                    else:
                                        p12 = (1-exp(r*(-1))) * q[u][ks2]
                                    if y2 == ks:
                                        p21 = exp(r*(-1)) + (1-exp(r*(-1))) * q[u][ks]
                                    else:
                                        p21 = (1-exp(r*(-1))) * q[u][ks]
                                    if y2 == ks2:
                                        p22 = exp(r*(-1)) + (1-exp(r*(-1))) * q[u][ks2]
                                    else:
                                        p22 = (1-exp(r*(-1))) * q[u][ks2]
                                    if seqs[u][v:v+1] == "0":
                                        beta[v][w2][y2] = beta[v][w2][y2]+p[w2][v][0]*p[y2][v][0]*beta[v-1][ks][ks2]*0.5*(p11*p22+p21*p12)
                                    #if allele in sample u, locus v is the second allele
                                    elif seqs[u][v:v+1] == "2":
                                        beta[v][w2][y2] = beta[v][w2][y2]+p[w2][v][1]*p[y2][v][1]*beta[v-1][ks][ks2]*0.5*(p11*p22+p21*p12)
                                    #if there is one of each at the locus, allele could be either, assume first?
                                    elif seqs[u][v:v+1] == "1":
                                        beta[v][w2][y2] = beta[v][w2][y2]+p[w2][v][0]*p[y2][v][1]*beta[v-1][ks][ks2]*0.5*(p11*p22+p21*p12)
                                    #if unkown, average?
                                    else:
                                        beta[v][w2][y2] = beta[v][w2][y2]+((p[ks][v][1]+p[ks][v][0])/2)*beta[v-1][ks][ks2]*((p[ks2][v][1]+p[ks2][v][0])/2)*0.5*(p11*p22+p12*p21)

            #backwards caluculation for the probabilities to sample z from
            #sample z at the same time since value used for one locus is used in next iteration
            #final loci probabilities are just beta_locus,population
            new_likelihood = 0.0
            total_pi = 0.0
            #print beta[len(seqs[u])-1]
            for w3 in range(k):
                for w3w in range(k):
                    pis[u][len(seqs[u])-1][w3][w3w] = beta[len(seqs[u])-1][w3][w3w]
                    total_pi = total_pi + beta[len(seqs[u])-1][w3][w3w]
            for x in range(k):
                for xx in range(k):
                    pis[u][len(seqs[u])-1][x][xx] = pis[u][len(seqs[u])-1][x][xx] / total_pi
            #print pis[u][len(seqs[u])-1]
            #sample z for last loci
            ran1 = random.random()
            flag = 0
            for bb in range(k):
                if flag == 1: break
                for cc in range(k):
                    ran1 = ran1 - pis[u][len(seqs[u])-1][bb][cc]
                    if ran1 < 0:
                        z1temp[u][len(seqs[u])-1] = bb
                        z2temp[u][len(seqs[u])-1] = cc
                        new_likelihood = new_likelihood + np.log(pis[u][len(seqs[u])-1][bb][cc])
                        flag = 1
                        break
            #calculate probabilities for rest of loci
            for v2 in range(len(seqs[u]) - 2, -1, -1):
                total_pi = 0.0
                for w4 in range(k):
                    for w4w in range(k):
                        if z1temp[u][v2+1] == w4:#z1temp?
                            p11 = exp(r*(-1)) + (1-exp(r*(-1)))*q[u][w4]
                        else:
                            p11 = (1-exp(r*(-1))) * q[u][w4]
                        if z1temp[u][v2+1] == w4w:
                            p12 = exp(r*(-1)) + (1-exp(r*(-1)))*q[u][w4w]
                        else:
                            p12 = (1-exp(r*(-1))) * q[u][w4w]
                        if z2temp[u][v2+1] == w4:
                            p21 = exp(r*(-1)) + (1-exp(r*(-1)))*q[u][w4]
                        else:
                            p21 = (1-exp(r*(-1))) * q[u][w4]
                        if z2temp[u][v2+1] == w4w:
                            p22 = exp(r*(-1)) + (1-exp(r*(-1)))*q[u][w4w]
                        else:
                            p22 = (1-exp(r*(-1))) * q[u][w4w]
                        pis[u][v2][w4][w4w] = beta[v2][w4][w4w]*0.5*(p22*p11+p21*p12)
                        total_pi = total_pi + beta[v2][w4][w4w]*0.5*(p22*p11+p21*p12)
                for x2 in range(k):
                    for xx2 in range(k):
                        pis[u][v2][x2][xx2] = pis[u][v2][x2][xx2] / total_pi
                #print pis[u][v2]
                #sample z for rest of loci
                ran2 = random.random()
                flag = 0
                for gg in range(k):
                    if flag == 1: break
                    for ff in range(k):
                        ran2 = ran2 - pis[u][v2][gg][ff]
                        if ran2 < 0:
                            z1temp[u][v2] = gg
                            z2temp[u][v2] = ff
                            new_likelihood = new_likelihood + np.log(pis[u][v2][gg][ff])
                            flag = 1
                            break

        #if likelihood is iproved, update z
        if new_likelihood > likelihood:
            print new_likelihood
            likelihood = new_likelihood
            for dd in range(len(seqs)):
                for ee in range(len(seqs[dd])):
                    z1[dd][ee] = z1temp[dd][ee]
                    z2[dd][ee] = z2temp[dd][ee]

    #print beta
    #print "_____________"
    #print pis[0]
    #print pis[80]
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
            total2 = total2 + 1
        else:
            total1 = total1 + 1
    print "normal from pop1: " + str(total1)
    print "normal from pop2: " + str(total2)


if __name__ == '__main__':
    main()
