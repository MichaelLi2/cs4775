#!/usr/bin/env python

''' Find the maximum likelihood tree for four sequences on the Jukes-Cantor
    model using Felsenstein's algorithm.

Arguments:
    -f: fasta file containing the multiple alignment
        (default is apoe.fa)
Outputs:
    likelihoods for three topologies in Figure 1 of pset4

Example Usage:
    python HW4Q1a.py -f hw4test.fa
'''

import argparse
import numpy as np
from math import exp


class Node():
    ''' Initializes a node with given parameters.

    Arguments:
        name: name of node (only relevant for leaves)
        left: left child (Node)
        right: right child (Node)
        branch_length: length of branch that leads to this node (float)
        branch_id: id of branch that leads to this node (int)
        probs: probability of observed bases beneath this node
                [list of 4 probs for 'ACGT'] (initialized to None]
    '''
    def __init__(self, name, left, right, branch_length, branch_id):
        self.name = name
        self.left = left
        self.right = right
        self.branch_length = branch_length
        self.branch_id = branch_id
        self.probs = [None for _ in range(4)]


''' Reads data from ```filename``` in fasta format.

Arguments:
    filename: name of fasta file to read
Returns:
    sequences: dictionary of outputs (string (sequence id) -> sequence (string))
    size: length of each sequence
'''
def read_data(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        sequences = {}
        output = ""
        size = 0
        curr = ""
        flag = False
        for l in lines:
            if l[0] == ">":
                if (len(output) != 0):
                    sequences[curr] = output
                    size = len(output)
                    output = ""
                curr = l[2:].strip()
            else:
                output += l.strip()
        sequences[curr] = output
    return sequences, size


''' Evaluates P(b|a, t) under the Jukes-Cantor model

Arguments:
    b: descendant base (string)
    a: ancestral base (string)
    t: branch length (float)
    u: mutation rate (float, defaults to 1)
Returns:
    prob: float probability P(b|a, t)
'''
def jcm(b, a, t, u = 1.0):
    prob = 0.0
    if b == a:
        prob = (1 + 3*exp(-4 * u * t / 3))/4
    else:
        prob = (1 - exp(-4 * u * t / 3))/4

    return prob


''' Constructs the ordering of the post-order traversal of ```index```
    topology from the pset.
Arguments:
    index: which topology to use
Returns:
    list of Nodes corresponding to post-order traversal of the topology
    branch_probs: 6x4x4 matrices, indexed as:
                  branch_probs[branch_id][a][b] giving P(b | a, t_branch_id)
'''
def initialize_topology(index):
    branch_lengths = np.array(
        [[0.07517, 0.03059, 0.03161, 0.11761, 0.14289],
        [0.20843, 0.03397, 0.03497, 0.24952, 0.00000],
        [0.20843, 0.03397, 0.03497, 0.24952, 0.00000]], dtype = float)

    names = ['human', 'mouse', 'rat', 'dog']
    bases = 'ACGT'
    branches = [0, 1, 2, 3]
    leaves = [Node(s, None, None, bl, i) for (s, i, bl) in
                zip(names, branches, branch_lengths[index, :])]
    ordering = None
    branch_probs = [np.zeros((4,4), dtype = float) for _ in range(6)]
    # Note that branch 5 (or 6 in 1-index) is the branch of 0-length
    if (index == 0):
        hdp = Node(None, leaves[0], leaves[3], 0, 5)
        mrp = Node(None, leaves[1], leaves[2], branch_lengths[index, 4], 4)
        root = Node('root', hdp, mrp, None, None)
        ordering = [leaves[0], leaves[3], hdp, leaves[1], leaves[2], mrp, root]

    elif (index == 1):
        hmp = Node(None, leaves[0], leaves[1], 0, 5)
        rdp = Node(None, leaves[2], leaves[3], branch_lengths[index, 4], 4)
        root = Node('root', hmp, rdp, None, None)
        ordering = [leaves[0], leaves[1], hmp, leaves[2], leaves[3], rdp, root]

    else:
        mdp = Node(None, leaves[1], leaves[3], 0, 5)
        hrp = Node(None, leaves[0], leaves[2], branch_lengths[index, 4], 4)
        root = Node('root', mdp, hrp, None, None)
        ordering = [leaves[1], leaves[3], mdp, leaves[0], leaves[2], hrp, root]

    #Assign 6x4x4 branch_probs values: branch_probs[branch_id][ancestor_base][descendant_base] '''
    for i in range (0,5):
        for j in range (0,4):
            for k in range (0,4):
                branch_probs[i][j][k] = jcm(bases[k:k+1], bases[j:j+1], branch_lengths[index, i], 1.0)
    for j in range (0,4):
        for k in range (0,4):
            branch_probs[5][j][k] = jcm(bases[k:k+1], bases[j:j+1], 0.0, 1.0)

    return ordering, branch_probs


''' Computes the likelihood of the data given the topology specified by ordering

Arguments:
    data: sequence data (dict: name of sequence owner -> sequence)
    seqlen: length of sequences
    ordering: postorder traversal of our topology
    bp: branch probabilities for the given branches: 5x4x4 matrix indexed as
        branch_probs[branch_id][a][d] giving P(b | a, t_branch_id)
Returns:
    total_log_prob: log likelihood of the topology given the sequence data
'''
def likelihood(data, seqlen, ordering, bp):
    prob = [0.0]*(seqlen)
    for j in range(0, seqlen):
        probs_all = [[0.0]*(len(ordering)), [0.0]*(len(ordering)), [0.0]*(len(ordering)), [0.0]*(len(ordering))]
        for i in range(0, len(ordering)-1):
            if ordering[i].left == None and ordering[i].right == None:
                seq = data[ordering[i].name]
                if seq[j:j+1] == "A":
                    probs_all[0][i] = 1
                elif seq[j:j+1] == "T":
                    probs_all[1][i] = 1
                elif seq[j:j+1] == "C":
                    probs_all[2][i] = 1
                else:
                    probs_all[3][i] = 1
            else:
                for a in range (0,4):
                    acc1 = 0
                    acc2 = 0
                    for b in range (0,4):
                        acc1 = acc1 + probs_all[b][i-1] * bp[ordering[i-1].branch_id][a][b]
                        acc2 = acc2 + probs_all[b][i-2] * bp[ordering[i-2].branch_id][a][b]
                    probs_all[a][i] = acc1 * acc2
        #calculate root probs and total prob
        total_prob = 0
        for a in range (0,4):
            acc1 = 0
            acc2 = 0
            for b in range (0,4):
                acc1 = acc1 + probs_all[b][len(ordering)-2] * bp[ordering[len(ordering)-2].branch_id][a][b]
                acc2 = acc2 + probs_all[b][2] * bp[ordering[2].branch_id][a][b]
            probs_all[a][len(ordering)-1] = acc1 * acc2
            total_prob = total_prob + probs_all[a][len(ordering)-1]
        prob[j] = total_prob/4

    log_prob = 0.0
    for p in prob:
        log_prob = log_prob + np.log(p)

    return log_prob



def main():
    parser = argparse.ArgumentParser(
        description='Compute the likelihood of the three topologies in Figure 1 given the alignment.')
    parser.add_argument('-f', action="store", dest="f", type=str, default='apoe.fa')

    args = parser.parse_args()
    fasta_file = args.f

    data, seqlen = read_data(fasta_file)
    for i in range(3):
        ordering, probs = initialize_topology(i)
        the_likelihood = likelihood(data, seqlen, ordering, probs)
        print "Log likelihood for topology %d is %.2f" % (i+1, the_likelihood)


if __name__ == "__main__":
    main()
