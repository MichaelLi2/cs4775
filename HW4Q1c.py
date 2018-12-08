#!/usr/bin/env python

''' Converts reconstructed sequences of bases to corresponding amino acids and
    displays the alignment of amino acids between the human and ancestral
    sequences.
    NOTE: you can simply slot in your code for computing MAP distributions from
          1b

Arguments:
    -f: fasta file containing the multiple alignment
        (default is apoe.fa)
    -t: topology index (1-3), corresponding to the maximum likelihood topology
        (from 1a)
Outputs:
    alignment of amino acids between human and MAP root sequence

Example usage:
    python HW4Q1c.py -f HW4test.fa -t 1
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
                  branch_probs[branch_id][a][d] giving P(b | a, t_branch_id)
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


''' Computes maximum posterior distribution of bases at the root of the tree
    given the topology specified by ordering

Arguments:
    data: sequence data (dict: name of sequence owner -> sequence)
    seqlen: length of sequences
    ordering: postorder traversal of our topology
    bp: branch probabilities for the given branches: 5x4x4 matrix indexed as
        branch_probs[branch_id][a][d] giving P(b | a, t_branch_id)
Returns:
    output: maximum a posteriori (MAP) estimate of root sequence
'''
def map_estimate(data, seqlen, ordering, bp):
    map_seq = ""
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
                total_prob = 0
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
            total_prob = total_prob + acc1 * acc2
        #find sequence
        for c in range (0,4):
            probs_all[c][len(ordering)-1] = probs_all[c][len(ordering)-1] / total_prob

        if (probs_all[0][len(ordering)-1] >= probs_all[1][len(ordering)-1]
        and probs_all[0][len(ordering)-1] >= probs_all[2][len(ordering)-1]
        and probs_all[0][len(ordering)-1] >= probs_all[3][len(ordering)-1]):
            map_seq = map_seq + "A"
        elif (probs_all[1][len(ordering)-1] >= probs_all[0][len(ordering)-1]
        and probs_all[1][len(ordering)-1] >= probs_all[2][len(ordering)-1]
        and probs_all[1][len(ordering)-1] >= probs_all[3][len(ordering)-1]):
            map_seq = map_seq + "T"
        elif (probs_all[2][len(ordering)-1] >= probs_all[0][len(ordering)-1]
        and probs_all[2][len(ordering)-1] >= probs_all[1][len(ordering)-1]
        and probs_all[2][len(ordering)-1] >= probs_all[3][len(ordering)-1]):
            map_seq = map_seq + "C"
        else:
            map_seq = map_seq + "G"

    return map_seq


''' Translates DNA to amino acids for the data in the ```data``` dictionary and
    the MAP root sequence.

Arguments:
    data: dictionary of sequences (name -> sequence)
    map_output: a string with maximum a posteriori sequence at root
Returns:
    amino_data: dictionary with amino acid sequences
    amino_map: MAP root sequence translated into amino acids
'''
def translate(data, map_output):
    # NOTE: stop codon encoded as '$'
    amino = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"$", "TAG":"$",
    "TGT":"C", "TGC":"C", "TGA":"$", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

    amino_seq = ""
    for i in range (0, len(map_output)/3):
        amino_seq = amino_seq + amino[map_output[i*3:i*3+3]]

    amino_dict = data
    for n, s in data.items():
        aminos = ""
        for i in range (0, len(s)/3):
            aminos = aminos + amino[s[i*3:i*3+3]]
        amino_dict[n] = aminos

    return amino_dict, amino_seq


''' Outputs the MAP estimate and the data in a way that can easily be examined
    for comparison.

Arguments:
    data: dictionary of sequences (name -> sequence)
    map_output: a string containing MAP amino acid sequence at root
'''
def output_alignment(data, map_output, per_line = 70):
    iter = len(map_output) / per_line
    names = ['human: ']
    if (len(map_output) % per_line != 0): iter += 1

    for i in range(iter):
        if (i == iter - 1):
            print 'root:  ' + map_output[i * per_line:]
            for name in names: print name + data[name.strip()[:-1]][i * per_line:]
        else:
            print 'root:  ' + map_output[i * per_line: (i+1) * per_line]
            for name in names: print name + data[name.strip()[:-1]][i * per_line: (i+1) * per_line]
            print '\n'


def main():
    parser = argparse.ArgumentParser(
        description='Compare human and MAP ancestral amino acid sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='apoe.fa')
    parser.add_argument('-t', action="store", dest="t", type=int, required=True)
    args = parser.parse_args()
    fasta_file = args.f
    topology_index = args.t

    assert topology_index in [1, 2, 3], "Invalid topology index."
    topology_index -= 1

    data, seqlen = read_data(fasta_file)
    ordering, probs = initialize_topology(topology_index)
    map_output = map_estimate(data, seqlen, ordering, probs)
    amino_data, amino_map = translate(data, map_output)
    output_alignment(amino_data, amino_map)


if __name__ == "__main__":
    main()
