#!/usr/bin/env python

''' Outputs the Newick formatted tree after performing the neighbor-joining
    algorithm on an arbitrary number of species.

Arguments:
    -f: distances file (symmetric matrix with 0 on the diagonal)
        (default is dist10.txt)
Outputs:
    Newick formatted tree after neighbor-joining

Example usage:
    python HW4Q2a.py -f dist10.txt
'''

import argparse


''' Reads the input file of distances between the sequences

Arguments:
    distances_file: file name of distances between sequences
Returns:
    D: matrix of distances (map of maps)
    mapping: index to name mapping (dictionary)
'''
def read_data(distances_file):
    with open(distances_file, "rb") as f:
        lines = [l.strip().split() for l in f.readlines()]
        mapping = {i: s for i, s in enumerate(lines[0])}
        lines = [l[1:] for l in lines[1:]]
        D = {i: {} for i in range(len(lines))}
        for i, l in enumerate(lines):
            for j, sval in enumerate(l):
                D[i][j] = float(sval)
    return D, mapping


''' Performs the neighbor joining algorithm on a given set of sequences.

Arguments:
    D: map of maps, defining distances between the sequences
       (initially n x n, symmetric, 0 on the diagonal)
       (index -> index -> distance)
Returns:
    edges: edges formed using the neighbor-join algorithm
    D: updated matrix of distances (original n x n are the same, adds new
       distances)
'''
def neighbor_join(D):
    edgs = []
    n = len(D)

    while n > 2:
        d_new = []
        for row in range(len(D)): d_new += [[0]*len(D)]
        col_dists = [0.0]*len(D)

        for a in range (0, len(D)):
            col_dist = 0.0
            for b in range (0, len(D)):
                col_dist = col_dist + D[a][b]
                d_new[a][b] = (n-2) * D[a][b]
            col_dists[a] = col_dist
        #subtract total distances
        for a1 in range (0, len(D)):
            for b1 in range (0, len(D)):
                if d_new[a1][b1] != 0:
                    d_new[a1][b1] = d_new[a1][b1] - col_dists[a1] - col_dists[b1]

        #find smallest distances
        smal = 1000
        nods = (0,0)
        for i in range (1, len(D)):
            for j in range (0, i):
                if d_new[i][j] < smal:
                    smal = d_new[i][j]
                    nods = (i,j)
        #add edges
        Dij = D[nods[0]][nods[1]]
        edgs.append((nods[0],len(D), (Dij + (col_dists[nods[0]] - col_dists[nods[1]])/(n-2))/2))
        edgs.append((nods[1],len(D), (Dij + (col_dists[nods[1]] - col_dists[nods[0]])/(n-2))/2))
        #add new node to D and remove old ones
        for x in range (0, len(D)):
            if D[nods[0]][x] != 0 and D[nods[1]][x] != 0:
                D[x][len(D)] = (D[nods[0]][x] + D[nods[1]][x] - Dij) / 2
            else:
                D[x][len(D)] = 0
        D[len(D)] = {}
        for y in range (0, len(D)-1):
            D[len(D)-1][y] = D[y][len(D)-1]
        #print D[len(D)-1]
        D[len(D)-1][len(D)-1] = 0
        for z in range (0, len(D)):
            D[nods[0]][z] = 0.0
            D[nods[1]][z] = 0.0
            D[z][nods[0]] = 0.0
            D[z][nods[1]] = 0.0

        n = n-1

    #find final edge between 2 nodes
    for w in range (0, len(D)):
        if D[w][len(D)-1] != 0:
            edgs.append((w, len(D)-1, D[w][len(D)-1]))

    return edgs, D


''' Helper function for defining a tree data structure.
    First finds the edge to add a root node to and then generates binary tree.
    Root node should be at the midpoint of the longest branch.

Arguments:
    edges: edges in the tree
    num: number of nodes in tree
Returns:
    root: root node
    tree_map: map from nodes in the tree -> list of children (leaves have
              empty lists)
'''
def assemble_tree(edges, num):
    root = edges[len(edges)-1]
    tmap = []
    for row in range(num): tmap.append([])
    node_in_tree = ["no"] * num
    node_in_tree[root[0]] = "yes"
    node_in_tree[root[1]] = "yes"

    for i in range (len(edges)-2, -1, -1):
        if node_in_tree[edges[i][1]] == "yes":
            tmap[edges[i][1]].append((edges[i][0], edges[i][2]))
            node_in_tree[edges[i][0]] = "yes"
        elif node_in_tree[edges[i][0]] == "yes":
            tmap[edges[i][0]].append((edges[i][1], edges[i][2]))
            node_in_tree[edges[i][1]] = "yes"

    return root, tmap


def recurs(node_lst, mping, tree_mp):
    retur = ""
    for n in node_lst:
        if tree_mp[n[0]] == []:
            retur += str(mping[n[0]]) + ":" + str(n[1]) + ", "
        else:
            retur += "(" + recurs(tree_mp[n[0]], mping, tree_mp) + "):" + str(n[1]) + ", "
    retur = retur[:-2]
    return retur

''' Returns a string of the Newick tree format for the tree rooted at `root`.

Arguments:
    root: root of the tree (int)
    tree_map: map from node to list, describing each node's immediate children
    D: distance matrix of all nodes
    mapping: index to name mapping (dictionary)
Returns:
    output: rooted tree in Newick tree format (string)
'''
def generate_newick(root, tree_map, D, mapping):
    newick_tree = "("
    if tree_map[root[0]] == []:
        newick_tree += str(mapping[root[0]]) + ":" + str(root[2]/2) + ", "
    else:
        node_list = tree_map[root[0]]
        newick_tree += "(" + recurs(node_list, mapping, tree_map) + "):" + str(root[2]/2) + ", "

    if tree_map[root[1]] == []:
        newick_tree += str(mapping[root[1]]) + ":" + str(root[2]/2) + ", "
    else:
        node_list = tree_map[root[1]]
        newick_tree += "(" + recurs(node_list, mapping, tree_map) + "):" + str(root[2]/2)

    return newick_tree + ");"


def main():
    parser = argparse.ArgumentParser(
        description='Neighbor-joining algorithm on a set of n sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='dist10.txt')
    args = parser.parse_args()
    distances_file = args.f

    D, mapping = read_data(distances_file)
    edges, D = neighbor_join(D)
    #print D
    #print edges
    root, tree_map = assemble_tree(edges, len(D))
    #print tree_map
    print generate_newick(root, tree_map, D, mapping)


if __name__ == "__main__":
    main()
