import networkx as nx
import matplotlib.pyplot as plt
import random as rand
import csv
from collections import defaultdict
import numpy as np
from scipy import sparse
from scipy.spatial.distance import euclidean

## FILE THAT CONVERTS GRAPHS AND MATRICES TO CSV OR DIFF. DATA STRUCTS ##

def read_kfile(filename):
    # read in file of:

    # N
    # num_edges
    # adj_matrix

    with open(filename) as f:
         ls = f.read()

    n, numEdges, adjMatrix = ls.splitlines()

    n = int(n)
    numEdges = int(numEdges)

    adjMatrix = map(int, adjMatrix.split(" "))

    return n, numEdges, adjMatrix

def read_file(filename):
    # read in file of g500 edgelist format

    with open(filename) as f:
         ls = f.read()

    startVerts, endVerts, weights = ls.splitlines()

    startVerts = map(float, startVerts.split(" ")[:-1])
    endVerts = map(float, endVerts.split(" ")[:-1])
    weights = map(float, weights.split(" ")[:-1])

    edgelist = [startVerts,endVerts,weights]

    return edgelist

def matToCSR(matrix, n):
    # convert a matrix to CSR

    matrix = list(matrix)
    matrix = np.reshape(matrix, (n,n))
    csr_mat = sparse.csr_matrix(matrix)

    return csr_mat

def CSRtoMat(csr, n):
    # Convert CSR to 1D matrix

    matrix = np.reshape(csr.toarray(), (1, n*n))
    return matrix[0] # reshape = [*[0 0 ... 0 0]*]

def matTolist(matrix, n):
    # convert a matrix to adjacency list
    # matrix = generator object of some type
    matrix = list(matrix)

    adjList = defaultdict(list)

    for index, k in enumerate(matrix):
        
          i = index//n
          j = index%n

          if k != 0:
            adjList[i].append(j)

    return adjList

def CSRtoDict(csr):
    # Convert CSR graph to Dictionary

    csr = csr.toarray()
    D = defaultdict(list)

    for i, row in enumerate(csr):
        for j, num in enumerate(row):
            if num:
                D[i].append((j, float(csr[i][j])))

    return D

def CSRtoDict2(csr):
    # Convert CSR graph to Dictionary

    csr = csr.toarray()
    D = defaultdict(list)

    for i, row in enumerate(csr):
        for j, num in enumerate(row):
            if num:
                D[(i,j)] = float(csr[i][j])

    return D

def DicttoCSR(dictG):

    N = len(dictG.keys())
    E = len(dictG.values())

    G = sparse.dok_matrix((N,E))

    for k, st in dictG.items():
        for v in st:
            G[k, v[0]] = v[1]

    G = G.transpose().tocsr()

    return G

def dictToCSV(D):
    fields = ['node1', 'node2', 'weight']
    filename = 'src/csv/graph_generation_000.csv'

    with open(filename, 'w') as csvfile:
        # creating a csv dict writer object
        writer = csv.writer(csvfile)

        # writing headers (field names)
        writer.writerow(fields)

        # write rows - node1, node2, weight
        for k, v in D.items():
            for x in v:
                writer.writerow([k, x[0], x[1]])

def edgelistToCSV(edgelist):
    fields = ['node1', 'node2', 'weight'] # write to csv file
    filename = 'csv/edgelist_000.csv'

    with open(filename, 'w') as csvfile:
        # creating a csv dict writer object
        writer = csv.writer(csvfile)

        # writing headers (field names)
        writer.writerow(fields)

        # write rows - node1, node2, weight
        for x, n in enumerate(edgelist[0]):
            writer.writerow([int(edgelist[0][x]), int(edgelist[1][x]), edgelist[2][x]])

def read_g500_file(filename):
    G = sparse.load_npz(filename)
    return G

def read_csv_file(N, filename, sym=False):
    # read in csv file - take in number of vertices
    # sparse matrix = [i][j] = weight
    sep=','
    i = []
    j = []
    w = []
    with open(filename, mode ='r')as file:
        csvFile = csv.reader(file)
        headers = next(csvFile, None)
        for data in csvFile:
            i.append(int(data[0]))
            j.append(int(data[1]))
            w.append(float(data[2]))
            if not sym:
                j.append(int(data[0]))
                i.append(int(data[1]))
                w.append(float(data[2]))


    G = sparse.csr_matrix((w, (i, j)), shape=(N, N), dtype=float) 

    return G

def nxtoG(NX):

    G = sparse.csr_matrix((NX.number_of_nodes(), NX.number_of_nodes()), dtype=float).toarray()

    pos = nx.spring_layout(NX)

    for (a,b,_) in NX.edges(data=True):
        G[a][b] = float(1./euclidean(pos[a], pos[b]))
        G[b][a] = float(1./euclidean(pos[b], pos[a]))
 
    return sparse.csr_matrix(G)