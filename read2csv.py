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

def read_suitesparse_rb_to_csv(filename, csv_file):

    # read suitesparse [case: RB]
    '''with open(filename, 'r') as rb_file:
        rb_file.readline()
        rb_file.readline()
        info = rb_file.readline()
        rb_file.readline()
        ls = rb_file.read()

    N = int(info.split()[1])
    rows = int(info.split()[1])
    cols = int(info.split()[2])
    print(N)'''

    # read suitespare [case: MM]
    with open(filename, 'r') as mm_file:
        for _ in range(44): mm_file.readline()
        row,col,verts = mm_file.readline().split()
        coords = mm_file.read()

    print(row,col,verts)

    true_vert = 0

    coords = coords.split('\n')

    # write to csv
    with open(csv_file, 'w') as csvfile:
        # creating a csv dict writer object
        writer = csv.writer(csvfile)

        # writing headers (field names)
        writer.writerow(['node1','node2','weight'])

        # write rows - node1, node2, weight
        for coord in coords:
            if coord:
                #x,y,w = coord.split()
                x,y = coord.split(' ')
                # 1 based not 0 based index:
                if x != y: # no self loops for now
                    true_vert+=1
                x = int(x) - 1
                y = int(y) - 1
                #w = float(w) # FLOAT
                w=1
                writer.writerow([x,y,w])
    
    return row

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

def read_csv_file(N, filename):
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

    G = sparse.csr_matrix((w, (i, j)), shape=(N, N), dtype=float) 

    return G

def main():

    filename = 'suitesparse/shock-9.mtx'

    csv_file = 'csv/shock-9.csv'

    N = read_suitesparse_rb_to_csv(filename, csv_file)
    print(N)

    return

if __name__ == '__main__':
    main()