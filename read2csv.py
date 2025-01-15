import networkx as nx
import matplotlib.pyplot as plt
import random as rand
import csv
import sys
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
        for _ in range(20): mm_file.readline()
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
                x,y,w = coord.split()
                #x,y = coord.split(' ')
                # 1 based not 0 based index:
                if x != y: # no self loops for now
                    true_vert+=1
                x = int(x) - 1
                y = int(y) - 1
                w = float(w) # FLOAT
                #w=1
                writer.writerow([x,y,w])
    
    return row

def matmart_to_ncol(filename, ncol_filename):
        # read suitespare [case: MM]
    with open(filename, 'r') as mm_file:
        for _ in range(1): mm_file.readline()
        row,col,verts = mm_file.readline().split()
        coords = mm_file.read()

    print(row,col,verts)

    true_vert = 0

    coords = coords.split('\n')

    # write to csv
    with open(ncol_filename, 'w') as csvfile:
        # creating a csv dict writer object
        writer = csv.writer(csvfile, delimiter=" ")

        # writing headers (field names)
        # writer.writerow(['node1','node2','weight'])

        # write rows - node1, node2, weight
        for coord in coords:
            if coord:
                x,y,w = coord.split()
                #x,y = coord.split(' ')
                # 1 based not 0 based index:
                if x != y: # no self loops for now
                    true_vert+=1
                x = int(x) - 1
                y = int(y) - 1
                w = float(w) # FLOAT
                #w=1
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

def txt_2_data(filename, N):
    
    itercount = 0
    edgenumcount = 0
    missizecount = 0
    cvertcount = 0
    timecount=0
    spncount=0

    iter=0
    numedge=0
    vertnum=0
    spn=0
    time=0
    M1_len = 0

    prev_edge = -1
    prev_iter = 0

    edge_ratio = 0
    edge_ratiocount = 0

    maxcoarsesz = 0
    maxedge=0
    
    with open(filename, 'r') as data_file:
        for data in data_file:
            
            if data[0] == '!':
                vertnum += int(data[3:].split(',')[0])

                if maxedge < int(data[3:].split(',')[1]) and int(data[3:].split(',')[2])>0:
                    maxedge = int(data[3:].split(',')[1])

                if int(data[3:].split(',')[0]) == N:
                    prev_iter = 0
                    prev_edge = -1

                if prev_edge == -1:
                    prev_edge = int(data[3:].split(',')[1])
                else:
                    edge_ratio += float(int(data[3:].split(',')[1])/prev_edge)
                    prev_edge = int(data[3:].split(',')[1]) 
                    edge_ratiocount+=1

                if int(data[3:].split(',')[2]) > prev_iter:
                    prev_iter = int(data[3:].split(',')[2])
                
                cvertcount+=1
                
            elif data[0] == 'M':
                M1_len += int(data[:-1].split(" ")[2])
                missizecount+=1
            elif data[0] == 'S':
                if maxcoarsesz < int(data.split(" ")[1]):
                    maxcoarsesz = int(data.split(" ")[1])
                spn += int(data.split(" ")[1])
                spncount+=1

            elif data[0] == 'T':
                iter+=prev_iter
                itercount+=1
                prev_iter = 0

                time += float(data.split(" ")[1])
                timecount+=1
        
    print(f"Avg MIS Sz: {M1_len/missizecount}")
    print(f"Avg Coarsen Sz: {spn/spncount}")
    print(f"Max. Coarsen Sz: {maxcoarsesz}")
    #print(f"Avg Edge Ratio: {edge_ratio/edge_ratiocount}")
    print(f"Max Edge #: {maxedge}")
    print(f"Avg Iter #: {iter/itercount}")
    print(f"Avg CPU Time: {time/timecount}")
    

def main():

    filename = 'suitesparse/bcsstk17.mtx'

    csv_file = 'csv/bcsstk17.csv'

    ncol_file = 'csv/bcsstk17.ncol'

    #N = read_suitesparse_rb_to_csv(filename, csv_file)
    #N = matmart_to_ncol(filename, ncol_file)
    #print(N)

    N=sys.argv[1] # 85
    txt_2_data(sys.argv[2],N) # "output/cilk1_a85.txt",N)
    print("\n")

    return

if __name__ == '__main__':
    main()