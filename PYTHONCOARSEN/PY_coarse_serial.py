import networkx as nx
import numpy as np
import random as rand
import copy
import time

import PY_convert_graph_formats as convert

# Measuring Performance #
# compare to a baseline
# how it scales on certain problems
# speedups
# mean times to run algs
# size&iteration count for varying prob sizes using structured probs
# grids used to vary grid dimension and see quality of result
# algorithmic scaling
# portability/scalability - running on diff. archs.
# run time = time spent setup+solve+aggregate times
# 

IN = 0
OUT = np.iinfo(np.int32).max
UNDECIDED = 2

NODES = 16

const = int('0x2545F4914F6CDD1D', 16) # ignore #

######### FIND MIS ###########

def kokkos_mis(G): # 3. Build MIS of supernodes based off Kokkos paper #
    print('3')
    numVerts = len(G.toarray()[0]) # assuming nv = numVerts

    rowStatus = [i for i in range(numVerts)]
    colStatus = [OUT]*numVerts

    worklist1 = set(range(numVerts)) # set of UNDECIDED vertices

    worklist2 = set(range(numVerts)) # set of vertices that are not OUT

    iter = 0

    nvbits = randPriori(numVerts)

    while worklist1:
    
        rowStatus = refreshRowStatus(rowStatus, worklist1, iter, nvbits, numVerts) # 3.a. # parallel

        # sync

        colStatus = refreshColStatus(colStatus, rowStatus, worklist2, numVerts, G)  # 3.b.

        # sync
        rowStatus, colStatus = decideInOut(worklist1, rowStatus, colStatus, G) # 3.c.

        # sync

        worklist1 = set(v for v in worklist1 if rowStatus[v] != IN and rowStatus[v] != OUT) 

        worklist2 = set(v for v in worklist2 if colStatus[v] != OUT)

        iter+=1

    return  set(v for v,k in enumerate(rowStatus) if k == IN)

def refreshRowStatus(rowStatus, worklist1, iter, nvBits, numVerts): # 3.a. assign random 'statuses' to vertices
    
    # paper uses an xorshift* but I could not quite figure it out so I xor random values b/w 1 and number of vertices-1

    for i in worklist1: # parallel for
        rand.seed()

        p1  = rand.randint(1, numVerts-1) # randval1
        p2  = rand.randint(1, numVerts-1) # randval2

        priority = p1^p2 # xor to get priority
        
        new_status = (i + 1) | (priority << nvBits)

        if (new_status == OUT): new_status-=1 # 1 or max

        rowStatus[i] = new_status

    # TODO: correctly implement, 64-bit xorshift* #
    return rowStatus

def refreshColStatus(colStatus, rowStatus, worklist2, nv, G): # 3.b. figure out if neighbors of a given vertex are in or out 

    for i in worklist2: # parallel for
        neighbors = adj(G.toarray()[i], i)

        s = rowStatus[i]

        for neigh in neighbors:

            if neigh < nv and neigh != i:
                neigh_stat = rowStatus[neigh]

                if (neigh_stat < s): 
                    s = neigh_stat

        if (s == IN): 
            s = OUT

        colStatus[i] = s


    return colStatus

def decideInOut(worklist1, rowStatus, colStatus, G): # 3.c. compare neighbors and see which vertices can be IN or OUT
    #print(f'Before: row-- {rowStatus}\ncol--{colStatus}')
    nv = len(G.toarray()[0]) # nv = numverts

    for i in worklist1:
        s = rowStatus[i]

        if (s == IN or s == OUT): continue # return or continue?

        neighbors = adj(G.toarray()[i], i)

        neigh_out = False
        neigh_mismatchS = False

        for neigh in neighbors: # ** # # if a given vertex passed in is UNDECIDED, see if neighbors are IN or OUT
        
            if neigh >= nv: continue # if it's not a valid index, skip
            neigh_stat = colStatus[neigh]

            if neigh_stat == OUT:
                neigh_out = True
                break
            elif neigh_stat != s:
                neigh_mismatchS = True

        if neigh_out: # update col stats for all neighbors of i

            rowStatus[i] = OUT # if a neighbor is out, so is given vertex

        elif not neigh_mismatchS: # s is min among all neighbors
            rowStatus[i] = IN # if a given vertex has the minimum priority of its neighbors, it is IN

    return (rowStatus, colStatus)

###############################

def xorshift64star(a): # IGNORE #
    
    #/* Algorithm "xorshift*" from Marsaglia, https://en.wikipedia.org/wiki/Xorshift && https://rosettacode.org/wiki/Pseudo-random_numbers/Xorshift_star*/
    x = a

    x ^= x >> 12
    x ^= x << 25
    x ^= x >> 27

    seed = x

    # this ain't happening right now dog
    # https://github.com/kokkos/kokkos/blob/develop/algorithms/src/Kokkos_Random.hpp#L925

    return (x*const) >> 32

def adj(pass_G, v): # return the neighboring adjacent vertices 
    # if graph is undirected, only need to check row of CSR matrix for neighbors
    neighbors = [j for j, x in enumerate(pass_G) if x != 0 and j != v]

    # if directed, check row and column
     
    return neighbors # TODO: return NEIGHBORS of v

def randPriori(numVerts):
     
    i = numVerts + 1
    nvBits = 0
    while(i > 0):
        i = int(i>>1)
        nvBits+=1

    return nvBits

###### ASSIGN SUPERNODES #######

def firstPass(G, M1): # 4. Root + Neighbors = Supernode
    print('4')
    # First Pass - Root + Neighbors = Supernode
    for supernode, m in enumerate(M1): # parallel for

        # build aggregate from neighbors of M1
        labels[m] = supernode

        neighbors = adj(G.toarray()[m], m)
        
        for neigh in neighbors: 
            if labels[neigh] == -1: # potential race condition? does it matter?
                labels[neigh] = supernode

        # supernode+=1 # new supernode
    return

def secondPass(G, vertex_set): # 5. Unassigned Neighbors of RootNeighbors = Supernode
    print('5')
    # Second Pass - Unassigned Neighbors of RootNeighbors = Supernode
    for unasgn in vertex_set: # parallel?

        if labels[unasgn] != -1: continue

        else:
            unasgn_neighs = adj(G.toarray()[unasgn], unasgn)

            for nn in unasgn_neighs: #
                if labels[nn] != -1: 
                    # race condition? benign?
                    labels[unasgn] = labels[nn] # not sure if this would ruin parallelness
                    break

    #print(labels)
    return

def thirdPass(G, vertex_set, lm1): # 6. Still Unassigned? New Supernode & Reassign Neighbors
    print('6')
    # Third Pass - Still Unassigned? New Supernode & Reassignment of Neighbors

    supernode = lm1

    for spnde, unasgn in enumerate(vertex_set): 

        if labels[unasgn] != -1: continue
        else:
            labels[unasgn] = supernode

            unasgn_neighs = adj(G.toarray()[unasgn], unasgn)

            for u in unasgn_neighs: # might not want to completely reassign/overwrite other relations, don't know yet
                labels[u] = supernode

            supernode+=1 # shared var
            #print('3rd pass: ', supernode, labels)
    return

###### COARSEN + DEPENDENCIES ###

def build_dependencies(G, labels): # 7.  Construct the subgraph and supernode dependencies
    print('7')
    supernodes = dict() # subgraph
    dependency = dict() # submeta - which vertices from orig. Graph have edges across >1 supernode and the edge weights
    sv = dict() # which vertices fall under which supernodes ex, {0: {2,3,4}, ...}

    total_weight = 0

    G_dict = convert.CSRtoDict(G)

    # 0: [(1,8)], 1:[(0,8), (2,9)] supernode: [(neighbor, weight),...]
    # 0: {(0,1,8), (1,0,8)} supernode:{(g1, g2, weight), (g2, g3, weight)}

    for i, v in enumerate(labels): # go through the labels array and construct 
        # v = supernode index/vertex i is assigned to in new SubGraph

        try: sv[v].add(i)
        except KeyError: sv[v] = {i}

        for c in G_dict[i]: # NOT parallel
            t, w = c          # c = short coodrinate
            coord = (i, t, w) # index of labels = G vertex, target, weight

            try: 
                if (t, i, w) not in supernodes[v]: 
                    supernodes[v].add(coord) 
            except KeyError: supernodes[v] = {coord}

    if len(supernodes.keys()) <= 1: # if there's one supernode, condense the graph and return
        for v in supernodes.values():
            for sets in v:
                x,y,w = sets
                total_weight+=w
        tup = (-1, total_weight)             
        dependency[0] = {tup} # divide weight by 2 if graph coordinates are symmetric
        return supernodes, dependency
    
    else: # build dependencies and keep track of weights and relationships
        for spnodes in sv.keys():
            dependency[spnodes] = set()

        for spnode, set_coor in supernodes.items(): # for each spnode=supernode
            tw = 0
            for coor in set_coor: # for each edge per supernode
                s,t,w = coor
                                        
                if t not in sv[spnode]: # if s and t in different spnodes but the spnodes have an edge                
                    
                    if len(dependency[spnode]): 
                        temp_set = copy.deepcopy(dependency[spnode])
                        
                        for ct in temp_set:

                            target, targ_w = ct
                            if target != labels[t]: continue 
                            else:
                                dependency[spnode].discard((target, targ_w)) # ? #
                                dependency[target].discard((spnode, targ_w))

                                dependency[spnode].add((target, (w/2)+targ_w)) # divide by 2 to not double count an edge weight
                                dependency[target].add((spnode, (w/2)+targ_w))

                        
                    else:                            
                        dependency[spnode] = {(labels[t], w/2)}
                        dependency[labels[t]] = {(spnode, w/2)}


    #print('\nSupernodes: ', supernodes) # printed in different function
    #print('Dependency: ', dependency, '\n')

    return supernodes, dependency

def kokkos_coarsen(G): # 2. Kokkos Coarsen #
    print('2')
    N = len(G.toarray()[0])

    supernode = 0
    vertex_set = set(range(N)) # set of vertices to assign to supernodes

    M1 = kokkos_mis(G) # parallel

    print(f'MIS2 IN SET: {M1}')

    # sync

    global labels # ??
    labels = [-1]*N # array to keep track of the supernode parents to graph vertices

    # First Pass - Aggregate a root and its immediate neighbors into a supernode

    firstPass(G, M1) # 

    # sync

    # 4.5. Clean up set of vertices to assign to supernodes
    for m in M1: # not parallel
        vertex_set.discard(m)  # Remove vertices from the vertex set that are in M1
    for i, l in enumerate(labels): # not parallel
        if l != -1:
            vertex_set.discard(i) # Remove vertices from the vertex set that have an assigned supernode
        else: continue

    # Second Pass - # If a vertex is unassigned and within 2 nodes of a supernode, assign it to that supernode
                    # Neighbor to a Supernode Neighbor

    secondPass(G, vertex_set)

    # sync

    # 5.5. Clean up set of vertices to assign to supernodes
    for l in labels: # not parallel
        if l != -1:
            vertex_set.discard(l)
        else: continue

    lm1 = len(M1) # next supernode id to create

    # Third Pass - Vertex Still Unassigned? Create a new Supernode & Reassign its neighbors

    thirdPass(G, vertex_set, lm1)

    # sync

    sub_meta, sub_graph = build_dependencies(G, labels) # this might have to be serial

    return sub_graph, sub_meta

def recursion_fcn(G, gc_verts): # 1. Recursive function for coarsening #
    subgraph_size = len(G.toarray()[0]) if len(G.toarray()) else 0

    # Base Case
    if subgraph_size <= gc_verts:
        return 1

    print('1')

    # Recurse

    subgraph, submeta = kokkos_coarsen(G) # subgraph = new coarsened graph, submeta = meta data on the parent graph verticess that have edges across supernodes and the new edge weight

    subgraph_size = len(subgraph.keys())

    print(f'Graph Coarsened: {subgraph})')#\nMetaDependencies:')
    #for k in submeta.keys():
    #    print(f'{k}:{submeta[k]}\n')
    print(f'vertices coarsened in this round: {len(G.toarray()[0])} - {len(G.toarray()[0])-len(subgraph.keys())} = {len(subgraph.keys())}')

    subgraph_CSR = convert.DicttoCSR(subgraph) # issue w/ key to coordinate translation

    return recursion_fcn(subgraph_CSR, gc_verts)

def main():

    ## Small Test Graph Options ##

    N = 85

    #### Simple Graphs ####
    #G = convert.read_csv_file(N, 'csv/simple_graph_000.csv',True) # N=5
    #G = convert.read_csv_file(N, 'csv/kk_simpleEx.csv',True)      # N=6
    #G = convert.read_csv_file(N, 'csv/simple_graph_001.csv',True) # N=7

    #### SuiteSparse Graphs ####
    G = convert.read_csv_file(N, 'csv/ash85.csv', False) # N=85
    #G = convert.read_csv_file(N, 'csv/netz4504.csv', False)  # N=1961

    #### NetworkX Graphs ####
    #NX = nx.fast_gnp_random_graph(N, 0.1)
    #G = convert.nxtoG(NX)

    ## Recursive Coarsen ##

    #*** G is in CSR format -- scipy, sparse ****#

    gc_verts = 6 #int((N/2)-1) # coarsen graph to N/2-1 of original size

    start_time = time.time()
    exit = recursion_fcn(G, gc_verts) # Graph, empty dependency graph, goal # coarsened vertices
    end = time.time() - start_time
    print(exit, end)

    return

if __name__ == '__main__':
        main()