import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as rand
import copy

import PY_convert_graph_formats as convert
import mtsim

# Ran with a conda environment with  - numpy - matplotlib - jupyter

IN = 0
OUT = np.iinfo(np.int32).max
UNDECIDED = 2

NODES = 16

const = int('0x2545F4914F6CDD1D', 16)

# TODO: do all single vars get labelled as args[var_name]

def kokkos_mis(G): # G is a CSR?

    numVerts = len(G.toarray()[0]) # assuming nv = numVerts?

    rowStatus = [i for i in range(numVerts)]
    mtsim.mt_array_malloc(rowStatus, mtsim.mt_block_cyclic, [0, 2, NODES])
    #### mt_allocate ####

    #rowStatus = np.random.default_rng().choice(numVerts, size=numVerts, replace=False) # random priorities
    colStatus = [OUT]*numVerts
    mtsim.mt_array_malloc(colStatus, mtsim.mt_block_cyclic, [0, 2, NODES])
    #### mt_allocate ####

    worklist1 = set(range(numVerts))
    # mt_allocate #

    worklist2 = set(range(numVerts))
    # mt_allocate #

    iter = 0

    nvbits = randPriori(numVerts)

    while worklist1:

        rowStatus = refreshRowStatus(rowStatus, worklist1, iter, nvbits, numVerts) # parallel

        colStatus = refreshColStatus(colStatus, rowStatus, worklist2, numVerts, G) 

        rowStatus, colStatus = decideInOut(worklist1, rowStatus, colStatus, G)

        worklist1 = set(v for v in worklist1 if rowStatus[v] != IN and rowStatus[v] != OUT)
        # mt_write? reallocate? #
        worklist2 = set(v for v in worklist2 if colStatus[v] != OUT)
        # mt_write? reallocate? #

        iter+=1

    # mt_write? mt_alloc? #
    return  set(v for v,k in enumerate(rowStatus) if k == IN)

def refreshRowStatus(rowStatus, worklist1, iter, nvBits, numVerts):

    for i in worklist1: # parallel for
        rand.seed()
        #priority = (xorshift64star(iter)^xorshift64star(i))
        #rand.seed(iter)
        p1  = rand.randint(1, numVerts-1)
        #if i != iter : rand.seed(i)
        #else: rand.seed()
        p2  = rand.randint(1, numVerts-1)
        priority = p1^p2
        
        new_status = (i + 1) | (priority << nvBits)
        if (new_status == OUT): new_status-=1 # 1 or max
        mtsim.mt_array_write(rowStatus, i, new_status) # rowStatus[i] = new_status
        #### mt_write ####

    # TODO: correctly implement, 64-bit xorshift*

    return rowStatus

def refreshColStatus(colStatus, rowStatus, worklist2, nv, G): # there's a parallel version of this function but I don't quite get it yet

    for i in worklist2: #parallel for
        neighbors = adj(G.toarray()[i], i)
        #### mt_allocate ####

        s = mtsim.mt_array_read(rowStatus, i) # rowStatus[i]
        #### mt_read ####

        for neigh in neighbors:
            if neigh < nv and neigh != i:
                neigh_stat = mtsim.mt_array_read(rowStatus, neigh) # rowStatus[neigh]
                #### mt_read ####

                if (neigh_stat < s): s = neigh_stat

        if (s == IN): s = OUT

        mtsim.mt_array_write(colStatus, i, s) # colStatus[i] = s
        #### mt_write ####

    return colStatus

def decideInOut(worklist1, rowStatus, colStatus, G): # there's a parallel version of this function but I don't quite get it yet

    nv = len(G.toarray()[0])# if nv = numverts

    for i in worklist1:
        s = mtsim.mt_array_read(rowStatus, i) # rowStatus[i]
        #### mt_read ####

        if (s == IN or s == OUT): continue # return or continue?

        neighbors = adj(G.toarray()[i], i)
        # does allocating in adj work?
        #### mt_allocate ####

        neigh_out = False
        neigh_mismatchS = False

        for neigh in neighbors: #** # account for 
            if neigh >= nv: continue
            neigh_stat = mtsim.mt_array_read(colStatus, neigh) # colStatus[neigh]
            #### mt_read ####

            if neigh_stat == OUT:
                neigh_out = True
                break
            elif neigh_stat != s:
                neigh_mismatchS = True
            

        if neigh_out: # update col stats for all neighbors of i

            mtsim.mt_array_write(rowStatus, i, OUT) # rowStatus[i] = OUT
            #### mt_write ####

        elif not neigh_mismatchS: # s is min among all neighbors

            mtsim.mt_array_write(rowStatus, i, IN) # rowStatus[i] = IN
            #### mt_write ####

    return (rowStatus, colStatus)

####################################################################

def xorshift64star(a):
    
    #/* Algorithm "xorshift*" from Marsaglia, https://en.wikipedia.org/wiki/Xorshift && https://rosettacode.org/wiki/Pseudo-random_numbers/Xorshift_star*/
    x = a

    x ^= x >> 12
    x ^= x << 25
    x ^= x >> 27

    seed = x

    # this ain't happening right now dog
    # https://github.com/kokkos/kokkos/blob/develop/algorithms/src/Kokkos_Random.hpp#L925

    return (x*const) >> 32

def adj(pass_G, v):
    # if graph is undirected, only need to check row of CSR matrix for neighbors
    neighbors = [j for j, x in enumerate(pass_G) if x != 0 and j != v]
    mtsim.mt_array_malloc(neighbors, mtsim.mt_block_cyclic, [0, 2, NODES])
    #### mt_allocate ####

    # if directed, check row and column
     
    return neighbors # TODO: return NEIGHBORS of v

def randPriori(numVerts):
     
    i = numVerts + 1
    nvBits = 0
    while(i > 0):
        i = int(i>>1)
        nvBits+=1

    return nvBits

####################################################################

def build_dependencies(G, labels): # even if it's not parallel does it need mtsim stuff??

    supernodes = dict()
    dependency = dict() #?
    sv = dict()
    # mtsim.mt_multimapn_malloc #

    total_weight = 0

    G_dict = convert.CSRtoDict(G)
    # mtsim.mt_multimaps_malloc #
    # 0: [(1,8)], 1:[(0,8), (2,9)]

    # 0: {(0,1,8), (1,0,8)}

    for i, v in enumerate(labels): # parallel???? i don't think so

        try: sv[v].add(i)  # mt_mmap_write #
        except KeyError: sv[v] = {i}

        for c in G_dict[i]: # NOT parallel
            t, w = c
            coord = (i, t, w)

            try: 
                if (t, i, w) not in supernodes[v]:  # mt_mmap_read #
                    supernodes[v].add(coord)  # mt_mmap_write #
            except KeyError: supernodes[v] = {coord}

        # try: superneigh_dict[v] = set(adj(G, i)).union(superneigh_dict[v]) 
        # except KeyError: superneigh_dict[v] = set(adj(G, i))
 
    if len(supernodes.keys()) <= 1:
        for v in supernodes.values():
            for sets in v:
                x,y,w = sets
                total_weight+=w             # mt_mmap_write #
        dependency[0] = (-1, total_weight) # divide by 2 if graph coordinates are symmetric
        return supernodes, dependency
    
    else:
        for spnode, set_coor in supernodes.items(): # for each supernode
            tw = 0
            for coor in set_coor: # for each edge per supernode
                s,t,w = coor
                                        # mt_mmap_read #
                if t not in sv[spnode]: # if s and t in different spnodes but the spnodes have an edge                
                    
                    try: 
                        temp_set = copy.deepcopy(dependency[spnode])
                        # mt_write?
                        
                        for ct in temp_set:
                            #print('here2 ', ct, coor, spnode)
                            target, targ_w = ct
                            if target != labels[t]: continue # mt_mmap_read #
                            else:
                                dependency[spnode].discard((target, targ_w)) # ? #
                                dependency[target].discard((spnode, targ_w))

                                dependency[spnode].add((target, (w/2)+targ_w)) # mt_mmap_write #
                                dependency[target].add((spnode, (w/2)+targ_w))

                        
                    except KeyError:                            # mt_mmap_read #
                        dependency[spnode] = {(labels[t], w/2)} # mt_mmap_write #
                        dependency[labels[t]] = {(spnode, w/2)}

    print(supernodes)
    print(dependency)

    return supernodes, dependency

def kokkos_coarsen(G): # TODO: figure out how to keep metadata with assigning supernodes

    N = len(G.toarray()[0])

    labels = [-1]*N
    mtsim.mt_array_malloc(labels, mtsim.mt_block_cyclic, [0, 2, NODES])
    #### mt_allocate ####

    supernode = 0
    vertex_set = set(range(N))
    # mt_allocate #

    M1 = kokkos_mis(G) # parallel

    print(f'IN: {M1}')

    # TODO: Passes -> Fcns?
    # First Pass - Root + Neighbors = Supernode
    for supernode, m in enumerate(M1): # parallel for

        # build aggregate from neighbors of M1
        mtsim.mt_array_write(labels, m, supernode) #labels[m] = supernode
        #### mt_write ####

        neighbors = adj(G.toarray()[m], m)
        #### mt_alloc ####
        
        for neigh in neighbors: 
            if mtsim.mt_array_read(labels, neigh) == -1: #labels[neigh] == -1: # potential race condition? does it matter?
                #### mt_read ####

                mtsim.mt_array_write(labels, neigh, supernode) # labels[neigh] = supernode
                #### mt_write ####

        # supernode+=1 # new supernode

     # sync

    for m in M1: # not parallel
        vertex_set.discard(m)
        # realloc? #
    for l in labels: # not parallel
        if l != -1:
            vertex_set.discard(l)
            # realloc? #
        else: continue

    # Second Pass - Unassigned Neighbors of RootNeighbors = Supernode
    for unasgn in vertex_set: # parallel?

        if mtsim.mt_array_read(labels, unasgn) != -1: #labels[unasgn] != -1:
            #### mt_read ####
            continue

        else:
            unasgn_neighs = adj(G.toarray()[unasgn], unasgn)

            for nn in unasgn_neighs: #
                if mtsim.mt_array_read(labels, nn) != -1: # labels[nn] != -1: #
                    #### mt_read ####

                    # race condition? benign?
                    mtsim.mt_array_write(labels,unasgn, mtsim.mt_array_read(labels, nn)) #labels[unasgn] = labels[nn] # not sure if this would ruin parallelness
                    #### mt_write ####
                    break

    # sync

    for l in labels: # not parallel
        if l != -1:
            vertex_set.discard(l)
            # realloc? #
        else: continue

    lm1 = len(M1)

    # Third Pass - Still Unassigned? New Supernode & Reassignment of Neighbors
    for spnde, unasgn in enumerate(vertex_set): # I think this can be parallel....?

        supernode = spnde + lm1

        if mtsim.mt_array_read(labels, unasgn) != -1: # labels[unasgn] != -1: 
            #### mt_read ####
            continue
        else:
            mtsim.mt_array_write(labels, unasgn, supernode) # labels[unasgn] = supernode
            #### mt_write ####

            unasgn_neighs = adj(G.toarray()[unasgn], unasgn)

            for u in unasgn_neighs: # might not want to completely reassign/overwrite other relations, don't know yet
                mtsim.mt_array_write(labels, u, supernode) # labels[u] = supernode
                #### mt_write ####

            # supernode+=1

    sub_meta, sub_graph = build_dependencies(G, labels) # this might have to be serial

    print('SG: ', sub_graph)

    mtsim.mt_die()

    return sub_graph, sub_meta

'''
def main():

    #N = 5
    #G = convert.read_csv_file(N, 'csv/simple_graph_000.csv')

    N = 6
    G = convert.read_csv_file(N, 'csv/kk_simpleEx.csv')

    #N = 7
    #G = convert.read_csv_file(N, 'csv/simple_graph_001.csv')

    subgraph, submeta = kokkos_coarsen(G)
    print(f'GC: {subgraph}')



    print(f'verts coarsened in this round: {len(G.toarray()[0])} - {len(G.toarray()[0])-len(subgraph.keys())} = {len(subgraph.keys())}')

    return

if __name__ == '__main__':
        main()
#'''