import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as rand
import time

import PY_convert_graph_formats as convert

# Measuring Performance #
# costs of constructing coarse graphs
# local improvement alg cost
# speed of coarsening
# compare to other algoithms
# overall run time
# how much memory it takes to run
# partitioning metric: "gain" -- net redutcion in w of cut edges that would result from switching a vertex to a diff part. set

#### Works for NetworkX and custom input ####

# G - (1,2): (3, 4, 5) - source, target: edgew, sourcew, targetw
# weight of vertex?
def random_graph_gen(n, G): # simple graph generator - random number of nodes and edges

    # n == rand number of nodes (by hundreds for now)

    # rand number of edges
    randEdges = rand.randint(0, (n*(n-1))) # max edges = n*n-1

    print(f"Vertices: {n}")
    print(f"Edges: {randEdges}")

    # rand verts
    for edge in range(randEdges):
        addEdge = False

        randx = rand.randint(0, n-1)
        randy = rand.randint(0, n-1)
        vert = (randx, randy)

        randweight= rand.randint(1, int(n/4))

        while not addEdge: # while so true amount of edges is added
            if vert not in G and randx != randy: # vertice not added yet and not self vert
                G.add_edge(randx, randy, weight=randweight)
                addEdge = True
            else:         # if edge not added, try again
                randx = rand.randint(0, n-1)
                randy = rand.randint(0, n-1)
                vert = (randx, randy)
                addEdge = False
  
    graph = {}

    for x, nbrs in G.adj.items():
        for nbr in nbrs.keys():
            graph[(x, nbr)] = (G[x][nbr]['weight'],1,1)


    #l=[(n, nbrdict) for n, nbrdict in G.adjacency()]
    #print(l)     

    # Matplot draw graph            
    #fig, ax = plt.subplots()
    #pos = nx.random_layout(G)
    #nx.draw_networkx(G, pos=pos, ax=ax)
    #ax.set_title(f"DiGraph n={n}, num_edges={randEdges}")
    #plt.show()

    return graph

def maximal_matching(G): # NetworkX Greedy MM Function
    r"""Find a maximal matching in the graph.

    A matching is a subset of edges in which no node occurs more than once.
    A maximal matching cannot add more edges and still be a matching.

    Parameters
    ----------
    G : NetworkX graph
        Undirected graph

    Returns
    -------
    matching : set
        A maximal matching of the graph.

    Examples
    --------
    >>> G = nx.Graph([(1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5)])
    >>> sorted(nx.maximal_matching(G))
    [(1, 2), (3, 5)]

    Notes
    -----
    The algorithm greedily selects a maximal matching M of the graph G
    (i.e. no superset of M exists). It runs in $O(|E|)$ time.
    """
    matching = set()
    nodes = set()

    for edge in G.keys():
        # If the edge isn't covered, add it to the matching
        # then remove neighborhood of u and v from consideration.
        u, v = edge
        if u not in nodes and v not in nodes and u != v:
            matching.add(edge)
            nodes.update(edge)

    return matching

def draw_subgraph(SubGraph):
    # Matplot draw graph 
    Graph = nx.Graph(SubGraph.keys())    
        
    options = {
        "font_size": 8,
        "node_size": 300,
        "node_color": "white",
        "edgecolors": "black",
        "linewidths": 0.5,
        "width": 1,
    } 
    nx.draw_networkx(Graph, **options)
    ax = plt.gca()
    ax.margins(0.20)
    plt.axis("off")
    plt.show()

def construct_adjl(SubGraph):

    SadjL = {}

    for coord in SubGraph.keys():
        if coord[0]!= coord[1]:
            try: SadjL[coord[0]].add(coord[1])
            except KeyError: SadjL[coord[0]] = {coord[1]}

            try: SadjL[coord[1]].add(coord[0])
            except KeyError: SadjL[coord[1]] = {coord[0]}  

    return SadjL

def reconstruct(SubGraph, x, y, v, previous_child):#, vertex_weight): #### Re-Assgin Parent/Child Nodes ####
    num=0
    if x != -1 and y != -1:
        num=1
        #edge_weight = SubGraph[(previous_child, x)][0] + SubGraph[(previous_child, y)][0] # TODO: Double Check
        edge_weight = SubGraph[(previous_child, x)] + SubGraph[(previous_child, y)] # TODO: Double Check
        #child_weight = SubGraph[(previous_child, x)][1]
        del SubGraph[(previous_child, x)]
        del SubGraph[(previous_child, y)]

    elif x != -1 and y == -1:
        num=2
        edge_weight = SubGraph[(previous_child, x)]#[0]
        #child_weight = SubGraph[(previous_child, x)][1]
        del SubGraph[(previous_child, x)]

    elif x == -1 and y != -1:
        num=3
        edge_weight = SubGraph[(previous_child, y)]#[0]
        #child_weight = SubGraph[(previous_child, y)][1]
        del SubGraph[(previous_child, y)]

    SubGraph[(previous_child, v)] = edge_weight#(edge_weight,child_weight,vertex_weight)
    #print(f'{num}: {previous_child}, {v}')

    return SubGraph
def hendrickson_leland(G, adjL, n, goal_n, rounds):

    print(f'n: {n}, goal-n: {goal_n}, round: {rounds}')
    start = time.time()
    # Base Case
    if n <= goal_n:
        #print(f'Graph:\n{G}')
        #print(f'AdjList:\n{adjL}\n')
        print(f'n: {n}, goal-n: {goal_n}, round: {rounds}')
        #draw_subgraph(G)
        return 1
    
    # Matplot draw graph 
    #draw_subgraph(G)

    SubGraph = {}

    previous = {}

    # Find Max Edges
    max_edges = maximal_matching(G) # very slightly modified NetworkX algorithm

    for pair in max_edges:
        x = pair[0]
        y = pair[1]
        v = (x,y)

        # calc vertex weight
        # vertex_weight = G[(x,y)][1] + G[(x,y)][2]

        # neighbors of both vertices in pair
        neighbors = adjL[x].union(adjL[y]) 
        neighbors.discard(x)
        neighbors.discard(y)

        # associate x,y w/ supernode v=(x,y)
        previous[x] = v
        previous[y] = v 

        # for all the children
        for child in neighbors:

            m = 1 if child in adjL[x] else 0

            # if x and y are bot adjacent to child
            if m and child in adjL[y]:             

                try: G_xtuple = G[(x, child)]
                except KeyError: G_xtuple = G[(child, x)]
                try: G_ytuple = G[(y, child)]
                except KeyError: G_ytuple = G[(child, y)]

                edge_weight = G_xtuple + G_ytuple #G_xtuple[0] + G_ytuple[0]
                #child_weight = G_xtuple[2]

                # if child is from a previous supernode - delete any trace of it and put the supernode relation into subgraph
                if child in previous.keys() and (previous[child], x) in SubGraph and (previous[child], y) in SubGraph:
                        SubGraph = reconstruct(SubGraph, x, y, v, previous[child])#, vertex_weight)
                        continue
                elif child in previous.keys() and (previous[child], x) not in SubGraph and (previous[child], y) not in SubGraph:
                        continue # if it has the supernode relation has already been removed, skip
                
                elif child in previous.keys() and (previous[child], v) in SubGraph:  # Option for if one prevchild gets deleted and not the other -- edge weight needs to accum.

                    if (previous[child], y) in SubGraph: del SubGraph[(previous[child], y)]
                    if (previous[child], x) in SubGraph: del SubGraph[(previous[child], x)]

                    edge_weight = edge_weight + SubGraph[(previous[child], v)]#[0]
                    #child_weight = child_weight + SubGraph[(previous[child], v)][1]
                    SubGraph[(previous[child], v)] = edge_weight #(edge_weight,child_weight,vertex_weight)
                    #print(f'5: {previous_child}, {v}')

                    continue

            else: # if only child of 1 # in pair
                if m: 
                    try: G_xtuple = G[(x, child)]
                    except KeyError: G_xtuple = G[(child, x)]

                    edge_weight = G_xtuple#[0]
                    #child_weight = G_xtuple[2]

                     # if child is from a previous supernode - delete any trace of it and put the supernode relation into subgraph
                    if child in previous.keys() and (previous[child], x) in SubGraph:
                        SubGraph = reconstruct(SubGraph, x, -1, v, previous[child])#, vertex_weight)
                        continue
                    elif child in previous.keys() and (previous[child], x) not in SubGraph:
                        continue

                else: 
                    try: G_ytuple = G[(y, child)]
                    except KeyError: G_ytuple = G[(child, y)]

                    edge_weight = G_ytuple#[0]
                    #child_weight = G_ytuple[2]

                     # if child is from a previous supernode - delete any trace of it and put the supernode relation into subgraph
                    if child in previous.keys() and (previous[child], y) in SubGraph:
                        SubGraph = reconstruct(SubGraph, -1, y, v, previous[child])#, vertex_weight)
                        continue
                    elif child in previous.keys() and (previous[child], y) not in SubGraph:
                        continue
            
            SubGraph[(v, child)] = edge_weight#(edge_weight, vertex_weight, child_weight)
            #print(f'4: {v}, {child}')
        
 
    # Construct Subgraph Adjancency List
    SadjL = construct_adjl(SubGraph)

    print("time: ", time.time() - start)
    # Recursively call until graph is coarsened enough
    return hendrickson_leland(SubGraph, SadjL, len(SadjL.keys()), goal_n, rounds+1)


def main():
    '''
    #### Simple Graph ####
    #G = nx.Graph()
    #n = 100
    #G = random_graph_gen(n, G)
    #al = construct_adjl(G)

    #### Very Simple Graph Test #### OLD!!!!
    #G = {(0,1): (3,6,3),(0,3): (16,6,12),(0,4): (19,6,9),(1,0): (3,3,6),(1,2): (12,3,10),(2,1): (12,10,3),(2,4): (7,10,9),(2,5): (3,10,4),(3,0): (16,12,6),(3,4): (11,12,9),(4,0): (19,9,6),(4,2): (7,9,10),(4,3): (11,9,12),(4,5): (14,9,4),(5,2): (3,4,10),(5,4): (14,4,9)}
    #al = {0:{1,3,4}, 1: {0,2}, 2: {1,4,5}, 3: {0,4}, 4:{0,2,3,5}, 5: {2,4}}

    #n = len(al.keys())

    #Gr = nx.Graph(G.keys())
    '''

    N = 1961
    goal_n = int(N/2) - 1
    rounds = 1
    Gr = convert.read_csv_file(N, 'csv/netz4504.csv', False)  # N=1961
    #Gr = convert.read_csv_file(N, 'csv/shock-9.csv', False)  # N=36476
    # NEW: {(x,y): edge_weight}

    start = time.time()

    G = convert.CSRtoDict2(Gr)
    al = construct_adjl(G)

    print('Graph & adjl constructed', time.time()-start)
 
    start_time = time.time()
    exit = hendrickson_leland(G, al, N, goal_n, rounds)
    end = time.time() - start_time
    print(exit, end)


if __name__ == '__main__':
    main()
