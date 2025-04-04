#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <pthread.h> //pthread library
#include "fib.c"
#include "igraph/igraph.h"

// shared memory, multithread, 
// 1) run different sizes - move to rickbobby****
// 2) parallelize
// get some data on performance
// random #s, -- test against kokkos paper graph

// Performance branches:
// runtime-whole/alg, scalability,size/iter count for varying prob sizes, portability
// NN? :"spectral properties of the graph Laplacian", 

// cut graph down to only necessary verts from benchmark

int randPriori(int);
void refreshRowStatus(int*, int, int, int);
void refreshColStatus(int*, int*, int, int, const igraph_t*);
void decideInOut(int, int*, int*, igraph_t*);
void display_graph_info(igraph_t*);
int cleanup(igraph_t*);

int IN = 0;
int OUT = BUFSIZ*128;
int UNDECIDED = 2;
int ITERATIONS = 0;

int* kokkos_mis(igraph_t* G){ // 3. Build MIS of supernodes roughly based off Kokkos paper

    assert(igraph_vector_int_size(&G->os) == igraph_vcount(G)+1);
    assert(igraph_vector_int_size(&G->is) == igraph_vcount(G)+1);

    // size of G
    igraph_integer_t numVerts = igraph_vcount(G);

    int rowStatus[(int)numVerts]; // helps track status of all vertices
    int colStatus[(int)numVerts]; // helps track status of local minima amongst neighbors
    for (int i = 0; i < (int)numVerts; i++){
        rowStatus[i] = i; 
        colStatus[i] = OUT;
    }

    //printf("rowStat: ");for (int j = 0; j < numVerts; j++){printf("%d ", rowStatus[j]);}printf("\n");
    //printf("colStat: ");for (int i = 0; i < numVerts; i++){printf("%d ", colStatus[i]);} printf("\n");

    // Initialize worklist1 & worklist2

    int* worklist1 = (int*)calloc((int)numVerts,sizeof(int)); // keeps track of UNDECIDED vertices
    int* worklist2 = (int*)calloc((int)numVerts,sizeof(int)); // keeps track of OUT vertices

    int k=0;
    for(  k = 0; k < (int)numVerts; k++){ worklist1[k] = k;}
    for(  k = 0; k < (int)numVerts; k++){ worklist2[k] = k;}
    int wrk1sz = (int)numVerts;
    int wrk2sz = (int)numVerts;

    //for (k = 0; k < wrk1sz; k++){printf("%d ", worklist1[k]);}printf("\n");
    //for ( k = 0; k < wrk2sz; k++){printf("%d ", worklist2[k]);} printf("\n");

    int iter = 0;
    int nvbits = randPriori((int)numVerts); // 3.1 Give undecided verticies random priorities

    while (wrk1sz > 0){ // while there exist undecided vertices 

        for(int i = 0; i < wrk1sz; i++){
            refreshRowStatus(rowStatus, worklist1[i], nvbits, (int)numVerts); // 3.a. Give random priorities to the undecided vertices
        }

        for(int i = 0; i < wrk2sz; i++){
            refreshColStatus(colStatus, rowStatus, worklist2[i], (int)numVerts, G); // 3.b. Mark local minima
        }

        for(int i = 0; i < wrk1sz; i++){
            decideInOut(worklist1[i], rowStatus, colStatus, G); // 3.c. Look at vertex & neighbors statuses, mark new Ins/Outs
        }

        //// 3.d. Reduce worklist1 & worklist2 ////
        int num_undecided = 0; // size
        for (int i = 0; i < wrk1sz; i++){ // tradeoff of two seq. for loops running throuch w1 and w2, or 1 going through all verts, but can check row and colstat
            if ( (rowStatus[worklist1[i]] != IN) && (rowStatus[worklist1[i]] != OUT) ){ num_undecided++; }
        }
        int not_out = 0;
        for (int i = 0; i < wrk2sz; i++){
            if (colStatus[worklist2[i]] != OUT) { not_out++; }
        }
        //// /////

        // allocate
        int j=0; k=0;

        int* w1 = (int*)calloc((int)num_undecided,sizeof(int)); // keeps track of UNDECIDED vertices
        int* w2 = (int*)calloc((int)not_out,sizeof(int)); // keeps track of OUT vertices

        //printf("wksz1: %d, wksz2: %d\n", wrk1sz, wrk2sz);
        
        for (int i = 0; i < wrk1sz; i++){ // Tradeoff of two seq. for loops running throuch w1 and w2, or 1 going through all verts, but can check row and colstat
            if ( (rowStatus[worklist1[i]] != IN) && (rowStatus[worklist1[i]] != OUT) ){ w1[j] = worklist1[i]; j++; }
        }
        for (int i = 0; i < wrk2sz; i++){
            if (colStatus[worklist2[i]] != OUT) { w2[k] = worklist2[i]; k++; }
        }

        //// ////

        worklist1 = (int*)realloc(worklist1, num_undecided*sizeof(int));
        worklist2 = (int*)realloc(worklist2, not_out*sizeof(int));
        memcpy(worklist1, w1, num_undecided*sizeof(int));
        memcpy(worklist2, w2, not_out*sizeof(int));
        wrk1sz = num_undecided;
        wrk2sz = not_out;
        free(w1); free(w2); 
        //// //// //// ////

        iter += 1;

    }
    //// 3.e. Get a list of indicies where rowStatus values == IN ////
    // size
    int num_in = 0;
    for (int i = 0; i < (int)numVerts; i++){ // tradeoff of two seq. for loops  or 1 going through all verts, but can check row and colstat
        if (rowStatus[i] == IN) { 
            num_in = num_in+1; }
    }
    // alocate
    // over allocate, have a delimeter to find the length later
    int* in_vertices = (int*)calloc((num_in+1),sizeof(int));
    k=0;
    for (int i = 0; i < (int)numVerts; i++){ 
        if (rowStatus[i] == IN) { in_vertices[k] = i; k++; }
    }
    in_vertices[num_in] = -1; // add a delimeter

    int* M1 = (int*)calloc((num_in+1),sizeof(int)); 
    memcpy(M1, in_vertices, (num_in+1)*sizeof(int));

    free(in_vertices);
    free(worklist1);
    free(worklist2);

    return M1; // return pointer to 2MIS
}

/****************************************/

void refreshRowStatus(int* rowStatus, int i, int nvBits, int numVerts){ // 3.a. Give random priorities to the undecided vertices
    // rand seed?

    int p1 = rand() % ((numVerts-1) - 1 + 1) + 1;
    int p2 = rand() % ((numVerts-1) - 1 + 1) + 1;

    int priority = p1^p2;

    int new_status = (i + 1) | (priority << nvBits);

    if(new_status == OUT){ new_status--; }

    rowStatus[i] = new_status;

    return;

}

void refreshColStatus(int* colStatus, int* rowStatus, int i, int nv, const igraph_t* G){ // 3.b. Mark local minima
    
    // !!!! get neighbors !!!! method
    //// Get neighbors of vertex i ////

    assert(igraph_vector_int_size(&G->os) == igraph_vcount(G)+1);
    assert(igraph_vector_int_size(&G->is) == igraph_vcount(G)+1);

    igraph_vector_int_t neighbors;
    igraph_vector_int_init(&neighbors, 0);

    igraph_neighbors(G, &neighbors, i, IGRAPH_ALL);
    //igraph_vector_int_print(&neighbors);

    igraph_integer_t neighbors_size = igraph_vector_int_size(&neighbors);

    int s = rowStatus[i];

    int neigh_stat;

    for(igraph_integer_t j = 0; j < neighbors_size; j++){
        if( ((int)VECTOR(neighbors)[j] < nv) && ((int)VECTOR(neighbors)[j] != i) ){ // check against neighbors: do I have the lowest priority?
            neigh_stat = rowStatus[(int)VECTOR(neighbors)[j]];

            if(neigh_stat < s){
                s = neigh_stat;
            }
        }
    }
    
    if(s == IN){ // if my neighbor is in set, I cannot be
        s = OUT;
    }

    colStatus[i] = s; // mark the local minima amongst i and neighbors

    igraph_vector_int_destroy(&neighbors);

    return;
}

void decideInOut(int i, int* rowStatus, int* colStatus, igraph_t* G){ // 3.c. Look at vertex & neighbors statuses, mark new Ins/Outs

    igraph_integer_t nv = igraph_vcount(G);

    int s = rowStatus[i];

    if( (s == IN) || (s==OUT) ){ return; }

    //// Get neighbors of vertex i ////
    // !!!! get neighbors !!!! method

    igraph_vector_int_t neighbors;
    igraph_vector_int_init(&neighbors, 0);

    igraph_neighbors(G, &neighbors, i, IGRAPH_ALL);
    // igraph_vector_int_sort(&neighbors);
    //igraph_vector_print(&neighbors);

    igraph_integer_t neighbors_size = igraph_vector_int_size(&neighbors);

    //// //// //// ////

    bool neigh_out = false;
    bool neigh_mismatchS = false;
    int neigh_stat;

    for(igraph_integer_t j=0; j < neighbors_size; j++){  // if me and my neighbors have the same priorty, I am the local minima and can be in MIS
        if((int)VECTOR(neighbors)[j] >= (int)nv){ continue; } // else: mismatch
        
        neigh_stat = colStatus[(int)VECTOR(neighbors)[j]];

        if (neigh_stat == OUT){ // if neighbor is out, i is out
            neigh_out = true;
            break;
        } else if(neigh_stat != s){
            neigh_mismatchS = true;
        }
    }

    if (neigh_out){
        rowStatus[i] = OUT;                 // Locked in parallel

    } else if ( !(neigh_mismatchS) ){
        rowStatus[i] = IN;                  // Locked in parallel
    }

    igraph_vector_int_destroy(&neighbors);

    return;

}

/****************************************/

int randPriori(int numVerts){ // 3.1 Give undecided verticies random priorities

    // alg. checks vert. + neighbors to see who has lowest priority 
    // & they become part of the MIS

    int i = numVerts + 1;
    int nvBits = 0;
    while(i > 0){
        i = (int)(i>>1);
        nvBits+=1;
    }
    return nvBits;
}

/****************************************/

void firstPass(igraph_t* G, int supernode, int m, igraph_integer_t* labels){ // 4. First Pass

    labels[m] = (igraph_integer_t)supernode;

    //// Get neighbors of vertex i ////
    // !!!! get neighbors !!!! method
    //printf("FP\n");
    igraph_vector_int_t neighbors;
    igraph_vector_int_init(&neighbors, 0);

    igraph_neighbors(G, &neighbors, m, IGRAPH_ALL);
    // igraph_vector_int_sort(&neighbors);
    //igraph_vector_print(&neighbors);

    igraph_integer_t neighbors_size = igraph_vector_int_size(&neighbors);

    //// //// //// ////

    for(igraph_integer_t i=0; i < neighbors_size; i++){ // Label all MIS vertices and their neighbors to the corresponding supernode
        if (labels[(int)VECTOR(neighbors)[i]] == -1){
            labels[(int)VECTOR(neighbors)[i]] = (igraph_integer_t)supernode; 
        } 
    }

    igraph_vector_int_destroy(&neighbors);

    return;
}

void secondPass(igraph_t* G, int unasgn, igraph_integer_t* labels){ // 5. Second Pass
    
    if (labels[unasgn] != -1){ return; } // if vertex has label,skip
    else{

    //// Get neighbors of vertex i ////
    // !!!! get neighbors !!!! method

    igraph_vector_int_t unasgn_neighs;
    igraph_vector_int_init(&unasgn_neighs, 0);

    igraph_neighbors(G, &unasgn_neighs, unasgn, IGRAPH_ALL);
    // igraph_vector_int_sort(&unasgn_neighs);
    //igraph_vector_int_print(&unasgn_neighs);

    igraph_integer_t unasgn_neighs_size = igraph_vector_int_size(&unasgn_neighs);

    //// //// //// ////

        for(igraph_integer_t i=0; i < unasgn_neighs_size; i++){ // mark undecided vertices based on if their neighbors have labels
            if (labels[(int)VECTOR(unasgn_neighs)[i]] != -1){
                // lock in parallel
                labels[unasgn] = labels[(int)VECTOR(unasgn_neighs)[i]];
                // unlock
            }
        }

        igraph_vector_int_destroy(&unasgn_neighs);
    }

    return;
}

void thirdPass(igraph_t* G, int unasgn, igraph_integer_t* labels, int* supernode_ptr){ // 6. Third Pass

    if (labels[unasgn] != -1){ return; }
    else{
    //// Get neighbors of vertex i ////
    // !!!! get neighbors !!!! method
    igraph_vector_int_t unasgn_neighs;
    igraph_vector_int_init(&unasgn_neighs, 0);

    igraph_neighbors(G, &unasgn_neighs, unasgn, IGRAPH_ALL);
    // igraph_vector_int_sort(&unasgn_neighs);
    //igraph_vector_int_print(&unasgn_neighs);

    igraph_integer_t unasgn_neighs_size = igraph_vector_int_size(&unasgn_neighs);

    //// //// //// ////
        //lock
        (*supernode_ptr)++;
        igraph_integer_t local_spn = (igraph_integer_t)*supernode_ptr; // make a new supernode
        // unlock

        labels[unasgn] = local_spn; // assign neighbors of new supernode to group [overwrites potentially]
        for(igraph_integer_t i=0; i < unasgn_neighs_size; i++){ 
            labels[(int)VECTOR(unasgn_neighs)[i]] = local_spn; 
        }
    
        igraph_vector_int_destroy(&unasgn_neighs);

    }

    return;
}

/****************************************/
igraph_t* build_edges(igraph_t* G, igraph_integer_t* labels, int spn){

    assert(igraph_cattribute_has_attr(G, IGRAPH_ATTRIBUTE_EDGE, "weight") == 1);

    const char* w = "weight";
    
    //printf("vcount: %lld\n", igraph_vcount(G));

    igraph_integer_t edge_count = igraph_ecount(G);

    //printf("ecount: %lld\n", igraph_ecount(G));

    igraph_vector_t weights;
    igraph_vector_init(&weights, 0);

    igraph_vector_int_t el;
    igraph_vector_int_init(&el, 0);
    igraph_es_t eids;
    igraph_es_all(&eids, IGRAPH_EDGEORDER_ID);

    igraph_get_edgelist(G, &el, 0);
    igraph_edges(G, eids, &el);
    igraph_cattribute_EANV(G, "weight", eids, &weights);

    igraph_integer_t i, j;
    igraph_integer_t sg_eid;

    igraph_real_t value, G_weight;

    igraph_integer_t target_label, source_label;

    cleanup(G);

    static igraph_t subgraph;
    igraph_empty(&subgraph, spn+1, IGRAPH_UNDIRECTED); //-!-

    for (i = 0, j = 0; i < edge_count; i++, j += 2) {

        target_label = labels[VECTOR(el)[j]];
        source_label = labels[VECTOR(el)[j + 1]];
        G_weight = VECTOR(weights)[i];

        //printf("i:%lld j:%lld, tl=%lld, sl=%lld\n", i,j, target_label, source_label);

        igraph_get_eid(&subgraph, &sg_eid, target_label, source_label, IGRAPH_UNDIRECTED, false);

        if(sg_eid == -1 && target_label != source_label){ // edge is not in subgraph

            //printf("!!!! %lld, %lld %f!!!!\n", target_label, source_label, G_weight);
            igraph_add_edge(&subgraph, target_label, source_label); //-!-
            igraph_get_eid(&subgraph, &sg_eid, target_label, source_label, IGRAPH_UNDIRECTED, false);
            SETEAN(&subgraph, w, sg_eid, G_weight); // -!- assign G_weight to new edge weight

        } else if(sg_eid != -1 && target_label != source_label){ // edge is in subgraph

            // get edge weight // new edge weight = G_weight + old
            value = G_weight + EAN(&subgraph, w, sg_eid);
            SETEAN(&subgraph, w, sg_eid, value);
            //printf("!! %lld, %lld %f!!\n", target_label, source_label, value);

        }

    }
    // Do not need larger graph anymore, free the memory
    // cleanup G when done with it, remove subgraph cleanup here
    
    igraph_es_destroy(&eids);
    igraph_vector_int_destroy(&el);
    igraph_vector_destroy(&weights);
    
    return &subgraph;

}
/****************************************/

igraph_t* kokkos_coarsen(igraph_t* G, igraph_integer_t N){ // 2. Kokkos Coarsen

    // list of vertices to assign to supernodes
    int* vertex_set = NULL;
    vertex_set = (int*)calloc(N,sizeof(int));
    for (int i = 0; i < N; i++){ vertex_set[i] = i; }

    //printf("here9998, %lld %lld\n",  igraph_vector_int_size(&G->os), igraph_vcount(G)+1);
    assert(igraph_vector_int_size(&G->os) == igraph_vcount(G)+1);
    assert(igraph_vector_int_size(&G->is) == igraph_vcount(G)+1);

    // 3. Find the 2-Maximal Independent Set
    int* M1= NULL;
    M1 = kokkos_mis(G); // !!!!
    if(!M1){printf("error\n");}

    int M1_length=0;                            // get the length of the MIS = # of starting supernodes
    while(M1[M1_length] != -1){ M1_length+=1;}

    printf("M1 Length: %d\n", M1_length);

    // print out MIS
    // printf("M1: "); for(int i=0;i<M1_length;i++){ printf("%d ", M1[i]); } printf("\n");

    // Create array where indices = vertices & values = supernode vertices is grouped under
    igraph_integer_t labels[N];
    memset(labels, -1, N*sizeof(igraph_integer_t));

    // 4. First Pass
    int* firstPasspntr = M1;
    for(int i = 0; i < M1_length; i++){             // Label MIS vertices and immediate neighbors with appropriate supernode label
        firstPass(G, i, firstPasspntr[i], labels);  // !!!!
    }
    free(M1);

    // cleanup vertex_set and labels == remove the vertices with labels 
    int new_len = 0;
    for(int i=0; i < N; i++){ 
        if(labels[i] == -1){ 
            new_len++;
        } 
    }

    int nvs[new_len]; 
    int k=0;
    for(int i=0; i < N; i++){ 
        if(labels[i] == -1){ 
            nvs[k] = i; k++;
        } 
    }

    vertex_set = (int*)realloc(vertex_set, new_len*sizeof(int)); // resive set of unlabeled vertices
    memcpy(vertex_set, nvs, sizeof(nvs));

    // 5. Second Pass
	new_len = sizeof(nvs)/sizeof(nvs[0]); // secure the size of the new_len, might be unnecessary


    int* secondPasspntr = vertex_set;
    for(int j=0; j < new_len; j++){                 //  Mark unlabeled vertices who are neighbors of labeled vertices
        secondPass(G, secondPasspntr[j], labels);   // !!!!
    }

    // cleanup vertex_set = only unlabeled vertices in the set
    new_len = 0;
    for(int i=0; i < N; i++){ 
        if(labels[i] == -1){ new_len++;} 
    }

    int nvs2[new_len]; 
    k=0;
    for(int i=0; i < N; i++){ 
        if(labels[i] == -1){ 
            nvs2[k] = i; 
            k++;
        } 
    }

    vertex_set = (int*)realloc(vertex_set, new_len*sizeof(int));
    memcpy(vertex_set, nvs2, sizeof(nvs2));

    int supernode = M1_length; // number of elements in M1
    int* spn_ptr = &supernode;

	// 6. Third Pass
    if(new_len > 0){
        int* thirdPasspntr = vertex_set;
        for(int i=0; i < new_len; i++){ // create new supernodes or group still unlabeled vertices
            thirdPass(G, thirdPasspntr[i], labels, spn_ptr); // !!!!
        }
    }

    printf("SPN: %d\n", supernode);
    //printf("Labels: ");for(igraph_integer_t y=0; y < N; y++){ printf("%lld ", labels[y]); }
    printf("\n");

    //7. Build Edges
    igraph_t* subgraph = build_edges(G, labels, supernode); 
    
    if (new_len) { free(vertex_set); } // free vertex set memory
    
    return subgraph; // return graph pointer
}

void display_ratios(igraph_t* G){
    printf("\n!: %" IGRAPH_PRId  ",", igraph_vcount(G));
    printf("%" IGRAPH_PRId  ",", igraph_ecount(G));
    printf("%d\n", ITERATIONS);
    return;
}

int cleanup(igraph_t* G){
	printf("NCLEANUP: %" IGRAPH_PRId "\n", igraph_vcount(G));

    igraph_destroy(G);

    return 0;

}

int recursion_fcn(igraph_t* G, int goal_verts){ // 1. Base Case/Recursive Step

    // get graph size
    igraph_integer_t graph_size = igraph_vcount(G);

    printf("\nI: %d, %" IGRAPH_PRId  ", %d\n\n", ITERATIONS, graph_size, goal_verts);

    // Base Case 
    if (graph_size <= goal_verts){
        cleanup(G); // no cleanup? -- sigabrts if you try to clean up this last graph
        return 0;
    }

    // Recursion
    igraph_t* subgraph = kokkos_coarsen(G, graph_size); // pass in pointer to G AND G's size

    ITERATIONS+=1;

    display_ratios(subgraph);

    return recursion_fcn(subgraph, goal_verts);
}

int main(int argc, char *argv[]){
    
    // random seed
    // srand(time(NULL));
    srand(0);

    igraph_set_attribute_table(&igraph_cattribute_table);

    int N = 10974;         // number of vertices in data graph

    int MAXNAME = 1024;
    char filename[MAXNAME];
    memset(filename, '\0', MAXNAME*sizeof(char) );
    strcpy(filename, "csv/bcsstk17.ncol");              // csv file to read from
    // "csv/ash85.csv",   // False, #N=85
    // "csv/netz4504.csv" // False, #N=1961
    // "csv/gemat11.csv"  // False, #N=4929
    // "csv/bcspwr10.csv" // False, #N=5300
    // "csv/bcsstk17.csv" // False, #N=10974
    // "csv/shock-9.csv"  // False, #N=36476
    // "csv/ecology1.csv" // False, #N=1000000

    igraph_t g = ncol_graph(filename, N);

    igraph_t* graph_ptr = &g;

    display_ratios(graph_ptr);
    //display_graph_info(graph_ptr); // Optional: Display initial graph

    int goal_verts = (int)(N/6) - 1;        // set minimum goal vertices
    printf("goal_verts: %d\n", goal_verts);

    clock_t start = clock();
	 
    int exit1 = recursion_fcn(graph_ptr, goal_verts); // 1. Begin recursive coarsening

    clock_t end = clock();

    double myTime = (((double)end - start)/CLOCKS_PER_SEC);
    //printf("# of Iterations: %d\n", ITERATIONS);
    //printf("Exit %d \n", exit1);
    printf("\nT: %f\n", myTime);
    //printf("Coarsened in %f seconds.\n", myTime); //

    return 0;
}