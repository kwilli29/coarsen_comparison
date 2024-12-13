#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <cilk/cilk.h>
#include <time.h>
#include <math.h>
#include <pthread.h> //pthread library
#include "convert_types.c"

// shared memory, multithread, 
// 1) run different sizes - move to ricky-bobby****
// 2) parallelize
// get some data on performance
// random #s, -- test against kokkos paper graph

// Performance branches:
// runtime-whole/alg, scalability,size/iter count for varying prob sizes, portability
// NN? :"spectral properties of the graph Laplacian", 

// cut graph down to only necessary verts from benchmark

int randPriori(int);
int neigh_sz(struct Graph*, int, int);
void refreshRowStatus(int*, int, int, int);
void refreshColStatus(int*, int*, int, int, struct Graph*);
void decideInOut(int, int*, int*, struct Graph*);
void display_graph_info(struct Graph*);
int cleanup(struct Graph*);

// set random seed

int IN = 0;
int OUT = BUFSIZ*128;
int UNDECIDED = 2;
int ITERATIONS = 0;

pthread_mutex_t m0; //define the lock
pthread_mutex_t m2; //define the lock
pthread_mutex_t m3; //define the lock

int* kokkos_mis(struct Graph* G){ // 3. Build MIS of supernodes roughly based off Kokkos paper
    int numVerts = G->N;   // size of G

    int rowStatus[numVerts]; // helps track status of all vertices
    int colStatus[numVerts]; // helps track status of local minima amongst neighbors
    for (int i = 0; i < numVerts; i++){
        rowStatus[i] = i; 
        colStatus[i] = OUT;
    }

    //printf("rowStat: ");for (int j = 0; j < numVerts; j++){printf("%d ", rowStatus[j]);}printf("\n");
    //printf("colStat: ");for (int i = 0; i < numVerts; i++){printf("%d ", colStatus[i]);} printf("\n");

    // Initialize worklist1 & worklist2

    int* worklist1 = (int*)calloc(numVerts,sizeof(int)); // keeps track of UNDECIDED vertices
    int* worklist2 = (int*)calloc(numVerts,sizeof(int)); // keeps track of OUT vertices

    int k=0;
    for(  k = 0; k < numVerts; k++){ worklist1[k] = k;}
    for(  k = 0; k < numVerts; k++){ worklist2[k] = k;}
    int wrk1sz = numVerts;
    int wrk2sz = numVerts;

    //for (k = 0; k < wrk1sz; k++){printf("%d ", worklist1[k]);}printf("\n");
    //for ( k = 0; k < wrk2sz; k++){printf("%d ", worklist2[k]);} printf("\n");

    int iter = 0;
    int nvbits = randPriori(numVerts); // 3.1 Give undecided verticies random priorities

    // !! Potential edit of grain size right inside while loop !!

    while (wrk1sz > 0){  // while there exist undecided vertices 
       
        cilk_for(int i = 0; i < wrk1sz; i++){
            refreshRowStatus(rowStatus, worklist1[i], nvbits, numVerts);
        }
        
        cilk_sync; // sync? !!!!
        
        cilk_for(int i = 0; i < wrk2sz; i++){
            refreshColStatus(colStatus, rowStatus, worklist2[i], numVerts, G);
        } 

        cilk_sync; // sync?

        cilk_for(int i = 0; i < wrk1sz; i++){
            decideInOut(worklist1[i], rowStatus, colStatus, G);
        }
        cilk_sync; // sync

        //// 3.d. Reduce worklist1 & worklist2 ////
        // size
        int num_undecided = 0;
        for (int i = 0; i < wrk1sz; i++){ // tradeoff of two seq. for loops running throuch w1 and w2, or 1 going through all verts, but can check row and colstat
            if ( (rowStatus[worklist1[i]] != IN) && (rowStatus[worklist1[i]] != OUT) ){ num_undecided++; }
        }
        int not_out = 0;
        for (int i = 0; i < wrk2sz; i++){
            if (colStatus[worklist2[i]] != OUT) { not_out++; }
        }
        //// /////
        // allocate
        int w1[num_undecided]; int j=0;
        int w2[not_out];  k=0;

        for (int i = 0; i < wrk1sz; i++){ // Tradeoff of two seq. for loops running throuch w1 and w2, or 1 going through all verts, but can check row and colstat
            if ( (rowStatus[worklist1[i]] != IN) && (rowStatus[worklist1[i]] != OUT) ){ w1[j] = worklist1[i]; j++; }
        }
        for (int i = 0; i < wrk2sz; i++){
            if (colStatus[worklist2[i]] != OUT) { w2[k] = worklist2[i]; k++; }
        }

        //// ////

        worklist1 = (int*)realloc(worklist1, num_undecided*sizeof(int));
        worklist2 = (int*)realloc(worklist2, not_out*sizeof(int));
        memcpy(worklist1, w1, sizeof(w1));
        memcpy(worklist2, w2, sizeof(w2));
        wrk1sz = num_undecided;
        wrk2sz = not_out;
        //// //// //// ////

        iter += 1;

    }

    //// 3.e. Get a list of indicies where rowStatus values == IN ////
    // size
    int num_in = 0;
    for (int i = 0; i < numVerts; i++){ // tradeoff of two seq. for loops  or 1 going through all verts, but can check row and colstat
        if (rowStatus[i] == IN) { 
            num_in = num_in+1; }
    }
    // alocate
    int in_vertices[num_in+1];  // over allocate, have a delimeter to find the length later
    k=0;
    for (int i = 0; i < numVerts; i++){ 
        if (rowStatus[i] == IN) { in_vertices[k] = i; k++; }
    }
    in_vertices[num_in] = -1; // add a delimeter


    int* M1 = (int*)calloc((num_in+1),sizeof(int)); 
    memcpy(M1, in_vertices, sizeof(in_vertices));

    free(worklist1);
    free(worklist2);

    return M1; // return pointer to 2MIS
}
/****************************************/

void refreshRowStatus(int* rowStatus, int i, int nvBits, int numVerts){ // 3.a. Give random priorities to the undecided vertices

    int p1 = rand() % (numVerts-1);
    int p2 = rand() % (numVerts-1);

    int priority = p1^p2;

    int new_status = (i + 1) | (priority << nvBits);

    if(new_status == OUT){ new_status--; }

    rowStatus[i] = new_status;

    return;

}

void refreshColStatus(int* colStatus, int* rowStatus, int i, int nv, struct Graph* G){ // 3.b. Mark local minima

    //// Get neighbors of vertex i ////
    int start_i = binary_search(G, G->active_index, i);
    int neighbors_size = neigh_sz(G, i, start_i);
    int neighbors[neighbors_size];
    memset(neighbors, 0, neighbors_size*sizeof(int));

    neighbors_size = sizeof(neighbors)/sizeof(neighbors[0]);

    for(int j = 0; j < neighbors_size; j++){
			if((start_i + j) < G->active_index){
        		neighbors[j] = G->col_ptr[start_i+j];
			}
    }
    //// //// //// ////

    int s = rowStatus[i];

    int neigh_stat;

    for(int j = 0; j < neighbors_size; j++){
        if( (neighbors[j] < nv) && (neighbors[j] != i) ){ // check against neighbors: do I have the lowest priority?
            neigh_stat = rowStatus[neighbors[j]];

            if(neigh_stat < s){
                s = neigh_stat;
            }
        }
    }

    if(s == IN){
        s = OUT;
    }

    colStatus[i] = s;

    return;
}

void decideInOut(int i, int* rowStatus, int* colStatus, struct Graph* G){ // 3.c. Look at vertex & neighbors statuses, mark new Ins/Outs
    
    int nv = G->N; // get size of G

    int s = rowStatus[i];

    if( (s == IN) || (s==OUT) ){ return; }

    //// Get neighbors of vertex i ////
    int start_i = binary_search(G, G->active_index, i);
    int neighbors_size = neigh_sz(G, i, start_i);
    int neighbors[neighbors_size];
    memset(neighbors, 0, neighbors_size*sizeof(int));

    neighbors_size = sizeof(neighbors)/sizeof(neighbors[0]);

    for(int j = 0; j < neighbors_size; j++){
        if((start_i + j) < G->active_index){
            neighbors[j] = G->col_ptr[start_i+j];
            //printf("i: %d, neigh[j]: %d == s: %d, colStat: %d\n", i, neighbors[j], s, colStatus[neighbors[j]]);
        }
    }
    //// //// //// ////

    bool neigh_out = false;
    bool neigh_mismatchS = false;
    int neigh_stat;

    for(int j=0; j < neighbors_size; j++){ // if me and my neighbors have the same priorty, I am the local minima and can be in MIS
        if(neighbors[j] >= nv){ continue; } // else: mismatch

        neigh_stat = colStatus[neighbors[j]];

        if (neigh_stat == OUT){ // if neighbor is out, i is out
            neigh_out = true;
            break;
        } else if(neigh_stat != s){
            neigh_mismatchS = true;
        }
    }

    //pthread_mutex_t m; //define the lock
    

    if (neigh_out){
        pthread_mutex_lock(&m0); // lock
        rowStatus[i] = OUT;                 // !!!! Does this need to be locked??
        pthread_mutex_unlock(&m0);// unlock
    } else if ( !(neigh_mismatchS) ){
        pthread_mutex_lock(&m0);// lock
        rowStatus[i] = IN;
        pthread_mutex_unlock(&m0);// unlock
    }



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

int neigh_sz(struct Graph* G, int v, int start_i){ // Get size of neighbor
    
    int end_i = 1;

    while((start_i + end_i < G->active_index) && G->row_ptr[start_i+end_i] == G->row_ptr[start_i]){
        end_i++;
    }

    end_i = start_i + end_i - 1;

    return (end_i - start_i + 1);

}

int* adj(int neighbors[], int neigh_sz, int start_i, struct Graph* G){

    // neigh_sz = start ---- end

    for(int i = 0; i < neigh_sz; i++){
        neighbors[i] = G->col_ptr[start_i+i];
        //printf(" %d,", neighs[i]);
    }
    //printf("\n");

    return neighbors; // !!!!
}

/****************************************/

void firstPass(struct Graph* G, int supernode, int m, int* labels){ // 4. First Pass

    labels[m] = supernode;

    //// Get neighbors of vertex m////
    int start_i = binary_search(G, G->active_index, m);
    int neighbors_size = neigh_sz(G, m, start_i);
    int neighbors[neighbors_size];
    memset(neighbors, 0, neighbors_size*sizeof(int));

    neighbors_size = sizeof(neighbors)/sizeof(neighbors[0]);

    for(int i = 0; i < neighbors_size; i++){
		if((start_i+i)<G->active_index){
            neighbors[i] = G->col_ptr[start_i+i];
		}
    }
    //// //// //// ////

    for(int i=0; i < neighbors_size; i++){ // Label all MIS vertices and their neighbors to the corresponding supernode
        if (labels[neighbors[i]] == -1){
            labels[neighbors[i]] = supernode; 
        } 
    }

    return;
}

void secondPass(struct Graph* G, int unasgn, int* labels){ // 5. Second Pass
    
    if (labels[unasgn] != -1){ return; } // if vertex has label,skip
    else{
        //// Get neighbors of vertex unasgn ////
        int start_i = binary_search(G, G->active_index, unasgn);
        int unasgn_neighs_size = neigh_sz(G, unasgn, start_i);
        int unasgn_neighs[unasgn_neighs_size];
        memset(unasgn_neighs, 0, unasgn_neighs_size*sizeof(int));

        unasgn_neighs_size = sizeof(unasgn_neighs)/sizeof(unasgn_neighs[0]);


        for(int i = 0; i < unasgn_neighs_size; i++){
            if( (start_i + i) < G->active_index){	
                unasgn_neighs[i] = G->col_ptr[start_i+i];
            }
        }
        //// //// //// ////

        //pthread_mutex_t m; //define the lock

        for(int i=0; i < unasgn_neighs_size; i++){ // mark undecided vertices based on if their neighbors have labels
            if (labels[unasgn_neighs[i]] != -1){
                pthread_mutex_lock(&m2); // lock // shared? !!!!
                labels[unasgn] = labels[unasgn_neighs[i]];
                pthread_mutex_unlock(&m2); // unlock
            }
        }

    }
    
    return;
}

void thirdPass(struct Graph* G, int unasgn, int* labels, int* supernode_ptr){ // 6. Thid Pass

    if (labels[unasgn] != -1){ return; }
    else{
        //// Get neighbors of vertex unasgn ////
        int start_i = binary_search(G, G->active_index, unasgn);
        int unasgn_neighs_size = neigh_sz(G, unasgn, start_i);
        int unasgn_neighs[unasgn_neighs_size];
        memset(unasgn_neighs, 0, unasgn_neighs_size*sizeof(int));

        unasgn_neighs_size = sizeof(unasgn_neighs)/sizeof(unasgn_neighs[0]);

        for(int i = 0; i < unasgn_neighs_size; i++){
            if((start_i+i)<G->active_index){
                unasgn_neighs[i] = G->col_ptr[start_i+i];
            }
        }
        //// //// //// ////

        pthread_mutex_t m; //define the lock
        // so several threads don't think they have the same supernode label

        pthread_mutex_lock(&m3); // lock // shared? !!!!
        (*supernode_ptr)++;
        int local_spn = *supernode_ptr; // make a new supernode
        pthread_mutex_unlock(&m3); // unlock

        labels[unasgn] = local_spn; // assign neighbors of new supernode to group [overwrites potentially]
        for(int i=0; i < unasgn_neighs_size; i++){ 
            labels[unasgn_neighs[i]] = local_spn; 
        } //!!!! 
    
    }

    return;
}

/****************************************/

struct Graph* build_edges_1D(struct Graph* G, int* labels, int spn){
    //printf("Labels ");for(int k=0;k<G->N;k++){ printf("%d ", labels[k]);}printf("\n");
    //printf("\nspn: %d\n", spn);

    int index1 = 0; // actual number of edges
    int row_spn=-1;
    int col_spn=-1;

    for(int i=0; i < G->active_index; i++){ // potentially parallelizable, but would need locks

        row_spn =  labels[G->row_ptr[i]];  // edge (row[i],col[i]) have same label, else an edge needs to be marked
        col_spn =  labels[G->col_ptr[i]];
        
        if(row_spn != col_spn){
            index1+=1;
        }
            
    }
    // calloc Edges Matrix
    float* edge_matrix = (float*)calloc(index1*index1,sizeof(float));

    // how to accumulate value correctly

    int index = 0;
    float val = 0.0;

    for(int i=0; i < G->active_index; i++){ // potentially parallelizable, but would need locks

        row_spn =  labels[G->row_ptr[i]];  // edge (row[i],col[i]) have same label, else an edge needs to be marked
        col_spn =  labels[G->col_ptr[i]];
        val = G->value_ptr[i];
        
        if(row_spn != col_spn){
            // is edge already logged in subgraph

            // no  = add it
            if ( (row_spn<spn) && (col_spn<spn) && edge_matrix[(row_spn*index1)+col_spn]==0 && edge_matrix[(col_spn*index1)+row_spn]==0 )  {
                
                edge_matrix[(row_spn*index1)+col_spn] = val;
                edge_matrix[(col_spn*index1) +row_spn] = val;
                index+=2;
                //index++;
                
            // yes = keep track of edge weights b/w supernodes only
            } else {
                edge_matrix[(row_spn*index1)+col_spn] = edge_matrix[(row_spn*index1)+col_spn] + val;
                edge_matrix[(col_spn*index1) +row_spn] = edge_matrix[(col_spn*index1) +row_spn] + val;
            }
        }
            
    }

    // Convert the edge matrix to the subgraph csr

	cleanup(G); // Do not need larger graph anymore, free the memory

	static struct Graph subgraph; // create new subgraph
    subgraph.N = spn;
	subgraph.active_index = index;

	subgraph.row_ptr = (int*)calloc(index,sizeof(int));
	subgraph.col_ptr  = (int*)calloc(index,sizeof(int));
	subgraph.value_ptr = (float*)calloc(index,sizeof(float));

    int ind = 0;

    // mark row,col,val based off edge matrix

    for(int i=0; i < index1*index1; i++){ // can also make edge_matrix not symmetric and only update/go through half of it for this part
        //for(int j=0; j < spn; j++){
            
            if ( (edge_matrix[i] == 0) || (edge_matrix[i] == 0.0) || !(edge_matrix[i]) || ind >= index){ continue; }
            else{

                subgraph.row_ptr[ind] = i%index1;
                subgraph.col_ptr[ind] = (int)(i/index1);
                subgraph.value_ptr[ind] = edge_matrix[i];
                ind = ind + 1;
                
            }
            
        //}
    }
    free(edge_matrix);

    quickSort(subgraph.row_ptr, 0, index-1, subgraph.col_ptr, subgraph.value_ptr); // sort the CSR
    
    return (&subgraph); // return pointer to subgraph

}

struct Graph* build_edges(struct Graph* G, int* labels, int spn){ // 7. Build Edges

    //printf("Labels ");for(int k=0;k<G->N;k++){ printf("%d ", labels[k]);}printf("\n");
    //printf("\nspn: %d\n", spn);

    // 2D Edges Matrix
    float edge_matrix[spn][spn];
    for(int k=0; k < spn; k++){  memset(edge_matrix[k],0.0,spn*sizeof(float)); }

    int index = 0; // actual number of edges
    int row_spn=-1;
    int col_spn=-1;
    float val = 0.0;

    for(int i=0; i < G->N; i++){ // potentially parallelizable, but would need locks

        row_spn =  labels[G->row_ptr[i]];  // edge (row[i],col[i]) have same label, else an edge needs to be marked
        col_spn =  labels[G->col_ptr[i]];
        val = G->value_ptr[i];
        
        if(row_spn != col_spn){
            // is edge already logged in subgraph

            // no  = add it
            if ( (row_spn<spn) && (col_spn<spn) && (edge_matrix[row_spn][col_spn]==0) && (edge_matrix[col_spn][row_spn]==0) ) {
                
                edge_matrix[row_spn][col_spn] = val;
                edge_matrix[col_spn][row_spn] = val;
                index+=2;
                
            // yes = keep track of edge weights b/w supernodes only
            } else {
                edge_matrix[row_spn][col_spn] = edge_matrix[row_spn][col_spn] + val;
                edge_matrix[col_spn][row_spn] = edge_matrix[col_spn][row_spn] + val;
            }
        }
            
    }

    // Convert the edge matrix to the subgraph csr

	cleanup(G); // Do not need larger graph anymore, free the memory

	static struct Graph subgraph; // create new subgraph
    subgraph.N = spn;
	subgraph.active_index = index;

	subgraph.row_ptr = (int*)calloc(index,sizeof(int));
	subgraph.col_ptr  = (int*)calloc(index,sizeof(int));
	subgraph.value_ptr = (float*)calloc(index,sizeof(float));

    int ind = 0;

    // mark row,col,val based off edge matrix

    for(int i=0; i < spn; i++){ // can also make edge_matrix not symmetric and only update/go through half of it for this part
        for(int j=0; j < spn; j++){
            
            if ( (edge_matrix[i][j] == 0) || !(edge_matrix[i][j]) || ind >= index){ continue; }
            else{

                subgraph.row_ptr[ind] = i;
                subgraph.col_ptr[ind] = j;
                subgraph.value_ptr[ind] = edge_matrix[i][j];
                ind = ind + 1;
                
            }
            
        }
    }

    quickSort(subgraph.row_ptr, 0, index-1, subgraph.col_ptr, subgraph.value_ptr); // sort the CSR
    
    return (&subgraph); // return pointer to subgraph   
}

/****************************************/

struct Graph* kokkos_coarsen(struct Graph* G, int N){ // 2. Kokkos Coarsen

    // list of vertices to assign to supernodes
    int* vertex_set = NULL;
    vertex_set = (int*)calloc(N,sizeof(int));
    for (int i = 0; i < N; i++){ vertex_set[i] = i; }

    // 3. Find the 2-Maximal Independent Set
    int* M1= NULL;
    M1 = kokkos_mis(G);
    if(!M1){printf("error\n");}

    int M1_length=0;                            // get the length of the MIS = # of starting supernodes
    while(M1[M1_length] != -1){ M1_length+=1;}

    printf("M1 Length: %d,\n", M1_length);

    // print out MIS
    //printf("M1: "); for(int i=0;i<M1_length;i++){ printf("%d ", M1[i]); } printf("\n");
    
    // Create array where indices = vertices & values = supernode vertices is grouped under
    int labels[N];
    memset(labels, -1, N*sizeof(int));

    
    // 4. First Pass
    int* firstPasspntr = M1;
    cilk_for(int i = 0; i < M1_length; i++){
        firstPass(G, i, firstPasspntr[i], labels);
    }

    cilk_sync; // sync
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
    cilk_for(int j=0; j < new_len; j++){ //  Mark unlabeled vertices who are neighbors of labeled vertices
        secondPass(G, secondPasspntr[j], labels);
    }

    cilk_sync; // sync

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
    int* thirdPasspntr = vertex_set;
    cilk_for(int i=0; i < new_len; i++){ // create new supernodes or group still unlabeled vertices
        thirdPass(G, thirdPasspntr[i], labels, spn_ptr);
    }

    cilk_sync; // sync

    // 7. Build the subgraph edge relationships based on the original graph
    struct Graph* subgraph;
    if (N >= 16000){
        subgraph = build_edges_1D(G, labels, supernode);
    }else {
        subgraph = build_edges(G, labels, supernode);
    }

    if (new_len) { free(vertex_set); } // free vertex set memory
    
    return subgraph; // return graph pointer
}

void display_graph_info(struct Graph* G){
    printf("Graph info:\n");
    printf("N: %d\n# Edges: %d\n", G->N, G->active_index);
    printf("Rows: ");for(int i=0; i < G->active_index;i++){ printf("%d ", G->row_ptr[i]);} printf("\n");
    printf("Cols: ");for(int i=0; i < G->active_index;i++){ printf("%d ", G->col_ptr[i]);} printf("\n");
    printf("Vals: ");for(int i=0; i < G->active_index;i++){ printf("%lf ", G->value_ptr[i]);} printf("\n");
    return;
}

int cleanup(struct Graph* G){
    printf("NCLEANUP: %d\n", G->N);
    free(G->row_ptr);
    free(G->col_ptr);
    free(G->value_ptr);

    return 0;

}

int recursion_fcn(struct Graph* G, int goal_verts){ // 1. Base Case/Recursive Step

    // get graph size
    int graph_size = G->N;

    // Base Case 
    if (graph_size <= goal_verts){
        cleanup(G); // no cleanup? -- sigabrts if you try to clean up this last graph
        return 0;
    }

    // Recursion
    struct Graph* subgraph = kokkos_coarsen(G, graph_size); // pass in pointer to G AND G's size

    // print subgraph & number of vertices coarsened
    display_graph_info(subgraph);
    printf("Vertices coarsened: %d - %d = %d\n", graph_size, (graph_size - subgraph->N), subgraph->N);

    ITERATIONS+=1;

    return recursion_fcn(subgraph, goal_verts);
}

int main(){

    // random seed
    // srand(time(NULL));
    srand(0);

    if (pthread_mutex_init(&m0, NULL) != 0) { perror("Mutex initialization failed"); return 1;}
    if (pthread_mutex_init(&m2, NULL) != 0) { perror("Mutex initialization failed"); return 1;}
    if (pthread_mutex_init(&m3, NULL) != 0) { perror("Mutex initialization failed"); return 1;}

    int MAXNAME = 1024; // max length of a file name

    int N = 85;        // number of vertices in data graph
    bool g_sym = false;  // is the data symmetric

    // load G in from a csv file
    char filename[MAXNAME];
    memset( filename, '\0', MAXNAME*sizeof(char) );
    strcpy(filename, "csv/ash85.csv"); // csv file to read from
     // "csv/simple_graph_000.csv"  //True,  # N=5
    // "csv/kk_simpleEx.csv",       //True,  # N=6
    // "csv/simple_graph_001.csv",  //True,  # N=7
    // "csv/ash85.csv",         //False, # N=85
    // "csv/netz4504.csv" //False, # N=1961
    // "csv/gemat11.csv" // False, #N=4929
    // "csv/bcspwr10.csv" // False, # N=5300
    // "csv/bcsstk17.csv" // False, #N=10974
    // "csv/shock-9.csv" // False, #N=36476

    struct Graph* graph_ptr; // G = ptr to graph struct w/ csr data
    if (N < 100){
        graph_ptr = csv_to_graph_small(filename, N, g_sym); // smaller graph
    } else {
        graph_ptr = csv_to_graph_large(filename, N, g_sym); // bigger graph
    }

    //display_graph_info(graph_ptr); // Optional: Display initial graph

    int goal_verts = (int)(N/2) - 1; // set minimum goal vertices
    printf("goal_verts: %d\n", goal_verts);

    clock_t start = clock();

    int exit1 = recursion_fcn(graph_ptr, goal_verts); // 1. Begin recursive coarsening

    clock_t end = clock();

    double myTime = (((double)end - start)/CLOCKS_PER_SEC);
    printf("# of Iterations: %d\n", ITERATIONS);
    printf("Exit %d \n", exit1);
    printf("Coarsened in %f seconds.\n", myTime); //


    pthread_mutex_destroy(&m0);
    pthread_mutex_destroy(&m2);
    pthread_mutex_destroy(&m3);

    return 0;
} 