#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <cilk/cilk.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h> //pthread library
#include "convert_types.c"

int randPriori(int);
int neigh_sz(struct Graph*, int, int);
void refreshRowStatus(int*, int, int, int);
void refreshColStatus(int*, int*, int, int, struct Graph*);
void decideInOut(int, int*, int*, struct Graph*);
void display_graph_info(struct Graph*);
void cleanup(struct Graph*);

// set random seed

int IN = 0;
int OUT = BUFSIZ;
int UNDECIDED = 2;

int kokkos_mis(struct Graph* G, int* M1){ // 3. Build MIS of supernodes roughly based off Kokkos paper
    int numVerts = G->N; // size of G

    int rowStatus[numVerts];
    int colStatus[numVerts];
    for (int i = 0; i < numVerts; i++){
        rowStatus[i] = i; 
        colStatus[i] = OUT;
    }

    for (int i = 0; i < numVerts; i++){printf("%d ", rowStatus[i]);}printf("\n");
    for (int i = 0; i < numVerts; i++){printf("%d ", colStatus[i]);} printf("\n");

    // worklist1 & worklist2
    // NEEEDS RIGHT DATA STURCTURE!!!!

    int* worklist1 = (int*)calloc(numVerts,sizeof(int));
    int* worklist2 = (int*)calloc(numVerts,sizeof(int));

    //int worklist1[numVerts];
    //int worklist2[numVerts];
    int k=0;
    for(  k = 0; k < numVerts; k++){ worklist1[k] = k;}
    for(  k = 0; k < numVerts; k++){ worklist2[k] = k;}
    int wrk1sz = numVerts;
    int wrk2sz = numVerts;

    for (k = 0; k < wrk1sz; k++){printf("%d ", worklist1[k]);}printf("\n");
    for ( k = 0; k < wrk2sz; k++){printf("%d ", worklist2[k]);} printf("\n");

    int iter = 0;
    int nvbits = randPriori(numVerts);

    // potential edit of grain size right inside while loop

    while (wrk1sz > 0){ //(sizeof(worklist1)){
       
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

        // reduce worklist1 & worklist2
        int num_undecided = 0;
        for (int i = 0; i < wrk1sz; i++){ // tradeoff of two seq. for loops running throuch w1 and w2, or 1 going through all verts, but can check row and colstat
            if ( (rowStatus[worklist1[i]] != IN) && (rowStatus[worklist1[i]] != OUT) ){ num_undecided++; }
        }
        int not_out = 0;
        for (int i = 0; i < wrk2sz; i++){
            if (colStatus[worklist2[i]] != OUT) { not_out++; }
        }

        int w1[num_undecided]; int j=0;
        int w2[not_out];  k=0;

        for (int i = 0; i < wrk1sz; i++){ // tradeoff of two seq. for loops running throuch w1 and w2, or 1 going through all verts, but can check row and colstat
            if ( (rowStatus[worklist1[i]] != IN) && (rowStatus[worklist1[i]] != OUT) ){ w1[j] = worklist1[i]; j++; }
        }
        for (int i = 0; i < wrk2sz; i++){
            if (colStatus[worklist2[i]] != OUT) { w2[k] = worklist2[i]; k++; }
        }

        worklist1 = (int*)realloc(worklist1, num_undecided*sizeof(int));
        worklist2 = (int*)realloc(worklist2, not_out*sizeof(int));
        memcpy(worklist1, w1, sizeof(w1));
        memcpy(worklist2, w2, sizeof(w2));
        wrk1sz = num_undecided;
        wrk2sz = not_out;

        iter += 1;

    }

    // get a list of indicies where rowStatus values == IN

    int num_in = 0;
    for (int i = 0; i < numVerts; i++){ // tradeoff of two seq. for loops running throuch w1 and w2, or 1 going through all verts, but can check row and colstat
        if (rowStatus[i] == IN) { 
            num_in = num_in+1; }
    }

    int in_vertices[num_in];  k=0;
    for (int i = 0; i < numVerts; i++){ // tradeoff of two seq. for loops running throuch w1 and w2, or 1 going through all verts, but can check row and colstat
        if (rowStatus[i] == IN) { in_vertices[k] = i; k++; }
    }

    M1 = &(in_vertices[0]);

    free(worklist1);
    free(worklist2);

    return sizeof(in_vertices)/sizeof(in_vertices[0]);
}
/****************************************/

void refreshRowStatus(int* rowStatus, int i, int nvBits, int numVerts){

    int p1 = rand() % (numVerts-1);
    int p2 = rand() % (numVerts-1);

    int priority = p1^p2;

    int new_status = (i + 1) | (priority << nvBits);

    if(new_status == OUT){ new_status--; }

    rowStatus[i] = new_status;

    return;

}

void refreshColStatus(int* colStatus, int* rowStatus, int i, int nv, struct Graph* G){

    printf("%d\n", i);
    int start_i = binary_search(G, G->N, i);
    int neighbors_size = neigh_sz(G, i, start_i);
    int neighbors[neighbors_size];
    memset(neighbors, 0, neighbors_size*sizeof(int));

    neighbors_size = sizeof(neighbors)/sizeof(neighbors[0]);

    for(int j = 0; j < neighbors_size; j++){
        neighbors[j] = G->col_ptr[start_i+j];
    }

    int s = rowStatus[i];

    int neigh_stat;

    for(int j = 0; j < neighbors_size; j++){
        if( (neighbors[j] < nv) && (neighbors[j] != i) ){
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

void decideInOut(int i, int* rowStatus, int* colStatus, struct Graph* G){
    

    int nv = G->N; // get size of G !!!!

    int s = rowStatus[i];

    if( (s == IN) || (s==OUT) ){ return; }

    int start_i = binary_search(G, G->N, i);
    int neighbors_size = neigh_sz(G, i, start_i);
    int neighbors[neighbors_size];
    memset(neighbors, 0, neighbors_size*sizeof(int));

    neighbors_size = sizeof(neighbors)/sizeof(neighbors[0]);

    for(int j = 0; j < neighbors_size; j++){
        neighbors[j] = G->col_ptr[start_i+j];
        printf("i: %d, neigh[j]: %d == s: %d, colStat: %d\n", i, neighbors[j], s, colStatus[neighbors[j]]);
    }

    bool neigh_out = false;
    bool neigh_mismatchS = false;
    int neigh_stat;

    for(int j=0; j < neighbors_size; j++){

        if(neighbors[j] >= nv){ continue; }
        neigh_stat = colStatus[neighbors[j]];

        if (neigh_stat == OUT){
            neigh_out = true;
            break;
        } else if(neigh_stat != s){
            neigh_mismatchS = true;
        }
    }

    pthread_mutex_t m; //define the lock

    if (neigh_out){
        pthread_mutex_lock(&m); // lock
        rowStatus[i] = OUT;                 // !!!! Does this need to be locked??
        pthread_mutex_unlock(&m);// unlock
    } else if ( !(neigh_mismatchS) ){
        pthread_mutex_lock(&m);// lock
        rowStatus[i] = IN;
        pthread_mutex_unlock(&m);// unlock
    }


    return;

}

/****************************************/

int randPriori(int numVerts){
    int i = numVerts + 1;
    int nvBits = 0;
    while(i > 0){
        i = (int)(i>>1);
        nvBits+=1;
    }
    return nvBits;
}

int neigh_sz(struct Graph* G, int v, int start_i){
    int end_i = 1;

    // could do another round of binary, but will just be linear here
    while(G->row_ptr[start_i+end_i] == G->row_ptr[start_i] && (start_i + end_i < G->N)){
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

    printf("neigh2: ");
    for(int i = 0; i < neigh_sz; i++){
        printf("%d,", neighbors[i]);
    }
    printf("\n");

    return neighbors; // !!!!
}

/****************************************/

/****************************************/

void kokkos_coarsen(struct Graph* G, int N){ // 2. Kokkos Coarsen

    int vertex_set[N]; // list of vertices to assign to supernodes
    for (int i = 0; i < N; i++){ vertex_set[i] = i; }

    int* M1;
    int M1_length= kokkos_mis(G, M1);

    // works?
    printf("M1: "); for(int i=0;i<M1_length;i++){ printf("%d ", M1[i]); } printf("\n");
    

    //if (G->row_ptr){cleanup(G);}

    return;
}

void display_graph_info(struct Graph* G){
    printf("Graph info:\n");
    printf("N: %d\n# Edges: %d\n", G->N, G->active_index);
    printf("Rows: ");for(int i=0; i < G->active_index;i++){ printf("%d ", G->row_ptr[i]);} printf("\n");
    printf("Cols: ");for(int i=0; i < G->active_index;i++){ printf("%d ", G->col_ptr[i]);} printf("\n");
    printf("Vals: ");for(int i=0; i < G->active_index;i++){ printf("%d ", G->value_ptr[i]);} printf("\n");
    return;
}

void cleanup(struct Graph* G){

    free(G->row_ptr);
    free(G->col_ptr);
    free(G->value_ptr);

    return;

}

void recursion_fcn(struct Graph* G, int goal_verts){ // double check pointer passing

    // get graph size
    int graph_size = G->N;

    // Base Case 
    if (graph_size <= goal_verts){
        return;
    }

    // Recursion
    kokkos_coarsen(G, graph_size); // pass in pointer to G

    //if (G->row_ptr) {cleanup(G);}

    return;
}

int main(){
    int MAXNAME = 1024;
    struct timeval start,end;

    int N = 6;

    // load G in from a csv file
    // "csv/simple_graph_000.csv"  //True,  # N=5
    char filename[MAXNAME];
    memset( filename, 0, MAXNAME*sizeof(char) );
    strcpy(filename, "csv/kk_simpleEx.csv");  //True,  # N=6
    // "csv/simple_graph_001.csv", //True,  # N=7
    // "csv/west0067.csv",         //False, # N=67

    struct Graph* graph_ptr = csv_to_graph(filename, N);
    // G = ptr to graph struct w/ csr data

    int goal_verts = (int)(N/2) - 1;
    printf("goal_verts: %d\n", goal_verts);

    gettimeofday(&start, NULL);

    recursion_fcn(graph_ptr, goal_verts);
    int exit = 0;
    gettimeofday(&end, NULL);

    double myTime = (end.tv_sec+ (double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);

    printf("Exit %d, Coarsened in %lf seconds.\n", exit, myTime);

    if (graph_ptr->row_ptr) { cleanup(graph_ptr); }

    return 0;
} 