#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <cilk/cilk.h>
#include <assert.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h> //pthread library
#include <string.h>
#include "zconvert_types.c"

// CSR
// neighbors DS
// some sort of tree?

// Performance branches:
// runtime-whole/alg, scalability,size/iter count for varying prob sizes, portability
// NN? :"spectral properties of the graph Laplacian", 

// cut graph down to only necessary verts from benchmark

int randPriori(int);
void refreshRowStatus(int*, int, int, int);
void refreshColStatus(int*, int*, int, int, const struct CSR*);
void decideInOut(int, int*, int*, struct CSR*);
// void display_graph_info(struct CSR*);
struct CSR* display_ratio(struct CSR* G);
int cleanup(struct CSR*);

int IN = 0;
int OUT = BUFSIZ*128;
int UNDECIDED = 2;
int ITERATIONS = 0;

pthread_mutex_t m0; //define the lock
pthread_mutex_t m2; //define the lock
pthread_mutex_t m3; //define the lock
pthread_mutex_t m4; //define the lock

int* kokkos_mis(struct CSR* G){ // 3. Build MIS of supernodes roughly based off Kokkos paper

    int numVerts = G->N;   // size of G

    //if (numVerts <= 1000000){
    //    callocstats = false;
    //     int rowStatus[(int)numVerts]; // helps track status of all vertices
    //   int colStatus[(int)numVerts]; // helps track status of local minima amongst neighbors
    //} 
    struct timeval start_time,end_time;
    gettimeofday(&start_time, NULL);int* rowStatus = (int*)calloc(numVerts,sizeof(int));gettimeofday(&end_time, NULL);
    printf("CALLOC RowStatus: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
    gettimeofday(&start_time, NULL);int* colStatus = (int*)calloc(numVerts,sizeof(int));gettimeofday(&end_time, NULL);
    printf("CALLOC ColStatus: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));

    for (int i = 0; i < numVerts; i++){
        rowStatus[i] = i; 
        colStatus[i] = OUT;
    }

    //printf("rowStat: \n");//for (int j = 0; j < numVerts; j++){printf("%d ", rowStatus[j]);}printf("\n"); //printf("colStat: ");for (int i = 0; i < numVerts; i++){printf("%d ", colStatus[i]);} printf("\n");

    // Initialize worklist1 & worklist2
    gettimeofday(&start_time, NULL);int* worklist1 = (int*)calloc(numVerts,sizeof(int));gettimeofday(&end_time, NULL); // keeps track of UNDECIDED vertices
    printf("CALLOC Wrkl1: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
    gettimeofday(&start_time, NULL);int* worklist2 = (int*)calloc(numVerts,sizeof(int));gettimeofday(&end_time, NULL); // keeps track of OUT vertices
    printf("CALLOC Wrkl2: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));

    int k=0;
    for(  k = 0; k < numVerts; k++){ worklist1[k] = k;}
    for(  k = 0; k < numVerts; k++){ worklist2[k] = k;}
    int wrk1sz = numVerts;
    int wrk2sz = numVerts;

    //for (k = 0; k < wrk1sz; k++){printf("%d ", worklist1[k]);}printf("\n"); //for (k = 0; k < wrk2sz; k++){printf("%d ", worklist2[k]);} printf("\n");

    int iter = 0;
    int nvbits = randPriori(numVerts); // 3.1 Give undecided verticies random priorities

    // !! Potential edit of grain size right inside while loop !!

    while (wrk1sz > 0){  // while there exist undecided vertices 
       
        cilk_for(int i = 0; i < wrk1sz; i++){
            refreshRowStatus(rowStatus, worklist1[i], nvbits, numVerts);
        }
        
        cilk_sync; // sync
        
        cilk_for(int i = 0; i < wrk2sz; i++){
            refreshColStatus(colStatus, rowStatus, worklist2[i], numVerts, G);
        } 

        cilk_sync; // sync

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
        // allocate //
        int j=0; k=0;
        gettimeofday(&start_time, NULL);int* w1 = (int*)calloc((int)num_undecided,sizeof(int)); gettimeofday(&end_time, NULL); // keeps track of UNDECIDED vertices
        printf("CALLOC-Re W1: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
        gettimeofday(&start_time, NULL); int* w2 = (int*)calloc((int)not_out,sizeof(int)); gettimeofday(&end_time, NULL); // keeps track of OUT vertices
        printf("CALLOC-Re W2: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));

        for (int i = 0; i < wrk1sz; i++){
            if ( (rowStatus[worklist1[i]] != IN) && (rowStatus[worklist1[i]] != OUT) ){ w1[j] = worklist1[i]; j++; }
        }
        for (int i = 0; i < wrk2sz; i++){
            if (colStatus[worklist2[i]] != OUT) { w2[k] = worklist2[i]; k++; }
        }
        // copy memory to the newly shortened lists //
        // !!!! Compact worklists with parallel prefix sums
        gettimeofday(&start_time, NULL);worklist1 = (int*)realloc(worklist1, num_undecided*sizeof(int));gettimeofday(&end_time, NULL);
        printf("C-REALLOC-Re W1: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
        gettimeofday(&start_time, NULL);worklist2 = (int*)realloc(worklist2, not_out*sizeof(int));gettimeofday(&end_time, NULL);
        printf("C-REALLOC-Re W2: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
        memcpy(worklist1, w1, num_undecided*sizeof(int));
        memcpy(worklist2, w2, not_out*sizeof(int));
        wrk1sz = num_undecided;
        wrk2sz = not_out;
        free(w1); free(w2); 
        //// //// //// ////

        iter += 1;

    }

    //// 3.e. Get a list of indicies where rowStatus values == IN ////
    // size //
    int num_in = 0;
    for (int i = 0; i < numVerts; i++){ 
        if (rowStatus[i] == IN) { 
            num_in = num_in+1; }
    }
    // alocate // 
    
    gettimeofday(&start_time, NULL);
    int* in_vertices = (int*)calloc((num_in+1),sizeof(int));gettimeofday(&end_time, NULL);  // over allocate, have a delimeter to find the length later
    printf("CALLOC in-vertices: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
    k=0;
    for (int i = 0; i < numVerts; i++){ 
        if (rowStatus[i] == IN) { in_vertices[k] = i; k++; }
    }
    in_vertices[num_in] = -1; // add a delimeter to later know length

    gettimeofday(&start_time, NULL);int* M1 = (int*)calloc((num_in+1),sizeof(int));gettimeofday(&end_time, NULL);
    printf("CALLOC MIS: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
    memcpy(M1, in_vertices, (num_in+1)*sizeof(int));

    free(in_vertices);
    free(worklist1);
    free(worklist2);
    free(rowStatus);
    free(colStatus);

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

void refreshColStatus(int* colStatus, int* rowStatus, int i, int nv, const struct CSR* G){ // 3.b. Mark local minima

    //// Get neighbors of vertex i ////

    // is i at the end of the vertex array or not

    int tot=0; 
    if(i!=(nv-1)){ tot=G->index_ptr[i+1];} 
    else{  tot=G->num_edges;} 

    //// //// //// ////

    int s = rowStatus[i];

    int neigh_stat;

    for(int j = G->index_ptr[i]; j < tot; j++){
        if( (G->col_ptr[j].coln  < nv) && (G->col_ptr[j].coln  != i) ){ // check against neighbors: do I have the lowest priority?
            neigh_stat = rowStatus[G->col_ptr[j].coln];

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

void decideInOut(int i, int* rowStatus, int* colStatus, struct CSR* G){ // 3.c. Look at vertex & neighbors statuses, mark new Ins/Outs
    
    int nv = G->N;

    int s = rowStatus[i];

    if( (s == IN) || (s==OUT) ){ return; }

    //// Get neighbors of vertex i ////

    int tot=0; 
    if(i!=(G->N-1)){ tot=G->index_ptr[i+1];} 
    else{  tot=G->num_edges;} 

    //// //// //// ////

    bool neigh_out = false;
    bool neigh_mismatchS = false;
    int neigh_stat;

    for(int j=G->index_ptr[i]; j < tot; j++){ // if me and my neighbors have the same priorty, I am the local minima and can be in MIS
        if( G->col_ptr[j].coln >= nv){ continue; } // else: mismatch

        neigh_stat = colStatus[G->col_ptr[j].coln];

        if (neigh_stat == OUT){ // if neighbor is out, i is out
            neigh_out = true;
            break;
        } else if(neigh_stat != s){
            neigh_mismatchS = true;
        }
    }
    

    if (neigh_out){
        pthread_mutex_lock(&m0); // lock
        rowStatus[i] = OUT;                 
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

/****************************************/

void firstPass(struct CSR* G, int supernode, int m, int* labels){ // 4. First Pass

    labels[m] = (int)supernode;

    //// Get neighbors of vertex i ////
    int tot=0; 
    if(m!=(G->N-1)){ tot=G->index_ptr[m+1];} 
    else{  tot=G->num_edges;} 
    //// //// //// ////

    for(int i=G->index_ptr[m]; i < tot; i++){ // Label all MIS vertices and their neighbors to the corresponding supernode
        if (labels[G->col_ptr[i].coln] == -1){
            labels[G->col_ptr[i].coln] = (int)supernode; 
        } 
    }

    return;
}

void secondPass(struct CSR* G, int unasgn, int* labels){ // 5. Second Pass - neighbors of MIS neighbors [2-MIS]
    
    if (labels[unasgn] != -1){ return; } // if vertex has label,skip
    else{
        //// Get neighbors of vertex i ////
        int tot=0; 
        if(unasgn!=(G->N-1)){ tot=G->index_ptr[unasgn+1];} 
        else{  tot=G->num_edges;}
        //// //// //// ////

        for(int i=G->index_ptr[unasgn]; i < tot; i++){ // mark undecided vertices based on if their neighbors have labels
            if (labels[G->col_ptr[i].coln] != -1){
                pthread_mutex_lock(&m2); // lock
                labels[unasgn] = labels[G->col_ptr[i].coln];
                pthread_mutex_unlock(&m2); // unlock
            }
        }
    }
    
    return;
}

void thirdPass(struct CSR* G, int unasgn, int* labels, int* supernode_ptr){ // 6. Thid Pass

    if (labels[unasgn] != -1){ return; }
    else{
        //// Get neighbors of vertex i ////
        int tot=0; 
        if(unasgn!=(G->N-1)){ tot=G->index_ptr[unasgn+1];} 
        else{  tot=G->num_edges;}
        //// //// //// ////

        pthread_mutex_lock(&m3); // lock
        (*supernode_ptr)++;
        int local_spn = (int)*supernode_ptr; // make a new supernode
        pthread_mutex_unlock(&m3); // unlock

        labels[unasgn] = local_spn; // assign neighbors of new supernode to group [overwrites potentially]
        for(int i=G->index_ptr[unasgn]; i < tot; i++){ 
            labels[G->col_ptr[i].coln] = local_spn; 
        } 
    }

    return;
}

/****************************************/
// potential parallel option:
// scan G and count edges in new subgraph
// have a storage vector of edges
// put edges in a vector in parallel?
// add_edges??

bool mark_edges(int r, int c, int N, bool** bool_pntr_array){ // !!!! I think idea would work but it is more just like a half matrix style marking 1s or 0s !!!!
    bool marked = false;
    // if they're not equal and the particular neighbor relation has already been marked
    // is there a way to build this graph without needing to know if an edge has been tracked

    // if not a bool type, could potentially store information of what index the edge will be in the col structure, 
    // but it's harder to parallelize

    if(c < r){           // 0,1 == 1?
        if (bool_pntr_array[r][c]){
            marked = true;
            return marked;
        }
    }else if(r<c){         // else 1,0 == 1?
        if (bool_pntr_array[c][r]){
            marked = true;
            return marked;
        }
    }

    marked = false;

    return marked;
}

struct CSR* build_edges(struct CSR* G, int* labels, int spn){ // 7. Build Edges

    //printf("Labels ");for(int k=0;k<G->N;k++){ printf("%d ", labels[k]);}printf("\n");
    //printf("\nspn: %d\n", spn);

    int row_spn=-1;
    int col_spn=-1;
    float val = 0.0;

    struct timeval start_time,end_time;
    gettimeofday(&start_time, NULL);int* numNeighbors = (int*)calloc((int)spn+1,sizeof(int));gettimeofday(&end_time, NULL);
    printf("CALLOC BE #Neigh: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
    gettimeofday(&start_time, NULL);bool** bool_pntr_array = (bool**)calloc(spn+1, sizeof(bool *));gettimeofday(&end_time, NULL);
    printf("CALLOC Boolptr: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));

    for(int p=1; p < spn+1; p++){
        //bool* w = calloc(p, sizeof(bool));
        bool w[p]; memset(w, 0, p*sizeof(bool));
        //bool_pntr_array[p] = (bool*)w;
        bool_pntr_array[p] = &w[0];
    }
    
    int ncnt = 0;
    int total_edges = 0;
    bool marked = false;
    for(int i=0; i < G->num_edges; i++){ // potentially parallelizable, but would need locks

        if(ncnt < (G->N-1)){
            if(i == G->index_ptr[ncnt+1]){ 
                ncnt++;
                
            }
        }else{}

        //    0 1 ....
        // n: 0 4
        // i: 0 1 2 | 3 4 .. n edges

        row_spn =  labels[ncnt];  // edge (row[i],col[i]) have same label, else an edge needs to be marked
        col_spn =  labels[G->col_ptr[i].coln];
        val = G->col_ptr[i].w;

        if(row_spn != col_spn){
            marked = false; // marked will probably need to be locked somehow !!
            marked = mark_edges(row_spn,col_spn,spn+1,bool_pntr_array);
            
            //printf("(%d,%d), marked:%d \t", row_spn, col_spn,marked);
            if(!marked){

                if(row_spn < col_spn ){ bool_pntr_array[col_spn][row_spn] = 1; }
                if(col_spn < row_spn){ bool_pntr_array[row_spn][col_spn] = 1; }

                numNeighbors[row_spn] += 1;
                numNeighbors[col_spn] += 1;
                total_edges+=2;
                marked=true;

            }
            //printf("rcnt:%d ccnt:%d\n",numNeighbors[row_spn],numNeighbors[col_spn] );

        }
            
    }//printf("\n");

    //cilk_sync;

    // index array
    int k = 0;

    int tns = total_edges;
    
    for(k=spn; k>=0; k--){
        tns = tns - numNeighbors[k];
        numNeighbors[k] = tns;
        //printf("%d: %d\t", k, numNeighbors[k]);
    }//printf("\ntotal #s: %d\n\n", total_edges);

    
    gettimeofday(&start_time, NULL);
    struct coltup* col = (struct coltup*)calloc(total_edges,sizeof(struct coltup));
    gettimeofday(&end_time, NULL);
    printf("CALLOC BE COL: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
    // 0 = [1,2,3]
    for(int k=0; k < total_edges; k++){
        col[k].coln =  -1; // 
    }

    /******************/
    int count = 0; int next_neigh = 0; // set up something where a thread can do row or col first
    for(int ptr=spn; ptr >= 0; ptr--){
        for(int i=0; i < ptr; i++){
            
            if(bool_pntr_array[ptr][i] == 1){
                
                count = 0;
                if (ptr < (spn)){ next_neigh = numNeighbors[ptr+1]; }
                else{ next_neigh = total_edges; }
                // 1) what is the index of ptr in col --> numNeigh
                while(col[numNeighbors[ptr]+count].coln != -1 && count < (next_neigh - numNeighbors[ptr]) && (next_neigh - numNeighbors[ptr])>0){
                    count++;
                }
                // 2) which index in col is it assigned to
                if((next_neigh - numNeighbors[ptr])>0){ col[numNeighbors[ptr]+count].coln = i; }


                count = 0;
                if (i < (spn)){ next_neigh = numNeighbors[i+1];}
                else{ next_neigh = total_edges; }
                while(col[numNeighbors[i]+count].coln != -1 && count < (next_neigh - numNeighbors[i])&& (next_neigh - numNeighbors[i])>0){ count++; }
                if((next_neigh - numNeighbors[i])>0){ col[numNeighbors[i]+count].coln = ptr;}

            }
        }
    }
    //cilk_sync;

    //for(int p=1;p<spn+1; p++){
    //    free((bool*)bool_pntr_array[p]);
    //}
    free(bool_pntr_array);

    /******************/ // ???? ????
    ncnt=0;
    for(int i=0; i < G->num_edges; i++){
        // so now, all of the edges should be recorded in col now, and just the weights need to be tracked
        if(ncnt < (G->N-1)){
            if(i == G->index_ptr[ncnt+1]){ 
                ncnt++;
                
            }
        }else{}

        row_spn =  labels[ncnt];     
        col_spn =  labels[G->col_ptr[i].coln];
        val = G->col_ptr[i].w;

        if(row_spn != col_spn){

            count = 0;
            if (row_spn < (spn)){ next_neigh = numNeighbors[row_spn+1]; }
            else{ next_neigh = total_edges; }
            while(col[numNeighbors[row_spn]+count].coln != col_spn && count < (next_neigh - numNeighbors[row_spn]-1) && (next_neigh - numNeighbors[row_spn])>0){
                count++;
            }
            
            if((next_neigh - numNeighbors[row_spn])>0){ 
                //printf("r:%d, c:%d, v:%f, colrspc:%f\t",row_spn,col_spn,val, col[numNeighbors[row_spn]+count].w);
                col[numNeighbors[row_spn]+count].w = col[numNeighbors[row_spn]+count].w+val; 
                //printf("cnt:%d  r:%d, c:%d,v:%f,colrspc:%f\n",count,row_spn,col_spn,val,col[numNeighbors[row_spn]+count].w);
            } 
            //printf("r:%d, c:%d,nnr:%d,#n:%d,cnt:%d,colrspc:%d\n",row_spn,col_spn,numNeighbors[row_spn],(next_neigh - numNeighbors[row_spn]),count,col[numNeighbors[row_spn]+count].coln);

            count=0;
            if (col_spn < (spn)){ next_neigh = numNeighbors[col_spn+1]; }
            else{ next_neigh = total_edges; }
            // col[numNeighbor[c]+offset].coln && col[numNeighbor[r]+offset].coln
            while(col[numNeighbors[col_spn]+count].coln != row_spn && count < (next_neigh - numNeighbors[col_spn]-1) &&(next_neigh - numNeighbors[col_spn])>0){
                count++;
            }
            if((next_neigh - numNeighbors[col_spn])>0){ 
                //printf("c:%d, r:%d, v:%f, colrspc:%f\t",col_spn,row_spn,val, col[numNeighbors[col_spn]+count].w);
                col[numNeighbors[col_spn]+count].w =col[numNeighbors[col_spn]+count].w + val; 
                //printf("cnt:%d  c:%d, r:%d,v:%f,colrspc:%f\n",count,col_spn,row_spn,val,col[numNeighbors[col_spn]+count].w);
            } 

        }

    }
    
    //cilk_sync;

    cleanup(G); // Do not need larger graph anymore, free the memory

    int tott=0; 
// */

    static struct CSR subgraph;
    subgraph.N = spn+1;
    subgraph.num_edges = total_edges;
    subgraph.index_ptr = numNeighbors;
    subgraph.col_ptr = col;

    /*int nct =0;
    for(int i=0; i < subgraph.num_edges; i++){
        printf("(%d, %d): %f\n", subgraph.index_ptr[nct], subgraph.col_ptr[i].coln, subgraph.col_ptr[i].w); 
        if(nct < (subgraph.N-1)){ if(i == subgraph.index_ptr[nct+1]){ nct++; }
        }else{}
    }*/

    return &subgraph;  
}

/****************************************/

struct CSR* kokkos_coarsen(struct CSR* G, int N){ // 2. Kokkos Coarsen
    struct timeval start_time,end_time;
    // list of vertices to assign to supernodes - unassigned vertices 
    int* vertex_set = NULL;
    gettimeofday(&start_time, NULL);vertex_set = (int*)calloc(N,sizeof(int));gettimeofday(&end_time, NULL);
    printf("CALLOC unassgin.v: %lf\n", (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000));
    for (int i = 0; i < N; i++){ vertex_set[i] = i; }

    // 3. Find the 2-Maximal Independent Set
    int* M1= NULL;
    
    gettimeofday(&start_time, NULL);
    M1 = kokkos_mis(G);
    gettimeofday(&end_time, NULL);
    double B1 = (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000);
    printf("1B1: %lf\n", B1);
    
    if(!M1){printf("error\n");}

    int M1_length=0;                            // get the length of the MIS = # of starting supernodes
    while(M1[M1_length] != -1){ M1_length+=1;}

    printf("M1 Length: %d\n", M1_length);

    // print out MIS
    // printf("M1: "); for(int i=0;i<M1_length;i++){ printf("%d ", M1[i]); } printf("\n");
    
    // Create array where indices = vertices & values = supernode vertices is grouped under
    gettimeofday(&start_time, NULL);
    //int labels[N];
    struct timeval st,et;
    gettimeofday(&st, NULL);int* labels = (int*)calloc(N,sizeof(int));gettimeofday(&et, NULL);
    printf("CALLOC labelArr: %lf\n", (et.tv_sec+ (double)et.tv_usec/1000000) - (st.tv_sec+(double)st.tv_usec/1000000));
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
    for(int i=0; i < N; i++){ // new size of vertex set
        if(labels[i] == -1){ 
            new_len++;
        } 
    }

    int nvs[new_len];           // why do I keep resizing vertex_set and just use nvs?
    int k=0;
    for(int i=0; i < N; i++){ 
        if(labels[i] == -1){ 
            nvs[k] = i; k++;
        } 
    }

    gettimeofday(&st, NULL);vertex_set = (int*)realloc(vertex_set, new_len*sizeof(int)); gettimeofday(&et, NULL); // resize set of unlabeled vertices
    printf("C-REALLOC unasgn.v1: %lf\n", (et.tv_sec+ (double)et.tv_usec/1000000) - (st.tv_sec+(double)st.tv_usec/1000000));
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

    gettimeofday(&st, NULL);vertex_set = (int*)realloc(vertex_set, new_len*sizeof(int));gettimeofday(&et, NULL);
    printf("C-REALLOC unasgn.v2: %lf\n", (et.tv_sec+ (double)et.tv_usec/1000000) - (st.tv_sec+(double)st.tv_usec/1000000));
    memcpy(vertex_set, nvs2, sizeof(nvs2));

    int supernode = M1_length; // number of elements in M1
    int* spn_ptr = &supernode;

    // 6. Third Pass
    if(new_len > 0){
        int* thirdPasspntr = vertex_set;
        cilk_for(int i=0; i < new_len; i++){ // create new supernodes or group still unlabeled vertices
            thirdPass(G, thirdPasspntr[i], labels, spn_ptr);
        }
    }
    gettimeofday(&end_time, NULL);
    double B2 = (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000);
    printf("2B2: %lf\n", B2);

    cilk_sync; // sync

    printf("SPN: %d\n", supernode);
    //printf("Labels: ");for(igraph_integer_t y=0; y < N; y++){
    //    printf("%lld ", labels[y]); }
    printf("\n");

    // 7. Build the subgraph edge relationships based on the original graph

    gettimeofday(&start_time, NULL);
    //start = clock();
    struct CSR* subgraph = build_edges(G, labels, supernode);
    //if(N <= 000000){ subgraph = build_edges(G, labels, supernode); }
    //if(N > 1000000){  subgraph = build_edges2(G, labels, supernode);}
    gettimeofday(&end_time, NULL);
    double B3 = (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000);
    //end = clock();
    //double B3 = (((double)end - start)/CLOCKS_PER_SEC);
    printf("3B3: %lf\n", B3);

    free(labels);

    if (new_len) { free(vertex_set); } // free vertex set memory

    return subgraph; // return graph pointer
}

struct CSR* display_ratio(struct CSR* G){
    printf("\n!: %d, %d, %d\n\n", G->N, G->num_edges, ITERATIONS);

    return G;
}

int cleanup(struct CSR* G){
	printf("NCLEANUP: %d\n", G->N);

    free(G->col_ptr); // col
    free(G->index_ptr); // indexptr
    // free(); // CSR pointer? 

    return 0;

}

int recursion_fcn(struct CSR* G, int goal_verts){ // 1. Base Case/Recursive Step

    // get graph size
    int graph_size = G->N;
    int edge_num = G->num_edges;
    printf("\nI: %d, %d, %d\n\n", ITERATIONS, graph_size, goal_verts);

    // Base Case 
    if (graph_size <= goal_verts || edge_num <= 0){
        cleanup(G);
        return 0;
    }

    // Recursion
    struct CSR* subgraph = kokkos_coarsen(G, graph_size); // pass in pointer to G AND G's size

    ITERATIONS+=1;

    subgraph = display_ratio(subgraph);

    return recursion_fcn(subgraph, goal_verts);
}

int main(int argc, char *argv[]){

    // random seed
    // srand(time(NULL));
    srand(0);

    if (pthread_mutex_init(&m0, NULL) != 0) { perror("Mutex initialization failed"); return 1;}
    if (pthread_mutex_init(&m2, NULL) != 0) { perror("Mutex initialization failed"); return 1;}
    if (pthread_mutex_init(&m3, NULL) != 0) { perror("Mutex initialization failed"); return 1;}
    if (pthread_mutex_init(&m4, NULL) != 0) { perror("Mutex initialization failed"); return 1;}

    int N = atoi(argv[1]);        // number of vertices in data graph

    // load G in from a ncol file
    int MAXNAME = 1024; // max length of a file name
    char filename[MAXNAME];
    memset( filename, '\0', MAXNAME*sizeof(char) );
    strcpy(filename, argv[2]); // ncol file to read from
    // "csv/ash85.csv",   // #N=85
    // "csv/netz4504.csv" // #N=1961
    // "csv/gemat11.csv"  // #N=4929
    // "csv/bcspwr10.csv" // #N=5300
    // "csv/bcsstk17.csv" // #N=10974
    // "csv/shock-9.csv"  // #N=36476
    // "csv/ecology1.csv" // #N=1000000

    bool sym = false;

    //struct CSR g = csv_to_hash(filename, N, sym);
    //struct CSR* graph_ptr = &g;
    struct CSR* graph_ptr = csv_to_hash(filename, N, sym);

    graph_ptr = display_ratio(graph_ptr);
    int ec = graph_ptr->num_edges;

    // To Write: display graph function

    int goal_verts = (int)(N/6) - 1;//  (int)N-1; //set minimum goal vertices
    printf("goal_verts: %d\n", goal_verts);

    struct timeval start_time,end_time;
    gettimeofday(&start_time, NULL);
    int exit1 = recursion_fcn(graph_ptr, goal_verts); // 1. Begin recursive coarsening
    gettimeofday(&end_time, NULL);

    double realTime = (end_time.tv_sec+ (double)end_time.tv_usec/1000000) - (start_time.tv_sec+(double)start_time.tv_usec/1000000);
    //printf("# of Iterations: %d\n", ITERATIONS);
    //printf("Exit %d \n", exit1);
    printf("RT: %lf", realTime);
    //printf("Coarsened in %f seconds.\n", myTime); //

    printf("\n!: %d,", N);
    printf("%d,", ec);
    printf("%d\n", 0);

    pthread_mutex_destroy(&m0);
    pthread_mutex_destroy(&m2);
    pthread_mutex_destroy(&m3);
    pthread_mutex_destroy(&m4);

    return 0;
} 