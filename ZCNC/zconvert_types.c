#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>

struct coltup{
    int coln;
    float w;    // weight
};
struct CSR{
    int num_edges; // number of unique edges
    int N; // # verts
    int* index_ptr;
    struct coltup* col_ptr;
    //float* value_ptr; // float
};

// 1,2,w
// 3,4,w

// 0,1,w
// 1,2,w
// 1,4,w
// 0,2,w
// 0,3,w
// 1,5,w

// CSV -> hash table? -> CSR
// 


// if sparse == # edges << max. edges
// CSV -> intermediate processing -> CSR
// s: 0 3
// t: 1 2 3 2 4 5 
// w: w w w w w w


void quickSort(int [], int, int, int [], float []);

struct CSR* csv_to_hash(char input_file[], int N, bool sym){

    static struct CSR graph;
    graph.N = N;

    int MAXSZ = 1024;       // max_num = max size of csv row
    //int max_matsz = N*N;

    // define hash table

    // how to setup hash table:
    // array 0 to N
    // at each index, a linked list for the neighborhood
    // size of a pntr vs. size of an int
    // will the hash table of the 

    // Read CSV and get # neighbors per N
    // int num_Neighbors[N];
    int* num_Neighbors = (int*)calloc(N,sizeof(int));
    // memset(num_Neighbors, 0, N*sizeof(int));

    int arr = 0;

    char csv_row[MAXSZ];
    char* token;

    int total_numbers = 0;

    FILE* stream = fopen(input_file, "r");
    fgets(csv_row, MAXSZ, stream); // get headers

    int csvarr0=0; int csvarr1=0;

    while (fgets(csv_row, MAXSZ, stream))
    {
        token = strtok(csv_row, ",");
        arr = 0;
        
        while(token != NULL) { // node1, node2, weight
            //printf("Token: %s\t", token);

            if (arr == 0){
                csvarr0 = atoi(token);
                num_Neighbors[atoi(token)]+=1;
                total_numbers++;
            } else if(arr == 1){
                csvarr1 = atoi(token);
                if(!sym && (csvarr0!=csvarr1)){
                    num_Neighbors[atoi(token)]+=1;
                    total_numbers++;}
            } else if (arr == 2){                
            }

            token = strtok(NULL, ",");
            arr++;
        }

    }

	fclose(stream);

    graph.num_edges = total_numbers;
    graph.index_ptr = num_Neighbors; // this may go poorly

    /************************/

    // index array
    int tns = total_numbers;
    for(int k=N-1; k>=0; k--){
        tns = tns - num_Neighbors[k];
        num_Neighbors[k] = tns;
        //if (k < 11 || k > (N-12)) {printf("%d: %d\t", k, num_Neighbors[k]);}
    }printf("\ntotal #s: %d\n", total_numbers);

    // struct 

    // reset csvrow?
    // char csv_row[MAXSZ];
    // return NULL;

    // init Table
    // int col[total_numbers];
    // int* col = (int*)calloc(total_numbers,sizeof(int));
    struct coltup* col = (struct coltup*)calloc(total_numbers,sizeof(struct coltup));
    // 0 = [1,2,3]
    for(int k=0; k < total_numbers; k++){
        col[k].coln =  -1; // (int*)calloc(num_Neighbors[k],sizeof(int));
    }

    int csvarr[3];
    stream = fopen(input_file, "r");
    fgets(csv_row, MAXSZ, stream); // get headers
    while (fgets(csv_row, MAXSZ, stream)) {

        token = strtok(csv_row, ",");
        arr = 0;
        
        while(token != NULL) { // node1, node2, weight
            //printf("Token: %s\n", token);

            if (arr == 0){
                // add to table
                csvarr[0] = atoi(token);

            } else if(arr == 1){
                // add to table
                // (0,1,w)
                // table[1] = [0, -1, -1];
                
                csvarr[1] = atoi(token);

                // printf("(%d, %d)\n", csvarr[0], csvarr[1]);

            } else if (arr == 2){  
                csvarr[2] = atoi(token);   

                if (col[num_Neighbors[csvarr[1]]].coln == -1){
                    col[num_Neighbors[csvarr[1]]].coln  = csvarr[0];
                    col[num_Neighbors[csvarr[1]]].w  = csvarr[2];
                } 

                else if(col[num_Neighbors[csvarr[1]]].coln != -1){

                    int cnt=0; int tot=0; bool flip = true;
                    if(csvarr[1]!=(N-1)){cnt = 0; tot=num_Neighbors[csvarr[1]+1] - num_Neighbors[csvarr[1]];} 
                    else{ cnt = 0; tot=total_numbers-num_Neighbors[csvarr[1]];} 

                    while( ((int)cnt < (int)tot) && (flip) ){ // i!=-1 and index
                        cnt++;
                        //if(atoi(token) == 1){printf("tkn:%s, col:%d, cnt:%d, tot:%d, nn:%d, csvarr[0]:%d, csvarr[1]:%d\n",token, col[num_Neighbors[atoi(token)]+cnt], cnt,tot, num_Neighbors[atoi(token)],csvarr[0],csvarr[1]);}
                        //printf("tkn:%s, col:%d, cnt:%d, tot:%d, nn:%d\n",token, col[num_Neighbors[atoi(token)]+cnt], cnt,tot, num_Neighbors[atoi(token)]);
                        if (col[num_Neighbors[csvarr[1]]+cnt].coln == -1){
                            col[num_Neighbors[csvarr[1]]+cnt].coln  = csvarr[0]; // mark the array
                            col[num_Neighbors[csvarr[1]]+cnt].w  = csvarr[2];
                            // printf("col:%d, nn:%d, csvarr[1]:%d, cnt:%d, tot:%d, csvarr[0]:%d \n", col[num_Neighbors[atoi(token)]+cnt], num_Neighbors[atoi(token)],csvarr[1], cnt,tot,csvarr[0]);
                            flip = false;
                        } 

                    }
                }

                if(!sym && (csvarr[0]!=csvarr[1]) ){
                    if(col[num_Neighbors[csvarr[0]]].coln == -1){
                        col[num_Neighbors[csvarr[0]]].coln = csvarr[1];
                        col[num_Neighbors[csvarr[0]]].w = csvarr[2];
                    }
                    else if(col[num_Neighbors[csvarr[0]]].coln != -1){
                        int cnt1=0; int tot1=0; bool flip1=true; 
                        if(csvarr[0]!=(N-1)){cnt1=0; tot1=num_Neighbors[csvarr[0]+1] - num_Neighbors[csvarr[0]];}  else{ cnt1=0; tot1=total_numbers-num_Neighbors[csvarr[0]];} 
                        while( ((int)cnt1 < (int)tot1) && ( flip1) ){ // i!=-1 and index
                            cnt1++; 
                            if (col[num_Neighbors[csvarr[0]]+cnt1].coln == -1){ 
                                col[num_Neighbors[csvarr[0]]+cnt1].coln  = csvarr[1];
                                col[num_Neighbors[csvarr[0]]+cnt1].w  = csvarr[2]; 
                                flip1=false;
                            } 
                        }
                    }
                }          
            }

            token = strtok(NULL, ",");
            arr++;
        }

    }

	fclose(stream);

    // col ptr
    int tott=0;
    /*printf("zcv Neighbor Lists: \n");
    for(int k=0; k < 2; k++){
      if (k < 11 || k > (N-12)) {
        printf("%d: [ ",k);

        if(k < (N-1)){ tott=num_Neighbors[k+1] - num_Neighbors[k];} 
        else{ tott = total_numbers-num_Neighbors[k];} 
        
        for(int l=num_Neighbors[k]; l < (num_Neighbors[k]+tott); l++){
            printf(" (%d, %f) ", col[l].coln, col[l].w);
        }
        printf("]\n"); }
    }*/

    graph.col_ptr = col;

    // for now - delete later
    //free(col);
    //free(num_Neighbors);

    return &graph; // &graph;
}

/********************************************************************/
// implement Quick Sort Algorithm
// Built upon the quickort algorithm at: https://www.geeksforgeeks.org/quick-sort-in-c/

void swap1(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}
void swap2(float* a, float* b) {
    float temp = *a;
    *a = *b;
    *b = temp;
}
int partition(int arr[], int low, int high, int link1[], float link2[]) {

    // Initialize pivot to be the first element
    int p = arr[low];
    int i = low;
    int j = high;

    while (i < j) {

        // Find the first element greater than
        // the pivot (from starting)
        while (arr[i] <= p && i <= high - 1) {
            i++;
        }

        // Find the first element smaller than
        // the pivot (from last)
        while (arr[j] > p && j >= low + 1) {
            j--;
        }
        if (i < j) {
            swap1(&arr[i], &arr[j]);
            swap1(&link1[i], &link1[j]);
            swap2(&link2[i], &link2[j]);
        }
    }
    swap1(&arr[low], &arr[j]);
    swap1(&link1[low], &link1[j]);
    swap2(&link2[low], &link2[j]);
    
    return j;
}
void quickSort(int arr[], int low, int high, int link1[], float link2[]) {
    if (low < high) {

        // call partition function to find Partition Index
        int pi = partition(arr, low, high, link1, link2);

        // Recursively call quickSort() for left and right
        // half based on Partition Index
        quickSort(arr, low, pi - 1, link1, link2);
        quickSort(arr, pi + 1, high, link1, link2);
    }
}

/********************************************************************/


/********************************************************************/
/*int main(){

    char file[1024] = "csv/ecology1.csv";
    int N =  1000000; //1971282; //
    bool sym = true;

    struct CSR* grptr = csv_to_hash(file, N, sym);

    //int n1 =  neigh_sz2(grptr, 0);
    //int n2 =  neigh_sz2(grptr, 5487);
    //int n3 =  neigh_sz2(grptr, 10973);
    //printf("\n");
    //for(int j=grptr->row_inds[10973]; j<(n3+grptr->row_inds[10973]-1);j++){
    //    printf("%d ", grptr->col_ptr[j]);
    //}printf("\n");

    return 0;
} */


/* int neigh_sz2(struct Graph* G, int v){
    int neighbor_size = 0;

    // 0 3 
    // n0 n0 n0 n1 n1 ....

    if (v >= (G->N-1)){
        neighbor_size = G->active_index - G->row_inds[v]; //
        
    } else {
        neighbor_size = G->row_inds[v+1] - G->row_inds[v];
    }

    int neighbors[neighbor_size];
    memset(neighbors, 0, neighbor_size*sizeof(int));

    for(int j=G->row_inds[v]; j<(neighbor_size+G->row_inds[v]);j++){
        //printf("%d ", G->col_ptr[j]);
        neighbors[j-G->row_inds[v]] = G->col_ptr[j];
    }

    printf("Neighbors of %d: ", v);
    for(int j = 0; j < neighbor_size; j++){
        printf("%d ", neighbors[j]);
    }printf("\n");

    return neighbor_size;
}*/
