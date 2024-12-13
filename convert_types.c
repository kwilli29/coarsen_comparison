#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>

struct Graph{
    int active_index; // number of unique edges
    int N; // subgraph
    int* row_ptr;
    int* col_ptr;
    float* value_ptr; // change to float!!!!

    // keep a vertex list? 
    // submeta if ever relevant
};

void quickSort(int [], int, int, int [], float []);

struct Graph* csv_to_graph_small(char input_file[], int N, bool sym){
    
    int MAXSZ = 1024;  // max_num = max size of a csv row
    int max_matsz = N*N;
    
    int max_row[max_matsz]; // aboslute max of a graph size
    int max_col[max_matsz];
    float max_value[max_matsz];

    int arr = 0;
    int index = 0;

    FILE* stream = fopen(input_file, "r");

    char csv_row[MAXSZ];
    char* token;

    fgets(csv_row, MAXSZ, stream); // get headers

    while (fgets(csv_row, MAXSZ, stream))
    {
        token = strtok(csv_row, ",");
        arr = 0;
        
        while(token != NULL) { // node1, node2, weight
            //printf("Token: %s\n", token);

            if (arr == 0){
                max_row[index] = atoi(token); // atoi bad?
            } else if(arr == 1){
                max_col[index] = atoi(token);
            } else if (arr == 2){
                max_value[index] = atof(token);
                for(int k=0; k < strlen(token); k++){
                    max_value[index] = (float)(token[k] - '0');
                    break;
                }
                
            }

            if ( (arr >= 2) && (sym == false) ){ // if data is not symmetric, make csr symmetric
                index++;
                max_row[index] = max_col[index-1];
                max_col[index] = max_row[index-1];
                max_value[index] = max_value[index-1];
            }

            token = strtok(NULL, ",");
            arr++;
        }

        index++;

        if(max_row[index-1] == max_col[index-1]){ 
            index--; 
        } // don't count self loops for now

    }

    //index--;

    static struct Graph ret_graph; // return graph
    ret_graph.N = N;
    ret_graph.active_index = index;
    ret_graph.row_ptr = (int*)calloc(index,sizeof(int));
    //ret_graph.row_ptr = (int*)calloc(N,sizeof(int));
    ret_graph.col_ptr = (int*)calloc(index,sizeof(int));
    ret_graph.value_ptr = (float*)calloc(index,sizeof(float)); // floats


    for(int i=0; i<index;i++){
        ret_graph.row_ptr[i] = max_row[i]; 
        ret_graph.col_ptr[i] = max_col[i]; 
        ret_graph.value_ptr[i] = max_value[i];
    }

    quickSort(ret_graph.row_ptr, 0, index-1, ret_graph.col_ptr, ret_graph.value_ptr);

	fclose(stream);

    return (&ret_graph);
}

struct Graph* csv_to_graph_large(char input_file[], int N, bool sym){
    
    int MAXSZ = 1024;  // max_num = max size of csv row
    int max_matsz = N*N;

    int* max_row = (int*)calloc(max_matsz,sizeof(int));
    int* max_col = (int*)calloc(max_matsz,sizeof(int));
    float* max_value = (float*)calloc(max_matsz,sizeof(float)); // floats

    int arr = 0;
    int index = 0;

    FILE* stream = fopen(input_file, "r");

    char csv_row[MAXSZ];
    char* token;

    fgets(csv_row, MAXSZ, stream); // get headers

    while (fgets(csv_row, MAXSZ, stream))
    {
        token = strtok(csv_row, ",");
        arr = 0;
        
        while(token != NULL) { // node1, node2, weight
            //printf("Token: %s\n", token);

            if (arr == 0){
                max_row[index] = atoi(token); // atoi bad?
            } else if(arr == 1){
                max_col[index] = atoi(token);
            } else if (arr == 2){
                max_value[index] = atof(token);
                /*for(int k=0; k < strlen(token); k++){
                    max_value[index] = (float)(token[k] - '0');
                    printf("token %lf\n", max_value[index]);
                    break;
                }*/
                
            }

            if ( (arr >= 2) && (sym == false) ){ // if data is not symmetric, make csr symmetric
                index++;
                max_row[index] = max_col[index-1];
                max_col[index] = max_row[index-1];
                max_value[index] = max_value[index-1];
            }

            token = strtok(NULL, ",");
            arr++;
        }

        index++;

        if(max_row[index-1] == max_col[index-1]){ 
            index--; 
        } // don't count self loops for now

    }

    //index--;

    static struct Graph ret_graph; // return graph
    ret_graph.N = N;
    ret_graph.active_index = index;
    ret_graph.row_ptr  = (int*)calloc(index,sizeof(int));
    //ret_graph.row_ptr = (int*)calloc(N,sizeof(int));
    ret_graph.col_ptr = (int*)calloc(index,sizeof(int));
    ret_graph.value_ptr = (float*)calloc(index,sizeof(float)); // floats


    for(int i=0; i<index;i++){
        ret_graph.row_ptr[i] = max_row[i]; ret_graph.col_ptr[i] = max_col[i]; ret_graph.value_ptr[i] = max_value[i];
    }

    quickSort(ret_graph.row_ptr, 0, index-1, ret_graph.col_ptr, ret_graph.value_ptr);

    /*int vert_track = row_verts[0];  // 0 0 0 1 1 1 = 0, 3, 6 ...
    int in = 0;
    ret_graph.row_ptr[0] = in;
    in++;
    for(int i = 0; i < index; i++){
        if (vert_track != row_verts[i]){
            ret_graph.row_ptr[in] = i;
            in++;
            vert_track = row_verts[i];
        }
    }
    */

    free(max_row);
    free(max_col);
    free(max_value);

	 fclose(stream);

    return (&ret_graph);
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

int binary_search(struct Graph* G, int NUM_EDGE, int v){ // https://stackoverflow.com/questions/13197552/using-binary-search-with-sorted-array-with-duplicates
 //   int start = 0;
  //  int end = NUM_EDGE-1;
    int ret_index = -1;

	for(int k=0; k< G->active_index;k++){ // linear search
		if (G->row_ptr[k] == v){
			ret_index = k;
            break;
		}
	}

    return ret_index; 

	// ACTING UP //
	/*
    while(start <= end){
        int middle = (int)((end-start)/2) + start;
        if(G->row_ptr[middle] > v){         // curr vert too big
            end = middle - 1;
        } else if (G->row_ptr[middle] == v){ // curr vert equal, must find earliest occurence
            ret_index = middle;
            break;
        } else if (G->row_ptr[middle] < v){ // curr vert too small
            start = middle + 1;
        }
    }

   
	 bool loopw=true;
    while(loopw){
		if ((ret_index > 0) && (G->row_ptr[ret_index]) == v) {
			loopw = false;
			break;
		}else if(ret_index >0){
        ret_index--;
		}
		else{
			ret_index = (G->N)-1;
		}
    	
	 }
    ret_index++;

    return ret_index; */
}


/********************************************************************/

/* 
int main(){

    char file[1024] = "csv/west0067.csv";
    int N = 67;
    struct Graph* grptr = csv_to_graph(file, N);

    return 0;
}
*/