/*
    Test the result stored in "data_output" against a serial implementation.

    -----
    Compiling:
    Include "Lab4_IO.c" to compile. Set the macro "LAB4_EXTEND" defined in the "Lab4_IO.c" file to include the extended functions
    $ gcc serialtester.c Lab4_IO.c -o serialtester -lm 

    -----
    Return values:
    0      result is correct
    1      result is wrong
    2      problem size does not match
    253    no "data_output" file
    254    no "data_input" file
*/
#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab4_IO.h"
#include "timer.h"
#include <mpi.h>

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

#define THRESHOLD 0.0001

int main (int argc, char* argv[]){
    struct node *nodehead;
    int nodecount;
    int *num_in_links, *num_out_links;
    double *r, *r_pre;
    int i, j;
    double damp_const;
    int iterationcount = 0;

    double cst_addapted_threshold;
    double error;

    double time_start, time_end;

 	

    // Parallelization Initialization
    int npes, myrank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    get_node_stat(&nodecount, &num_in_links, &num_out_links);
    node_init(&nodehead, num_in_links, num_out_links, 0, nodecount);

    double* piece_buf = malloc(nodecount/npes * sizeof(double));

    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));

    // i = myrank; i < ((myrank * nodecount/npes) + nodecount/npes)
    for (i=0; i < nodecount; i++ ){
        r[i] = 1.0 / nodecount;
    }
    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;
    // CORE CALCULATION
    int proper_index;
	GET_TIME(time_start);
    do{
        ++iterationcount;
        vec_cp(r, r_pre, nodecount);
        for (i = myrank * nodecount/npes; i < (myrank*(nodecount/npes) + nodecount/npes); ++i){
            proper_index = i - (myrank * nodecount/npes);
            piece_buf[proper_index] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                piece_buf[proper_index] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            piece_buf[proper_index] *= DAMPING_FACTOR;
            piece_buf[proper_index] += damp_const;

        }
        MPI_Allgather(piece_buf, nodecount/npes, MPI_DOUBLE, r, nodecount/npes, MPI_DOUBLE, MPI_COMM_WORLD);
//        if(myrank == 0){
//            error = rel_error(r, r_pre, nodecount) ;
//        }
//        MPI_Bcast(&error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }while(rel_error(r, r_pre, nodecount) >= EPSILON);
    //printf("Program converges at %d th iteration.\n", iterationcount);

   




    // post processing
   // if(myrank == 0){
        GET_TIME(time_end);
        Lab4_saveoutput(r, nodecount, (time_end - time_start));
   // }
    free(r);
    //free(piece_buf);
    free(r_pre);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    node_destroy(nodehead, nodecount);
    free(num_in_links); free(num_out_links);
}
