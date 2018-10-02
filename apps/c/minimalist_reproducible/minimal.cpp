
//
// test program for new OPlus2 development
//

//
// standard headers
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// global constants

double alpha;

// jac header file

//
// OP header file
//

#include "op_seq.h"

//
// kernel routines for parallel loops
//
#include "increment_log.h"

// Error tolerance in checking correctness

#define TOLERANCE 1e-12

// define problem size

#define NN 6
#define NITER 2




int main(int argc, char **argv) {
  // OP initialisation
  op_init(argc, argv, 5);

  // timer
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  const char* file = "data.h5";
  // declare sets, pointers, and datasets

  op_set nodes = op_decl_set_hdf5(file, "nodes");
  op_set edges = op_decl_set_hdf5(file, "edges");

  op_map ppedge = op_decl_map_hdf5(edges, nodes, 2, file, "ppedge");

  op_dat p_A = op_decl_dat_hdf5(edges, 1, "double", file, "p_A");
  op_dat p_r = op_decl_dat_hdf5(nodes, 1, "double", file, "p_r");

  alpha = 1.0f;
  op_decl_const(1, "double", &alpha);

  op_diagnostic_output();

  // initialise timers for total execution wall time
  op_timers(&cpu_t1, &wall_t1);

  // main iteration loop

  double u_sum, u_max, beta = 1.0f;
  
  op_partition("PARMETIS", "KWAY", edges, ppedge, p_r);
  create_reversed_mapping();
  
  op_par_loop(increment_log,"increment_log",edges,
        op_arg_dat(p_A,-1,OP_ID,1, "double", OP_READ),
        op_arg_dat(p_r,0,ppedge,1,"double",OP_RW),
        op_arg_dat(p_r,1,ppedge,1,"double",OP_RW));
  
  
  op_fetch_data_hdf5_file(p_r,"repr_kimenet.h5");
  
  op_exit();
  
  
  
}
