//
// auto-generated by op2.py
//

//user function
#include "../res_calc2.h"
#include "op_lib_mpi.h"

// host stub function
void op_par_loop_res_calc2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(1);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: res_calc2\n");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

#ifdef COMM_AVOID
  int exec_size = 0;
  for(int l = 0; l < 1; l++){
    exec_size += OP_aug_import_exec_lists[l][set->index]->size;
  }
  set_size = set->size + exec_size;
#endif

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map0idx;
      int map1idx;
      int map2idx;
      int map3idx;
      map0idx = arg0.map_data[n * arg0.map->dim + 0];
      map1idx = arg0.map_data[n * arg0.map->dim + 1];
      map2idx = arg2.map_data[n * arg2.map->dim + 0];
      map3idx = arg2.map_data[n * arg2.map->dim + 1];


      res_calc2(
        &((double*)arg0.data)[1 * map0idx],
        &((double*)arg0.data)[1 * map1idx],
        &((double*)arg2.data)[2 * map2idx],
        &((double*)arg2.data)[2 * map3idx]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }

#ifndef COMM_AVOID
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);
#endif
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[1].name      = name;
  OP_kernels[1].count    += 1;
  OP_kernels[1].time     += wall_t2 - wall_t1;
  OP_kernels[1].transfer += (float)set->size * arg0.size;
  OP_kernels[1].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[1].transfer += (float)set->size * arg0.map->dim * 4.0f;
  OP_kernels[1].transfer += (float)set->size * arg2.map->dim * 4.0f;
}
