//
// auto-generated by op2.py
//

//user function
inline void res_calc2(double *res1, double *res2, 
                      double *node1,
                      double *node2) {
  for(int i = 0; i < 2; i++){
      // r1[i] =  res1[i];
      // r2[i] =  res2[i];
    
      node1[i] += res1[0] + 1;
      node2[i] += res2[0] - 1;
  }
}#ifdef VECTORIZE
//user function -- modified for vectorisation
#if defined __clang__ || defined __GNUC__
__attribute__((always_inline))
#endif
inline void res_calc2_vec( double res1[][SIMD_VEC], double res2[][SIMD_VEC], double node1[][SIMD_VEC], double node2[][SIMD_VEC], int idx ) {
  for(int i = 0; i < 2; i++){



      node1[i][idx] = res1[0][idx] + 1;
      node2[i][idx] = res2[0][idx] - 1;
  }

}
#endif

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
  //create aligned pointers for dats
  ALIGNED_double const double * __restrict__ ptr0 = (double *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr1 = (double *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr2 = (double *) arg2.data;
  DECLARE_PTR_ALIGNED(ptr2,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr3 = (double *) arg3.data;
  DECLARE_PTR_ALIGNED(ptr3,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(1);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: res_calc2\n");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      if ((n+SIMD_VEC >= set->core_size) && (n+SIMD_VEC-set->core_size < SIMD_VEC)) {
        op_mpi_wait_all(nargs, args);
      }
      ALIGNED_double double dat0[1][SIMD_VEC];
      ALIGNED_double double dat1[1][SIMD_VEC];
      ALIGNED_double double dat2[2][SIMD_VEC];
      ALIGNED_double double dat3[2][SIMD_VEC];
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx0_1 = 1 * arg0.map_data[(n+i) * arg0.map->dim + 0];
        int idx1_1 = 1 * arg0.map_data[(n+i) * arg0.map->dim + 1];

        dat0[0][i] = (ptr0)[idx0_1 + 0];

        dat1[0][i] = (ptr1)[idx1_1 + 0];

        dat2[0][i] = 0.0;
        dat2[1][i] = 0.0;

        dat3[0][i] = 0.0;
        dat3[1][i] = 0.0;

      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        res_calc2_vec(
          dat0,
          dat1,
          dat2,
          dat3,
          i);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx2_2 = 2 * arg2.map_data[(n+i) * arg2.map->dim + 0];
        int idx3_2 = 2 * arg2.map_data[(n+i) * arg2.map->dim + 1];

        (ptr2)[idx2_2 + 0] += dat2[0][i];
        (ptr2)[idx2_2 + 1] += dat2[1][i];

        (ptr3)[idx3_2 + 0] += dat3[0][i];
        (ptr3)[idx3_2 + 1] += dat3[1][i];

      }
    }

    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
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
        &(ptr0)[1 * map0idx],
        &(ptr1)[1 * map1idx],
        &(ptr2)[2 * map2idx],
        &(ptr3)[2 * map3idx]);
    }
  }

  if (exec_size == 0 || exec_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

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
