//
// auto-generated by op2.py
//

//user function
int opDat0_dirichlet_stride_OP2CONSTANT;
int opDat0_dirichlet_stride_OP2HOST=-1;
//user function
//#pragma acc routine
inline void dirichlet( double *res) { *res = 0.0; }

// host stub function
void op_par_loop_dirichlet(char const *name, op_set set,
  op_arg arg0){

  int nargs = 1;
  op_arg args[1];

  args[0] = arg0;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(1);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[1].name      = name;
  OP_kernels[1].count    += 1;

  int  ninds   = 1;
  int  inds[1] = {0};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: dirichlet\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_1
    int part_size = OP_PART_SIZE_1;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set->size >0) {

    if ((OP_kernels[1].count==1) || (opDat0_dirichlet_stride_OP2HOST != getSetSizeFromOpArg(&arg0))) {
      opDat0_dirichlet_stride_OP2HOST = getSetSizeFromOpArg(&arg0);
      opDat0_dirichlet_stride_OP2CONSTANT = opDat0_dirichlet_stride_OP2HOST;
    }

    //Set up typed device pointers for OpenACC
    int *map0 = arg0.map_data_d;

    double *data0 = (double *)arg0.data_d;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;
    int set_size1 = set->size + set->exec_size;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      #pragma acc parallel loop independent deviceptr(col_reord,map0,data0)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map0idx = map0[n + set_size1 * 0];


        dirichlet(
          &data0[1 * map0idx]);
      }

    }
    OP_kernels[1].transfer  += Plan->transfer;
    OP_kernels[1].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[1].time     += wall_t2 - wall_t1;
}
