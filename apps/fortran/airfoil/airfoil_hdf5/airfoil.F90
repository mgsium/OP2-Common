program airfoil
    use op2_fortran_hdf5_declarations
    use op2_fortran_reference
    use op2_fortran_rt_support

    use airfoil_constants
    use airfoil_kernels

    use, intrinsic :: iso_c_binding

    implicit none

    character(*), parameter :: file_name = "new_grid.h5"

    integer(4), parameter :: niter = 1000
    integer(4) :: iter, i, k, stat

    integer(4) :: nnode, ncell, nbedge, nedge
    real(8), dimension(2) :: rms

    type(op_set) :: nodes, edges, bedges, cells
    type(op_map) :: pedge, pecell, pcell, pbedge, pbecell
    type(op_dat) :: p_bound, p_x, p_q, p_qold, p_adt, p_res

    real(kind=c_double) :: start_time, end_time

    real(8) :: diff

    call op_init_base(0, 0)

    print *, "Declaring OP2 sets"
    call op_decl_set_hdf5(nnode, nodes, file_name, "nodes")
    call op_decl_set_hdf5(nedge, edges, file_name, "edges")
    call op_decl_set_hdf5(nbedge, bedges, file_name, "bedges")
    call op_decl_set_hdf5(ncell, cells, file_name, "cells")

    print *, "Declaring OP2 maps"
    call op_decl_map_hdf5(edges, nodes, 2, pedge, file_name, "pedge", stat)
    call op_decl_map_hdf5(edges, cells, 2, pecell, file_name, "pecell", stat)
    call op_decl_map_hdf5(bedges, nodes, 2, pbedge, file_name, "pbedge", stat)
    call op_decl_map_hdf5(bedges, cells, 1, pbecell, file_name, "pbecell", stat)
    call op_decl_map_hdf5(cells, nodes, 4, pcell, file_name, "pcell", stat)

    print *, "Declaring OP2 data"
    call op_decl_dat_hdf5(bedges, 1, p_bound, "integer(4)", file_name, "p_bound", stat)
    call op_decl_dat_hdf5(nodes, 2, p_x, "real(8)", file_name, "p_x", stat)
    call op_decl_dat_hdf5(cells, 4, p_q, "real(8)", file_name, "p_q", stat)
    call op_decl_dat_hdf5(cells, 4, p_qold, "real(8)", file_name, "p_qold", stat)
    call op_decl_dat_hdf5(cells, 1, p_adt, "real(8)", file_name, "p_adt", stat)
    call op_decl_dat_hdf5(cells, 4, p_res, "real(8)", file_name, "p_res", stat)

    print *, "Declaring OP2 constants"
    call op_decl_const(gam, 1, "real(8)")
    call op_decl_const(gm1, 1, "real(8)")
    call op_decl_const(cfl, 1, "real(8)")
    call op_decl_const(eps, 1, "real(8)")
    call op_decl_const(mach, 1, "real(8)")
    call op_decl_const(alpha, 1, "real(8)")
    call op_decl_const(qinf, 4, "real(8)")

    call op_partition("PTSCOTCH", "KWAY", edges, pecell, p_x)
    call op_timers(start_time)

    do iter = 1, niter
        call op_par_loop_2(save_soln, cells, &
            op_arg_dat(p_q,    -1, OP_ID, 4, "real(8)", OP_READ), &
            op_arg_dat(p_qold, -1, OP_ID, 4, "real(8)", OP_WRITE))

        do k = 1, 2
            call op_par_loop_6(adt_calc, cells, &
                op_arg_dat(p_x,    1, pcell, 2, "real(8)", OP_READ), &
                op_arg_dat(p_x,    2, pcell, 2, "real(8)", OP_READ), &
                op_arg_dat(p_x,    3, pcell, 2, "real(8)", OP_READ), &
                op_arg_dat(p_x,    4, pcell, 2, "real(8)", OP_READ), &
                op_arg_dat(p_q,   -1, OP_ID, 4, "real(8)", OP_READ), &
                op_arg_dat(p_adt, -1, OP_ID, 1, "real(8)", OP_WRITE))

            call op_par_loop_8(res_calc, edges, &
                op_arg_dat(p_x,    1, pedge,  2, "real(8)", OP_READ), &
                op_arg_dat(p_x,    2, pedge,  2, "real(8)", OP_READ), &
                op_arg_dat(p_q,    1, pecell, 4, "real(8)", OP_READ), &
                op_arg_dat(p_q,    2, pecell, 4, "real(8)", OP_READ), &
                op_arg_dat(p_adt,  1, pecell, 1, "real(8)", OP_READ), &
                op_arg_dat(p_adt,  2, pecell, 1, "real(8)", OP_READ), &
                op_arg_dat(p_res,  1, pecell, 4, "real(8)", OP_INC),  &
                op_arg_dat(p_res,  2, pecell, 4, "real(8)", OP_INC))

            call op_par_loop_6(bres_calc, bedges, &
                op_arg_dat(p_x,      1, pbedge,  2, "real(8)",    OP_READ), &
                op_arg_dat(p_x,      2, pbedge,  2, "real(8)",    OP_READ), &
                op_arg_dat(p_q,      1, pbecell, 4, "real(8)",    OP_READ), &
                op_arg_dat(p_adt,    1, pbecell, 1, "real(8)",    OP_READ), &
                op_arg_dat(p_res,    1, pbecell, 4, "real(8)",    OP_INC),  &
                op_arg_dat(p_bound, -1, OP_ID,   1, "integer(4)", OP_READ))

            rms = 0.0
            call op_par_loop_5(update, cells, &
                op_arg_dat(p_qold, -1, OP_ID, 4, "real(8)", OP_READ),  &
                op_arg_dat(p_q,    -1, OP_ID, 4, "real(8)", OP_WRITE), &
                op_arg_dat(p_res,  -1, OP_ID, 4, "real(8)", OP_RW),    &
                op_arg_dat(p_adt,  -1, OP_ID, 1, "real(8)", OP_READ),  &
                op_arg_gbl(rms,     2,           "real(8)", OP_INC))
        end do

        rms(2) = sqrt(rms(2) / real(ncell))

        if (op_is_root() .eq. 1 .and. mod(iter, 100) == 0) then
            print *, iter, rms(2)
        end if
    end do

    call op_timers(end_time)
    call op_timing_output()

    print *
    print *, "Time =", end_time - start_time, "seconds"

    if (niter == 1000 .and. ncell == 720000) then
        diff = abs((100.0_8 * (rms(2) / 0.0001060114637578_8)) - 100.0_8)

        print *
        write (*, "(A, I0, A, E16.7, A)") " Test problem with ", ncell , &
            " cells is within ", diff, "% of the expected solution"

        if(diff < 0.00001_8) THEN
            print *, "Test PASSED"
        else
            print *, "Test FAILED"
        end if
    end if

    call op_exit()
end program airfoil
