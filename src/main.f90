Program  firm_finance

    use global_variables_mod
    use functions_mod
    use grids_mod
    use LinInterpModule
    use omp_lib
    use solve_values_ss_mod
    use policy_to_mm
    use calibration_estimation_mod

    implicit none
    integer :: i
    integer :: ik, id, iz

    !CALL OMP_SET_NUM_THREADS(threadnum)

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    !print *, 'Hello from process ', my_rank

    ! Conversion between (ik,id,iz) and iloc
    i = 1
    do iz = 1, nz
        do id = 1, nd
            do ik = 1, nk
                i_k(i) = ik   
                i_d(i) = id
                i_z(i) = iz
                kdz_i(ik,id,iz) = i
                i = i + 1
            end do 
        end do
    end do

    nii = int(dble(ni-1.0d0)/dble(nproc)) + 1
    if (my_rank == 0) then
        print *, nproc, ' processes alltogether'
    end if
    itop = my_rank*nii + 1
    iend = min((my_rank+1)*nii, ni)
    print *, 'Process ', my_rank, ' will solve nodes ', itop, ':', iend
    allocate(firmv_loc(nii), firmv_agg(nii*nproc))
    allocate(pol_k_loc(nii), pol_k_agg(nii*nproc))
    allocate(pol_d_loc(nii), pol_d_agg(nii*nproc))
    allocate(pol_e_loc(nii), pol_e_agg(nii*nproc))
    if (my_rank == idmaster .and. scr_id .ne. 6) open(scr_id, file='results_screen.txt')

    call Main_ss_cali

    deallocate(firmv_loc, firmv_agg, pol_k_loc, pol_k_agg, &
        pol_d_loc, pol_d_agg, pol_e_loc, pol_e_agg)

    if (my_rank == idmaster .and. scr_id .ne. 6) close(scr_id)
    call MPI_Finalize(ierr)

end  Program firm_finance
    
    
    
  
    
    
    
