module solve_values_ss_mod  !to solve all value functions at the steady state

    use global_variables_mod
    use functions_mod
    use grids_mod
    use LinInterpModule
    use omp_lib

    implicit none

contains

    subroutine solve_v(last_iter)  
        ! solve for the firm's value function; 
        ! here is just updating  the Bellman eq by one iteration, 
        ! not getting the final convergence yet. 

        integer:: ik, id, iz,ikp, idp,inx12(2), izp
        real(8):: k,dd,zlevel,kp,ddp,equity,tempconti, temp,nextfv
        real (8), ALLOCATABLE ::  tempmatrix(:,:),  tempmatrix2(:,:,:)
        integer :: iloc
        integer :: last_iter

        allocate( tempmatrix(nkp,ndp), tempmatrix2(nkp,ndp,nz))

        if (my_rank == idmaster) write(scr_id, *), 'starting to interpolate'

        !*** interpolate first, this is used for all ik, id, iz
        !$omp parallel do private(izp, ikp, idp, kp, ddp,  nextfv)
        do ikp=1,nkp
            do idp=1,ndp
                kp = vec_kp(ikp);
                ddp = vec_dp(idp); !ddp is for debt in the next period
                do izp=1,nz  
                    nextfv = lininterp2(kp, ddp, vec_k, &
                        vec_d, firmv(:,:,izp)) ! relatively slow
                    tempmatrix2(ikp,idp,izp) = nextfv    
                enddo
            enddo
        enddo  
        !$omp end parallel do
        if (my_rank == idmaster) write(scr_id, *), 'end interpolate';

        !$omp parallel do private(ik, id, iz, ikp, idp, inx12, izp, k, dd, zlevel, kp, ddp, equity, tempconti, tempmatrix, temp, nextfv)
        do iloc = itop, iend
            !do iz = 1, nz
            !do id = 1, nd
            !do ik = 1, nk
            ik = i_k(iloc)
            id = i_d(iloc)
            iz = i_z(iloc)
            k = vec_k(ik);
            dd = vec_d(id); !dd is for debt
            zlevel = vec_z(iz);
            do ikp = 1, nkp
                do idp = 1,ndp
                    kp = vec_kp(ikp);
                    ddp = vec_dp(idp); !ddp is for debt in the next period
                    equity = (1.0-tax_c)*k**theta*zlevel*Az &
                        + delta_k*k*tax_c &
                        + ddp &
                        - dd*(1.0+rrate*(1.0-tax_c)) &
                        - (kp-(1.0-delta_k)*k) &
                        - k_adjustmentcost(k,kp) &
                        + eta*dd*indicatorfunc(dd)
                    tempconti = 0.0;
                    do izp = 1, nz    
                        tempconti = tempconti + &
                            Tran_z(iz,izp)*tempmatrix2(ikp,idp,izp);
                    enddo
                    tempmatrix(ikp,idp) = equity + equitycost(equity) + &
                        tempconti/(1+rrate);
                enddo
            enddo

            temp=maxval(tempmatrix);
            !firmv_new(ik,id,iz)=temp; 
            firmv_loc(iloc-itop+1) = temp
            inx12=maxloc(tempmatrix);
            ikp=inx12(1);
            idp=inx12(2);
            !pol_k(ik,id,iz)=vec_kp(ikp);
            pol_k_loc(iloc-itop+1) = vec_kp(ikp)
            !pol_debtp(ik,id,iz)=vec_dp(idp);
            pol_d_loc(iloc) = vec_dp(idp)
            kp=vec_kp(ikp);
            ddp=vec_dp(idp); !ddp is for debt in the next period

            equity = (1.0-tax_c)*k**theta*zlevel + delta_k*k*tax_c + ddp &
                - dd*(1.0+rrate*(1.0-tax_c)) - (kp-(1-delta_k)*k) - &
                k_adjustmentcost(k,kp) + eta*ddp*indicatorfunc(ddp)

            !pol_equity(ik,id,iz)=equity;
            pol_e_loc(iloc-itop+1) = equity
            !enddo
            !enddo
            !enddo
        enddo
        !$omp end parallel do

        call MPI_ALLGATHER(firmv_loc(1), nii, MPI_DOUBLE_PRECISION, &
            firmv_new(1,1,1), nii, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        if (last_iter == 1) then
            call MPI_ALLGATHER(pol_k_loc(1), nii, MPI_DOUBLE_PRECISION, &
                pol_k(1,1,1), nii, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_ALLGATHER(pol_d_loc(1), nii, MPI_DOUBLE_PRECISION, &
                pol_debtp(1,1,1), nii, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
            call MPI_ALLGATHER(pol_e_loc(1), nii, MPI_DOUBLE_PRECISION, &
                pol_equity(1,1,1), nii, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        end if

        !call MPI_ALLGATHER(firmv_loc(1), nii, MPI_DOUBLE_PRECISION, &
            !firmv_agg(1), nii, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        !call MPI_ALLGATHER(pol_k_loc(1), nii, MPI_DOUBLE_PRECISION, &
            !pol_k_agg(1), nii, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        !call MPI_ALLGATHER(pol_d_loc(1), nii, MPI_DOUBLE_PRECISION, &
            !pol_d_agg(1), nii, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
        !call MPI_ALLGATHER(pol_e_loc(1), nii, MPI_DOUBLE_PRECISION, &
            !pol_e_agg(1), nii, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

        !firmv_new = reshape(firmv_agg, (/nk, nd, nz/))
        !pol_k = reshape(pol_k_agg, (/nk, nd, nz/))
        !pol_debtp = reshape(pol_d_agg, (/nk, nd, nz/))
        !pol_equity = reshape(pol_e_agg, (/nk, nd, nz/))

        deallocate( tempmatrix, tempmatrix2 )

    end subroutine solve_v  

    subroutine solve_v2 
        integer:: count1
        real(8) :: temp_metric1, temp0, temp1, temp_metric2
        real(8) :: time_start, time_finish

        if (my_rank == idmaster) write(scr_id, *), ' to  find  VFI '

        call make_grids

        !=============================== 
        ! Solve the optimization problem (Steady State)    

        count1 = 0  ! used to control for the loop times 
        temp_metric1=10.0

        firmv =1./(rrate); 

        firmv_new=firmv;

        pol_k=0.0;  
        pol_debtp=0.0;  
        pol_equity=0.0;  

        if (my_rank == idmaster) print *, 'NumIter = ', NumIter
        find_values: do
            if (my_rank == idmaster) write(scr_id, *), 'count1,  Start VF interation ',  count1, temp_metric1
            firmv=firmv_new 
            time_start = MPI_Wtime()
            call solve_v(0) 
            time_finish = MPI_Wtime()
            if (my_rank == idmaster) write(scr_id, '(a, f20.1, a)'), 'elapsed time: ', time_finish - time_start, ' seconds'
            ! firmv is not updated there, but  firmv_new  is the updated VF
            temp_metric1=maxval(abs(firmv_new -firmv))
            if( temp_metric1 < error .or. count1 > NumIter) then 
                exit find_values
            end if
            count1 = count1 + 1      
        end do find_values
        call solve_v(1) 
        if (my_rank == idmaster) call  value_v_print
    end subroutine solve_v2 

    subroutine value_v_print
        integer::i,ierror, ik,id, iz

        open(unit=5,file='v.txt',status='replace',iostat=ierror,POSITION ='REWIND')

        do ik=1,nk
            do id=1,nd
                do iz=1,nz
                    write(5,52) firmv_new(ik,id,iz), pol_k(ik,id,iz), pol_debtp(ik,id,iz),  pol_equity(ik,id,iz)
                end do       
            end do          
        end do       

        5 format(1500a)
        52    format(4f14.6 )     

        Close (unit=5)  

    end subroutine value_v_print 

    subroutine value_v_read
        integer::i,ierror, ik,id, iz

        open(unit=5,file='v.txt',status='OLD',iostat=ierror)

        do ik=1,nk
            do id=1,nd
                do iz=1,nz
                    read(5,52) firmv_new(ik,id,iz), pol_k(ik,id,iz), pol_debtp(ik,id,iz),  pol_equity(ik,id,iz)
                end do       
            end do          
        end do       

        5 format(1500a)
        52    format(4f14.6 )     

        Close (unit=5)  

    end subroutine value_v_read

end module
