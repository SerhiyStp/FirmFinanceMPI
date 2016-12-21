module grids_mod

    ! INCLUDE 'link_fnl_shared.h'
    !!include  'mkl.fi'

    use global_variables_mod
    use TauchenMod
    use LinInterpModule
    use functions_mod

    !use Numerical_Libraries
    ! USE IMSL_LIBRARIES


    implicit none

contains

    subroutine make_grids

        implicit none

        real (8) :: temp, vec_z_0(nz),vec_zb_0(nz_big), mean,kmax,  plimit
        integer :: i,j,tempi,ierror,ipart,itemp

        !==============================
        !     ! productivity process
        !==============================

        do i=1,nz
            temp=-4.0*sig_e+8.0*sig_e/(nz-1)*(i-1);
            vec_z_0(i)=temp; 
            vec_z(i)=exp(temp);    ! this is at levels, should be careful when used
        enddo  

        mean=0.0d0
        ! Get Tauchen transition matrix
        temp=sig_e/(dsqrt(1.0d0-pho_z**2.0d0));
        call TauchenTr(mean, temp, pho_z, vec_z_0, Tran_z);

        !now for a finer spaced z:  vec_zb(nz_big),Tran_zb(nz_big,nz_big),Tran_zb_cum(nz_big,nz_big)  ! vec_z(nz) is the levels, not the logs of z
        do i=1,nz_big
            temp=-4.0*sig_e+8.0*sig_e/(nz_big-1)*(i-1);
            vec_zb_0(i)=temp; 
            vec_zb(i)=exp(temp);    ! this is at levels, should be careful when used
        enddo  
        mean=0.0d0
        ! Get Tauchen transition matrix
        temp=sig_e/(dsqrt(1.0d0-pho_z**2.0d0));
        call TauchenTr(mean, temp, pho_z, vec_zb_0, Tran_zb);

        ! to obtain the cumulative transition matrix
        Tran_zb_cum=0.0;
        do i=1,nz_big
            do j=1,nz_big
                tempi=max(j-1,1);
                Tran_zb_cum(i,j)= Tran_zb_cum(i,tempi)+ Tran_zb(i,j);
            enddo
        enddo

        !==============================
        !    the discrete space for k
        !==============================

        !kmax is set such that (1-tax_c)*zmax*kmax**theta =  kmax * delta_k

        temp=(1-tax_c)*exp(4.0*0.4)*Az/delta_k;
        kmax=temp**(1/(1.0-theta));

        vec_k(nk)=kmax;
        itemp=nk/5;

        do j=1,itemp
            i=nk-j;
            vec_k(i)=vec_k(i+1)*(1-delta_k)**(10.0*0.25);
        enddo

        !print*, vec_k(nk-itemp);
        vec_k(1)=vec_k(nk-itemp)*(1-delta_k)**(2.0*0.4*(nk-itemp));
        temp=vec_k(nk-itemp);

        do i=1,nk-itemp
            vec_k(i)=vec_k(1)+(temp-vec_k(1))/(nk-itemp-1)*(i-1);
        enddo

        ipart=nkp/5*4;
        temp=vec_k(nk-itemp);

        do i=1,ipart
            vec_kp(i)=vec_k(1)+(temp-vec_k(1))/(ipart-1)*(i-1);
        enddo

        do i=ipart+1,nkp  
            vec_kp(i)=temp+(vec_k(nk)-temp)/(nkp-ipart)*(i-ipart);
        enddo

        !==============================
        !    the discrete space for debt p
        !==============================
        temp=theta*Az/(rrate+delta_k);

        kss=temp**(1/(1.0-theta));
        plimit=0.9*kss;  !taken from Deangelo et al. paper

        do i=1,nd
            vec_d(i)=-0.5*plimit+(1.5*plimit)/(nd-1)*(i-1);
        enddo

        !more space for dp
        do i=1,ndp
            vec_dp(i)=-0.5*plimit+(1.5*plimit)/(ndp-1)*(i-1);
        enddo

        !print*, 'kss, kmax', kss, kmax; !kss is about 200

        !==============================
        ! Save and print out 
        !==============================

        if (my_rank == idmaster) then
            open(unit=1,file='tauchen_grid.txt',status='replace',iostat=ierror)
            write(1,*) '**********************************************************'
            write(1,*) '************              Grids            ***************'
            write(1,*) 'z-y ; levels (here not taking logs)'
            do i=1,nz
                write(1,*) vec_z(i)
            end do
            close(unit=1)
            open(unit=1,file='k_space.txt',status='replace',iostat=ierror)
            ! write(1,*) 'k space'
            do i=1,nk
                write(1,*)  vec_k(i)
            end do
            close(unit=1)

            write(1,*) 'kp space'
            do i=1,nkp
                write(1,*)  vec_kp(i)
            end do

            open(unit=1,file='debt_space.txt',status='replace',iostat=ierror)
            !write(1,*) 'debt_space'
            do i=1,nd
                write(1,*) vec_d(i)
            end do
            close(unit=1)
        end if

        1   format(20f20.5)        
        11  format(2000a)

    end subroutine make_grids

end module
