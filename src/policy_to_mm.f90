module policy_to_mm

    use global_variables_mod
    use LinInterpModule
    use functions_mod
    use solve_values_ss_mod
    use TauchenMod
    use grids_mod

    use omp_lib 



    implicit none

contains

    !*******************************************************************************
    !*******************************************************************************
    !******************      Generate steady-state distribution from value functions and decision rules    ****************
    !*******************************************************************************
    !*******************************************************************************   




    subroutine  moment_iterate_ss

        implicit none

        integer :: t, time

        !  initilizing the index and real-vectors for the iteration of dist 

        m_index=(nz_big+1)/2; 
        m_index_new=(nz_big+1)/2;      !iz0, integers

        m_r_index_new=vec_zb((nz_big+1)/2);
        m_r_index_new(:,1)=vec_k(10);
        m_r_index_new(:,2)=vec_d(5);
        m_r_index=m_r_index_new;      !for real moments, to be reported; order: k,p, 

        if (my_rank == idmaster) write(scr_id, *), 'Now to reach ss'   
        do t=1, TimeLast
            ! print*,t; 
            m_index=m_index_new; !for integer index
            m_r_index=m_r_index_new;       !for real variables index

            call one_step_dist_update;
            ! ! for given m_index, m_r,   pol_k,  pol_debtp, pol_equity firmv; update and generate new  m_index_new m_r_new in the next period
        end do ! for many times to reach a ss

        if (my_rank == idmaster) write(scr_id, *), 'Now to keep some data, afer reaching ss'      

        !now keep a short panel; get the agg statistics
        do time=1,keep_period+1
            ! print*, time;
            m_index=m_index_new; !for integer index
            m_r_index=m_r_index_new;       !for real variables index

            call one_step_dist_update
            !updating: get the micro records
            m_state_t(:,time,:)=m_r_index;  !  for individual states
            mm_t(:,time,:)=m_data;     ! for individual mm
        enddo

        call panel_mm_calculate    

    end  subroutine  moment_iterate_ss 

    subroutine  one_step_dist_update  
        ! for given m_index, m_r,    pol_k,  pol_debtp, pol_equity firmv
        !to get one_step_dist_update: using policies, to update the whole dist of state variables, for tomorrow
        ! the results are in m_index_new    m_r_new
        implicit none

        real(8) ::  temp, testprob(nz_big),zlevel, k, dd, equity, kp, fv,ddp
        integer ::  inx1(1),iz,i

        !$omp parallel do private (temp, testprob,zlevel, k, dd, equity, kp, fv,ddp, inx1,iz, i)   

        do i=1,Nindi
            iz=m_index(i);
            zlevel=vec_zb(iz);

            k=m_r_index(i,1);
            dd=m_r_index(i,2);


            call random_number (temp)  

            temp=min((Tran_zb_cum(iz,nz_big)-1.0e-7),temp);
            testprob=Tran_zb_cum(iz,:)-temp;
            inx1=minloc(testprob,MASK = testprob .GE. 0.0);

            m_index_new(i)=inx1(1);  ! tomorrow's z'
            m_r_index_new(i,3)=vec_zb(inx1(1));

            !to use policies:

            kp=lininterp3(k, dd,  zlevel, vec_k, vec_d, vec_z, pol_k);
            kp=min(vec_k(nk)-error,kp);
            kp=max(vec_k(1)+error,kp);

            fv=lininterp3(k, dd,  zlevel, vec_k, vec_d, vec_z, firmv_new);
            ddp=lininterp3(k, dd,  zlevel, vec_k, vec_d, vec_z, pol_debtp);

            ddp=min(vec_d(nd)-error,ddp);
            ddp=max(vec_d(1)+error,ddp);

            equity = (1.0-tax_c)*k**theta*zlevel*Az + delta_k*k*tax_c+ddp &
                - dd*(1.0+rrate*(1.0-tax_c)) - (kp-(1.0-delta_k)*k) &
                - k_adjustmentcost(k,kp) + eta*dd*indicatorfunc(dd)

            !the ddp choice is highly nonlinear!!!

            m_r_index_new(i,1)=kp;
            m_r_index_new(i,2)=ddp; 

            !now compute needed moments
            m_data(i,1)=equity/k;
            m_data(i,2)=(kp-(1.0-delta_k)*k)/k; !investment
            m_data(i,3)=(fv+dd)/k; !firm value Q
            m_data(i,4)=k**theta*zlevel*Az/k; !firm income
            m_data(i,5) = max(dd, 0.0d0)/k; !leverage
            m_data(i,6)=max(-dd,0.0+error)/k; ! cash; adding error, just in case that no firms with cash, so problematic for regressions and stats
            m_data(i,6)=max(-dd,0.0+error)/k; ! cash; adding error, just in case that no firms with cash, so problematic for regressions and stats
            m_data(i,7)=k;
            m_data(i,8)=k**theta*zlevel*Az; 
            ! print*,'m_data',  kp,ddp
        end do
        !$omp end parallel do

    end  subroutine  one_step_dist_update       

    subroutine sum_stat(x, stats)
        integer :: n
        real(8) :: x(:)
        real(8) :: stats(2)
        integer :: i
        real(8) :: xsum, x2sum

        n = size(x)
        xsum = 0.0d0
        x2sum = 0.0d0
        do i = 1, n
        xsum = xsum + x(i)
        x2sum = x2sum + x(i)**2.0d0
        end do
        stats(1) = xsum/dble(n)
        stats(2) = x2sum/dble(n) - stats(1)**2d0
        stats(2) = stats(2)*dble(n)/dble(n-1)
    end subroutine sum_stat

    subroutine lin_reg(x, y, a, b, sum_err_sq)
        real(8) :: x(:), y(:) 
        real(8) :: a, b, sum_err_sq
        integer :: i, n
        real(8) :: x_sum, y_sum, xy_sum, x2_sum
        real(8), allocatable:: errs(:)

        n = size(x)
        allocate(errs(n))
        x_sum = 0.0d0
        y_sum = 0.0d0
        xy_sum = 0.0d0
        x2_sum = 0.0d0
        do i = 1, n
        x_sum = x_sum + x(i)
        y_sum = y_sum + y(i)
        xy_sum = xy_sum + x(i)*y(i)
        x2_sum = x2_sum + x(i)*x(i)
        end do
        b = ((xy_sum)/dble(n) - x_sum/dble(n)*y_sum/dble(n))
        b = b/(x2_sum/dble(n)-(x_sum/dble(n))**2d0)
        a = y_sum/dble(n) - b*x_sum/dble(n)
        errs = y - a - b*x
        sum_err_sq = 0.0d0
        do i = 1, n
        sum_err_sq = sum_err_sq - errs(i)*errs(i)
        end do

        deallocate(errs)
    end subroutine lin_reg

    subroutine  panel_mm_calculate
        ! given the short panel, get the agg statistics
        implicit none

        real(8) :: STAT(2), X(keep_period+1,1)
        real(8) :: sum_err_sq
        real(8) ::  alpha,  STAT2(12)

        integer ::  i, j,IPRINT

        real (8), ALLOCATABLE ::  mm_stat_i(:,:,:),xdata(:), xdata_temp(:,:), ydata(:), ydata_temp(:,:)     

        allocate( mm_stat_i(Nindi,8,2), xdata(keep_period*Nindi), &
            xdata_temp(keep_period*Nindi,1), ydata(keep_period*Nindi), &
            ydata_temp(keep_period*Nindi,1))

        if (my_rank == idmaster) write(scr_id, *), 'Now calculate moments';   

        mm_stat_i=-99.0; 
        mm_stat=-99.0;

        do i=1,Nindi
            do j=1,8
                X(1:(keep_period+1),1)=mm_t(i,1:(keep_period+1),j);
                ! 'the time series data for a given i'; 
                !CALL UVSTA(X, STAT, IPRINT = 0);
                call sum_stat(X(:,1), stat);
                mm_stat_i(i,j,1)=STAT(1); ! for the mean
                mm_stat_i(i,j,2)=STAT(2); ! for the variance
            end do
        end do

        if (my_rank == idmaster) write(scr_id, *), 'Now compute the averages across i'   

        do j=1,8
            do i=1,2
                ! print*,  j, i;  
                mm_stat(j,i)=sum( mm_stat_i(:,j,i))/(Nindi); 
                !print*, mm_stat(j,i), j, i;
            end do
        end do

        if (my_rank == idmaster) write(scr_id, *), 'Now OLS'  

        xdata_temp=reshape( mm_t(:,1:keep_period,4),  (/ keep_period*Nindi, 1 /) );
        ydata_temp=reshape( mm_t(:,2:(keep_period+1),4),  (/ keep_period*Nindi, 1 /)  );
        xdata=xdata_temp(1:keep_period*Nindi,1);
        ydata=ydata_temp(1:keep_period*Nindi,1);

        !CALL RLINE (xdata, ydata, alpha, pho_income, STAT=STAT2)
        !ols regression
        call lin_reg(xdata, ydata, &
            alpha, pho_income, &
            sum_err_sq)

        !pho is the auto correlation;  
        sig_income_error=STAT2(11)/(keep_period*Nindi);
        deallocate(  mm_stat_i,xdata, xdata_temp, ydata, ydata_temp)
    end  subroutine     panel_mm_calculate   

    subroutine mm_to_print
        implicit none
        integer::i,ierror

        !mm_stat_i(Nindi,8,2)
        ! for the mean and var of: 
        !      m_data(i,1)=equity/k;                         m_data(i,2)=(kp-(1.0-delta_k)*k)/k; !investment          m_data(i,3)=(fv+dd)/k; !firm value Q
        ! m_data(i,4)=k**theta*zlevel*Az/k; !firm income      m_data(i,5)=max(dd,0.0)/k; !leverage                   m_data(i,6)=max(-dd,0.0)/k; ! cash
        ! means of k, mean of y  

        data_stat=0.0;
        data_stat(1,1)=0.0407;
        data_stat(2,1)= 0.1868 ;  
        data_stat(3,1)= 3.67 ;  
        data_stat(4,1)= 0.1915 ;  
        data_stat(5,1)= 0.2393 ;  
        data_stat(6,1)= 0.1399 ;  
        data_stat(1,2)=0.0079;
        data_stat(2,2)=0.0385;    
        data_stat(5,2)=0.0118;


        if (my_rank == idmaster) then
            !open(unit=13, file='final_mm.txt',status='OLD',iostat=ierror, position='Append')
            open(11, file='final_mm.txt')
            write(11,*) 'Indionsyncratic Volatility: std', sig_e
            WRITE(*,*)
            WRITE(*,*)

            write(11,*) '************ Data moment and Model moment***************'  

            write(11,'(9A14)') 'Data:Mean', 'equity/k', 'I/K', 'Tobin Q', 'income/K', &
                'Debt/k', 'cash/k', 'k', 'y'      
            write(11,"(a14,8f14.5)"), 'Data:Mean',  data_stat(:,1)
            write(11,"(a14,8f14.5)"), 'Model:Mean',  mm_stat(:,1)   ! the means for 8 var
            WRITE(*,*)
            WRITE(*,*)

            write(11,'(9A14)'), 'Data:Variance', 'equity/k', 'I/K', 'Tobin Q', 'income/K',&
                'Debt/k', 'cash/k', 'k', 'y'     
            write(11,"(a14,8f14.5)"), 'Data:Variance',  data_stat(:,2) 
            write(11,"(a14,8f14.5)"), 'Model:Variance',  mm_stat(:,2)   ! the means for 8 var
            WRITE(*,*)
            WRITE(*,*)     
            write(11,*) 'pho_income/k       ', 'variance of innovation to income/k'
            write(11,*) 0.6635, 0.0048 
            write(11,*)  pho_income, sig_income_error        
            !write(11,*) '************Means, Variances, pho, errors***************'    
            !write(11,1) sig_e
            !write(11,1) mm_stat(:,1) ! the means for six var
            !write(11,1) mm_stat(:,2)! the variances for six var
            !write(11,1) pho_income, sig_income_error
            close(11)
        end if

        !51         format(A14,  8f14.5)
        !52         format(9A14)
        !! 53   format(2A14)
        !1   format(1f10.5)  
    end subroutine mm_to_print



end  module policy_to_mm
