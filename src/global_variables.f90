module global_variables_mod

    implicit none

    include 'mpif.h'

    integer :: my_rank
    integer :: nproc
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: scr_id = 6   ! Where to direct output:
                                       ! - to screen if scr_id = 6
                                       ! - to a file 'results_screen.txt' otherwise
    integer, parameter :: idmaster = 0


    !save

    !real(8),parameter::error = 1.0e-5
    real(8),parameter::error = 1.0e-3
    real(8),parameter::big_num = 100000000.

    integer,parameter::threadnum =36 ! thread num 
    !integer,parameter:: NumIter=1500; ! for the maximal times of   value function iteration 
    integer,parameter:: NumIter=2 ! for the maximal times of   value function iteration 
    integer,parameter::inside_vf_times= 20 ! inside VF iteration times
    integer,parameter::inside_howard_times= 10 ! inside Howard iteration times


    integer,parameter::TimeLast=1000      !!! assume at this period, reaches steady state
    !integer,parameter::Nindi=50000    !!! numbe of individual firms to simulate
    integer,parameter::Nindi=1000    !!! numbe of individual firms to simulate
    integer,parameter::keep_period=50;  !after long periods of simulation, keep a short panel



    integer,parameter::   nk=70; ! number of discrete points in the k space
    integer,parameter::   nd=30; ! number of discrete points in the debt d space

    integer,parameter::   nkp=1000;    ! number of choice points in the k space
    integer,parameter::   ndp=100;  ! number of choice points in the debt d space

    integer,parameter::  nz=9;       ! number of discrete points in the log(z)  space
    integer,parameter::  nz_big=151; ! number of discrete points in the log(z)  space: for simulation purpose

    integer, parameter :: ni = nk*nd*nz
    integer :: nii
    integer :: itop
    integer :: iend
    integer :: i_k(ni), i_d(ni), i_z(ni)
    integer :: kdz_i(nk,nd,nz)

    real(8), allocatable :: firmv_loc(:), firmv_agg(:)
    real(8), allocatable :: pol_k_loc(:), pol_k_agg(:)
    real(8), allocatable :: pol_d_loc(:), pol_d_agg(:)
    real(8), allocatable :: pol_e_loc(:), pol_e_agg(:)

    !vectors
    real(8)::kss, vec_k(nk), vec_d(nd), vec_z(nz), Tran_z(nz,nz), &
        vec_zb(nz_big), Tran_zb(nz_big,nz_big), &
        Tran_zb_cum(nz_big,nz_big)  ! vec_z(nz) is the levels, not the logs of z
    real(8):: vec_kp(nkp), vec_dp(ndp) ! 


    !! Value functions  
    real(8),dimension(nk,nd,nz):: firmv, firmv_new

    !! Policy function
    real(8),dimension(nk,nd,nz):: pol_k,  pol_debtp, pol_equity







    !!! not to calibrate parameters
    real(8),parameter::      delta_k=0.15;  ! annual depreciation rate
    real(8),parameter::      rrate=0.015;  ! annual  interest rate
    real(8),parameter::      tax_c=0.35;  ! the effective corporate tax rate



    real(8),parameter::   lamda1=0.1615  ! equity issuance cost: linear part
    real(8),parameter::   lamda2=0.0041  ! equity issuance cost: quadratic part

    real(8),parameter::   eta=0.0077     ! the cost of holding cash, per unit of cash


    !production
    real(8),parameter::   theta=0.788  ! decreasing return to scale for k in the production

    !ajustment cost for k
    real(8),parameter::   gamma_k=0.0034 !fixed adjustment cost, per unit of k
    real(8),parameter::   a_cost_k=0.1519 ! quadratic part adjustment cost

    !productivity process
    real(8),parameter::   pho_z=0.7280  ! persistence for log(z)

    !real(8),parameter::   sig_e=0.2843  ! std  for the innovations in log(z)

    real(8)::   sig_e   !=0.2843  ! std  for the innovations in log(z)


    real(8),parameter::   Az=1.0/3.0;   !/4.0  ! adjust the level of techology, so that the capital space needed is not so large




    !!!! for interpolations
    !integer,parameter::     KXORD=2;
    ! integer,parameter::     KYORD=2; ! orders of spline
    !   integer,parameter::     KZORD=2; ! orders of spline
    !   
    !   
    ! integer,parameter::     NXDATA=nk;
    ! integer,parameter::     NYDATA=nd;
    ! integer,parameter::     NZDATA=nz;
    ! 
    !
    ! integer,parameter::  NXKNOT=NXDATA+KXORD; 
    !  integer,parameter:: NYKNOT=NYDATA+KYORD;
    !  integer,parameter:: NZKNOT=NZDATA+KZORD;
    !  
    !real(8):: XDATA(NXDATA), YDATA(NYDATA), ZDATA(NZDATA), FDATA(NXDATA,NYDATA,NZDATA), XKNOT(NXKNOT), YKNOT(NYKNOT),ZKNOT(NZKNOT), &
    !    BSCOEF_e(NXDATA,NYDATA,NZDATA),   BSCOEF_kp(NXDATA,NYDATA,NZDATA),  BSCOEF_fv(NXDATA,NYDATA,NZDATA), BSCOEF_dp(NXDATA,NYDATA,NZDATA)


    real(8):: tempmatrix3_fv(nkp,ndp,nz_big), tempmatrix3_kp(nkp,ndp,nz_big), tempmatrix3_dp(nkp,ndp,nz_big)



    ! for simulations

    real(8):: m_r_index(Nindi,3), m_r_index_new(Nindi,3),  m_data(Nindi,8), &
        m_state_t(Nindi,keep_period+1,3),  mm_t(Nindi,keep_period+1,8)

    real(8)::avg_I_K, avg_income_K,avg_equity_K,avg_d_K,avg_cash_K,avg_Q

    real(8)::  mm_stat(8,2),pho_income, sig_income_error

    !mm_stat_i(Nindi,8,2)
    ! for the mean and var of: 
    !      m_data(i,1)=equity/k;                         m_data(i,2)=(kp-(1.0-delta_k)*k)/k; !investment          m_data(i,3)=(fv+dd)/k; !firm value Q
    ! m_data(i,4)=k**theta*zlevel*Az/k; !firm income      m_data(i,5)=max(dd,0.0)/k; !leverage                   m_data(i,6)=max(-dd,0.0)/k; ! cash
    ! means of k, mean of y

    real(8)::  data_stat(8,2)    


    integer::  m_index(Nindi), m_index_new(Nindi)             
    !        for m_r, m_r_new: will be iterating for many times;   mindex  mshort  will serve as storages

end module
