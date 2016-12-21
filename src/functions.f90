module functions_mod

  use global_variables_mod
  
  implicit none
  
  contains
    
    !========================================
    ! Define I matrix
    !========================================
   
function eye_mat(nn)

    implicit none


   integer,intent(in) :: nn
	
  real(8), dimension(nn,nn):: eye_mat
   integer :: i,j

do i=1,nn
do j=1,nn
if (i==j) then
eye_mat(i,j)=1.0;
else
eye_mat(i,j)=0.0;
endif


enddo
enddo


end function

  
function indicatorfunc(x)
    implicit none
    
real(8),intent(in) :: x
real(8)            :: indicatorfunc

if ( x .LE. 0.0) then
indicatorfunc=1.0;
else
indicatorfunc=0.0;
endif


end function
  
  
function indicatorfunc2(x)
    implicit none
    
real(8),intent(in) :: x
real(8)            :: indicatorfunc2

if ( abs(x) .GE. error) then
indicatorfunc2=1.0;
else
indicatorfunc2=0.0;
endif


end function
    
  
  
function equitycost(e)
    implicit none
    
real(8),intent(in) :: e

real(8)            ::equitycost

equitycost=indicatorfunc(e)*(lamda1*e-0.50*lamda2*e**2.0) 

end function



function equitycost_fd(e) !firs order derivative
    implicit none
    
real(8),intent(in) :: e

real(8)            ::equitycost_fd

equitycost_fd=indicatorfunc(e)*(lamda1-lamda2*e) 

end function

    
    
function k_adjustmentcost(k,kp)

    implicit none
    
real(8),intent(in) :: k, kp



real(8)            ::k_adjustmentcost, invest


invest=kp-(1.0-delta_k)*k;
k_adjustmentcost=indicatorfunc2(invest)*k*gamma_k+0.5*a_cost_k*(invest/k)**2.0*k

	



  end function
  
function adj_fd_k(k,kp) !firs order derivative

    implicit none
    
real(8),intent(in) :: k, kp



real(8)            ::adj_fd_k, invest


invest=kp-(1-delta_k)*k;
adj_fd_k=indicatorfunc2(invest)*(gamma_k-0.5*a_cost_k*(invest/k)**2.0-(1-delta_k)*a_cost_k*(invest/k)    ); 


end function
  
 function adj_fd_kp(k,kp) !firs order derivative

    implicit none
    
real(8),intent(in) :: k, kp



real(8)            ::adj_fd_kp, invest


invest=kp-(1-delta_k)*k;
adj_fd_kp=indicatorfunc2(invest)*( a_cost_k*(invest/k)    );

	



  end function       
    end module
    