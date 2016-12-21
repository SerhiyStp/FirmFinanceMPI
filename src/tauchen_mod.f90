!	----------------------------------------------------------------------
!	File name: Tauchen.f90
!
!	This routine calculates the transition probability matrix of a Markov 
!	process. It is used to approximate the AR(1) process of shock. And the
!	calculated transition matrix is used for numerical integration when 
!	evaluating Emax functions.
!	----------------------------------------------------------------------	



module TauchenMod

!INCLUDE 'link_fnl_shared.h'
!use Numerical_Libraries


implicit none


contains

subroutine TauchenTr(mean, sigy, lambday, ygrid, prob)
! inputs are all statistics about y; ygrid should be already taken logs
!should be the case that evenly spaced!


	implicit none

	integer :: node, i, j
	real(8), intent(in) :: mean, sigy, lambday, ygrid(:)  ! the mean should be the unconditional mean of the vector ygrid(:)
	real(8), intent(out) :: prob(size(ygrid),size(ygrid))
	real(8) :: sig, ymax, ymin, ystep, cond_mean
	real(8) :: cdf2, cdf1

	node = size(ygrid)
	
	sig = sigy*dsqrt(1.0d0-lambday**2.0d0)
    ! sigy is the sd of y: the earnings process;  
    !sig is the sd for the error term
    

	ymax = ygrid(node)
	ymin = ygrid(1)
	ystep = (ymax - ymin)/(node-1)  !should be the case that evenly spaced!

	!	calculate transition probability matrix

	do i = 1, node
		cond_mean = (1.0d0-lambday)*mean + lambday*ygrid(i)
		call normal_01_cdf((ygrid(1) - cond_mean + &
		    0.5d0*ystep)/sig, cdf2)
		!prob(i,1) = dnordf((ygrid(1) - cond_mean + 0.5d0*ystep)/sig)
        prob(i,1) = cdf2
		do j = 2, node-1
            call normal_01_cdf((ygrid(j) - cond_mean + &
                0.5d0*ystep)/sig, cdf2)
            call normal_01_cdf((ygrid(j) - cond_mean - &
                0.5d0*ystep)/sig, cdf1)
			!prob(i,j) = dnordf((ygrid(j) - cond_mean + 0.5d0*ystep)/sig)	&
					  !- dnordf((ygrid(j) - cond_mean - 0.5d0*ystep)/sig)
            prob(i,j) = cdf2 - cdf1
		end do
        call normal_01_cdf((ygrid(node) - cond_mean - &
            0.5d0*ystep)/sig, cdf1)
		!prob(i,node) = 1.0d0 - dnordf((ygrid(node) - cond_mean - 0.5d0*ystep)/sig)
        prob(i,node) = 1.0d0 - cdf1
	end do

end subroutine TauchenTr


subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
  real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
  real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
  real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
  real ( kind = 8 ), parameter :: b1 = 3.8052D-08
  real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
  real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
  real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
  real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
  real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ) q
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end subroutine normal_01_cdf

end module TauchenMod
