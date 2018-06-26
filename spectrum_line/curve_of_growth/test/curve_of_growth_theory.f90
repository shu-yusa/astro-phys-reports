!====================================================================!
!     Residual flux in a spectrum line  is gained by solving equation!
!     of Radiative transfar using various apporoximations.           !
!     For some cases, it can be written in terms of Voigt function.  !
!     Since Voigt function is defined by integral which cannot be    !
!     calculated analytically, we have to integrate it numerically   !
!     to plot residual flux.  In this program, Voigt function and    !
!     residual flux are calculated. Simpson's formula is used to     !
!     carry out the integration. The results are written in files,   !
!     and we can plot them by using, say, gnuplot.                   !  
!====================================================================!
      program  Spectrum_line  
!--------------------------------------------------------------------!
!     Definitions of variables and parameters.                       !
!     a.b-----Range of integration of Voigt function.                ! 
!     c-------One of the arguments of Voigt function 'H(c,v)'.       !
!     hh------Sampling interval of Voigt function and residual flux. !
!     v0,v1---Residual flux 'R' is calculated in the range [v0,v1].  !
!     v(i)----One of the arguments of Voigt function. Residual flux  !
!             is a function of v, and we'll plot 'R' to 'v'.         !
!     R(k,j,i)--Residual flux.                                       !
!     beta------beta=beta0*H(c,v)                                    !
!--------------------------------------------------------------------!
      implicit none
      integer :: i, j, n_beta
      real(8), parameter :: A0 = 0.5d0
      real(8), parameter :: c = 1.0d-3
      real(8), allocatable, dimension(:) :: beta0, W, fac
      real(8), parameter :: logb_min = -2.0d0
      real(8), parameter :: logb_max = 6.0d0
      real(8), parameter :: d_logb = 0.05d0
      real(8), parameter :: dv = 0.02d0
      
      n_beta = nint((logb_max - logb_min) / d_logb)
      allocate(beta0(n_beta+1), W(n_beta+1), fac(n_beta+1))

!$    call omp_set_num_threads(8)
!$OMP parallel default(none) shared(n_beta,fac,beta0,W) &
!$    private(i,j)
!$OMP do
      do i=1, n_beta
        beta0(i) = 10.0d0 ** (- 2.0d0 + d_logb * dble(i-1))
        fac(i) = 2.0d0 * A0 * beta0(i)
      enddo
!$OMP end do

!$OMP do
      do j=1, 80
        call equ_width(c, fac(j), beta0(j), 0.01d0, 3.5d0, W(j))
 
      end do
!$OMP end do

!$OMP do
      do j=81, 100
        call equ_width(c, fac(j), beta0(j), 0.01d0, 15.0d0, W(j))

      enddo 
!$OMP end do

!$OMP do
      do j=101, 161
        call equ_width(c, fac(j), beta0(j), 0.05d0, 120.0d0, W(j))
      enddo
!$OMP end do
!$OMP end parallel

      open(unit=7, file='Voigt.dat3')
      do j=1, n_beta
        write(7,*)  log10(beta0(j)), log10(W(j))
      end do
      close(7)
      write(6,*) 'Completed!'
 
      stop
      end program
!======================================================================!
      function f(a, v, y) result(g)
!----------------------------------------------------------------------!
!     Integrand of Voigt-function.                                     !
!----------------------------------------------------------------------!
      implicit none
      real(8), intent(in) :: a, v, y
      real(8), parameter :: PI = 3.141592653589793d0
      real(8) :: g
      
      g = a / PI * exp(- y * y) / ((v - y) ** 2 + a * a)
      
      return
      end function
!====================================================================!
      subroutine Simpson2(v,c,H)
!--------------------------------------------------------------------!
!     'a' and 'b' are the lower limit and upper limit of integration !
!     respectively.                                                  !
!     'c' is the first argument of f(c,v,y).                         !
!--------------------------------------------------------------------!
      implicit none
      integer :: i, n
      real(8), intent(in) :: c, v
      real(8), intent(out) :: H
      real(8), external :: f
      real(8), parameter :: hh = 0.002d0
      real(8) :: a, b
      real(8) :: S1, S2
      
      a = - 3.0d0
      b = 3.0d0
      n = int((b - a) / hh)  
      b = a + dble(n) * hh

      S1 = 0.0d0
      S2 = 0.0d0
      do i=2, n-2, 2
        S1 = S1 + f(c,v,a+dble(i)*hh)
        S2 = S2 + f(c,v,a+dble(i-1)*hh)
      end do
      S2 = S2 + f(c,v,a+dble(n-1)*hh)
     
      
      H = (f(c,v,a) + 2.0d0*S1 + 4.0d0*S2 + f(c,v,b)) * hh / 3.0d0

      return
      end subroutine
!======================================================================!
      subroutine equ_width(c, fac, beta0, dv, vmax, W)
      implicit none
      integer :: i, n
      real(8), intent(in) :: c, fac, beta0, dv, vmax
      real(8), intent(out) :: W
      real(8) :: W1, W2, v, H

      n = nint(vmax/dv)
      v = 0.0d0
      W1 = 0.0d0 ; W2 = 0.0d0
      call Simpson2(v,c,H)
      W = fac * H / (1.0d0 + beta0 * H)
      do i=2, n-2, 2 
        v = v + dv ! v = dble(i-1) * hh1
        call Simpson2(v,c,H)
        W1 = W1 + fac * H / (1.0d0 + beta0 * H)
        v = v + dv ! v = dble(i) * hh1
        call Simpson2(v,c,H)
        W2 = W2 + fac * H / (1.0d0 + beta0 * H)
      end do
      v = v + dv
      call Simpson2(v,c,H)
      W1 = W1 + fac * H / (1.0d0 + beta0 * H)
      v = v + dv
      call Simpson2(v,c,H)
      W = W + fac * H / (1.0d0 + beta0 * H)
      W = (W + 2.0d0 * W1 + 4.0d0 * W2) * dv / 3.0d0

      return
      end subroutine
