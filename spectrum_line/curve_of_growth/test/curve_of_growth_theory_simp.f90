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
      integer :: i, j
      integer, parameter :: n1 = 5000
      integer, parameter :: n2 = 10000
      integer, parameter :: n3 = 19000
      integer, parameter :: nn = 160
      real(8), parameter :: A0 = 0.5d0
      real(8), parameter :: c = 1.0d-3
      real(8), parameter :: v0 = 0.0d0
      real(8), parameter :: v11 = 5.0d1
      real(8), parameter :: v12 = 1.0d2
      real(8), parameter :: v13 = 1.9d2
      real(8), external :: R
      real(8) :: hh1, hh2, hh3, h1
      real(8) :: v, H
      real(8), dimension(200) :: beta0, W
      real(8) :: fac

      hh1 = (v11 - v0) / dble(n1)
      hh2 = (v12 - v0) / dble(n2)
      hh3 = (v13 - v0) / dble(n3)

      open(unit=10, file='Voigt.dat3')
      
      h1 = 5.0d-2

      do i=1, nn+1
        beta0(i) = 10.0d0 ** (- 2.0d0 + h1 * dble(i-1))
      enddo

      do j=1, 100
        v = 0.0d0
        fac = 2.0d0 * A0 * beta0(j)
!       call trap(c, v, H)
        call Simpson2(v,c,H)
        W(j) = 0.5d0 * fac * H / (1.0d0 + beta0(j) * H)
        do i=2, n1-1
          v = v + hh1
!         call trap(c, v, H)
          call Simpson2(v,c,H)
          W(j) = W(j) + fac * H / (1.0d0 + beta0(j) * H)
        end do
        v = v + hh1
!       call trap(c, v, H)
        call Simpson2(v,c,H)
        W(j) = W(j) + 0.5d0 + fac * H / (1.0d0 + beta0(j) * H)
      
        W(j) = W(j) * hh1
 
        write(10,*)  log10(beta0(j)), log10(W(j))
        write(6,*) 'log W =' , log10(W(j))
        write(6,*) j,'回目終了'
      end do

      do j=101, 120
        v = 0.0d0
        fac = 2.0d0 * A0 * beta0(j)
!       call trap(c, v, H)
        call Simpson2(v,c,H)
        W(j) = 0.5d0 * fac * H / (1.0d0 + beta0(j) * H)
        do i=1, n2
          v = v + hh2
!         call trap(c, v, H)
          call Simpson2(v,c,H)
          W(j) = W(j) + fac * H / (1.0d0 + beta0(j) * H)
        end do
        v = v + hh2
!       call trap(c, v, H)
        call Simpson2(v,c,H)
        W(j) = W(j) + 0.5d0 + fac * H / (1.0d0 + beta0(j) * H)
      
        W(j) = W(j) * hh2

        write(10,*) log10(beta0(j)),log10(W(j))
        write(6,*) 'log W =' , log10(W(j))
        write(6,*) j,'回目終了'
      enddo 

      do j=121, 161
        v = 0.0d0
        fac = 2.0d0 * A0 * beta0(j)
        call Simpson2(v,c,H)
!       call trap(c, v, H)
        W(j) = 0.5d0 * fac * H / (1.0d0 + beta0(j) * H)
        do i=1, n3
          v = v + hh3
!         call trap(c, v, H)
          call Simpson2(v,c,H)
          W(j) = W(j) + fac * H / (1.0d0 + beta0(j) * H)
        enddo
        v = v + hh3
!       call trap(c, v, H)
        call Simpson2(v,c,H)
        W(j) = W(j) + 0.5d0 + fac * H / (1.0d0 + beta0(j) * H)
      
        W(j) = W(j) * hh3
        write(10,*) log10(beta0(j)), log10(W(j))
        write(6,*) j,'回目終了'
      enddo

      close(10)
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
      real(8), parameter :: eps = 1.0d-8
      real(8), parameter :: hh = 0.001d0
      real(8) :: a, b
      real(8) :: S1, S2
      
      a = 0.0d0
      b = 0.0d0                
      a = - 3.0d0
      b = 3.0d0

!    Determine the upper limit of integration.
!     do 
!       if(f(c,v,b) < eps) then  
!         exit
!       else
!         b = b + hh 
!       endif
!     end do
!     
!    Determine the lower limit of integration.
!     do
!       if(f(c,v,a) < eps) then 
!         exit
!       else 
!         a = a - hh
!       endif 
!     end do

      n = int((b - a) / hh)  
      b = a + dble(n) * hh    ! Redifine 'b' so that 'a' and 'b' 
                              ! to be consistent with n.
      S1 = 0.0d0     
      do i = 2, n-2, 2
        S1 = S1 + f(c,v,a+dble(i)*hh)
      end do
     
      S2 = 0.0d0
      do i=1, n-1, 2
        S2 = S2 + f(c,v,a+dble(i)*hh)
      end do
      
      H = (f(c,v,a) + 2.0d0*S1 + 4.0d0*S2 + f(c,v,b) ) * hh / 3.0d0

      return
      end subroutine
!======================================================================!
      subroutine trap(c, v, H)
      implicit none
      integer :: i, imesh
      real(8), intent(in) :: c, v
      real(8), intent(out) :: H
      real(8), external :: f
      real(8), parameter :: dy = 0.002d0
      real(8) :: ymin, ymax, y
      
      ymin = - 3.0d0!min(v,0.0d0) - 10.0d0
      ymax = 3.0d0!max(v,0.0d0) + 10.0d0
      imesh = nint((ymax - ymin) / dy)
      ymax = ymin + dble(imesh) * dy

      H = 0.5d0 * (f(c,v,ymin) + f(c,v,ymax))
      y = ymin + dy
      do i=2, imesh-1
        H = H + f(c,v,y)
        y = y + dy
      end do
      H = H * dy
     
      return
      end subroutine
