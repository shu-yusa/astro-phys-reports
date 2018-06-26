!====================================================================!
!     Residual flux in a spectrum line is gained by solving the      !
!     equation of Radiative transfar using various apporoximations.  !
!     For some cases, it can be written in terms of Voigt function.  !
!     Since Voigt function is defined by integral which cannot be    !
!     calculated analytically, we have to integrate it numerically   !
!     to plot residual flux.                                         !
!     In this program, Voigt function and residual flux are          !
!     calculated. Simpson's formula is used to carry out the         !
!     integration. The results are written in files, and we can plot !
!     them by using, say, gnuplot.                                   !
!====================================================================!
      program  Spectrum_line  
!--------------------------------------------------------------------!
!     a.b-------Range of integration of Voigt function.
!     c---------One of the arguments of Voigt function 'H(c,v)'.
!     hh--------Sampling interval of Voigt function and residual flux.
!     v0,v1-----Residual flux 'R' is calculated in the range [v0,v1].
!     v(i)------One of the arguments of Voigt function. Residual flux 
!     is a function of v, and we'll plot 'R' to 'v'.
!     R(k,j,i)--Residual flux.     
!     beta------beta=beta0*H(c,v)
!--------------------------------------------------------------------!
      real*8 a, b, c, hh
      real*8 v0, v1, v(10000), H(10000)
      real*8 R, Br(3), beta(5,10000), beta0(5)
      external R
      integer i, j, k, n, m
      parameter(c=1.0d-3, v0=0.0d0, v1=1.2d1, n=600)
      character FNAME(15)*40
      data FNAME /'residual_flux_B1_beta_0_abs.dat', & 
     &            'residual_flux_B1_beta_1_abs.dat', &
     &            'residual_flux_B1_beta_2_abs.dat', &
     &            'residual_flux_B1_beta_3_abs.dat', &
     &            'residual_flux_B1_beta_4_abs.dat', &
     &            'residual_flux_B2_beta_0_abs.dat', &
     &            'residual_flux_B2_beta_1_abs.dat', &
     &            'residual_flux_B2_beta_2_abs.dat', &
     &            'residual_flux_B2_beta_3_abs.dat', &
     &            'residual_flux_B2_beta_4_abs.dat', &
     &            'residual_flux_B3_beta_0_abs.dat', &
     &            'residual_flux_B3_beta_1_abs.dat', &
     &            'residual_flux_B3_beta_2_abs.dat', &
     &            'residual_flux_B3_beta_3_abs.dat', &
     &            'residual_flux_B3_beta_4_abs.dat'/
      data Br /1.0d0,2.0d0,3.0d0/
      hh = (v1 - v0) / dble(n)
      
      v(1)=v0
!---- We calcutate for cases beta0=1,10,100,1000,10^4.--------

      do j=1, 5
          beta0(j) = 1.0d0 * 10 ** (j - 1)
      enddo
      
      do j=1, 15
          open(unit=10*j, file=FNAME(j))
      enddo       

      do k=1, 3
        do j=1, 5
           do i=1, n
              v(i+1) = v(i) + hh
              call Simpson2(a, b, m, v, c, i, H)
              beta(j,i) = beta0(j) * H(i)
              write((5*(k-1)+j)*10,*) v(i), R(beta(j,i), Br(k))
           enddo
           write(6,*) 5 * (k - 1) + j, '-th integration has completed.'
       enddo
      enddo
      
      
      do i=1,15    
         close(10*i)
      enddo

      write(6,*) 'Completed!'
      write(6,*) 'lower limit was', a
      write(6,*) 'upper limit was', b
 
      stop
      end
!=========================================================================
!     Integrand of Voigt-function.
!=========================================================================
      real*8 function f(a, v, y)
      implicit none
      real*8 :: a, v, y, pi
      parameter(PI=3.14159265389793d0)
      
      f = a / pi * exp(-y**2) / ((v - y) ** 2 + a ** 2)
      
      return
      end
!====================================================================!
!    Residual flux
!====================================================================!
      real*8 function R(beta,Br)
      implicit none
      real*8 :: beta, Br, S

      R = (sqrt(3.0d0) + Br / (1.0d0 + beta)) / (dsqrt(3.0d0) + Br)
      return
      end
!====================================================================!
!     'a' and 'b' are the lower limit and upper limit of integration !
!     respectively
!     'c' is the first argument of f(c,v,y).
!====================================================================!
      subroutine Simpson2(a, b, m, v, c, j, H)
      implicit none
      real*8 :: a, b, c, H(10000)
      real*8 v(10000), hh, S1, S2, vv
      integer i, j, n, m
      real*8 :: eps
      real*8, external :: f
      parameter(eps=1.0d-10,hh=1.0d-4) 
      
      a = 0.0d0
      b = 0.0d0                
!----Determine the upper limit of integration.----
      do 
         if (f(c,v(j),b) <= eps) then
            exit
         else
            b = b + hh 
         end if
      end do

!----Determine the lower limit of integration.----
      do 
         if (f(c,v(j),a) <= eps) then 
            exit
         else 
            a = a - hh
         endif 
      end do

      n = int((b - a) / hh)  
      b = a + n * hh              !Redifine 'b' so that 'a'
                                  !and 'b' to be consistent with n.
      
      S1 = 0.0d0     
      do i=2, n-2, 2
          S1 = S1 + f(c,v(j),a+i*hh)
      end do
      
      S2 = 0.0d0
      do i=1, n-1, 2
        S2 = S2 + f(c,v(j),a+i*hh)
      end do
      
      H(j) = (f(c,v(j),a) + 2.0*S1 + 4.0*S2 + f(c,v(j),b)) * hh / 3.0d0
      m = n
      return
      end
