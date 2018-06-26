!======================================================================!
!     Title  : transonic_flow.f90                                      !
!     Author : Yusa Shusaku                                            !
!     Data   : 2008-2-18-Mon                                           !
!                                                                      !
!     A Program which solves the following equation with respect to    ! 
!     y by bisection method:                                           !
!          y^4+6*y^2-16*(1+AA*x^2)*y+9 = 0.                            !
!     We assign some value for 'AA' and find the solutions as a        !
!     function of x.                                                   ! 
!     This equation is obtained by solving differential equations      !
!     of hydrodynamics for the problem of Laval nozzle.                !
!     In solving the differential equations, we have assumed that      ! 
!     the flow is static and adiabatic, and that the gas is an ideal   ! 
!     gas.'y' is a Much number, x is a position and 'A' determines     !
!     the shape of the Laval nozzle.                                   !
!----------------------------------------------------------------------!
      module Com_var
      implicit none
      integer, parameter :: kk=4
      real*8, parameter :: xmin=-2.0d0, xmid=0.0d0, xmax=2.0d0
      real*8, parameter :: h=0.1d0, hh=0.01d0, AA=0.15d0
      end module
!======================================================================!
      program transonic_flow
!----------------------------------------------------------------------!
!     This is a main program.                                          !
!     hh....sampling interval                                          !
!     h.....An interval from which we seek the solution.               !
!     ax,bx,cx....We solve the equation in the range [ax,cx],but we    !
!                 divide the region at 'bx'. In the range [ax,bx],     ! 
!                 the solution is always smaller than 1.0, while in    !
!                 the range [bx,cx], the solution is always bigger     !
!                 than 1.0.                                            !
!     ay,by,cy....The equation has solution between [ax,cx].           !
!     y1,y2.......Starting points of bisection method. These values    !
!                 are found by subroutine 'prep'.                      !
!------------------------------------1----------------------------------!
      use Com_var, only : kk, xmin, xmid, xmax, h, hh, AA
      implicit none
      real*8 :: y1, y2, x(10000), y(10000)
      real*8 :: ay, by, cy, x0(4), y0(4), z(10000)
      integer :: i, j, m1, m2, k
      character :: FNAME(6)*30 
      real*8, external :: f, df
      data FNAME/'trans_crit.dat','trans_subsonic.dat',           &
                  'trans_front.dat','trans_crit2.dat',            &
                  'trans_supersonic.dat','trans_back.dat'/
      data x0 /0.0d0,0.0d0,-1.0d0,0.0d0/
      data y0 /1.0d0,0.5d0, 1.0d0,1.6d0/
      parameter(ay=0.0d0,by=1.0d0,cy=15.0d0)

      do k=1,3
         open(unit=100*k, file=FNAME(k), form='formatted')
         open(unit=101*k, file=FNAME(k+3), form='formatted')
      enddo

      do k=1,kk
         m1 = int((xmid - xmin) / hh)
         m2 = int((xmax - xmin) / hh)
!------Checking of parameter x0 and y0-------------------------------
        if((x0(k) /= 0.0d0).and.(y0(k) /= 1.0d0)) then
           write(6,*) 'Format of parameters (x0,y0) is wrong.'
           stop
        endif
        if(x0(k) > 0.0d0) then
           write(6,*)  'Assign x0 a negative value.'
           stop
        endif
!---------------------------------------------------------------------
!      open(unit=100*k,err=99,status='old',
!     *     file=FNAME(k),form='formatted')
!      
!      open(unit=101*k,err=99,status='old',
!     *     file=FNAME(k+kk),form='formatted')
      
      x(1) = xmin
!*****The case that x0 is equal to 0.***************************      
      if(x0(k) == 0.0d0) then     
!------The case that the solution passes through a critical point.---
        if(y0(k) == 1.0d0) then                                  
          do i=1 ,m1+1
             call prep(ay,by,x0(k),y0(k),x(i),y1,y2)
             call bisec(x0(k),y0(k),x(i),y(i),y1,y2)
             x(i+1) = x(i) + hh
          enddo
          do i=m1+2, m2+1
             call prep(by, cy, x0(k), y0(k), x(i), y1,y2)
             call bisec(x0(k),y0(k),x(i),y(i),y1,y2)
             x(i+1) = x(i) + hh
          enddo
          do i=1, m2+1
             write(100*k,*) x(i), y(i)
             write(101*k,*) x(i), y(m2-i+2)
          enddo
!-------Not passing through a critical point.-------------------------------
        elseif(y0(k) < 1.0d0) then
          do i=1, m2+1
             call prep(ay,by,x0(k),y0(k),x(i),y1,y2)
             call bisec(x0(k),y0(k),x(i),y(i),y1,y2)
             write(100*k,*) x(i),y(i)
!            call prep(by,cy,x0(k),y0(k),x(i),f,y1,y2)
!           call bisec(x0(k),y0(k),x,z,y1,y2,f,i)
!           write(101*k,*) x(i),z(i)
!            x(i+1)=x(i)+hh
          enddo
        else 
          do i=1, m2+1
             call prep(by,cy,x0(k),y0(k),x(i),y1,y2)
             call bisec(x0(k),y0(K),x(i),y(i),y1,y2)
             x(i+1) = x(i) + hh
             write(101*(k-2),*) x(i),y(i)
          enddo
        endif
!******The case that x0 is not equal to 0.*************************
      else 
           x(1) = xmin
           m2 = int((x0(k) - xmin) / hh + 0.1d0)   !kiriage
           do i=1,m2
             call prep(ay,by,x0(k),y0(k),x(i),y1,y2)
             call bisec(x0(k),y0(k),x(i),y(i),y1,y2)
             call prep(by,cy,x0(k),y0(k),x(i),y1,y2)
             call bisec(x0(k),y0(k),x(i),z(i),y1,y2)
             x(2*m2+2-i) = x(i)
             y(2*m2+2-i) = z(i)
             x(i+1) = x(i) + hh
           enddo
!-----Only for x=x0, we use Newton method.-------------------------
             call Newton(x0(k),y0(k),x(m2+1),y(m2+1))
             x(m2+1) = x(m2) + hh
           do i=1, 2*m2+1
             write(100*k,*)   x(i), y(i)
             write(101*k,*) - x(i), y(i)
           enddo
!--------------------------------------------------------------------
       endif
             
      enddo              ! enddo for loop of k.
      
      do k=1,3
        close (100*k)
        close (101*k)
      enddo

      stop
      end program
!====================================================================!
      subroutine  prep(a, b, x0, y0, x, y1, y2)
!--------------------------------------------------------------------!
!     A subroutine which finds starting points of bisection method   !
!     in the range [a,b]. These points are such that f(y1)*f(y2) is  !
!     negative.                                                      !
!     Input  : a, b, h, AA, x, f and i                               !
!     Output : y1 and y2                                             ! 
!--------------------------------------------------------------------!
      use Com_var, only : h
      implicit none
      integer :: j, n
      real*8 :: a, b, x
      real*8 :: y, y1, y2, x0, y0
      real*8, external :: f

      n = int((b - a) / h)
      y = a
      do j=1, 2*n
         if (f(x,y,x0,y0)*f(x,y+h,x0,y0) < 0.0d0) exit
         y = y + 0.5d0 * h
      end do

      if (f(x,y,x0,y0)*f(x,y+h,x0,y0) > 0.0d0) then 
        write(6,*) 'We could not find starting points.'
        write(6,*) 'x(i) =',x
        write(6,*) 'f(y) =',f(x,y,x0,y0)
        write(6,*) 'f(y+h) =',f(x,y+h,x0,y0)
        stop
      endif

      y1 = y
      y2 = y + h

      return
      end subroutine
!====================================================================!
      subroutine bisec(x0, y0, x, y, y1, y2)
!--------------------------------------------------------------------!
!     A subroutine which executes bisection method.                  !
!     Input  :  AA, x, y1, y2, f, i                                  !
!     Output :  y                                                    !
!--------------------------------------------------------------------!
      implicit none
      integer :: i, j
      real*8 :: y1, y2, x, y, x0, y0
      real*8,parameter :: eps=1.0d-9
      real*8, external :: f

      j = 1
      if( abs(f(x,y1,x0,y0)) < eps) then 
         y = y1
         return
      else if( abs(f(x,y2,x0,y0)) < eps) then
         y = y2
         return
      endif

      do while(j < 100)
         y = (y1 + y2) * 0.5d0

         if (abs(f(x,y,x0,y0)) < eps) exit

         if ((f(x,y1,x0,y0)*(f(x,y,x0,y0)) > 0.0d0)) & 
          then
            y1 = y
         else 
            y2 = y
         endif

         j = j + 1
         
         if (j > 100) then
            write(6,*)  'There may be something wrong in using    &
     &                   bisection method.'
            write(6,*) y1, y2
            stop
         endif
       
      end do
      end subroutine
!====================================================================!
      function f(x, y, x0, y0)  result(P)
!--------------------------------------------------------------------!
!     Definition of function.                                        !
!--------------------------------------------------------------------!
      use Com_var, only : AA
      implicit none
      real*8, intent(in) :: x, y, x0, y0
      real*8 :: P

     P = ((y*y + 6.0d0)*y   &
        -(y0*y0+3.0d0)**2/y0*(1.0d0+AA*x*x)/(1.0d0+AA*x0*x0))*y+9.0d0
  
      end function
!====================================================================!
      function df(x, y, x0, y0)  result(dP)
!--------------------------------------------------------------------!
!     Definition of derivative function of f.                        !
!--------------------------------------------------------------------!
      use Com_var, only : AA
      implicit none
      real*8, intent(in) :: x, y, x0, y0
      real*8 :: dP

      dP = (y * y + 3.0d0) * y    &
          - (y0 * y0 + 3.0d0)**2/y0*(1.0d0+AA*x*x)/(1.0d0+AA*x0*x0)

      end function
!====================================================================!
      subroutine Newton(x0, y0, x, y)                         
!--------------------------------------------------------------------!
!     subroutine of Mewton method.                                   !
!     Input  : AA,x0,y0,f,df,x,i                                     !
!     Output : y                                                     !
!--------------------------------------------------------------------!
      implicit none
      integer :: j
      real*8  :: x, x0, y0, y
      real*8, external :: f, df

      y = 1.0d0 

      do j=1, 500
         y = y - f(x,y,x0,y0) / df(x,y,x0,y0)
         if(abs(f(x,y,x0,y0)) < 1.0d-9) exit
      end do
      
      if (j > 500) then 
        write(6,*) 'There may be something wrong with Newton method.'
      end if

      end subroutine
