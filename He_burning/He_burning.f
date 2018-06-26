c=====================================================================
c     A program which solves simultaneous differential equation 
c             dy_j/dx=f_j(x,y_i)(1<=i,j<=4)
c     by Runge-Kutta method.
c=====================================================================
c     In helium burning stage of a star, helium converts to carbon by
c     triple alpha process. Subsequentry, carbon captures alpha to
c     convert oxigen and oxigen captures alpha to convert neon.
c     Although this chain does not end and neon becomes heavier element,
c     we can approximately terminate the reaction up to neon.
c     This approximation is valid as long as the temperature of the
c     star is not so high. In such a situation, the abundance of
c     elements changes according to the four ordinary differential
c     equation. In this program, we solve these four differential
c     equation by Runge-Kutta method. Initial condition is that there
c     is only helium at the beginnig of burning and if helium run out,
c     the program stops. We give the density and temperature of the
c     star and assume that they are constant during the burning. 
c=====================================================================
c-----Definitions of various characters and initial value-------------   
      implicit none
      character FNAME(4)*20
      real*8 s(4,20000),x(20000)
      real*8 x0,x1,d,rho,NA,T9,yr,eps2
      integer n,ni,l,i,j,nfile
      parameter(nfile=4,yr=3.15d7,eps2=1.0d-5)
      parameter(rho=1.0d4,T9=0.15d0,NA=6.0225d23)
      data n,x0,x1,ni/4, 0.0d0, 5.0d0,10/  
      data FNAME /'abundance_He.dat','abundance_C.dat',
     *            'abundance_O.dat','abundance_Ne.dat'/ 
c=====================================================================
c     FNAME(j)----filename of data:(x,y_j)
c     x0,x1-------integration interval[x0,x1]
c     s(j,i)------y_j(i) or y_j(x(i)) 
c     Initial value is given at x0 (x2<x0<x1)
c     "n" is the number of equation.
c     "d" is sampling interval.This is not constant. Its value varies
c     depending on the behavior of the y(i)'s.
c     "nfile" is the number of pairs(x,y(i)).
c     ni...Arrange the sampling interval to d/ni in computation of integration.
c     l....Total number of pair of data.
c     eps2...If the value of X4 become smaller than eps2, then integration stops.   
c     rho...Density.
c     T9...Temperature.
c     yr...1 year in second.
c==========================================================================

c-----Write down initial values and conditions.--------------
      write(*,*) ' *** Numerical Analysis of ODE ***'
      write(*,*) '      (With Runge-Kutta Method)'
      write(*,*)
      write(*,240) n,x0
      write(6,*) 'd  = Sampling interval'
      write(6,*) 'h  = Division interval'
      write(6,*) 'ni = Division number'
      write(6,*) 'sum = sum of abundance'
      l=1
      x(1)=x0
      call INITIAL(s)
      write(6,*) 
      write(*,*) 'Initial Condition'
      write(*,210) x(1),(s(i,1),i=1,n)

c-----Solve the differential equation by Runge-Kutta method.---
      call DRUNGE(l,x,n,s,d,ni)
c-----Write down solution in the file.---------------------------------

c      write(6,*)'Input the output filename.'   !Input the filenames.
c      do j=1,nfile
c      write(6,230) j,'th filename'
c      read(5,*) FNAME(j)
c      end do
      
      do j=1,nfile
        open(100*j,file=FNAME(j))               !Open files.
      end do
      write(6,*) 
      write(6,*) 'Final results--abundance of end products'
      write(6,250) 'At  t =',x(l),'*10^5 yr' 
      write(6,*) 
      write(6,*) ' He = ',s(1,l)
      write(6,*) ' C  = ',s(2,l)
      write(6,*) ' O  = ',s(3,l)
      write(6,*) ' Ne = ',s(4,l)
      write(6,*) 'sum = ',s(1,l)+s(2,l)+s(3,l)+s(4,l)
      
      do i=1,l
 20     do j=1,nfile
          write(100*j,220) x(i),s(j,i)            !Write down solutions in the files.
        end do
      enddo

      do j=1,nfile                        !Close files.
        close(100*j)
      end do
      
c-----Setting of format.--------------------------------------------------         
 240   format(1h ,'number of equations',i2
     */' lower limit=',f5.1)
 210   format(1h ,2x,'x=',f8.4,4x,'y=(',f8.5,','
     *,f8.5,',',f8.5,',',f8.5,')')
 220   format(f10.5,' ',f14.12)
c 230   format(' ',i1,a)
 250  format(' ',a,f7.2,a)
      end
c------------------------------------------------






c-----Insert initial values.---------------------------------
      subroutine INITIAL(s)
      real*8 s(4,20000)
      s(1,1)=1.0d0
      s(2,1)=0.0d0
      s(3,1)=0.0d0
      s(4,1)=0.0d0
      return
      end
c------------------------------------------------






c
c-----Differential equations---------------------------------
c     'tmsc' is the timescale. We should assign it an appropiate 
c     value for computation.
c-----------------------------------------------------------
      subroutine EQATION(x,y,f)
      real*8 x,y(4),f(4),NA,yr,tmsc
      real*8 rho,T9,lam1,lam2,lam3
      external lam1,lam2,lam3
      parameter(NA=6.0225d23,rho=1.0d4,T9=0.15d0)
      parameter(tmsc=1.0d6,yr=3.15d7)
      
      f(1) = (-3.0/1.6d1*rho**2*NA**2*lam1(rho,T9)*y(1)**3
     *       -rho*NA/1.2d1*lam2(rho,T9)*y(1)*y(2)
     *       -rho*NA/1.6d1*lam3(rho,T9)*y(1)*y(3))*yr*tmsc 
      f(2) = (3.0d0/1.6d1*rho**2*NA**2*lam1(rho,T9)*y(1)**3
     *       -2.5d-1*rho*NA*lam2(rho,T9)*y(1)*y(2))*yr*tmsc
      f(3) = (1.0d0/3.0d0*rho*NA*lam2(rho,T9)*y(1)*y(2)
     *       -2.5d-1*rho*NA*lam3(rho,T9)*y(1)*y(3))*yr*tmsc
      f(4) = (1.25d0*rho*NA*lam3(rho,T9)*y(1)*y(3))*yr*tmsc
      
      return
      end
c------------------------------------------------







c-----Driver routine of the Runge-Kutta-Method.-----
c==========================================================================
c     "x(i)" is independent variable x.  s(j,i) represents y_j.
c     eps3 determines sampling interval. It is nothing to do with precision 
c     of the calculation. 
c     Precision is determined mainly by eps.
c     eps2 gives the condition of when to terminate the calculation.
c==========================================================================
      subroutine DRUNGE(l,x,n,s,d,ni)   !Outputs are x,s,l.
      real*8 s(4,20000),y(4),eps,eps2,eps3,f(4)
      real*8 x(20000),h,d,hh(4),dd(4)
      integer ni
      parameter(eps=1.0d-4,eps2=1.0d-4,eps3=1.0d-2)

      y(1)=s(1,1)
      y(2)=s(2,1)
      y(3)=s(3,1)
      y(4)=s(4,1)
 10   x(l+1)=x(l)
c============================================================================================   
c    In the below, we determine the intervals between sampling points and 
c    interval of integration.
c-------------------------------------------------------------------------
      call EQATION(x(l),y,f)
  
      do i=1,n
         hh(i) = dabs(y(i)/(f(i)+1.0d-10))*eps                         
         dd(i) = dabs(1.0d0/(f(i)+1.0d-10))*eps3
         if(hh(i).lt.5.0d-9) then
            hh(i)=5.0d-9
         elseif(hh(i).gt.1.0d-1) then
            hh(i)=1.0d-1
         endif

         if(dd(i).gt.0.1d0) then
            dd(i) = 0.1d0 
         elseif(dd(i).lt.5.0d-4) then
            dd(i) = 5.0d-4
         endif
      enddo

        h = dmin1(hh(1),hh(2),hh(3),hh(4))
        d = dmin1(dd(1),dd(2),dd(3),dd(4))
        ni =int(d/h)
        h=d/real(ni)
        write(6,40) l,'th:',' d=',d,', h=',h,', ni=',ni,
     *              ', sum=',y(1)+y(2)+y(3)+y(4)
 40     format(' ',i7,a,a,f8.5,a,f12.10,a,i6,a,f16.13)
c===============================================================================================
        do i=1,ni
           call RUNGE(h,n,x(l+1),y)   !Caluculate the quantity y_j at x(l+1),that is, s(j,l+1). 
        end do
        l=l+1
        s(1,l)=y(1)
        s(2,l)=y(2)
        s(3,l)=y(3)
        s(4,l)=y(4)
        if(y(1).lt.eps2) goto 30  !If the abundance of He becomes smaller than eps2,stop the culculation.  
      goto 10
 30   return
      end
c------------------------------------------------








c-----Runge-Kutta-Method-------------------------
      subroutine RUNGE(h,n,x,y)                         !Outputs are x,y.
      real*8 f(4),y(4),s1(4),s2(4),s3(4),s4(4),ysav(4)
      real*8 x,xsav,h,hh(4),eps
      parameter(eps=1.0d-5)


      xsav=x
      do  j=1,n
            ysav(j)=y(j)
      enddo
      call EQATION(x,y,f)     !Insert the present value of x and y into the f. 

c-----Calculation of s1.-----------------------------------------------
      
      do j=1,n
            s1(j)=h*f(j)
      enddo
c-----Caliculation of s2---------------------------------------------
      x=xsav+0.5d0*h
      do j=1,n
            y(j)=ysav(j)+0.5d0*s1(j)
      enddo

      call EQATION(x,y,f)    !Insert the present value of x and y into the f. 
      
      do j=1,n
            s2(j)=h*f(j)
c------------------------------------------------------------------   

            y(j)=ysav(j)       !Recover the original value of y(j) to use it in the next step.
      enddo

      x=xsav+0.5d0*h           !This sentence isn't necessary. It is inserted only to confirm the value of x.

c-----Calculation of s3---------------------------------------------      
      do j=1,n
            y(j)=ysav(j)+0.5d0*s2(j)
      enddo

      call EQATION(x,y,f)

      do j=1,n
            s3(j)=h*f(j)
c-----Calculation of s4-------------------------------------------------    
            y(j)=ysav(j)
      enddo
         x=xsav+h               !This sentence is necessary.
      do j=1,n     
            y(j)=ysav(j)+s3(j)
      enddo
      call EQATION(x,y,f)
      do j=1,n
            s4(j)=h*f(j)
c---------------------------------------------------------------------------
            y(j)=ysav(j)
      enddo

c-----Final procesure of constructing y(x+h) from y(x).----------------------      
      do j=1,n
            y(j)=ysav(j)+(s1(j)+2.0d0*s2(j)+2.0d0*s3(j)+s4(j))/6.0d0
      enddo
      return
      end
c------------------------------------------------------------------------------

c=====Definition of functions used in the above program.==================================================

c-----lam1 represents the function 'lambda' of triple alpha reaction.------------------------------------
      real*8 function lam1(rho,T9)
      implicit none
      real*8 rho,T9
      
      lam1 = 1.3d-56/T9**3*dexp(-4.4027d0/T9)

      return
      end
c----------------------------------------------------------------------------

c-----lam2 represents the function 'lambda' of the reaction of alpha and 12C----------------
      real*8 function lam2(rho,T9)
      implicit none
      real*8 rho,T9

      lam2 = 1.7d-16*T9**(-2.0d0/3.0d0)*
c     *       (1.0d0+0.0489d0*T9**(-2.0d0/3.0d0))**2*
     *       dexp(-32.12d0*T9**(-1.0d0/3.0d0)-(0.286d0*T9)**2)

      return 
      end
c----------------------------------------------------------------------------------------

c-----lam3 represents the function 'lambda' of the reaction of alpha and 16O----------------
      real*8 function lam3(rho,T9)
      implicit none
      real*8 rho,T9
      
      lam3 = 1.6d-14*T9**(-2.0d0/3.0d0)*
     *       dexp(-39.757d0*T9**(-1.0d0/3.0d0)-(0.631d0*T9)**2)

      return
      end
c----------------------------------------------------------------------------------------------
