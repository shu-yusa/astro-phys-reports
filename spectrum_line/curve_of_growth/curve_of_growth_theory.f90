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
      real*8 a,b,c,hh1,hh2,hh3,A0,h1,h2,s
      real*8 v0,v11,v12,v13
      real*8 v(19000),H(19000),g(19000,200)
      real*8 R,Br(3),beta(5,19000),beta0(19000)
      real*8 S1,S2,W(200)
      external R
      integer i,j,k,n1,n2,n3,m,nn
      parameter(A0=0.5d0,c=1.0d-3,v0=0.0d0)
      parameter(v11=5.0d1,n1=5000)
      parameter(v12=1.0d2,n2=10000)
      parameter(v13=1.9d2,n3=19000)
      parameter(nn=160)
      character FNAME(1)*40
      data FNAME /'Voigt.dat'/
      data Br /1.0d0,2.0d0,3.0d0/

      hh1=(v11-v0)/dble(n1)
      hh2=(v12-v0)/dble(n2)
      hh3=(v13-v0)/dble(n3)
      v(1)=v0

      open(unit=10, file=FNAME(1))
      
      h1=5.0d-2

      do i=1,nn+1
        beta0(i) = 10.0d0**(-2.0d0+h1*(i-1))
!      write(6,*) dlog10(beta0(i))
      enddo

      do j=1,100
       do i=1,n1
           v(i+1)=v(i)+hh1
           call Simpson2(a,b,m,v,c,i,H)
           g(i,j) = 2*A0*beta0(j)*H(i)/(1+beta0(j)*H(i))
!          write(10,*) v(i),g(i,j)
       enddo
        write(6,*) j,'回目終了'
      
      S1=0.0d0
      S2=0.0d0
      
      do k=2,n1-2,2
        S1=S1+g(k,j)
      enddo
       
      do k=1,n1-1,2
        S2=S2+g(k,j)
      enddo
 
      W(j)=(g(1,j)+2.0d0*S1+4.0d0*S2+g(n1,j))*hh1/3.0d0
        write(10,*)  dlog10(beta0(j)),dlog10(W(j))
        write(6,*) 'log W =' ,dlog10(W(j))
      enddo

      do j=101,120
       do i=1,n2
         v(i+1)=v(i)+hh2
         call Simpson2(a,b,m,v,c,i,H)
         g(i,j) = 2*A0*beta0(j)*H(i)/(1+beta0(j)*H(i))
       enddo
        write(6,*) j,'回目終了'

       S1=0.0d0
       S2=0.0d0
       
       do i=1,n2-2,2
         S1=S1+g(i,j)
       enddo

       do i=1,n2-1,2
         S2=S2+g(i,j)
       enddo

       W(j)=(g(1,j)+2.0d0*S1+4.0d0*S2+g(n2,j))*hh2/3.0d0
         write(10,*) dlog10(beta0(j)),dlog10(W(j))
         write(6,*) 'log W =' ,dlog10(W(j))
       enddo 

      do j=121,161
         do i=1,n3
         v(i+1)=v(i)+hh3
         call Simpson2(a,b,m,v,c,i,H)
         g(i,j) = 2*A0*beta0(j)*H(i)/(1+beta0(j)*H(i))
         enddo
         write(6,*) j,'回目終了'

         S1=0.0d0
         S2=0.0d0

         do i=1,n3-2,2
            S1=S1+g(i,j)
         enddo
  
          do i=1,n3-1,2
            S2=S2+g(i,j)
          enddo

            W(j)=(g(1,j)+2.0d0*S1+4.0d0*S2+g(n3,j))*hh3/3.0d0
            write(10,*) dlog10(beta0(j)),dlog10(W(j))
      enddo

      close(10)
      write(6,*) 'Completed!'
 
      goto 101 
 99   write(6,*) 'An error in status of file.'
 101  continue     
      stop
      end

!====================================================================!
      real*8 function f(a,v,y)
!--------------------------------------------------------------------!
!     Integrand of Voigt-function.                                   !
!--------------------------------------------------------------------!
      implicit none
      real*8 a,v,y,pi
      parameter(pi=3.141592653)
      
      f=a/pi*dexp(-y**2)/((v-y)**2+a**2)
      
      return
      end

!====================================================================!
      real*8 function R(beta,Br)
!--------------------------------------------------------------------!
!     Residual flux                                                  !
!--------------------------------------------------------------------!
      implicit none
      real*8 beta,Br,S

      S=dsqrt(3.0d0/(1.0d0+beta))      
      R=(S+Br/(1.0d0+beta))/(S+1.5d0)*(dsqrt(3.0d0)+1.5d0)   &
           /(dsqrt(3.0d0)+Br)
      return
      end

!====================================================================!
      subroutine Simpson2(a,b,m,v,c,j,H)
!--------------------------------------------------------------------!
!     'a' and 'b' are the lower limit and upper limit of integration !
!     respectively.                                                  !
!     'c' is the first argument of f(c,v,y).                         !
!--------------------------------------------------------------------!
      implicit none
      real*8 a,b,c,H(10000)
      real*8 v(10000),hh,S1,S2,vv
      integer i,j,n,m
      real*8 f,eps
      external f
      parameter(eps=1.0d-8,hh=1.0d-4) 
      
      a=0.0d0
      b=0.0d0                

!    Determine the upper limit of integration.
      do 
        if(f(c,v(j),b).lt.eps) then  
          exit
        else
          b=b+hh 
        endif
      end do
      
!    Determine the lower limit of integration.
      do
        if(f(c,v(j),a).lt.eps) then 
          exit
        else 
          a=a-hh
        endif 
      end do

      n=int((b-a)/hh)  
      b=a+n*hh                ! Redifine 'b' so that 'a' and 'b' 
                              ! to be consistent with n.
      S1=0.0d0     
      do i=2,n-2,2
        S1=S1+f(c,v(j),a+i*hh)
      end do
     
      S2=0.0d0
      do i=1,n-1,2
        S2=S2+f(c,v(j),a+i*hh)
      end do
      
      H(j)=( f(c,v(j),a) + 2.0d0*S1 + 4.0d0*S2 +    &
             f(c,v(j),b) )*hh/3.0d0
      m=n
      return
      end
     
