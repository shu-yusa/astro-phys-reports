!======================================================================!
!     Title  : Exp_Int.f90                                             !
!     Author : Yusa Shusaku                                            !
!     Date   : ???                                                     !
!     Last modified : 2008-8-11-Mon                                    !
!                                                                      !
!     シンプソンの公式によりExponential Integral--E_n(x)と呼ばれる     !
!     関数を数値的に求めるプログラム。E_n(x)は積分で定義されている     !
!     ため、各xに対して積分をすることで値が決まる。                    !
!     以下のプログラムでは区間[x0,x1]におけるE_n(x)の値を求める。      !
!     nはExponential Integralの次数を表していて,今回はn=1(2)からn=kk   !
!     まで求める。                                                     !
!     Sympsonのsubrotineが二つあるがsubroutine Sympsonの方は           !
!     積分区間が固定,subrotine Sympson2の方は積分区間を被積分関数の    !
!     値により変わるようにしている。                                   !
!     後者のほうが効率がよく精度もよい.                                !
!     このプログラムの計算は時間がかかる。やや精度を犠牲にして短縮     !
!     するには                                                         !
!     1.subroutine Sympson2 内のepsを1.0d-6~7ぐらいにする。            !
!     2.main program のparameter内のnをn=3000~1500ぐらいにする。       !
!     などをすればよい。                                               !
!     またこのプログラムは恒星物理学のレポートも兼ねている。           !
!======================================================================!
      module  Com_var_1
      implicit none
      integer, parameter :: Maxn=10000, kmax=6
      real*8, parameter :: x0=0.0d0, x1=3.0d0, dx=0.01d0
      real*8, parameter :: a=1.0d0
      end module
!======================================================================!
      program  Exponential_Integral
!----------------------------------------------------------------------!
!     a,b----シンプソンの公式を用いて積分する区間[a,b] for subrotine   !
!            Simpson                                                   !
!     c------subroutine Simpson2が用いる積分区間[a,c]                  !
!            この積分区間は被積分関数の値で決まる。こちらを用いると    !
!            効率よく計算できる。                                      !
!     m------シンプソンの公式で積分する区間の分割数。                  !
!     x------Exponential Integralの引数                                !
!     S(k,i)------k次のExponential Integralのx(i)における値            !
!     FNAME(k)----k次のExponential Integralが入る。組(x(i),S(k,i))     !
!     x1-----Exponential Integralの大きいほうの端点。                  !
!     n------Exponential Integralの方の区間の分割数。                  !
!     h------Exponential Integralのサンプリング幅。x1とnから決まる。   !
!     T------恒星物理学Ⅰのレポート用。フラックスのエラー               !
!     FLUX---恒星物理学Ⅰのレポート用。                                 !
!----------------------------------------------------------------------!
      use Com_var_1, only : Maxn, x0, x1, kmax, dx
      implicit none
      integer :: i, k, n
      real*8 :: x(Maxn)
      real*8, allocatable :: S(:,:)
      character(len=20) :: FNAME(6), FLUX
      character(len=22), parameter :: FM1='(a,f5.1,a,f5.1,a)'
      character(len=20), parameter :: FM2='(a,f8.5)'
      data FNAME /'Exp_Int1.dat','Exp_Int2.dat','Exp_Int3.dat',  &
     &            'Exp_Int4.dat','Exp_Int5.dat','Exp_Int6.dat'/
      data FLUX /'Error_of_flux.dat'/
      
      write(6,FM1) 'Range      : [',x0,',',x1,']'
      write(6,FM2) 'Step width :', dx
      
      x(1) = x0
      n = 1
      do 
        x(n+1) = x(n) + dx
        if (x(n+1) > x1) then 
           exit
        end if
        n = n + 1
        if (n == Maxn) stop 'Too small array size.'
      end do
      allocate(S(2:kmax, n))
!$OMP parallel do private(i) num_threads(5) 
      do k=2, kmax
         open(unit=k+6, file=FNAME(k))   
         do i=1, n
            call GausLagu(10, k, x(i), S(k,i))
            write(k+6,*) x(i), S(k,i)
         enddo
         write(6,*) k-1,'th integration has been done.' 
         close(k+6)
      end do
!$OMP end parallel do
!      S = f(0.0d0,m) * SQRTPI
!
!      do n=2, MAXny
!          call Draw_Line
!          S0 = S
!          call GausLagu(n, m, x, S)
!          DS = abs(S - S0)
!          
!          if(DS < epsa + epsr * (abs(S) + abs(S0))) exit
!          
!          if(n == MAXn) then
!              write(*,*) 'Integration did not converge.'
!              exit
!          endif
!      enddo
      
! Computation for Report
      open(unit=kmax+7,file=FLUX) 
      do i=1, n
         write(kmax+7,*) x(i), S(3,i) - 1.5d0 * S(4,i)
      end do
      close(kmax+7)
      
      stop
      end program
!====================================================================!
      function f(k,t,x)  result(g)
!--------------------------------------------------------------------!
!     Definition of integrand of Exponential Integral of order 'k'.  !
!     't' is an integration variable.                                !
!--------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: k
      real*8, intent(in)  :: t, x
      real*8 :: g
      
      g = exp(t - x * (t + 1.0d0)) * (t + 1.0d0) ** (-k)
      
      return
      end function
!====================================================================!
      module Com_var
!--------------------------------------------------------------------!
!     Definition of parameters common to several program units.      !
!--------------------------------------------------------------------!
        integer :: MAXn
        real*8  :: epsa, epsr, PI, SQRTPI
        parameter(MAXn=150)
        parameter(epsr=1.0d-15, epsa=1.0d-200)
        parameter(PI=3.1415926535897932d0)
        parameter(SQRTPI = sqrt(PI))
      end module 
!====================================================================!
      subroutine GausLagu(n, m, x, S)
!--------------------------------------------------------------------!
!     Calculation of integration by using Gauss-Laguerre formula.    !
!     If logical variable 'useDKA' is 'T' and 'n' is smaller than 12,!
!     we use DKA method for 'n < 12' and bisection method for        !
!     'n > 13'.                                                      !
!     If 'useDKA' is 'F', we use only bisection method throughout    !
!     this program.                                                  !
!     Input  : n                                                     !
!     Output : S                                                     !
!--------------------------------------------------------------------!
      use Com_var, only : MAXn
      implicit none
      integer, intent(in) :: n, m
      real*8, intent(in) :: x
      real*8, intent(out) :: S
      integer :: j
      real*8  :: xj, wj, t(Maxn)
      real*8, external :: f
      character(len=35) :: f100,f200,f300
      parameter(f100='("| ",f19.15," | ",1pd22.15," |")')
      parameter(f200='("|",3x,a,i2,12x," | ",23x,"|")')
      parameter(f300='("| ",5x,a,6x,"| ",7x,a,9x,"|")')

!      write(6,f200) 'n =', n
!      write(6,f300) 'Abscissas','Weights'
      
      call zeros(n, t)

      S = 0.0d0
      do j=n, 1, -1
          call weight(n, t(j), wj)
!          write(6,f100) t(j), wj
          S = S + wj * f(m,t(j),x)
      end do

      return
      end subroutine 
!====================================================================!
      subroutine  Pre_Bisec(n, x, xLEF, xRI)
!--------------------------------------------------------------------!
!     A subroutine which seeks for startig points of bisection       !
!     method.                                                        !
!     The zeros of Laguerre polynomial of order 'n' are such that    !
!         0 <= x(i) <= 4*n - 3                                       !
!     Therefore we seek for the starting points in this reagion.     !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in)    :: n
      real*8, intent(inout) :: x
      real*8, intent(out)   :: xLEF, xRI 
      real*8  :: dL, dx=0.05d0
      real*8,external :: L

      do
          if (L(n, x, dL) * L(n, x + dx, dL) < 0.0d0) then
              exit
          else 
              x = x + dx
          end if

          if (x > dble(4 * n - 3)) then
              write(6,*) 'ERROR'
              stop
          end if
      end do

      xLEF = x
      xRI  = x + dx
     
      return
      end  subroutine
!====================================================================!
      subroutine  Bisection(n, x)
!--------------------------------------------------------------------!
!     A subroutine that seeks for zeros of Hermite polynomial by     !
!     using bisection method. Startig points of this method are      !
!     computed from subroutine 'pre_bisec'.                          !
!--------------------------------------------------------------------!
      use Com_var, only : epsr, epsa
      implicit none
      integer,intent(in)  :: n
      real*8, intent(out) :: x
      integer :: i
      real*8  :: x0, xLEF, xRI
      real*8  :: sgn, dL
      real*8,external :: L

      call Pre_Bisec(n, x, xLEF, xRI)
      
      x = 0.0d0
      do
          x0 = x
          x  = 0.5d0 * (xLEF + xRI)

          !if (abs(x - x0) < epsr * (abs(x) + abs(x0))) exit
          if (abs(x - x0) == 0.0d0) exit
          
          sgn = sign(1.0d0, L(n,xLEF,dL)) * sign(1.0d0, L(n,x,dL))
          if (sgn > 0.0d0) then
              xLEF = x
          else
              xRI  = x
          end if
      end do
      
      return
      end subroutine
!====================================================================!
      subroutine  zeros(n, x)
!--------------------------------------------------------------------!
!     Caluculation of zeros of Laguerre polynomial by calling        !
!     subroutine 'bisection'. The array of zeros 'x(i)' are already  !
!     arranged after finding all zeros.                              !
!     Since we seek for zeros succssesively, 'y' is the zero of      !
!     Laguerre polynomial after each end of do-loop. Therfore, we    !
!     must shift 'y' slightly forward. Otherwise, the bisection      !
!     subroutine will not work well.                                 !
!--------------------------------------------------------------------!
      use Com_var, only : MAXn
      implicit none
      integer,intent(in) :: n
      real*8,intent(out) :: x(MAXn)
      integer :: i
      real*8  :: y, dy=0.01d0

      y = 0.0d0
      do i=n, 1, -1
          y = y + dy
          call Bisection(n, y)
          x(i) = y
      end do

      return
      end subroutine
!====================================================================!
      subroutine  weight(n, xj, wj)
!--------------------------------------------------------------------!
!     Calculation of a weight factor 'wj' from a zero point 'xj'.    !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in)  :: n
      real*8, intent(in)  :: xj
      real*8, intent(out) :: wj
      integer :: k
      real*8  :: dL
      real*8, external :: L
      
      wj = xj / (dble(n) * L(n-1,xj,dL)) ** 2

      return
      end subroutine
!====================================================================!
      real*8 function L(n, x, dL)
!--------------------------------------------------------------------!
!     Definition of Laguerre polynomial.                             !
!     To obtain the polynomial of order n, we use the recursion      !
!     relation.                                                      !
!     Derivative function of the Laguerre polynomial is also         !
!     calculated.                                                    !
!--------------------------------------------------------------------!
      implicit none
      integer,intent(in) :: n
      real*8,intent(in)  :: x
      real*8  :: L0, L1, dL
      integer :: j

      L1 = 1.0d0
      L  = - x + 1.0d0
      
      do j=2, n
          L0 = L1
          L1 = L
          L = (- (x - 2.0d0 * dble(j) + 1.0d0) * L1 &
            &  - (dble(j) - 1.0d0) * L0) / dble(j)
      enddo
       
      dL = dble(n) * (L - L1) / x

      return
      end  function
!====================================================================!
      subroutine pict
!--------------------------------------------------------------------!
!     This subroutine draws a picture...                             !
!--------------------------------------------------------------------!
      call Draw_Line
      write(6,*)
      write(6,*) '|*****************************************|'
      write(6,*) '| The function we calculated.             |'
      write(6,*) '|                                         |'
      write(6,*) '|     inf                                 |'
      write(6,*) '|     /                                   |'
      write(6,*) '|     ]             x^m                   |'
      write(6,*) '| S = | dx exp(-x) ------ =  1            |'
      write(6,*) '|     ]              m!                   |'
      write(6,*) '|     /                                   |'
      write(6,*) '|     0                                   |'
      write(6,*) '|                                         |'
      write(6,*) '*******************************************'

      return 
      end  subroutine
!====================================================================!
      subroutine  Draw_Line
      
      write(6,*) '-----------------------&
                 &-----------------------'

      return 
      end subroutine
!====================================================================!
      subroutine  Draw_Title
      implicit none
      character :: fm*10
      parameter(fm='(8x,a)')

      write(6,*) 
      write(6,fm) '**************************'
      write(6,fm) '  Gauss-Laguerre formula  '
      write(6,fm) '**************************'

      return 
      end subroutine
!====================================================================!
      subroutine  Draw_Result(m, n, S, DS)
      implicit none 
      integer, intent(in) :: n, m
      real*8, intent(in) :: S, DS
      character :: f100*20, f50*15, f10*20
      parameter(f10 ='(1x,a,2x,i3)')
      parameter(f50 ='(1x,a,f20.15)')
      parameter(f100='(1x,a,2x,1pd22.15)')

      write(6,*) 
      write(6,f10)  'm     =',m
      write(6,f10)  'n     =',n
      write(6,f50)  'S     =',S
      write(6,f100) 'Error =',DS
      write(6,*)  

      return
      end subroutine  
