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
      module  Com_var
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
      use Com_var, only : Maxn, x0, x1, kmax, dx
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
      
      open(unit=7, file=FNAME(2))   
      do i=1, n
         call Simpson2(x(i), 2, S(2,i))
         write(7,*) x(i), S(2,i)
      enddo
      write(6,*) 1,'th integration has been done.' 
      close(7)

!$OMP parallel do private(i,k) num_threads(4) 
      do k=3, kmax
         open(unit=k+7, file=FNAME(k))   
         do i=1, n
            call Simpson2(x(i), k, S(k,i))
            write(k+7,*) x(i), S(k,i)
         enddo
         write(6,*) k-1,'th integration has been done.' 
         close(k+7)
      end do
!$OMP end parallel do
      
! Computation for Report
      open(7,file=FLUX) 
      do i=1, n
         write(7,*) x(i), S(3,i) - 1.5d0 * S(4,i)
      end do
      close(7)
      
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
      
      g= exp(-x*t) * t ** (-k)
      
      return
      end function
!====================================================================!
      subroutine Simpson2(x, k, S)
!--------------------------------------------------------------------!
!     A subroutine which executes integration by using Simpson's     !
!     formula. In this suroutine, we also compute the range of       !
!     integration.                                                   !
!--------------------------------------------------------------------!
      use Com_var, only : a
      implicit none
      integer :: i, n
      integer, intent(in) :: k
      real*8, intent(out) :: S
      real*8, intent(in) :: x
      real*8, external :: f
      real*8, parameter :: h=0.01d0, eps=1.0d-10
      real*8 :: S1, S2, b

      b = a
      do
         if(f(k,b,x) < eps) then
           exit
         else 
           b = b + h
         endif
      end do
      
      n = int( (b - a) / h )

      S1=0.0d0     
      do i=2, n-2, 2
        S1 = S1 + f(k,a+dble(i)*h,x)
      end do
      
      S2=0.0d0
      do i=1, n-1, 2
        S2 = S2 + f(k,a+dble(i)*h,x)
      end do
      
      S = ( f(k,a,x) + 2.0d0*S1 + 4.0d0*S2 + f(k,b,x) ) * h / 3.0d0

      return
      end subroutine
     
