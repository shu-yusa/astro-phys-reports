!====================================================================!
!     シンプソンの公式によりExponential Integral--E_n(x)と呼ばれる   !
!     関数を数値的に求めるプログラム。E_n(x)は積分で定義されている   !
!     ため、各xに対して積分をすることで値が決まる。                  !
!     以下のプログラムでは区間[x0,x1]におけるE_n(x)の値を求める。    !
!     nはExponential Integralの次数を表していて,今回はn=1(2)からn=kk !
!     まで求める。                                                   !
!     Sympsonのsubrotineが二つあるがsubroutine Sympsonの方は         !
!     積分区間が固定,subrotine Sympson2の方は積分区間を被積分関数の  !
!     値により変わるようにしている。                                 !
!     後者のほうが効率がよく精度もよい.                              !
!     このプログラムの計算は時間がかかる。やや精度を犠牲にして短縮   !
!     するには                                                       !
!     1.subroutine Sympson2 内のepsを1.0d-6~7ぐらいにする。          !
!     2.main program のparameter内のnをn=3000~1500ぐらいにする。     !
!     などをすればよい。                                             !
!     またこのプログラムは恒星物理学のレポートも兼ねている。         !
!====================================================================!
      program  Exponential_Integral
!--------------------------------------------------------------------!
!     a,b----シンプソンの公式を用いて積分する区間[a,b] for subrotine !
!            Simpson                                                 !
!     c------subroutine Simpson2が用いる積分区間[a,c]                !
!            この積分区間は被積分関数の値で決まる。こちらを用いると  !
!            効率よく計算できる。                                    !
!     m------シンプソンの公式で積分する区間の分割数。                !
!     x------Exponential Integralの引数                              !
!     S(k,i)------k次のExponential Integralのx(i)における値          !
!     FNAME(k)----k次のExponential Integralが入る。組(x(i),S(k,i))   !
!     x1-----Exponential Integralの大きいほうの端点。                !
!     n------Exponential Integralの方の区間の分割数。                !
!     h------Exponential Integralのサンプリング幅。x1とnから決まる。 !
!     T------恒星物理学Ⅰのレポート用。フラックスのエラー             !
!     FLUX---恒星物理学Ⅰのレポート用。                               !
!--------------------------------------------------------------------!
      implicit none
      integer i,k,kk,n,m
      real*8 a,h,x(10000),S(4,10000) 
      real*8 x0,x1,c,T(10000)
      parameter(a=1.0d0,m=5000,n=6000,x0=0.0d0,x1=3.0d0,kk=6)
      character FNAME(6)*20,FLUX*20
      data FNAME /'Exp_Int1.dat','Exp_Int2.dat','Exp_Int3.dat',  &
                  'Exp_Int4.dat','Exp_Int5.dat','Exp_Int6.dat'/
      data FLUX /'Error_of_flux.dat'/
      
      h=x1/real(n)           !Exponential Integralのサンプリング幅      
      x(1)=x0                !Exponential Integralの小さい方の端点

! 種々のパラメータの画面への書き出し
      
      write(6,100) '積分値を求める区間 : [',x0,',',x1,']'
      write(6,200) '区間の刻み幅 =',h  
      write(6,300)   '計算するExponential Integralの数 :',kk-1
!
      do k=2,kk
         open(unit=k*10,err=500,status='old',  &
              file=FNAME(k),form='formatted')   
      end do
      
      
      do k=2,kk            !2次から4次までのExponential Integralを計算
         do i=1,n          !xの値を色々サンプリング
            x(i+1)=x(i)+h                     !xの集合も作る
            call Simpson2(a,c,m,x,k,i,S)      !シンプソンの公式で計算
            write(k*10,*) x(i),S(k,i)         !ファイルへの書き込み
         enddo
         write(6,400) k-1,'つ目の積分が終わりました。' 
      end do
      
!  恒星物理学Ⅰのレポート用の計算 

      open(unit=90,err=500,status='old',file=FLUX,form='formatted') 
      
      do i=1,n
        T(i)=S(3,i)-1.5d0*S(4,i)
        write(90,*) x(i),T(i)
      end do
      close(90)
!      
      do k=2,kk
         close(k*10)         !ファイルを閉じる。
      end do
      
 100  format(' ',a,f5.2,a,f5.2,a)
 200  format(' ',a,f7.4)
 300  format(' ',a,i2)
 400  format(' ',i1,a) 
      goto 600
 500  write(6,*) 'fileにエラー'
 600  stop
      end
      
!====================================================================!
      subroutine Simpson(a,b,n,x,k,j,S)
!--------------------------------------------------------------------!
!     シンプソンの公式により積分を行う副プログラム                   !
!     a,b-----積分区間[a,b]。入力変数。bはすべての(k,x)に共通の定数  !
!             としている。                                           !
!     n-------積分区間[a,b]を分割する数。入力変数。                  !
!     k,x,j---被積分関数に関する変数。被積分関数f(k,t,x)はkとx(j)を  !
!             与えると一つのtの関数として決まる。入力変数。          !
!     S-------積分結果。出力データ。                                 !
!     S1,S2---作業変数。                                             !
!     h-------分割間隔。この副プログラムでは一定間隔にしている。     !
!--------------------------------------------------------------------!
      real*8 a,b,S(4,10000),x(10000),h,S1,S2
      integer i,j,n
      real*8 f
      external f
      
      h=(b-a)/dble(n)         !積分区間の分割の幅(=0.01)
      
      S1=0.0d0     
      do i=2,n-2,2
        S1=S1+f(k,a+i*h,x(j))
      end do
      
      S2=0.0
      do i=1,n-1,2
        S2=S2+f(k,a+i*h,x(j))
      end do
      
      S(k,j)=( f(k,a,x(j)) + 2.0*S1 + 4.0*S2 + f(k,b,x(j)) )*h/3.0d0
      
      return
      end
   
!====================================================================!
      real*8 function f(k,t,x)
!--------------------------------------------------------------------!
!     Exponential Integralの被積分関数の定義                         !
!     k次のExponential Integral,E_k(x)の被積分関数f(k,t,x)の定義     !
!     tは積分変数を表す。f=t^(-k)*e^(-xt)                            !
!--------------------------------------------------------------------!
      implicit none
      real*8 t,x
      integer k
      
      f=1.0d0/(t**k)*exp(-x*t)
      
      return
      end
      
!====================================================================!
      subroutine Simpson2(a,b,m,x,k,j,S)
!--------------------------------------------------------------------!
!     シンプソンの公式を用いて積分をする副プログラム                 !
!     subroutine Simpsonと違い、積分範囲を(k,x)に依存して決めている。!
!     すなわち被積分関数の値がepsの値より小さくなるまで積分を行う。  !
!     プログラムではまず積分範囲の上限をループを用いて決めている。   !
!--------------------------------------------------------------------!
      real*8 a,b,S(4,10000),x(10000),h,S1,S2
      integer i,j,n,m
      real*8 f,eps,c
      external f
      parameter(eps=1.0d-10,h=1.0d-2)

      b=a                                ! b=aから適当な上限を探してゆく
       
      do                                 ! 積分の上限を決めるためのループ
         if(f(k,b,x(j)).lt.eps) then     ! 被積分関数が十分小さいか
           exit                          ! どうかの判定
         else 
           b=b+h                         ! 上限の値を少しずつ大きくしていく
         endif
      end do
      
      n=int((b-a)/h)                     ! 積分範囲が決まったので分割間隔を
                                         ! 一定にするように分割数を決定
      S1=0.0d0     
      do i=2,n-2,2
        S1=S1+f(k,a+i*h,x(j))
      end do
      
      S2=0.0d0
      do i=1,n-1,2
        S2=S2+f(k,a+i*h,x(j))
      end do
      
      S(k,j)=( f(k,a,x(j)) + 2.0d0*S1 + 4.0d0*S2 +  &
               f(k,b,x(j)) )*h/3.0d0
      n=m

      return
      end
     
