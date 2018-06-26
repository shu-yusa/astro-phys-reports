!====================================================================!
!     Runge-Kutta法によりLane-Emden方程式を解くプログラム。          !
!     プログラムを実行することにより5個のファイルができる。          !
!     ただしポリトロープ指数は整数のものに限っている。               !
!====================================================================!
      program main
!--------------------------------------------------------------------!
!     FNAME(j)----積分結果(x,y_j)を入れるファイル名                  !
!     a-----------ファイルを作るかどうかを決めるため                 !
!     x0,x1-------積分区間[x0,x1]                                    !
!     s(j,i)------y_j(i) or y_j(x(i)) のこと                         !
!     初期値を与える点はx0 (x2<x0<x1)                                !
!     nは方程式の数                                                  !
!     kkは解くLane-Emden方程式の最大ポリトロープ指数(<=4)            !
!     kはポリトロープ指数として使う。                                !
!     dはサンプリングの間隔                                          !
!     nfileは出力したい(x,y(i))の数(つくるファイルの数)              !
!     niは積分の計算のときの分割数をd/niとする。                     !
!====================================================================!
      implicit none
      character :: FNAME(0:4)*20, a
      real*8    :: s(0:10,4,100000),x(100000)
      real*8    :: x0,x1(0:10),d
      integer   :: n,ni,i,k,kk,l
      parameter(kk=4)
      data n,x0,d,ni/2, 1.0d-6, 2.0d-4, 10/  !定義域は0以上のところ
      data FNAME/'Lane_Emden0.dat','Lane_Emden1.dat', &
                 'Lane_Emden2.dat','Lane_Emden3.dat', &
                 'Lane_Emden4.dat'/
!-----初期値および条件の書き出し------
      write(*,*) ' *** Numerical Analysis of ODE ***'
      write(*,*) '      (With Runge-Kutta Method)'
      write(*,*)
      write(*,200) n,x0,d,ni
      l=1
      x(1)=x0
      call INITIAL(kk,s,x0)
      
      write(*,*) 'Initial condition'
      do k=0,kk
      write(*,210) 'n=',k,x(1),s(k,1,1)
      end do
      write(*,*)
!-----Solve differential equations vy Runge-Kutta method----
      do k=0,kk                       !Polytrope index from n=0 to n=kk.
        call DRUNGE(k,x,n,s,d,ni,x1)
        write(6,240) 'n=',k,' th node is:',x1(k)   !Display first nodes.
      end do
!-------------------------------------------
      write(6,*) 
      write(6,*) 'Do you make data files?'
      write(6,*) 'If you want to, type "y"'
      write(6,*) 'If you do not want to, type "n"'
      read(5,*) a
      if (a.eq.'n') goto 250
!----常微分方程式の解のファイルへの書き出し--------------------------------
      
      do k=0,kk                      ! ファイルを開く。
        open(k+11,file=FNAME(k))     ! k+11としているのは装置番号に
                                     ! 5,6がこないようにするため。          
      end do
      
      do k=0,kk
        i=1
 20   write(k+11,220) x(i),s(k,1,i)    ! ファイルに書き込む。
        if(s(k,1,i).ge.0.0d0) then     ! 方程式の解は零点を通る
                                       ! ところまで必要
            i=i+1
            goto 20
        end if
      end do
     
      do k=0,kk                        !ファイルを閉じる。
      close(k+11)
      end do
      
!-----format----------------------------------------------------------
 200   format(1h ,'number of equation',i2/' upper limit=',f5.1 &
      /' sampling interval=',f7.4  &
      /' division number=',i2)
 210   format(1h ,a,i1,2x,'x=',f7.4,4x,'y=(',f8.5,')')
 220   format(f8.5,' ',f8.5)
 240   format(' ',a,i1,a,f9.5) 
 250  stop
      end
!====================================================================!
      subroutine INITIAL(kk,s,x0)
!--------------------------------------------------------------------!
!     解く方程式は1/xの因子が入っているためx=0からは積分できない。   !
!     そこで原点から少しだけずらした点から積分する。                 !
!     原点付近の展開式を用いて初期値を代入する。                     !
!--------------------------------------------------------------------!
      real*8 s(0:10,4,100000),x0
      integer k,kk
      do k=0,kk
      s(k,1,1)=(k/1.2d2*x0*x0-1.0d0/6.0d0)*x0*x0+1.0d0
      s(k,2,1)=(k/3.0d1*x0*x0-1.0d0/3.0d0)*x0
      s(k,3,1)=1.0d0
      s(k,4,1)=1.0d0
      end do
      return
      end
!-----differential equations --------
      subroutine EQATION(k,x,y,f)
      real*8 x,y(4),f(4)
      integer k
      f(1)=y(2)
      f(2)=-y(1)**k-2.0d0/x*y(2)
      f(3)=0.0d0
      f(4)=0.0d0
      return
      end
!-----Driver routines of Runge-Kutta-Method-----
!===================================================================
!     x(i)は独立変数x。s(j,i)は求めたい関数y_j
!     kはポリトロープ指数。
!===================================================================
      subroutine DRUNGE(k,x,n,s,d,ni,x1)     !出力はx,s,l
      real*8 s(0:10,4,100000),y(4)
      real*8 x(200000),h,d,x1(0:10)
      integer ni,k,l
      l=1
      h=d/real(ni)
      y(1)=s(k,1,1)
      y(2)=s(k,2,1)
      y(3)=s(k,3,1)
      y(4)=s(k,4,1)
 10   x(l+1)=x(l)
        do i=1,ni
           call RUNGE(k,h,n,x(l+1),y)   !x(l+1)におけるy_j,つまりs(j,l+1)になる量を計算。
        end do
        l=l+1
        s(k,1,l)=y(1)
        s(k,2,l)=y(2)
        s(k,3,l)=y(3)
        s(k,4,l)=y(4)
        if(s(k,1,l).le.0.0d0)  goto 30  !解がゼロ点まで行ったら終わり。
      goto 10
 30   x1(k)=x(l)       !最初のゼロ点(node)をチェック
      return
      end
!------------------------------------------------
!-----Runge-Kutta-Method-------------------------
      subroutine RUNGE(k,h,n,x,y)                         !出力はx,y
      real*8 f(4),y(4),s1(4),s2(4),s3(4),s4(4),ysav(4)
      real*8 x,xsav,h
      integer k
      xsav=x
      do 10 j=1,n
            ysav(j)=y(j)
 10   continue
      call EQATION(k,x,y,f)
      do 20 j=1,n
            s1(j)=h*f(j)
 20   continue
      x=xsav+0.5d0*h
      do 30 j=1,n
            y(j)=ysav(j)+0.5d0*s1(j)
 30   continue
      call EQATION(k,x,y,f)
      do 40 j=1,n
            s2(j)=h*f(j)
            y(j)=ysav(j)
 40   continue
      x=xsav+0.5d0*h
      do 50 j=1,n
            y(j)=ysav(j)+0.5d0*s2(j)
 50   continue
      call EQATION(k,x,y,f)
      do 60 j=1,n
            s3(j)=h*f(j)
            y(j)=ysav(j)
 60   continue
      x=xsav+h
      do 70 j=1,n
            y(j)=ysav(j)+s3(j)
 70   continue
      call EQATION(k,x,y,f)
      do 80 j=1,n
            s4(j)=h*f(j)
            y(j)=ysav(j)
 80   continue
      do 90 j=1,n
            y(j)=ysav(j)+(s1(j)+2.0d0*s2(j)+2.0d0*s3(j)+s4(j))/6.0d0
 90   continue
      return
      end
!------------------------------------------------
