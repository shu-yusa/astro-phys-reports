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
      module global_const
      implicit none
!$    integer :: thread = 8
      real(8), parameter :: A0 = 0.5d0
      real(8), parameter :: c = 1.0d-3
      end module
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
      integer :: j
      real(8), external :: kai2
      real(8), parameter :: A0 = 0.5d0
      real(8), parameter :: c = 1.0d-3
      real(8), dimension(85) :: x, y, beta0, fac, W
      real(8) :: pms(9)
      
      pms(1) = 3.85d0
      pms(2) = 5.35d0
      pms(3) = - 2.0d0
      pms(4) = 5.0d0
      pms(5) = 5.1d0
      pms(6) = 5.35d0
      pms(7) = 0.8d0
      pms(8) = 5.0d0
      pms(9) = 0.92d0
      call read_data(pms, x, y)
      call calc_th(85, x, W)

      call OROF_method(85, 9, pms)
      write(6,*) kai2(85, y, W)

      open(unit=7, file='theory.dat')
      do j=1, 85
        write(7,*)  x(j), log10(W(j))
      end do
      close(7)
 
      stop
      end program
!======================================================================!
      subroutine calc_th(ndata, x, W)
      use global_const, only : A0, c
!$    use global_const, only : thread
      implicit none
      integer, intent(in) :: ndata
      integer :: j
      real(8), intent(in) :: x(ndata)
      real(8), intent(out) :: W(ndata)
      real(8), dimension(ndata) :: beta0, fac

      beta0 = 10.0d0 ** x
      fac = 2.0d0 * A0 * beta0

!$    call omp_set_num_threads(thread)
!$OMP parallel default(none) shared(ndata,fac,beta0,W) private(j)
!$OMP do
      do j=1, ndata
        if (beta0(j) <= 100.0d0) then
          call equ_width(c, fac(j), beta0(j), 0.01d0, 3.5d0, W(j))
        else if (beta0(j) <= 1.0d3) then
          call equ_width(c, fac(j), beta0(j), 0.01d0, 15.0d0, W(j))
        else
          call equ_width(c, fac(j), beta0(j), 0.05d0, 120.0d0, W(j))
        end if
      end do
!$OMP end do
!$OMP end parallel

      return
      end subroutine
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
!======================================================================!
      subroutine read_data(pms, x, y)
      implicit none
      real(8), intent(in) :: pms(9)
      real(8), intent(out), dimension(85) :: x, y
      real(8), parameter :: m = 0.001d0
      real(8), dimension(85) :: W, lambda, loggf, kai, c, e
      real(8) :: theta

!-----FeI-------------------------------------------------------
      lambda(1) = 5166.29d0
      lambda(2) = 5247.06d0
      lambda(3) = 5254.96d0
      lambda(4) = 5110.41d0
      lambda(5) = 5168.90d0
      lambda(6) = 5225.53d0
     
      W(1) = 115.0d0 * m
      W(2) =  59.0d0 * m 
      W(3) =  92.0d0 * m
      W(4) = 126.0d0 * m
      W(5) = 114.0d0 * m
      W(6) =  68.0d0 * m

      loggf(1) = - 3.68d0
      loggf(2) = - 4.50d0
      loggf(3) = - 4.23d0
      loggf(4) = - 3.34d0
      loggf(5) = - 3.49d0
      loggf(6) = - 4.26d0
      
      kai(1:6) = 0.05d0
!---FeI-Multiplet15------------------
      lambda(7)  = 5328.05d0
      lambda(8)  = 5405.78d0
      lambda(9)  = 5397.15d0
      lambda(10) = 5429.70d0
      lambda(11) = 5446.92d0
      lambda(12) = 5455.61d0

      W(7)  = 375.0d0 * m
      W(8)  = 266.0d0 * m
      W(9)  = 239.0d0 * m
      W(10) = 285.0d0 * m
      W(11) = 238.0d0 * m
      W(12) = 219.0d0 * m

      loggf(7)  = - 1.43d0
      loggf(8)  = - 1.78d0
      loggf(9)  = - 1.85d0
      loggf(10) = - 1.76d0
      loggf(11) = - 1.86d0
      loggf(12) = - 2.01d0

      kai(7:12) = 0.97d0

!----FeI-Muliplet66---------------------
      lambda(13) = 5145.10d0
      lambda(14) = 5131.48d0
      lambda(15) = 5098.70d0
      lambda(16) = 5079.23d0
      lambda(17) = 5250.65d0
      lambda(18) = 5198.71d0

      W(13) =  44.0d0 * m
      W(14) =  72.0d0 * m
      W(15) = 102.0d0 * m
      W(16) = 100.0d0 * m
      W(17) = 104.0d0 * m
      W(18) =  87.0d0 * m

      loggf(13) = -2.40d0
      loggf(14) = -1.92d0
      loggf(15) = -1.40d0
      loggf(16) = -1.45d0
      loggf(17) = -1.52d0
      loggf(18) = -1.50d0

      kai(13:18) = 2.20d0

!----FeI-Multiplet687------------------
      lambda(19) = 4966.10d0
      lambda(20) = 4946.39d0
      lambda(21) = 4882.15d0
      lambda(22) = 4963.65d0
      lambda(23) = 4875.90d0
      lambda(24) = 4855.68d0
      lambda(25) = 4843.16d0
      lambda(26) = 4838.52d0
      lambda(27) = 5039.26d0
      lambda(28) = 5002.80d0
      lambda(29) = 4950.11d0
      lambda(30) = 4907.74d0

      W(19) = 114.0d0 * m
      W(20) = 113.0d0 * m
      W(21) =  70.0d0 * m
      W(22) =  48.0d0 * m
      W(23) =  55.0d0 * m
      W(24) =  60.0d0 * m
      W(25) =  67.0d0 * m
      W(26) =  51.0d0 * m
      W(27) =  73.0d0 * m
      W(28) =  85.0d0 * m
      W(29) =  76.0d0 * m
      W(30) =  61.0d0 * m

      loggf(19) = - 0.30d0
      loggf(20) = - 0.74d0
      loggf(21) = - 1.10d0
      loggf(22) = - 1.21d0
      loggf(23) = - 1.39d0
      loggf(24) = - 1.33d0
      loggf(25) = - 1.30d0
      loggf(26) = - 1.39d0
      loggf(27) = - 0.89d0
      loggf(28) = - 1.03d0
      loggf(29) = - 1.08d0
      loggf(30) = - 1.33d0

      kai(19:30) = 3.40d0

!---TiI-Multiplet3--------------------
      lambda(31) = 5460.51d0
      lambda(32) = 5426.26d0
      lambda(33) = 5490.84d0

      W(31) = 8.5d0 * m
      W(32) = 5.5d0 * m
      W(33) = 2.5d0 * m

      loggf(31) = - 2.36d0
      loggf(32) = - 2.48d0
      loggf(33) = - 2.67d0
 
      kai(31:33) = 0.03d0

!---TiI-Multiplet4-----------------------------------
      lambda(34) = 5210.39d0
      lambda(35) = 5192.97d0
      lambda(36) = 5173.74d0
      lambda(37) = 5219.71d0
      lambda(38) = 5152.20d0
      lambda(39) = 5147.48d0

      W(34) = 86.0d0 * m
      W(35) = 80.0d0 * m
      W(36) = 67.0d0 * m
      W(37) = 25.0d0 * m
      W(38) = 38.0d0 * m
      W(39) = 36.0d0 * m

      loggf(34) = - 0.90d0
      loggf(35) = - 0.96d0
      loggf(36) = - 1.06d0
      loggf(37) = - 1.90d0
      loggf(38) = - 1.73d0
      loggf(39) = - 1.71d0

      kai(34:39) = 0.03d0
!---TiI-Multiplet---------------------------------
      lambda(40) = 4981.73d0
      lambda(41) = 4991.07d0
      lambda(42) = 4999.50d0
      lambda(43) = 5016.16d0
      lambda(44) = 5020.03d0
      lambda(45) = 5022.87d0
      lambda(46) = 5024.84d0
      lambda(47) = 5045.40d0
      lambda(48) = 5043.58d0
      lambda(49) = 5040.64d0

      W(40) = 112.0d0 * m
      W(41) = 102.0d0 * m
      W(42) = 104.0d0 * m
      W(43) =  60.0d0 * m
      W(44) =  86.0d0 * m
      W(45) =  72.0d0 * m
      W(46) =  62.0d0 * m
      W(47) =  10.0d0 * m
      W(48) =  14.0d0 * m
      W(49) =  16.0d0 * m

      loggf(40) =   0.57d0
      loggf(41) =   0.45d0
      loggf(42) =   0.38d0
      loggf(43) = - 0.44d0
      loggf(44) = - 0.29d0
      loggf(45) = - 0.30d0
      loggf(46) = - 0.47d0
      loggf(47) = - 1.49d0
      loggf(48) = - 1.30d0
      loggf(49) = - 1.37d0

      kai(40:49) = 0.82d0

!---TiI-Multiplet183---------------------------
      lambda(50) = 5224.30d0
      lambda(51) = 5224.56d0
      lambda(52) = 5223.62d0
      lambda(53) = 5222.69d0
      lambda(54) = 5263.48d0
      lambda(55) = 5247.29d0
      lambda(56) = 5186.33d0
      lambda(57) = 5194.04d0
      lambda(58) = 5201.10d0
      lambda(59) = 5207.85d0

      W(50) = 36.0d0 * m
      W(51) = 68.0d0 * m 
      W(52) = 11.0d0 * m
      W(53) = 23.0d0 * m
      W(54) = 13.0d0 * m
      W(55) = 10.0d0 * m
      W(56) =  7.0d0 * m
      W(57) = 10.0d0 * m
      W(58) = 11.0d0 * m
      W(59) =  8.0d0 * m

      loggf(50) =   0.42d0
      loggf(51) = - 0.03d0
      loggf(52) = - 0.09d0
      loggf(53) = - 0.05d0
      loggf(54) = - 0.27d0
      loggf(55) = - 0.15d0
      loggf(56) = - 0.36d0
      loggf(57) = - 0.08d0
      loggf(58) = - 0.22d0
      loggf(59) = - 0.16d0

      kai(50:59) = 2.09d0

!---FeII-Multiplet 37-------------------------------
      lambda(60) = 4472.93d0
      lambda(61) = 4489.19d0
      lambda(62) = 4491.41d0
      lambda(63) = 4555.89d0
      lambda(64) = 4582.84d0
      lambda(65) = 4520.23d0
      lambda(66) = 4534.17d0

      W(60) = 39.0d0 * m
      W(61) = 61.0d0 * m
      W(62) = 66.0d0 * m
      W(63) = 77.0d0 * m
      W(64) = 49.0d0 * m
      W(65) = 69.0d0 * m
      W(66) = 53.0d0 * m

      loggf(60) = - 4.76d0
      loggf(61) = - 3.69d0
      loggf(62) = - 2.78d0
      loggf(63) = - 2.42d0
      loggf(64) = - 3.17d0
      loggf(65) = - 3.20d0
      loggf(66) = - 3.33d0

      kai(60:66) = 2.83d0

!---FeII-Multiplet 38---------------------
      lambda(67) = 4508.28d0
      lambda(68) = 4576.34d0
      lambda(69) = 4620.51d0
      lambda(70) = 4541.52d0

      W(67) = 74.0d0 * m
      W(68) = 56.0d0 * m
      W(69) = 47.0d0 * m
      W(70) = 58.0d0 * m

      loggf(67) = - 2.42d0
      loggf(68) = - 2.89d0
      loggf(69) = - 3.13d0
      loggf(70) = - 2.93d0

      kai(67:70) = 2.82d0

!---TiII-Multiplet 18,19,20-------------------
      lambda(71) = 4469.16d0
      lambda(72) = 4493.53d0
      lambda(73) = 4395.03d0
      lambda(74) = 4443.80d0
      lambda(75) = 4450.49d0
       
      W(71) =  49.0d0 * m
      W(72) =  26.0d0 * m
      W(73) = 135.0d0 * m
      W(74) = 124.0d0 * m
      W(75) =  79.0d0 * m

      loggf(71) = - 3.00d0
      loggf(72) = - 3.62d0
      loggf(73) = - 0.55d0
      loggf(74) = - 0.73d0
      loggf(75) = - 1.56d0
  
      kai(71:75) = 1.08d0

!---TiII-Multiplet 50,51,60---------------
      lambda(76) = 4533.97d0
      lambda(77) = 4563.76d0
      lambda(78) = 4568.30d0
      lambda(79) = 4589.96d0
      lambda(80) = 4432.09d0
      lambda(81) = 4399.77d0
      lambda(82) = 4394.06d0
      lambda(83) = 4418.34d0
      
      W(76) = 109.0d0 * m
      W(77) = 120.0d0 * m
      W(78) =  25.0d0 * m
      W(79) =  70.0d0 * m
      W(80) =  15.0d0 * m
      W(81) = 115.0d0 * m
      W(82) =  72.0d0 * m
      W(83) =  70.0d0 * m

      loggf(76) = - 0.70d0
      loggf(77) = - 0.86d0
      loggf(78) = - 3.45d0
      loggf(79) = - 1.72d0
      loggf(80) = - 3.55d0
      loggf(81) = - 1.32d0
      loggf(82) = - 1.71d0
      loggf(83) = - 2.24d0
  
      kai(76:83) = 1.22d0

!---TiII-Multiplet 105-----------------------------
      lambda(84) = 4386.85d0
      
      W(84) = 59.0d0 * m

      loggf(84) = - 0.79d0

      kai(84) = 2.59d0

!---TiII Multiplet 115---------------------------------
      lambda(85) = 4488.32d0
      W(85)      = 45.0d0 * m
      loggf(85)  = - 0.62d0
      kai(85)    = 3.10d0
!---------------------------------------------------------
!     i=1~30-----FeI
!     i=31~59----TiI
!     i=60~70----FeII
!     i=71~85----TiII
!---------------------------------------------------------
      c(1:30)  = pms(1)
      c(31:59) = pms(3)
      c(60:70) = pms(5)
      c(71:85) = pms(7)

      e(1:30)  = pms(2) 
      e(31:59) = pms(4)
      e(60:70) = pms(6)
      e(71:85) = pms(8)

      theta = pms(9)

!-----------------------------------------------
      
      x = loggf + log10(lambda) - kai * theta + c
      y = log10(W / lambda) + e
      
      call sort(x,y,85,1,85)

      return
      end subroutine
!====================================================================!
      subroutine sort(a, b, N, p, q)
      integer :: i,j,k,N,p,q
      real(8) :: a(N),amin,b(N),s1,s2
      
      do i=p, q-1
        amin = a(i)
        k = i
        do j=i+1, q
          if (a(j) < amin) then
            amin = a(j)
            k = j
          end if
        enddo
        s1 = a(i)
        s2 = b(i)
        a(i) = amin
        b(i) = b(k)
        a(k) = s1
        b(k) = s2
      enddo

      return
      end subroutine
!======================================================================!
!     program main
!----------------------------------------------------------------------!
!     Main program.                                                    !
!----------------------------------------------------------------------!
!     use Com_var, only : Lc, Apro, Atar, pn, Maxdata
!     use fitting_parameter
!     use exp_data
!     use Disp_pm
!     implicit none
!     integer :: n
!     real(8), dimension(Maxdata) :: th, ds, er, Ruth
!     real(8), external :: kai2
!     real(8) :: sigma(0:Lc), pms(pn), dkai, z, tmp
!     complex(8) :: Sn(0:Lc)
!     type(fitp) :: pm
!     type(xdata), allocatable :: xd(:)
!     character(len=30), parameter :: FM='(1x,a,f12.3)'

!     pm = fitp(-10.0d0, -10.0d0, 1.0d0, 1.0d0, 6.0d0, 6.0d0)
!     pms(1) = pm%V0  ;  pms(2) = pm%W0
!     pms(3) = pm%a1  ;  pms(4) = pm%a2
!     pms(5) = pm%R1  ;  pms(6) = pm%R2

!     call Read_Data(n, th, ds, er, Ruth, 'experiment.dat')
!     allocate(xd(n))
!     xd(1:n)%th = th(1:n)  ;  xd(1:n)%ds   = ds(1:n)
!     xd(1:n)%er = er(1:n)  ;  xd(1:n)%Ruth = Ruth(1:n)

!     write(6,*)
!     write(6,*) '** Before fitting **'
!     write(6,FM) 'kai square =', kai2(pms, n, xd)
!     call  Disp_parameters(pms)
!     call  Phase_Shift(Sn, sigma, pms)
!     call  Gnuplot(sigma, Sn, 7, 'gnubefore', '12C_on_28Si_bf.eps')
!     call  OROF_method(n, xd, pms)
!     write(6,*)
!     write(6,*) '** After  fitting **'
!     write(6,FM) 'kai square =', kai2(pms, n, xd)
!     call  Disp_parameters(pms)
!     call  Phase_Shift(Sn, sigma, pms)
!     call  Disp_diff_CS(sigma, Sn)
!     call  Gnuplot(sigma, Sn, 8, 'gnuopt', '12C_on_28Si.eps')

!     stop
!     end program
!======================================================================!
      subroutine OROF_method(ndata, npms, pms)
!----------------------------------------------------------------------!
!     A subroutine which executes Oak Ridge and Oxfort method.         !
!     We compute increment of parameters by subroutine 'Solve_dp' and  !
!     compute difference of kai squares with new and old parameters.   !
!     If the difference becomes sufficiently small, we terminate the   !
!     iteration.                                                       !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: npms, ndata
      integer, parameter :: MaxIter = 1000
      integer :: i, k
      real(8), intent(inout) :: pms(npms)
      real(8), external :: kai2
      real(8), parameter :: epsr = 1.0d-8
      real(8) :: kai0, dkai
      real(8), dimension(npms) :: dp, pms0, x, th, xd
      character(len=11), parameter :: FM='(1x,a,i5)'

      do i=1, MaxIter
        pms0 = pms
        call read_data(pms0, x, xd)
        call calc_th(ndata, x, th)
        kai0 = kai2(ndata, th, xd)
        call solve_dp(ndata, npms, dp, pms)
        call read_data(pms, x, xd)
        call calc_th(ndata, x, th)
        dkai = (abs(kai2(ndata, th, xd) - kai0)) / kai0
        if (dkai < epsr) then
          write(6,FM) 'Iteration :',i
          exit
        end if
      end do
      if (i==MaxIter .or. kai0>10.0d0) print *,  'We failed to fitting.'

      return
      end subroutine
!======================================================================!
      subroutine  Solve_dp(ndata, npms, dp, pms)
!----------------------------------------------------------------------!
!     This subroutine returns the increment of fitting parameters      !
!     'dp(pn)'.                                                        !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: ndata, npms
      integer :: i
      real(8), intent(in)  :: pms(npms)
      real(8), intent(out) :: dp(npms)
      real(8), dimension(ndata) :: xd, th
      real(8) :: dF(npms), ddF(npms,npms), phi, z

      call  Deriv(ndata, npms, pms, dF, ddF)
      call  Conj_Grad(npms, ddF, -dF, dp)
      phi = maxval(dp / pms)
      if (phi < 0.1d0) then
         z = 1.0d0
      else if (phi < 0.3d0) then
         z = 0.5d0
      else if (phi < 0.5d0) then
         z = 0.2d0
      else if (phi < 1.0d0) then
         z = 0.1d0
      else if (phi >= 1.0d0) then
         z = 0.05d0 / phi
      end if
      call  Parab_Approx(npms, ndata, dp, z, pms)

      return
      end subroutine
!======================================================================!
      subroutine  Deriv(ndata, npms, pms0, dF, ddF)
!----------------------------------------------------------------------!
!     This subroutine computes the first and second derivatives of     !
!     kai square(F) as a function of fitting parameters.               !
!     We increment the parameters by one percent and compute the       !
!     derivatives of the cross secion.                                 !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: ndata, npms
      integer :: i, j, k
      real(8), intent(in) :: pms0(npms)
      real(8), intent(out) :: dF(npms), ddF(npms,npms)
      real(8), external :: dsigma, dsc
      real(8) :: pms(npms), S, rf(npms,ndata)
      real(8), dimension(ndata) :: xd, xd0, x, kai, th, th0

      ddF = 0.0d0

      pms = pms0
      do i=1, npms
        if (i /= 1) pms(i-1) = pms0(i-1)
        pms(i) = 1.01d0 * pms0(i)
        call read_data(pms0, x, xd0)
        call calc_th(ndata, x, th0)
        call read_data(pms, x, xd)
        do k=1, ndata
          kai(k) = (xd0(k) - th(k)) / (th(k) * 0.1d0) ** 2
          rf(i,k) = xd(k) - xd0(k)
        end do
        dF(i) = 2.0d2 * dot_product(kai, rf(i,1:ndata)) / pms0(i)
        do j=1, i
          S = 0.0d0
          do k=1, ndata
            S = S + rf(i,k) * rf(j, k) / (th(k) * 0.1d0) ** 2
          end do
          ddF(i,j) = 2.0d4 * S / (pms0(i) * pms0(j))
        end do
      end do
      ddF = ddF + transpose(ddF)

      return
      end subroutine
!====================================================================!
      subroutine  Parab_Approx(npms, ndata, dp, z, p0)
!--------------------------------------------------------------------!
!     This subroutine computes a factor 'z' which is multipled to    !
!     the increment of fitting parameters 'dp(pn)' and returns new   !
!     fitting parameters. The factor 'z' is computed according to    !
!     the parabolic approximation.                                   !
!--------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: npms, ndata
      real(8), intent(inout) :: p0(npms), z
      real(8), intent(in) :: dp(npms)
      real(8), dimension(npms) :: pa, pb
      real(8), external :: kai2
      real(8), dimension(ndata) :: x, xd, th
      real(8) :: y0, ya, yb, z0, za, zb, zm, ym, w, concave

      call read_data(p0, x, xd)
      call calc_th(ndata, x, th)
      y0 = kai2(ndata, th, xd)
      pa = p0 + z * dp
      call read_data(p0, x, xd)
      call calc_th(ndata, x, th)
      ya = kai2(ndata, th, xd)
      if (ya <= y0) then
        za = z
        zb = 2.0d0 * z
      else
        za = 0.5d0 * z
        zb = z
      end if
      pa = p0 + za * dp
      call read_data(pa, x, xd)
      call calc_th(ndata, x, th)
      ya = kai2(ndata, th, xd)
      pb = p0 + zb * dp   
      call read_data(pb, x, xd)
      call calc_th(ndata, x, th)
      yb = kai2(ndata, th, xd) 
      concave = (yb - y0) * za + (y0 - ya) * zb
      if (concave > 0.0d0) then
        w  = 3.0d0 * y0 - 4.0d0 * ya + yb
        zm = 0.5d0 * za * w / (y0 - 2.0d0 * ya + yb)
        if (zm <= - zb) then
          z = - zb
        else if (zm <= zb) then
          z = zm
        else if (zm > zb) then
          z = 2.0d0 * zb
        end if
      else
        if (yb < y0) then
          z = 3.0d0 * za
        else if (yb >= y0) then
          z = - za
        end if
      end if
      p0 = p0 + z * dp

      return
      end subroutine
!======================================================================!
      subroutine  Conj_Grad(N, a, b, x)
!----------------------------------------------------------------------!
!     A subroutine which solves coupled linear equation by conjugate-  !
!     gradient method.                                                 !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: N
      integer :: i, Step
      real(8), intent(in) :: a(N,N), b(N)
      real(8), intent(inout) :: x(N)
      real(8), dimension(N) :: r, p, q
      real(8) :: C1, C2, alpha, beta

      call  Guess_x(N, x)
      r = b - Matmul(a, x)
      p = r
      C2 = Dot_Product(x, b + r)

      Step = 0
      do 
        Step = Step + 1
        C1 = C2
        q = Matmul(a, p)
        alpha = Dot_product(p, r) / Dot_Product(p, q)
        x = x + alpha * p
        if (Mod(Step, 3) == 0) then
            r = b - Matmul(a, x)
        else
            r = r - alpha * q
        end if
        C2 = Dot_Product(x, b + r)
        if (C1 >= C2) exit
        beta = - Dot_Product(r, q) / Dot_Product(p, q)
        p = r + beta * p
      end do
          
      return
      end subroutine
!======================================================================!
      subroutine  Guess_x(N, x)
!----------------------------------------------------------------------!
!     A subroutine which gives the initial guess of the solutions of   !
!     coupled linear equation.                                         !
!----------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: N
      integer :: i
      real(8), intent(out) :: x(N)

      x = 1.0d0

      return
      end subroutine
!====================================================================!
      function kai2(n, th, xd) result(kai)
      implicit none
      integer, intent(in) :: n
      integer :: i
      real(8), intent(in), dimension(n) :: th, xd
      real(8) :: kai

      kai = 0.0d0
      do i=1, n
        kai = kai + ((th(i) - xd(i)) / (th(i) * 0.1d0)) ** 2
      end do
      kai = kai / dble(n)

      return
      end function
!======================================================================!
!     subroutine  Show_Matrix(N, a)
!     implicit none
!     integer, intent(in) :: N
!     integer :: i, j
!     real(8), intent(in) :: a(N,N)
!     character(len=30) :: FM

!     FM = '(x,a,' // CHAR(48+N) // 'es9.1,a)'
!     do i=1, N
!         write(6,FM) ' |', (a(i,j), j=1, N), ' |'
!     end do
!     write(6,*)

!     return
!     end subroutine
