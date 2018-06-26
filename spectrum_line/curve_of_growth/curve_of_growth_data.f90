      program curve_of_growth
      implicit none
      integer i
      real(8), parameter :: m = 0.001d0
      real(8), dimension(100) :: W, lambda, loggf, c, y, x, kai, e
      real(8) :: theta, d
      character(len=40) :: FNAME(5)
      data FNAME /'Equ_width_FeI.dat',   &
                  'Equ_width_TiI.dat',   &
                  'Equ_width_FeII.dat',  &
                  'Equ_width_TiII.dat',  &
                  'Equ_width_data.dat'/ 
          
      do i=1,5
        open(unit=100*i, file=FNAME(i))
      enddo
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
      
      do i=1, 6
        kai(i) = 0.05d0
      enddo
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

      do i=7, 12
        kai(i) = 0.97d0
      enddo

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

      do i=13, 18
        kai(i) = 2.20d0
      enddo

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

      do i=19, 30
        kai(i) = 3.40d0
      enddo

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
 
      do i=31, 33
        kai(i) = 0.03d0
      enddo

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

      do i=34, 39
        kai(i) = 0.03d0
      enddo
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

      do i=40, 49
        kai(i) = 0.82d0
      enddo

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

      do i=50, 59
        kai(i) = 2.09d0
      enddo

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

      do i=60, 66
        kai(i)=2.83d0
      enddo

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

      do i=67, 70
        kai(i) = 2.82d0
      enddo

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
  
      do i=71, 75
        kai(i) = 1.08d0
      enddo

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
  
      do i=76, 83
        kai(i) = 1.22d0
      enddo

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
      do i=1,30
        c(i) = 3.85d0
        e(i) = 5.35d0 
      enddo
      do i=31,59
        c(i) = - 2.00d0
        e(i) = 5.00d0
      enddo
      do i=60,70
        c(i) = 5.1d0
        e(i) = 5.35d0
      enddo
      do i=71,85
        c(i) = 0.8d0
        e(i) = 5.0d0
      enddo

      theta=0.92d0

!-----------------------------------------------
      
      do i=1,85
        x(i) = loggf(i) + log10(lambda(i)) - kai(i) * theta + c(i)
        y(i) = log10(W(i) / lambda(i)) + e(i)
      enddo
      
      call sort(x,y,85,1,30)
      
      do i=1, 30
        write(100,*) x(i),y(i)
      enddo
   
      call sort(x,y,85,31,59)
      
      do i=31,59
        write(200,*) x(i),y(i)
      enddo

      call sort(x,y,85,60,70)
      
      do i=60,70
        write(300,*) x(i),y(i)
      enddo

      call sort(x,y,85,71,85)
      
      do i=71,85
        write(400,*) x(i),y(i)
      enddo

      call sort(x,y,85,1,85)
      
      do i=1,85
        write(500,*) x(i),y(i)
      enddo

      close(100)

      stop
      end program

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




