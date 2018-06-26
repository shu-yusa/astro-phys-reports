!======================================================================!
!     Title  : Exp_Int.f90                                             !
!     Author : Yusa Shusaku                                            !
!     Last modification : 2010-1-12-Tue                                !
!                                                                      !
!     using library : mylib, pgplot                                    !
!                                                                      !
!======================================================================!
      program main
      use gaus_weights
      implicit none
      integer :: i, j, n, MESH
      integer :: ios, pgopen
      integer, parameter :: nmax = 5, m = 100  ! input
      real(8), parameter :: dx = 0.01d0        ! input
      real(8), parameter :: xmax = 3.0d0       ! input
      real(8), dimension(m) :: t, w
      real(8), allocatable :: x(:), En(:,:)
      character(len=30) :: FM

      MESH = nint(xmax / dx) + 1
      allocate(x(MESH), En(nmax,MESH))
      En = 0.0d0
      forall (j=1:MESH) x(j) = dble(j-1) * dx

      call weights_lagu(m, t(1:m), w(1:m))
      do n=1, nmax
        do j=1, MESH
          do i=1, m
            En(n,j) = En(n,j) + w(i) * f(n,x(j),t(i))
          end do
        end do
        En(n,:) = En(n,:) * exp(- x)
      end do

      ios = pgopen('/xserv')
      call pgsch(1.0)
      call pgenv(0.0, real(xmax), 0.0, 2.0, 0, 1)
      call pgscf(2)
      call pglab('x', 'En(x)', 'Exponential Integral')
      call pgslw(3)
      do i=1, nmax
        call pgsci(i+1)
        call pgtext(2.5, 1.9-0.12*real(i) ,'n='//char(48+i))
        call pgline(MESH, real(x), real(En(i,:)), 0, 1)
      end do
      call pgsci(1)
      call pgenv(0.0, real(xmax), -0.030001, 0.0, 0, 1)
      call pglab('x', 'dF(x)/F(x)', 'Error of flux')
      call pgsci(2)
      call pgline(MESH, real(x), real(En(3,:)-1.5d0*En(4,:)), 0, 1)
      call pgend

      open(7,file='ExpInt.dat')
      open(8,file='Error_of_flux.dat')
      write(FM,*) nmax
      FM = '(f8.4,'//trim(adjustl(FM))//'es13.4)'
      do j=1, MESH
        write(7,FM) x(j), En(:,j)
        write(8,*) x(j), En(3,j) - 1.5d0 * En(4,j)
      end do
      close(7)
      close(8)

      contains
!======================================================================!
      function f(n, x, t) result(g)
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: x, t
      real(8) :: g

      g = exp(t * (1.0d0 - x)) / (t + 1.0d0) ** n

      return
      end function
!======================================================================!
      end program
