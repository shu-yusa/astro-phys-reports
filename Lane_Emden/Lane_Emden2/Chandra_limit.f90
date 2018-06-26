!======================================================================!
      module const
      implicit none
      real(8), parameter :: mpl = 2.17645*1.0d-5
      real(8), parameter :: mp = 1.6749*1.0d-24
      real(8), parameter :: me = 9.1093826*1.0d-28
      real(8), parameter :: mpc2 = 938.272d0
      real(8), parameter :: hbarc = 197.327d-13
      real(8), parameter :: Msun = 1.989*1.0d33
      real(8), parameter :: Rsun = 6.955*1.0d10
      end module
!======================================================================!
      program main
      use const
      implicit none
      integer :: rgrid
      real(8) :: x, dz
      real(8) :: zc, q, Mass, Radius
      real(8), dimension(2) :: y, y0
      real(8), parameter :: h = 1.0d-2
      real(8), parameter :: epsa = 1.0d-3
      real(8), parameter :: PI = 3.141592653589793d0
      real(8) :: w
      character(len=30), parameter :: FM = '(1x,a,f8.5)'

      w = sqrt(3.0d0 * PI / 4.0d0)
      zc = 1.0001d0
      dz = 0.01d0
      open(7,file='test')
      do 
        if (zc < 1.2d0) then
          dz = 1.0d-4
        else 
          dz = 0.2d0
        end if
        x = 1.0d-10
        q = sqrt(1.0d0 - 1.0d0 / (zc * zc))
        call Init(q, x, y)
        do
          call runge_kutta(zc, x, y, h)
          if (y(1) < 1.0d0 / zc) exit
          x = x + h
        end do
        Mass = w * mpl ** 3 / (mp * mp) * (- x * x * y(2))
        Radius = w * hbarc / mpc2 * mpl / me * x / zc
        write(7,*) zc, x, y(1), Mass/Msun, Radius/Rsun
        write(6,*) zc, x, y(1), Mass/Msun, Radius/Rsun
        zc = zc + dz
        if (Radius/Rsun <= epsa) exit
      end do
      close(7)

      stop
      end program
!======================================================================!
      subroutine Init(q, x, y)
      implicit none
      real(8), intent(in) :: x, q
      real(8), intent(out) :: y(2)

      y(1) = 1.0d0 - q ** 3 /6.0d0 * x * x + (q * x) ** 4 / 40.0d0
      y(2) = - q ** 3 / 3.0d0 * x + 0.1d0 * q ** 4 * x ** 3

      return
      end subroutine
!======================================================================!
      subroutine runge_kutta(zc, x, y, h)
      implicit none
      real(8), intent(in) :: x, h, zc
      real(8), intent(inout) :: y(2)
      real(8), dimension(2) :: k1, k2, k3, k4
      real(8) :: h2

      h2 = 0.5d0 * h
      call fct(zc, x, y, k1)
      call fct(zc, x+h2, y+h2*k1, k2)
      call fct(zc, x+h2, y+h2*k2, k3)
      call fct(zc, x+h, y+h*k3, k4)
      y = y + h * (k1 + 2.0d0 * (k2 + k3) + k4) / 6.0d0

      return
      end subroutine
!======================================================================!
      subroutine fct(zc, x, y, f)
      implicit none
      real(8), intent(in) :: x, y(2), zc
      real(8), intent(out) :: f(2)

      f(1) = y(2)
      f(2) = - 2.0d0/x*y(2) - abs(y(1)*y(1) - 1.0d0/(zc*zc)) ** 1.5d0

      return
      end subroutine
!======================================================================!
