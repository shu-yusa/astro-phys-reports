!======================================================================!
      program main
      implicit none
      integer :: n
      real(8) :: xp, x
      real(8), dimension(2) :: y, y0
      real(8), parameter :: h = 1.0d-6
      character(len=30), parameter :: FM = '(1x,a,f8.5)'

      open(7,file='test')

      xp = 0.0d0

      x = 1.0d-10
      n = 1
      call Init(n, x, y)
      do
        if (x > xp) then
          write(7,*) x, y(1)
          xp = xp + 0.01d0
        end if
        y0 = y
        call runge_kutta(n, x, y, h)
        if (y(1) < 0.0d0) exit
        x = x + h
      end do
      x = x - h
      call zero(n, x, y0, h)

      write(6,*)  'Order :', n
      write(6,FM) 'The first zero point :', x

      close(7)

      stop
      end program
!======================================================================!
      subroutine zero(n, x, y, h)
      implicit none
      integer, intent(in) :: n
      real(8), intent(inout) :: x, y(2), h
      real(8), parameter :: eps = 1.0d-10
      real(8) :: x0, h2, d, y0(2)

      d = 0.5d0
      do 
        h2 = d * h
        y0 = y
        call runge_kutta(n, x, y, h2)
        if (abs(y(1)) < eps) then
          x = x + h2
          return
        else if (y(1) > 0.0d0) then
          x = x + h2
        else
          d = 0.5d0 * d
          y = y0
        end if
      end do

      return
      end subroutine
!======================================================================!
      subroutine Init(n, x, y)
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: x
      real(8), intent(out) :: y(2)

      y(1) = 1.0d0 - 1.0d0/6.0d0 * x * x + dble(n)/120.0d0 * x ** 4
      y(2) = - 1.0d0 / 3.0d0 * x + dble(n) / 30.0d0 * x ** 3

      return
      end subroutine
!======================================================================!
      subroutine runge_kutta(n, x, y, h)
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: x, h
      real(8), intent(inout) :: y(2)
      real(8), dimension(2) :: k1, k2, k3, k4
      real(8) :: h2

      h2 = 0.5d0 * h
      call fct(n, x, y, k1)
      call fct(n, x+h2, y+h2*k1, k2)
      call fct(n, x+h2, y+h2*k2, k3)
      call fct(n, x+h, y+h*k3, k4)
      y = y + h * (k1 + 2.0d0 * (k2 + k3) + k4) / 6.0d0

      return
      end subroutine
!======================================================================!
      subroutine fct(n, x, y, f)
      implicit none
      integer, intent(in) :: n
      real(8), intent(in) :: x, y(2)
      real(8), intent(out) :: f(2)

      f(1) = y(2)
      f(2) = - y(1) ** n - 2.0d0 / x * y(2)

      return
      end subroutine
!======================================================================!
