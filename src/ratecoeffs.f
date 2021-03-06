      module ratecoeffs
        use, intrinsic :: iso_fortran_env, only : REAL64, INT16, INT32, INT64
        private

        real(REAL64), allocatable, public :: k(:,:)

        real(REAL64) tstep, endt

        public :: calculate_ratecoeffs_evolution, ratecoeff
        public :: print_ratecoeffs
        public :: clean_up_ratecoeffs
      contains

        subroutine calculate_ratecoeffs_evolution()
          use eedf, only : Needfbins, Ntimes
          use physics, only : NCollProc
          implicit none

          integer it

          allocate(k(0:Ntimes,NCollProc))

          do it=0,Ntimes
            call calculate_ratecoeffs(it)
          end do

        end subroutine calculate_ratecoeffs_evolution


        subroutine print_ratecoeffs()
          use eedf, only : Ntimes
          use physics, only : dt, tfin
          use io, only : rate_f

          implicit none

          real(REAL64) t
          integer :: it

          open(rate_f, file='rate.dat')

          do it = 0, Ntimes
            t = min(dt*it,tfin)
            write(rate_f,*), t, k(it,:)*1e6 ! convert to cm^3/s
          end do
        end subroutine print_ratecoeffs

        subroutine calculate_ratecoeffs(it)
          use eedf, only : Needfbins, de
          use physics, only : NCollProc, cross_section, me, e

          implicit none

          integer, intent(in) :: it
          real(REAL64), parameter :: pi = 3.14159265358979323846264338327950288419716
          real(REAL64), save :: propconst

          integer ip, ie

          ! This puts k in m^3/s
          propconst = 1/sqrt(2*me)*de*sqrt(e)

          ! calculate rate coefficient for process i, and put it in k(i)
          do ip=1,int(NCollProc)
            ! integrate \int e*f(e)*cs(e)*de using the trapezoidal rule
            k(it,ip) = propconst * sum((/ (integrand(ip,it,ie)+integrand(ip,it,ie+1), ie=1, Needfbins-1) /))
          end do
        end subroutine calculate_ratecoeffs

        function ratecoeff(t, ip)
          use physics, only : dt, tfin
          use eedf, only : Ntimes
          implicit none
          
          real(REAL64), intent(in) :: t
          integer, intent(in) :: ip

          integer it
          real(REAL64) delta, ratecoeff, grad

          it = int(t/dt)
          delta = t - it*dt

          if (it .lt. Ntimes) then
            grad = (k(it+1,ip) - k(it,ip)) / (min(tfin, (it+1)*dt) - it*dt)
            ratecoeff = k(it,ip) + delta*grad
          else
            ratecoeff = k(it,ip)
          end if
        end function ratecoeff

        function integrand(ip, it, ie)
          use eedf, only : eedfbins, de
          use physics, only : cross_section

          implicit none

          integer, intent(in) :: ie, ip, it
          real(REAL64) :: integrand

          integrand = (ie*de)*eedfbins(ie,it)*cross_section(ie*de,ip)
                      ! e        f(e)               sigma(e)         
        end function integrand

        subroutine clean_up_ratecoeffs()
          implicit none
          deallocate(k)
        end subroutine clean_up_ratecoeffs

      end module ratecoeffs