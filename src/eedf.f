      module eedf
        use, intrinsic :: iso_fortran_env, only : REAL64, INT16, INT32, INT64, stdout=>output_unit
        private

        real(REAL64), allocatable :: eedfsum(:,:), totsum(:)
        real(REAL64), allocatable :: norm_factor(:), convergence_criteria(:)
        integer(INT32), allocatable, public :: total_e(:), e_below_10eV(:), e_e0(:)
        real(REAL64), allocatable, public :: eedfbins(:,:), aver_energy_eV(:), normalized_electr_distr(:,:)

        integer, parameter, public    :: Needfbins = int(5e2)
        integer, public :: Ntimes, recycled_e = 0
        real(REAL64), public :: de

        public :: init_eedfbins
        public :: calculate_totals
        public :: print_eedf
        public :: cleanup_histogram

      contains

        ! =======================================================

        subroutine init_eedfbins(e0)
          implicit none

          real(REAL64) e0

          allocate(eedfbins(Needfbins,0:Ntimes))
          de = e0/Needfbins
          eedfbins = 0
          allocate(aver_energy_eV(0:Ntimes))
          aver_energy_eV = 0
          allocate(total_e(0:Ntimes))
          total_e = 0
          allocate(e_below_10eV(0:Ntimes))
          e_below_10eV = 0
          allocate(e_e0(0:Ntimes))
          e_e0 = 0
          ! Normalized (integral=1) electron energy distribution function
          allocate(normalized_electr_distr(Needfbins,0:Ntimes))
          normalized_electr_distr = 0
          allocate(convergence_criteria(0:Ntimes))
          convergence_criteria = 0

        end subroutine init_eedfbins

        ! =======================================================

        subroutine calculate_totals(tidx)
          use io, only : simulation_mode
          implicit none

          integer ie, tidx

          ! normalize entire eedf to 1
          eedfbins(:,tidx) = eedfbins(:,tidx)/sum(eedfbins(:,tidx))

          ! Normalized (integral=1) electron energy distribution function
          normalized_electr_distr(:,tidx) = eedfbins(:,tidx)

          if (simulation_mode .eq. 0) then
            ! Since normalized_electr_distr is normalized to 1, criteria for
            ! stationary case is sum of absolute values of difference of EEDF.
            convergence_criteria(tidx) = &
                 sum(abs(normalized_electr_distr(:,tidx)-normalized_electr_distr(:,tidx-1)))
            write(stdout,'((A),(E15.2))') "Convergence criteria:", convergence_criteria(tidx)
          end if

          ! normalize to f(e) according to paper
          do ie=1,Needfbins
            eedfbins(ie,tidx) = normalize(ie*de, eedfbins(ie,tidx))
          end do
        end subroutine calculate_totals

        ! =======================================================

        function normalize(e, N)
          use physics, only : me
          implicit none

          real(REAL64) e, N, normalize

          real(REAL64), parameter :: pi = 3.14159265358979323846264338327950288419716
          real(REAL64), save :: K = 1 !/(4*pi)*(me/2)**(3/2)

          normalize = N * K * 1/(de*sqrt(e))

        end function normalize

        ! =======================================================

        subroutine print_eedf()
          use physics, only : dt, tfin
          use io, only : eedf_f, totals_f, norm_e_distr_f, convergence_f, simulation_mode
          implicit none

          integer it, i
          ! variables for output
          real(REAL64) t, e, p

          open(eedf_f, FILE="eedf.dat")
          do it=0, Ntimes
            do i=1, Needfbins
              t = min(dt*it, tfin)
              e = i*de
              p = eedfbins(i,it)
              write(eedf_f,'(3(E15.8))') t, e, p ! time, energy, probability
            end do
          end do
          close(eedf_f)

          open(norm_e_distr_f, FILE="norm-e-distr.dat")
          do it=0, Ntimes
            do i=1, Needfbins
              t = min(dt*it, tfin)
              e = i*de
              p = normalized_electr_distr(i,it)
              write(norm_e_distr_f,'(3(E15.8))') t, e, p ! time, energy, probability
            end do
          end do
          close(norm_e_distr_f)

          ! Prints to file time, total energy (eV) and number of electrons for each time step
          open(totals_f, FILE="totals.dat")
          do it=0, Ntimes
            t = min(dt*it, tfin)
            write(totals_f,'(2(E15.5),3(I10))') t, aver_energy_eV(it), total_e(it), e_below_10eV(it), e_e0(it)
          end do
          close(totals_f)

          open(convergence_f, FILE="convergence.dat")
          if (simulation_mode .eq. 0) then
            do it=0, Ntimes
              t = min(dt*it, tfin)
              write(convergence_f,'(2(E15.5))') t, convergence_criteria(it)
            end do
          end if
          close(convergence_f)

        end subroutine print_eedf

        ! =======================================================

        subroutine cleanup_histogram()
          implicit none
          deallocate(eedfbins)
          deallocate(aver_energy_eV)
          deallocate(total_e)
          deallocate(e_below_10eV)
          deallocate(e_e0)
          deallocate(normalized_electr_distr)
          deallocate(convergence_criteria)
        end subroutine

        ! =======================================================

      end module eedf
