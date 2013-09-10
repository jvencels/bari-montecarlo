      module particles
        use, intrinsic :: iso_fortran_env, only : REAL64, INT16, INT32, INT64, stdout=>output_unit
        use physics, only : particle_type, new_particles, first_particle, last_particle, handle_collision
   
        ! Some parameters for simulation
        integer tidx
        real(REAL64) e0
        integer(INT64) i
        integer(INT64) N_init_particles ! Number of initial particles (needed for stationary case)

       public :: simulate

      contains

        ! =======================================================

        subroutine simulate()
          use eedf, only : Ntimes, calculate_totals
          use physics, only : dt, tfin, new_particles
          use fparser, only : evalf
          use io, only : simulation_mode
          implicit none

          real(REAL64) t
          integer(REAL64) particles_to_add

          ! Initialize simulation
          t = 0.0_REAL64
          tidx = 1

          ! Initialize particles for time=0 (equal to the constant in
          ! particle initialization expression)
          N_init_particles = int(evalf(1, (/ .0_REAL64 /)))
          if (N_init_particles .gt. 0) then
            ! Initialize new particles with fparser
            call new_particles(N_init_particles, e0)
          else
            write(stdout,'(A)') "Error: At the beginning must be initialized at least 1 particle!"
            stop
          end if

          ! Collect initial state
          call collect(0)

          ! Run simulation
          do while (tidx .le. Ntimes .and. t .le. tfin)
            
            if (simulation_mode .eq. 1) then
              particles_to_add = int(evalf(1, (/ t+dt /)) - evalf(1, (/ t /)))
              if (particles_to_add .lt. 0) then
              write(stdout,'(A,E10.3,A)') "Error: Derivative of the particle initialization function at time ", &
                   t, " has negative value (can't delete particles)"
                exit
              else
                ! Initialize new particles with fparser
                call new_particles(particles_to_add, e0)
              end if
            end if

            ! propagate for dt
            call propagate(dt)
            ! collect snapshot histogram
            call collect(tidx)
            call calculate_totals(tidx)
            ! update current time and index
            t = t+dt
            
            if (simulation_mode .eq. 1) write(stdout, '(A,E10.3)') "Time :", t
            tidx = tidx + 1

          enddo
          ! if tfin/dt is not even, propagate to tfin
          if (t .lt. tfin .and. (abs(tfin-t) .gt. dt*1e-6)) then
            call propagate(tfin-t)
            ! collect final value
            call collect(tidx-1)
            call calculate_totals(tidx)
          end if

          call del_all_particles()

        end subroutine simulate

        ! =======================================================

        subroutine del_all_particles()
          implicit none

          type(particle_type), pointer :: tmp_p, next_p

          tmp_p => first_particle

          do
            if (.not. associated(tmp_p%next_particle)) then
              deallocate(tmp_p)
              exit
            end if
            next_p => tmp_p%next_particle
            deallocate(tmp_p)
            tmp_p => next_p
          enddo
        end subroutine

        ! =======================================================

        subroutine propagate(dt)
          implicit none

          type(particle_type), pointer :: tmp_particle
          real(REAL64) dt, dtp

          tmp_particle => first_particle
          do
            dtp = dt
            do while(tmp_particle%tcs .le. dtp .and. tmp_particle%tcs .gt. 0)
              dtp = dtp - tmp_particle%tcs
              call handle_collision(tmp_particle)
            enddo
            tmp_particle%tcs = tmp_particle%tcs - dtp
            if(.not.associated(tmp_particle%next_particle))then
              exit
            end if
            tmp_particle => tmp_particle%next_particle
          end do

        end subroutine propagate

        ! =======================================================

        subroutine collect(tidx)
          use eedf, only : eedfbins, Needfbins, de, total_e, aver_energy_eV, e_below_10eV, e_e0, recycled_e
          use io, only : simulation_mode
          implicit none

          type(particle_type), pointer :: tmp_particle
          integer idx, tidx

          tmp_particle => first_particle
          do while (associated(tmp_particle))
            idx = int( tmp_particle%E_eV / de )
            if (idx .eq. 0) then
              idx = 1
            end if
            eedfbins(idx, tidx) = eedfbins(idx, tidx) + 1
            total_e(tidx) = total_e(tidx) + 1
            ! Later energy will be divided by "total_e"
            aver_energy_eV(tidx) = aver_energy_eV(tidx) + tmp_particle%E_eV

            if (idx .eq. Needfbins) then
              e_e0(tidx) = e_e0(tidx) + 1
            end if

            if (simulation_mode .eq. 0 .and. idx .le. 5) then
              e_below_10eV(tidx) = e_below_10eV(tidx) + 1
              tmp_particle%E_eV = e0
            end if

            tmp_particle => tmp_particle%next_particle
          enddo
          aver_energy_eV(tidx) = aver_energy_eV(tidx) / total_e(tidx)

          if (simulation_mode .eq. 0 .and. total_e(tidx) .GE. 2*N_init_particles) then
            call delete_odd_particles()
          end if

        end subroutine collect

        ! =======================================================

        subroutine delete_odd_particles()
          implicit none

          type(particle_type), pointer :: odd_p, even_p
          even_p => first_particle%next_particle

          deallocate(first_particle)
          first_particle => even_p

          do
            if (.not. associated(even_p%next_particle)) then
              ! odd_p doesn't exist
              even_p%next_particle => null()
              last_particle => even_p
              exit
            end if
            ! odd_p exist
            odd_p => even_p%next_particle
            if (.not. associated(odd_p%next_particle)) then
              ! even_p doesn't exist
              deallocate(odd_p)
              even_p%next_particle => null()
              last_particle => even_p
              exit
            end if
            even_p%next_particle => odd_p%next_particle
            even_p => odd_p%next_particle
            deallocate(odd_p)
          end do

        end subroutine delete_odd_particles

        ! =======================================================

      end module particles
