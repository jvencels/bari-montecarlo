      module physics
        use, intrinsic :: iso_fortran_env, only : REAL64, INT16, INT32, INT64, stdout=>output_unit
        
        private
        !use interpolation

        ! nitrogen gas pressure and number density
        real(REAL64), public :: p, n
        ! number of interpolation points and collision processes
        integer(INT64), public :: Ncs, NCollProc
        ! cross section interpolations and boundary values
        ! the lower boundary is also the energy loss
        real(REAL64), allocatable, public :: cs(:,:), cs_min(:), cs_max(:)
        ! number of electrons created in processes
        integer, allocatable, public :: products(:)
        ! total collision frequency; constant through simulation
        real(REAL64), public :: total_collision_frequency

        real(REAL64), public :: me, e
        parameter(me=9.11E-31,e = 1.602176E-19)
        real(REAL64), public :: dt, tfin

        ! Type for making linked lists. Contains information about
        ! particular particle and address of the next particle.
        ! Code from A. Stock, "An Introduction to Fortran Pointer Techniques"
        type, public :: particle_type
          type(particle_type), pointer :: next_particle => null()
          real(REAL64)            :: E_eV = 0  ! Particle energy [eV]
          real(REAL64)            :: tcs = 0 ! Time till the next collision [s]
        end type particle_type

        type(particle_type), pointer, public :: first_particle, last_particle

        public :: init_physics
        public :: cross_section
        public :: collision_time
        public :: collision_frequency
        public :: handle_collision
        public :: velocity
        public :: new_particles

      contains

        ! =======================================================

        subroutine init_physics()
          implicit none

          integer(INT64) ie

          real(REAL64) :: cfreqs(Ncs)
          real(REAL64), parameter   :: T = 298.0
          real(REAL64), parameter   :: kB = 1.3806488e-23
          real(REAL64), parameter :: Torr2Pa = 101325/760
          

          ! calculate nitrogen number density from gas pressure
          n = p * Torr2Pa / (kB*T)

          ! calculate total collision frequency for each energy
          cfreqs = (/ (velocity(cs(ie,1))*n*sum(cs(ie,2:)), ie=1, Ncs) /)

          ! add the null collisions, by taking the maximum
          total_collision_frequency = maxval(cfreqs)
        end subroutine init_physics

        ! =======================================================

        function collision_time()
          use random

          implicit none
          real(REAL64) :: collision_time, r

          r = random_real()
          collision_time = -log(r) / total_collision_frequency
        end function collision_time

        ! =======================================================

        function collision_frequency(eV, process)
          implicit none

          integer process
          real(REAL64) eV, collision_frequency

          collision_frequency = velocity(eV)*n*cross_section(eV, process)
        end function collision_frequency

        ! =======================================================

        function cross_section(eV, process)
          implicit none

          ! Input argument and return value
          real(REAL64) eV, cross_section
          integer process

          ! Intermediates
          real(REAL64) dx
          integer idx

          if (eV .lt. cs_min(process) .or. eV .gt. cs_max(process)) then
            cross_section = 0
          else
            dx = cs(2,1)-cs(1,1)
            idx = 1 + int((eV-cs(1,1))/dx)
            cross_section = cs(idx,1+process)
          endif
        end function cross_section

        ! =======================================================

        subroutine handle_collision(tmp_particle)
          use random, only : random_real
          ! select a collision process, 

          implicit none
          ! other variables
          real(REAL64) r, cfsum,  E_2nd_particle
          integer ip
          type(particle_type), pointer :: tmp_particle

          ! select a collision process
          r = random_real()
          cfsum = 0.0
          ip = 0
          selectproc: do 
            ip = ip+1
            if (ip .gt. NCollProc) then 
              ! null collision!
              ! just give the particle a new collision time and exit
              !print *, "# null collision at process ", ip
              tmp_particle%tcs = collision_time()
              return
            end if

            cfsum = cfsum + collision_frequency(tmp_particle%E_eV, ip) / total_collision_frequency

            if (cfsum .gt. r) then
              ! this the process we want
              ! ip already has correct index
              exit selectproc
            end if
          end do selectproc

          ! subtract energy loss from available energy
          tmp_particle%E_eV = tmp_particle%E_eV - cs_min(ip)

          ! will a new electron be created?
          if (products(ip) .eq. 1) then

            E_2nd_particle = secondary_energy(tmp_particle%E_eV)
! CHECK how it works
            ! create a new electron which shares the available energy
            call new_particles(1_INT64,E_2nd_particle)

            ! we must also subtract this energy from the primary's
            tmp_particle%E_eV = tmp_particle%E_eV - E_2nd_particle
          end if

          ! update collision time
          tmp_particle%tcs = collision_time()
        end subroutine handle_collision

        ! =======================================================

        function secondary_energy(ea)
          use random

          implicit none
          real(REAL64) r, ea, secondary_energy, E
          parameter(E=13)
          r = random_real()
          secondary_energy = E*tan(atan(ea/(2*E))*r)    
        end function secondary_energy

        ! =======================================================

        function velocity(eV)
          ! Calculates the velocity in m/s for an electron with kinetic
          ! energy eV
          
          implicit none
          real(REAL64) eV, velocity

          velocity = sqrt(2*eV*e/me)
        end function velocity

        ! =======================================================

        function energy(v)
          ! Calculates the kinetic energy in eV for an electron which 
          ! moves at a velocity v (in m/s)
          !
          ! v is a 3-array that represents three coordinates
          ! (the coordinate system is not important as long as
          ! |v| = sqrt(v(1)^2+v(2)^2+v(3)^2)
          
          implicit none
          real(REAL64) v(3), energy
          
          energy = .5*me*(v(1)*v(1) + v(2)*v(2) + v(3)*v(3)) / e
        end function energy

        ! =======================================================

          subroutine new_particles(N,particle_E)
            implicit none

            ! Note: This function is here (not in the "particles" module) because
            ! it uses function collision_time(), which is initialized in module
            ! "physics". Make firstly compiles "particles" and then "physics"

            real(REAL64) particle_E
            integer(INT64) i, N
            type(particle_type), pointer :: new_p

            if (N .le. 0) return

            if (.not. associated(first_particle)) then
              allocate(first_particle)
              new_p => first_particle
              new_p%E_eV = particle_E
              new_p%tcs = collision_time()
            else
              allocate(last_particle%next_particle)
              new_p => last_particle%next_particle
              new_p%E_eV = particle_E
              new_p%tcs = collision_time()
            end if

            do i=2, N
              allocate(new_p%next_particle)
              new_p => new_p%next_particle
              new_p%E_eV = particle_E
              new_p%tcs = collision_time()
            end do

            new_p%next_particle => null()
            last_particle => new_p
            
          end subroutine new_particles

          ! =======================================================

      end module physics
