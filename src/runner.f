      program run_simulation
        use, intrinsic :: iso_fortran_env, only : REAL64, INT16, INT32, INT64, stdout=>output_unit

        use io, only : NDataFiles, read_program_input, clean_up_io, paramformat, padr, simulation_mode
        use physics, only : p, init_physics, NCollProc, dt, tfin, n, total_collision_frequency
        use interpolation, only : nri, init_interpolation, interpolate, clean_up_interp
        use random, only : seed_rand_0
        use particles, only : simulate, e0
        use eedf, only : init_eedfbins, Ntimes, calculate_totals, print_eedf, recycled_e
        use ratecoeffs, only : calculate_ratecoeffs_evolution, print_ratecoeffs, clean_up_ratecoeffs
        use populations, only : calculate_pops, print_pops, clean_up_pops, init_pops, Npops, neexpr

        implicit none

        ! VARIABLE DECLARATIONS       
        integer(INT64) i
        integer :: clock_start, clock_end, clock_rate
        !

        ! INITIALIZATION

          ! start timer
          call system_clock(count_rate=clock_rate)
          call system_clock(count=clock_start)

          ! get input from stdin        
          call read_program_input(simulation_mode, tfin, dt, e0, p, Npops, neexpr)
          
          ! initialize cross-section interpolations and population ODE's
          call init_interpolation()

          NCollProc = NDataFiles

          do i=1,NCollProc
            call interpolate(i)
          end do
          
          call init_pops()

          ! clean up allocated memory
          call clean_up_io()

          ! allocate space for histogram, initialize to zero
          Ntimes = int(ceiling(tfin/dt))
          call init_eedfbins(e0)

          ! initialize physics module
          call init_physics()

          write(stdout, paramformat) padr("Initial energy", 30), e0, padr("eV", 5)
          write(stdout, paramformat) padr("End time", 30), tfin*1e9, padr("ns", 5)
          write(stdout, paramformat) padr("Time step", 30), dt*1e9, padr("ns", 5)
          write(stdout, paramformat) padr("Gas pressure", 30) , p, padr("Torr", 5)
          write(stdout, paramformat) padr("Gas density", 30) , n*1e-6, padr("cm^-3", 5)
          write(stdout, '((A),(A))') padr("Electron density expression", 30), neexpr

          ! seed pseudo-random number generator
          call seed_rand_0()

        ! SIMULATION

          ! Main simulation, @particles
          call simulate()

        ! POST PROCESSING

          call print_eedf()

          ! calculate rate coefficients
          call calculate_ratecoeffs_evolution()
          call print_ratecoeffs()

          call calculate_pops()
          call print_pops()

          call clean_up_ratecoeffs()
          call clean_up_pops()

          ! stop timer
          call system_clock(count=clock_end)
          write(stdout,paramformat) padr('It took', 30), (clock_end-clock_start)/real(clock_rate,REAL64), padr("s", 5)
        ! CLEAN UP
          !  deallocate(bins)g
          call clean_up_interp()
        end program run_simulation
