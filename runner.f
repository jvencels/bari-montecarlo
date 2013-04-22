      program interptest
        use precision

        use mpi, only : rnk, nproc, init_mpi, barrier, share_data, finalize_mpi, wtime
        use io, only : raw, e0raw, e1raw, fnames, NDataFiles, read_program_input, read_interpolation_data, clean_up_io
        use physics, only : p, eI, NCollProc, clean_up_physics
        use interpolation, only : interp, nri, interp_min, interp_max, init_interpolation, interpolate, clean_up_interp
        use random, only : seed_rand_0
        use single_particle, only : onepart, e0
        use histogram

        implicit none

910   format('# ', a, ' took ', Es10.3, ' s')

        ! VARIABLE DECLARATIONS       
          integer(lkind) Nruns, i
          real(rkind) tfin, dt
          integer :: Ntimes
          REAL(REAL64) :: t_init, t_sim, t_hist, t_all
        !

        ! INITIALIZATION

          ! initialize mpi
            call init_mpi()
            !
          !

          ! start init and all timers
            call barrier()
            if (rnk.eq.0) then
              t_init = wtime()
              t_all = wtime()
            end if
          !

          ! master thread: read indata and interpolate cross-sections
            if (rnk.eq.0) then
              ! get input from stdin        
              call read_program_input(Nruns, tfin, dt, e0, p, eI)
              ! initialize e vector of interpolation
              call init_interpolation(e0)

              ! interplolate each data series
              do i=1,NDataFiles
                call read_interpolation_data(fnames(i))
                call interpolate(i, e0raw(i), e1raw(i))
              end do

              call clean_up_io()

              NCollProc = NDataFiles

            end if
          !

          ! share indata and interpolations between all processes
          
            call share_data(Nruns,tfin,dt,e0,p,eI,nri,NCollProc,interp,interp_min,interp_max)
            !print *, "# interp_min/max: ", interp_min, interp_max
            !
          !

          ! allocate space for histogram, initialize to zero
            Ntimes = int(ceiling(tfin/dt))
            call init_bins(Ntimes)
          !

          ! seed pseudo-random number generator
            call seed_rand_0()
            !
          !

          ! stop init timer
            call barrier()
            if (rnk.eq.0) then
              t_init = wtime() - t_init
              write(*,910), 'Initialization', t_init
            end if
          !
        !

        ! SIMULATION

          ! start simulation timer
            call barrier()
            if (rnk.eq.0) then
              t_sim = wtime()
            end if
          !

          ! split work on all threads
            do i=1, Nruns/nproc
              call onepart(dt, tfin, Ntimes)
            end do
          !
          
          ! master thread: run remaining, if Nruns/nproc is not even
            if (rnk.eq.0) then
              do i=1, Nruns-nproc*(Nruns/nproc)
                call onepart(dt, tfin, Ntimes)
              end do
            end if
          !

          ! stop simulation timer
            call barrier()
            if (rnk.eq.0) then
              t_sim = wtime() - t_sim
              write(*,910), 'Simulation', t_sim
            end if
        !

        ! POST PROCESSING
          ! stop postprocessing timer
            call barrier()
            if (rnk.eq.0) then
              t_hist = wtime() - t_hist
            end if
          !

          ! summarize and output histogram
            call calculate_totals(Ntimes, dt, tfin, e0)
            !
          !

          ! stop postprocessing and all timers
            call barrier()
            if (rnk.eq.0) then
              t_hist = wtime() - t_hist
              t_all = wtime() - t_all
              write(*,910) 'Postprocessing', t_hist
              write(*,910) 'Everything', t_all
            end if
          !
        !

        ! CLEAN UP
          !  deallocate(bins)
          call clean_up_interp()
          call clean_up_physics()

          call finalize_mpi()
        !
      end program interptest