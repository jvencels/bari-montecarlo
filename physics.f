      function collision_time(eV)
        ! This function simply uses the nonrandom version, and then adds
        ! a random factor of -log(rand()) 
        
        implicit none
          ! function input parameters
          real*8 eV
          ! parameter: air molecular density
          real*8 n
          parameter(n=3e21)
          ! return value and auxiliary function
          real*8 collision_time, nonrandom_collision_time
          
          collision_time = -log(rand())*nonrandom_collision_time(n,eV)
      end
      
      function nonrandom_collision_time(n, eV)
        ! Calculates the collision time for a particle of energy eV,
        ! given a particle density of the medium with which it collides
        ! of n m^(-3).
        
        implicit none
          ! function input parameters
          real*8 n, eV
          ! function return value, and helper functions
          real*8 nonrandom_collision_time, cross_section, velocity
          ! intermediate quantity
          real*8 collision_frequency

          collision_frequency = n*velocity(eV)*cross_section(eV)
          nonrandom_collision_time = 1.0 / collision_frequency
      end
      
      function cross_section(eV)
        ! Calculates the collision cross-section for electrons in a
        ! nitrogen gas, based on a fit to data from 
        ! www.lxcat.laplace.univ-tlse.fr
        
        implicit none
          ! Input argument and return value
          real*8 eV, cross_section
        
          ! Parameters from curve fitting in Matlab
          real q0, q1, q2, q3, q4, l0, l1, E0
          parameter(q0=-2.2007449E-20,q1=2.855227E-21,q2=-4.6951205E-23)
          parameter(q3=3.031989421E-25,q4=-6.335964E-28)
          parameter(l0=1.21665367425970e-20, l1=1.86036008409668e-18)
          parameter(E0=70)

          if (eV.LT.E0) then
              cross_section=q0+q1*eV+q2*eV*eV+q3*eV*eV*eV+q4*eV*eV*eV*eV
          else 
              cross_section=l0+l1/eV
          endif
      end

      function velocity(eV)
        ! Calculates the velocity in m/s for an electron with kinetic
        ! energy eV
        
        implicit none
          real*8 eV, velocity
          real me, e
          parameter(me=9.11E-31,e = 1.602176E-19)

          velocity = sqrt(2*eV*e/me)
      end

      function energy(v)
        ! Calculates the kinetic energy in eV for an electron which 
        ! moves at a velocity v (in m/s)
        !
        ! v is a 3-array that represents three coordinates
        ! (the coordinate system is not important as long as
        ! |v| = sqrt(v(1)^2+v(2)^2+v(3)^2)
        
        implicit none
          real*8 v(3), energy
          real me, e
          parameter(me=9.11E-31,e = 1.602176E-19)
          
          energy = .5*me*(v(1)*v(1) + v(2)*v(2) + v(3)*v(3)) / e
      end