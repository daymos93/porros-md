
      subroutine wall93(inter_type)
        
! Computes fluid-wall interactions:
!
! r0: posiciones de las particulas
! z_space_wall: altura de la muestra (ancho del canal,D)
! a_wall: epsilon partícula pared (vector con dos compo porque contempla paredes distintas TOP y BOTTOM)
! sigma_wall: sigma de interacción pared=partícula

      use commons ; implicit none
      real (kind=8) :: inv_z
      integer, intent(in) :: inter_type
      logical, parameter :: debug_fw=.false.

!
      v_fluid_wall = 0.

      select case(inter_type)




      case(2)     !  (1/z)^9-(1/z)^3 interaction



          do i_part = 1,n_mon_tot

! *** Bottom wall interaction
              i_type = a_type(i_part)


              inv_z = 1./r0(3,i_part) 
              r_dummy = sigma_wall(i_type)*inv_z
              v_fluid_wall = v_fluid_wall + abs(a_wall(2,i_type))*r_dummy**9 - a_wall(2,i_type)*r_dummy**3    
              force(3,i_part) = force(3,i_part) +   9.*abs(a_wall(i_type))*(sigma_wall(i_type))**9*(inv_z)**10 
              force(3,i_part) = force(3,i_part)   - 3.*a_wall(i_type)*(sigma_wall(i_type))**3*(inv_z)**4                              

! ***  Top wall interaction (possibly DIFFERENT from bottom wall)


              inv_z = 1./(z_space_wall-r0(3,i_part))
              r_dummy = sigma_wall(i_type)*inv_z

              v_fluid_wall = v_fluid_wall + abs(a_wall(1,i_type))*r_dummy**9 - a_wall(1,i_type)*r_dummy**3        

!        note: different sign in top wall

              force(3,i_part) = force(3,i_part)    - 9.*abs(a_wall(i_type))*(sigma_wall(i_type))**9*(inv_z)**10  
              force(3,i_part) = force(3,i_part)    + 3*a_wall(i_type)*(sigma_wall(i_type))**3*(inv_z)**4

          end do

      case   default
          print *," sub FLUID_WALL: inter_type value must be 1: LJ or 2: int -1/z^3 + 1/z^9 3: drop 4: hard top and lower walls "
          print*, "Change it! Stopping here "
          stop
      end select


         end subroutine wall93
