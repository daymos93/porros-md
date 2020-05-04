subroutine thermal_walls(mode)
use commons
use ziggurat

implicit none
integer, intent(in) :: mode
real(kind=8) :: v_new(3)
real(kind=8) :: fac

! CURSO intro_sims
! * r0: global, posición de las partículas. Cambiar la variable de posición de sus programas
! * z_space_wall: altura en z de la muestra. Ancho del canal, D
! * thermal_skin: longitud para dimensionar la "zona pared". Valor usual: thermal_skin= 1.2 sigma 

    select case(mode)

    case(2) ! *** Thermal walls in top and bottom walls ***

!     print *,kb,thermal_skin,top_thermal_wall,bottom_thermal_wall !deb

    do i_part = 1 , n_mon_tot
 
! Top wall

!NOTE: Changes velocity if the particle is above interwall spacing -thermal skin  AND going towards the wall 

        if (r0(3,i_part)>(z_space_wall-thermal_skin) .and. (v(3,i_part)>0.0)) then 
            fac= sqrt(kb*top_thermal_wall*inv_mass(i_part))
            v_new(:) = fac*(/rnor(),rnor(),-sqrt(2.)*sqrt(-log( uni() ) )/)
            v(:,i_part) = v_new(:)
        endif

! Bottom wall

       if ( (r0(3,i_part)<thermal_skin) .and. (v(3,i_part)<0.0)) then
           fac= sqrt(kb*bottom_thermal_wall*inv_mass(i_part))
           v_new(:) = fac*(/rnor(),rnor(),sqrt(2.)*sqrt(-log( uni() ) )/)
           v(:,i_part) = v_new(:)
       end if
  end do

end select 

end subroutine thermal_walls
