

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine force() 
    
    use vars
    integer :: i, j, k
    real(kind=8)::  drx,dry,drz,rmod,fpx,fpy,fpz

     f(:,:) = 0.0
     ener = 0.0
     vir = 0.0

    do i=1,N+N_wall-1
        do j=i+1,N+N_wall

            drx = r(1,j)-r(1,i)
            drx = drx-l*int(2*drx/l)   
            dry = r(2,j)-r(2,i)
            dry = dry-l*int(2*dry/l) 
            drz = r(3,j)-r(3,i)
            drz = drz-l*int(2*drz/l) 

            rmod = sqrt(drx**2+dry**2+drz**2)

            if (rmod<=rc) then
				sigma = sigma_mat(ptype(i),ptype(j))
				epsi = epsi_mat(ptype(i),ptype(j))
		        vrc=4*epsi*((sigma/rc)**12-(sigma/rc)**6)
                ener = ener + 4*epsi*((sigma/rmod)**12-(sigma/rmod)**6) - vrc

                fpx = 4*epsi*((12*sigma**12*rmod**(-13))-(6*sigma**6*rmod**(-7)))*drx/rmod
                fpy = 4*epsi*((12*sigma**12*rmod**(-13))-(6*sigma**6*rmod**(-7)))*dry/rmod
                fpz = 4*epsi*((12*sigma**12*rmod**(-13))-(6*sigma**6*rmod**(-7)))*drz/rmod
                
                f(1,j) = f(1,j)+fpx
                f(1,i) = f(1,i)-fpx
                f(2,j) = f(2,j)+fpy
                f(2,i) = f(2,i)-fpy
                f(3,j) = f(3,j)+fpz
                f(3,i) = f(3,i)-fpz

!call dpd_force

				vir = vir + (fpx*drx + fpy*dry + fpz*drz)
				
           	end if
        
        end do
   
    end do


  ! dpd  
    do i=1,N1-1
        do j=i+1,N1

            drx = r(1,j)-r(1,i)
            drx = drx-l*int(2*drx/l)   
            dry = r(2,j)-r(2,i)
            dry = dry-l*int(2*dry/l) 
            drz = r(3,j)-r(3,i)
            drz = drz-l*int(2*drz/l) 

            rmod = sqrt(drx**2+dry**2+drz**2)

            if (rmod<=rc) then
                  !call dpd_force

				
           	end if
        
        end do
    end do




    return

    end subroutine force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine verlet_positions() 
    use vars 

    integer :: i,j,k
     
        do j=1,N1
        
           do k=1,3

             r(k,j) = r(k,j) + v(k,j)*dt + 0.5*f(k,j)*dt**2
             v(k,j) = v(k,j) + 0.5*f(k,j)*dt

    		
    		end do

           do k=1,2
           	if (r(k,j) > l) r(k,j) = r(k,j)-l
            if (r(k,j) < 0) r(k,j) = r(k,j)+l
           end do
        end do

    
    return
    end subroutine verlet_positions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
subroutine verlet_velocities() 
    use vars
    integer :: i,j,k

        do j=1,N1

            do k=1,3
              v(k,j) = v(k,j) + 0.5*f(k,j)*dt
            end do
        end do

        do j=N1+1,N+N_wall
			v(:,j) = 0.
		end do
    return
    end subroutine verlet_velocities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine measurements() 
    
    use vars
    integer :: i,j,k
    real(kind=8)::  vmod,kinetic
	       
	kinetic = 0.0
	Etotal = 0.0
	Kin = 0.0
	P = 0.0
	       
        do j=1,N
        
           	vmod = v(1,j)**2+v(2,j)**2+v(3,j)**2
           	kinetic = kinetic + 0.5*vmod

        end do

        Kin = kinetic/N

        T = 2*kinetic/(3*N-3)
        !T = Kin*0.67

        P = (N*T + 0.33*vir/N)/(L*L*L)
        
        Etotal = Kin + ener/N
    
    return
    end subroutine measurements

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lgv_force() 
    use ziggurat
    use vars
    real(kind=8) :: f_visc(3,N)
    integer :: i, j, k

!f_visc(:,:) = 0

do i=1,N
   
    f_visc(1,i) = -1*gama*v(1,i)
    f_visc(2,i) = -1*gama*v(2,i)
    f_visc(3,i) = -1*gama*v(3,i)

    f(1,i) = f(1,i) + f_visc(1,i) + sqrt(2*temp*gama/dt)*rnor()
    f(2,i) = f(2,i) + f_visc(2,i) + sqrt(2*temp*gama/dt)*rnor()
    f(3,i) = f(3,i) + f_visc(3,i) + sqrt(2*temp*gama/dt)*rnor()

end do

return
    end subroutine lgv_force


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine external_force()

    use vars
    implicit none
    integer :: i, j, k

    eforce(:) = 10.
    !eforce(1) = 4.0 
 
    do i = 1,N1
        f(1,i) = f(1,i) + eforce(i) 
    end do 

end subroutine external_force



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 subroutine wall93(inter_type)
        
! Computes fluid-wall interactions:
!
! r: posiciones de las particulas
! z_space_wall: altura de la muestra (ancho del canal,D)
! a_wall: epsilon partícula pared (vector con dos compo porque contempla paredes distintas TOP y BOTTOM)
! sigma_wall: sigma de interacción pared=partícula

      use vars ; implicit none
      real (kind=8) :: inv_z
      integer, intent(in) :: inter_type
      logical, parameter :: debug_fw=.false.

!
      v_fluid_wall = 0.

      select case(inter_type)




      case(2)     !  (1/z)^9-(1/z)^3 interaction



          do i_part = 1,N1

! *** Bottom wall interaction
              i_type = ptype(i_part)
            
            
              inv_z = 1./r(3,i_part) 
              r_dummy = sigma_wall(i_type)*inv_z
              v_fluid_wall = v_fluid_wall + abs(a_wall(2,i_type))*r_dummy**9 - a_wall(2,i_type)*r_dummy**3    
              f(3,i_part) = f(3,i_part) + 9.*abs(a_wall(2,i_type))*(sigma_wall(i_type))**9*(inv_z)**10 
              f(3,i_part) = f(3,i_part) - 3.*a_wall(2,i_type)*(sigma_wall(i_type))**3*(inv_z)**4                              

! ***  Top wall interaction (possibly DIFFERENT from bottom wall)


              inv_z = 1./(L-r(3,i_part))
              r_dummy = sigma_wall(i_type)*inv_z

              v_fluid_wall = v_fluid_wall + abs(a_wall(1,i_type))*r_dummy**9 - a_wall(1,i_type)*r_dummy**3        

!        note: different sign in top wall

              f(3,i_part) = f(3,i_part)    - 9.*abs(a_wall(1,i_type))*(sigma_wall(i_type))**9*(inv_z)**10  
              f(3,i_part) = f(3,i_part)    + 3*a_wall(1,i_type)*(sigma_wall(i_type))**3*(inv_z)**4

          end do

      case   default
          print *," sub FLUID_WALL: inter_type value must be 1: LJ or 2: int -1/z^3 + 1/z^9 3: drop 4: hard top and lower walls "
          print*, "Change it! Stopping here "
          stop
      end select


         end subroutine wall93

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!