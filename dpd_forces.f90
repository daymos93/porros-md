subroutine dpd_forces(inv_sqrt_r_2)
    use commons

use ziggurat, only: rnor,uni    

    implicit none
    real(kind=8) , intent(in) :: inv_sqrt_r_2
    real(kind=8) :: g_rand,r_versor(3),delta_v(3),rrc,vec_dummy(3)
    integer :: i,j 
!
!   -------- DPD forces calculation : Version inside the force caculation
!
! WARN: check-out if works with different masses as is.

!    *** Gaussian random distribution. sig^2 = 1 , <x> = 0

! CURSO Intro Sims 
! * Se llama desde adentro del doble loop de la rutina de fuerzas conservativas. Es decir a i_part,j_part=constante
!  variables: 
! 
! * r_cut_dpd_2: cutoff de DPD. Igual a Rc de LJ normalmente
! * sig está definida en otra parte del programa como: sig = sqrt( 2.*temp*friction(1)/dt ) 
! * i_part y j_part son globales y valen lo mismo que en el punto de la rutina de fuerzas en que son llamada
! * delta_r viene de forces


          g_rand =  rnor()      

         

! --- Weight functions

             w_d = (1. - sqrt(r_2/r_cut_dpd_2) )

          w_r = w_d
          w_d = w_d*w_d


! ----  Random force computation:

          r_versor(:) = delta_r(:)*inv_sqrt_r_2
         vec_dummy(:) = sig*w_r*g_rand*r_versor(:) 


      force(:,i_part) =  force(:,i_part) + vec_dummy(:)   
      force(:,j_part) =  force(:,j_part) - vec_dummy(:)   

     
! Dissipative force computation

      delta_v(:) = v(:,i_part) - v(:,j_part)

      r_dummy = delta_v(1)*r_versor(1) + delta_v(2)*r_versor(2) + delta_v(3)*r_versor(3)

      vec_dummy(:) = -1.*friction(1)*w_d*r_dummy*r_versor(:)
        

      force(:,i_part) =  force(:,i_part) + vec_dummy(:)   
      force(:,j_part) =  force(:,j_part) - vec_dummy(:)   



end subroutine dpd_forces
