module vars

    real(kind=8),allocatable:: r(:,:), v(:,:), f(:,:), eforce(:)
    integer,allocatable:: ptype(:)
    real(kind=8)             :: dt, L, Rc, ener, vrc, T, temp, vir
    real(kind=8)             :: Kin, Etotal, P
    integer                  :: N, N1, i_type, i_part, N_wall
    real(kind=8)             :: sigma, epsi, gama, sigma_mat(2,2), epsi_mat(2,2)
    real(kind=8)             :: r_dummy, v_fluid_wall, sigma_wall(2), a_wall(2,2)

end module vars
