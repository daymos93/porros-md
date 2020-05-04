program md
    use ziggurat !es modulo de fortran (coleccion de rutinas o variables)
    use vars

    implicit none !nada implicito, defino todas las variables
    logical :: es !definicion de variable logica (type::var_name)
    integer :: seed,i,j,k,Nt 


    ! Lee cantidad de part. (N) y lado del cubo (L) de input.dat

        open(unit = 10, file = 'input.dat', status = 'old')
        read(10,*)
        read(10,*) L, N, N1, Nt, dt, rc, temp, gama
        close(10)

        sigma_mat(1,1) = 1.
        sigma_mat(1,2) = 1.
		sigma_mat(2,1) = 1.
        sigma_mat(2,2) = 1.

        epsi_mat(1,1) = 1.
        epsi_mat(1,2) = 0.5
		epsi_mat(2,1) = 0.5
        epsi_mat(2,2) = 0.

        a_wall(:,:) = epsi_mat(:,:)
        sigma_wall(:) = sigma_mat(1,:)

        N_wall = int(2*0.4*L*L)
        allocate(r(3,N+N_wall),v(3,N+N_wall),f(3,N+N_wall), ptype(N+N_wall), eforce(N1))


    ![NO TOCAR] Inicializa generador de número random

        inquire(file='seed.dat',exist=es)!pregunta al SO si seed.dat existe y guarda la rta logica en es
        if(es) then
            open(unit=10,file='seed.dat',status='old')
            read(10,*) seed
            close(10)
            print *,"  * Leyendo semilla de archivo seed.dat"
        else
            seed = 24583490
        end if

        call zigset(seed)

    ![FIN NO TOCAR]   


    ! Verifica si existe o genere config inicial

        inquire(file='matrix.dat',exist=es)!pregunta al SO si matrix.dat existe y guarda la rta logica en es
        if(es) then
            open(unit=10,file='matrix.dat',status='old')
            read(10,*) r(:,:), v(:,:)
            close(10)
            print*, "  * Usando matriz OLD" 
        else
            do i=1,3
                do j=1,N1
                    ptype(j)=1
                    r(i,j) = uni()*L
                    v(i,j) = rnor()                    
                end do
			do j=N1+1,N
                    ptype(j)=2
                    r(i,j) = uni()*L
                    v(i,j) = 0.                   
            end do

            


            end do

            !!!!rugosidad de la pared con particulas typo dos

            

            do j=N+1,N+N_wall
                    ptype(j)=2
                    r(1,j) = uni()*L
                    r(2,j) = uni()*L
                    r(3,j) = nint(uni())*L
                    v(i,j) = 0.                   
            end do

        end if
 
    call force() 
 	call external_force ()
	call wall93 (2) 
    call measurements() 

    open(unit=13,file='movie.vtf',status='unknown')
    write(13,*) '# STRUCTURE BLOCK'
    write(13,*) '# define the default atom'
    write(13,'(a,i1,a)') 'atom 0:',N1-1,'      radius 0.5 name S' ! idigits
    write(13,'(a,i2,a,i2,a)') 'atom ',N1,':',N-1,'      radius 0.5 name O' ! idigits
    write(13,'(a,i2,a,i2,a)') 'atom ',N,':',N+N_wall-1,'      radius 0.5 name O' ! idigits
    !write(13,*) 'atom 10:20      radius 0.5 name O'



    write(13,*) 'timestep'
    write(13,*) 

    open(unit=10,file='evol.dat',status='unknown')
    write(10,*) '#	t,	U,	K,	E,  T,	P'
    write(10,*) 0.0, ener/N, Kin, Etotal, T, P

    print *, '  * Take a cup of coffee and come back next week :)'

    do j=1,Nt

        call verlet_positions() 

        call force() 
 		call external_force ()
	    call wall93 (2) 
        call verlet_velocities() 


         
      
        if (mod(j,500).eq.0) then
            call measurements()
            write(10,*) dt*j, ener/N, Kin, Etotal, T, P
            
        	do i=1,N+N_wall 
        		write(13,*) r(:,i)
        	end do
        	write(13,*) 'timestep'
        	write(13,*) 
        end if

    end do

    close(10)
    close(13)


    !Escribe la ultima config

        !open(unit=10,file='matrix.dat',status='unknown') 
        !write(10,*) r(:,:), v(:,:)
        !close(10)


    print*, "  ************** END **************"



    !! FIN FIN edicion
    !! 
    ![No TOCAR]
    ! Escribir la última semilla para continuar con la cadena de numeros aleatorios 

        open(unit=10,file='seed.dat',status='unknown')
        seed = shr3() 
        write(10,*) seed
        close(10)

    ![FIN no Tocar]        


    end program md

