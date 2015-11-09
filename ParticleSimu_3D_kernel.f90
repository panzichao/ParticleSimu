!dec$ if .not. defined(__ParticleSimu_3D_kernel_f90)
!dec$ define __ParticleSimu_3D_kernel_f90

module ParticleSimu_3D_kernel
    implicit none
    contains
        
        ! generate the ref radius for ellipsoid, and also the sphere itself
        ! ref_sphere_radius     :   [total_particle_num]
        ! container_volume      :   volume of the containder
        subroutine get_ref_sphere_radius(ref_sphere_radius,container_volume)
            use ParticleSimu_global_parameter
            implicit none
            real(kind=8),allocatable    ::  ref_sphere_radius(:)
            real(kind=8),intent(in)     ::  container_volume
            character(len=100)          ::  sphere_rad_file
            real(kind=8),allocatable    ::  particle_needed_vol(:), particle_real_vol(:)
            real(kind=8)                ::  last_particle_vol, r(1), particle_real_vol_in_current_sieve, low_value, high_value, temp
            integer(kind=4)             ::  total_particle_num, i, j, k
            
            if (allocated(particle_needed_vol)) deallocate(particle_needed_vol)
            if (allocated(particle_real_vol)) deallocate(particle_real_vol)
            allocate(particle_needed_vol(total_sieve_num-1))
            allocate(particle_real_vol(total_sieve_num-1))
            
            do i = 1, total_sieve_num-1
                particle_needed_vol(i) = (gradation(2,i+1) - gradation(2,i)) * particle_occupy_percentage * container_volume
            end do
            
            ! randomize
            call random_seed()
 
            sphere_rad_file='sphere_radius.tmp'   
            open(unit=100,file=sphere_rad_file,form='formatted',status='replace',action='write')
            particle_real_vol = 0.0
            total_particle_num = 0
            do i = total_sieve_num-1, 1, -1
                particle_real_vol_in_current_sieve = 0.0 
                do while (particle_real_vol(i) .lt. particle_needed_vol(i))
                    ! generate the circle
                    call random_number(r)
                    low_value = gradation(1,i) / 2.0
                    high_value = gradation(1,i+1) / 2.0
                    r(1) = r(1) * (high_value - low_value) + low_value
                    particle_real_vol_in_current_sieve = particle_real_vol_in_current_sieve + r(1) ** 3.0 * 3.14159 * 4.0 / 3.0
                    if (particle_real_vol_in_current_sieve .gt. particle_needed_vol(i)) then
                        last_particle_vol = particle_needed_vol(i) - particle_real_vol(i)
                        r(1) = (last_particle_vol * 3.0 / 4.0 / 3.14159) ** (1.0 / 3.0)
                    end if
                    write(100,'(e15.5)') r(1)
                    total_particle_num = total_particle_num + 1
                    particle_real_vol(i) = particle_real_vol(i) + r(1) ** 3.0 * 3.14159 * 4.0 / 3.0
                end do
            end do
            close(unit=100)

            ! read radius from file 
            if (allocated(ref_sphere_radius)) deallocate(ref_sphere_radius)
            allocate(ref_sphere_radius(total_particle_num))
            ref_sphere_radius = 0.0
            open(unit=100,file=sphere_rad_file,form='formatted',status='old',action='read')
            do i = 1, total_particle_num
                read(100,*) ref_sphere_radius(i)
            end do
            close(unit=100,status='delete')

            ! sort
            do i = 1, total_particle_num - 1
                k = i
                do j = i+1, total_particle_num
                    if (ref_sphere_radius(k) .lt. ref_sphere_radius(j)) k = j
                end do
                if (k .ne. i) then
                    temp = ref_sphere_radius(k)
                    ref_sphere_radius(k) = ref_sphere_radius(i)
                    ref_sphere_radius(i) = temp
                end if
            end do
        end subroutine

        ! generate spherical particles
        ! spherical_particle_info    size:   [total_particle_num,4]
        ! each row contains:     radius      x_pos       y_pos      z_pos
        subroutine sphere_without_periodic_bc(spherical_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  spherical_particle_info(:,:)
            real(kind=8),allocatable    ::  ref_sphere_radius(:)
            real(kind=8)                ::  r(1), low_value, high_value, distance, model_x_len, model_y_len, model_z_len
            integer(kind=4)             ::  total_particle_num, i, j, progress
            logical                     ::  is_overlapped

            model_x_len = container_size1
            model_y_len = container_size2
            model_z_len = container_size3

            call random_seed()

            call get_ref_sphere_radius(ref_sphere_radius,model_x_len*model_y_len*model_z_len)
            total_particle_num = size(ref_sphere_radius)

            if (allocated(spherical_particle_info)) deallocate(spherical_particle_info)
            allocate(spherical_particle_info(total_particle_num,4))
            spherical_particle_info = 0.0
            do i = 1, total_particle_num
                spherical_particle_info(i,1) = ref_sphere_radius(i)
            end do
            
            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            ! generate the position of circular particles
            i = 1
            do while (i .le. total_particle_num)
                is_overlapped = .false.

                call random_number(r)
                spherical_particle_info(i,2) = r(1) * (model_x_len - 2 * spherical_particle_info(i,1)) + spherical_particle_info(i,1)
                spherical_particle_info(i,2) = spherical_particle_info(i,2) - 0.5 * model_x_len

                call random_number(r)
                spherical_particle_info(i,3) = r(1) * (model_y_len - 2 * spherical_particle_info(i,1)) + spherical_particle_info(i,1)
                spherical_particle_info(i,3) = spherical_particle_info(i,3) - 0.5 * model_y_len

                call random_number(r)
                spherical_particle_info(i,4) = r(1) * (model_z_len - 2 * spherical_particle_info(i,1)) + spherical_particle_info(i,1)
                spherical_particle_info(i,4) = spherical_particle_info(i,4) - 0.5 * model_z_len

                do j = 1, i-1
                    distance = sqrt((spherical_particle_info(i,2) - spherical_particle_info(j,2)) ** 2.0 + &
                                 &  (spherical_particle_info(i,3) - spherical_particle_info(j,3)) ** 2.0 + &
                                 &  (spherical_particle_info(i,4) - spherical_particle_info(j,4)) ** 2.0)
                    if (distance .lt. (spherical_particle_info(i,1) + spherical_particle_info(j,1))) then
                        is_overlapped = .true.
                        exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    i = i + 1
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate spherical particles
        ! spherical_particle_info    size:   [total_particle_num,4]
        ! each row contains:     radius      x_pos       y_pos      z_pos
        ! for cylinderical container
        subroutine sphere_in_cylinder_container(spherical_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  spherical_particle_info(:,:)
            real(kind=8),allocatable    ::  ref_sphere_radius(:)
            real(kind=8)                ::  r(1), low_value, high_value, distance, container_radius, container_height, polar_angle, polar_radius
            integer(kind=4)             ::  total_particle_num, i, j, progress
            logical                     ::  is_overlapped, is_out_region

            container_radius = container_size1
            container_height = container_size2

            call random_seed()

            call get_ref_sphere_radius(ref_sphere_radius,container_radius**2.0*3.14159*container_height)
            total_particle_num = size(ref_sphere_radius)

            if (allocated(spherical_particle_info)) deallocate(spherical_particle_info)
            allocate(spherical_particle_info(total_particle_num,4))
            spherical_particle_info = 0.0
            do i = 1, total_particle_num
                spherical_particle_info(i,1) = ref_sphere_radius(i)
            end do
            
            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            ! generate the position of circular particles
            i = 1
            do while (i .le. total_particle_num)
                is_out_region = .false.
                is_overlapped = .false.

                ! fisrt generate position in polar coordinate
                call random_number(r)
                polar_angle = r(1) * 2.0 * 3.14159
                call random_number(r)
                low_value = -0.5 * container_height + spherical_particle_info(i,1)
                high_value = 0.5 * container_height - spherical_particle_info(i,1)
                spherical_particle_info(i,4) = r(1) * (high_value - low_value) + low_value
                call random_number(r)
                low_value = 0.0
                high_value = container_radius - spherical_particle_info(i,1)
                polar_radius = r(1) * (high_value - low_value) + low_value

                spherical_particle_info(i,2) = polar_radius * cos(polar_angle)
                spherical_particle_info(i,3) = polar_radius * sin(polar_angle)

                do j = 1, i-1
                    distance = sqrt((spherical_particle_info(i,2) - spherical_particle_info(j,2)) ** 2.0 + &
                                 &  (spherical_particle_info(i,3) - spherical_particle_info(j,3)) ** 2.0 + &
                                 &  (spherical_particle_info(i,4) - spherical_particle_info(j,4)) ** 2.0)
                    if (distance .lt. (spherical_particle_info(i,1) + spherical_particle_info(j,1))) then
                        is_overlapped = .true.
                        exit
                    end if
                end do

                if ((is_overlapped .eq. .false.) .and. (is_out_region .eq. .false.)) then
                    i = i + 1
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate spherical particles
        ! spherical_particle_info    size:   [total_particle_num,4]
        ! each row contains:     radius      x_pos       y_pos      z_pos
        subroutine sphere_with_periodic_bc(spherical_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  spherical_particle_info(:,:)
            real(kind=8),allocatable    ::  real_sphere_info(:,:), ghost_sphere_info(:,:)
            real(kind=8),allocatable    ::  ref_sphere_radius(:)
            real(kind=8)                ::  r(1), distance, ghost_sphere_pos(7,3), x_pos, y_pos, z_pos, model_x_len, model_y_len, model_z_len
            integer(kind=4)             ::  total_particle_num, i, j, k, flag, progress, ghost_sphere_status(7), total_real_sphere_num, total_ghost_sphere_num
            logical                     ::  is_overlapped, is_out_region

            model_x_len = container_size1
            model_y_len = container_size2
            model_z_len = container_size3

            call random_seed()

            call get_ref_sphere_radius(ref_sphere_radius,model_x_len*model_y_len*model_z_len)
            total_real_sphere_num = size(ref_sphere_radius)

            if (allocated(ghost_sphere_info)) deallocate(ghost_sphere_info)
            allocate(ghost_sphere_info(total_real_sphere_num*7,4))
            ghost_sphere_info = 0.0

            if (allocated(real_sphere_info)) deallocate(real_sphere_info)
            allocate(real_sphere_info(total_real_sphere_num,4))
            real_sphere_info = 0.0
            do i = 1, total_real_sphere_num
                real_sphere_info(i,1) = ref_sphere_radius(i)
            end do
            
            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a17)') 'Particle shape done.','(',total_real_sphere_num,'real particles).'

            ! generate the position of circular particles
            total_ghost_sphere_num = 0
            i = 1
            do while (i .le. total_real_sphere_num)
                is_overlapped = .false.

                call random_number(r)
                real_sphere_info(i,2) = r(1) * model_x_len - 0.5 * model_x_len

                call random_number(r)
                real_sphere_info(i,3) = r(1) * model_y_len - 0.5 * model_y_len

                call random_number(r)
                real_sphere_info(i,4) = r(1) * model_z_len - 0.5 * model_z_len

                ! check whether the current real sphere is intersected by the boundary of model
                ! if yes, generate ghost spheres
                is_out_region = .false.
                if ((-0.5 * model_x_len - real_sphere_info(i,2) ** 2.0 .lt. real_sphere_info(i,1) ** 2.0) .or. &
                  & ((0.5 * model_x_len - real_sphere_info(i,2)) ** 2.0 .lt. real_sphere_info(i,1) ** 2.0)) is_out_region = .true.

                if ((-0.5 * model_y_len - real_sphere_info(i,3) ** 2.0 .lt. real_sphere_info(i,1) ** 2.0) .or. &
                  & ((0.5 * model_y_len - real_sphere_info(i,3)) ** 2.0 .lt. real_sphere_info(i,1) ** 2.0)) is_out_region = .true.

                if ((-0.5 * model_z_len - real_sphere_info(i,4) ** 2.0 .lt. real_sphere_info(i,1) ** 2.0) .or. &
                  & ((0.5 * model_z_len - real_sphere_info(i,4)) ** 2.0 .lt. real_sphere_info(i,1) ** 2.0)) is_out_region = .true.

                ghost_sphere_status = 0
                if (is_out_region .eq. .true.) then
                    ! generate ghost particle based on the current real particle
                    ! and detect whether the ghost particle is intersected with 
                    ! the region
                    ! for sphere, seven possible ghost particles exist
                    if (real_sphere_info(i,2) .gt. 0.0) then
                        x_pos = real_sphere_info(i,2) - model_x_len
                    else
                        x_pos = real_sphere_info(i,2) + model_x_len
                    end if
                    if (real_sphere_info(i,3) .gt. 0.0) then
                        y_pos = real_sphere_info(i,3) - model_y_len
                    else
                        y_pos = real_sphere_info(i,3) + model_y_len
                    end if
                    if (real_sphere_info(i,4) .gt. 0.0) then
                        z_pos = real_sphere_info(i,4) - model_z_len
                    else
                        z_pos = real_sphere_info(i,4) + model_z_len
                    end if

                    ghost_sphere_pos(1,1) = x_pos
                    ghost_sphere_pos(1,2) = real_sphere_info(i,3)
                    ghost_sphere_pos(1,3) = real_sphere_info(i,4)

                    ghost_sphere_pos(2,1) = x_pos
                    ghost_sphere_pos(2,2) = y_pos
                    ghost_sphere_pos(2,3) = real_sphere_info(i,4)

                    ghost_sphere_pos(3,1) = x_pos
                    ghost_sphere_pos(3,2) = y_pos
                    ghost_sphere_pos(3,3) = z_pos

                    ghost_sphere_pos(4,1) = x_pos
                    ghost_sphere_pos(4,2) = real_sphere_info(i,3)
                    ghost_sphere_pos(4,3) = z_pos

                    ghost_sphere_pos(5,1) = real_sphere_info(i,2)
                    ghost_sphere_pos(5,2) = y_pos
                    ghost_sphere_pos(5,3) = real_sphere_info(i,4)

                    ghost_sphere_pos(6,1) = real_sphere_info(i,2)
                    ghost_sphere_pos(6,2) = y_pos
                    ghost_sphere_pos(6,3) = z_pos

                    ghost_sphere_pos(7,1) = real_sphere_info(i,2)
                    ghost_sphere_pos(7,2) = real_sphere_info(i,3)
                    ghost_sphere_pos(7,3) = z_pos

                    do j = 1, 7
                        flag = 0
                        if ((ghost_sphere_pos(j,1) .gt. -0.5 * model_x_len - real_sphere_info(i,1)) .and. (ghost_sphere_pos(j,1) .lt. 0.5 * model_x_len + real_sphere_info(i,1))) then
                            flag = flag + 1
                        end if
                        if ((ghost_sphere_pos(j,2) .gt. -0.5 * model_y_len - real_sphere_info(i,1)) .and. (ghost_sphere_pos(j,2) .lt. 0.5 * model_y_len + real_sphere_info(i,1))) then
                            flag = flag + 1
                        end if
                        if ((ghost_sphere_pos(j,3) .gt. -0.5 * model_z_len - real_sphere_info(i,1)) .and. (ghost_sphere_pos(j,3) .lt. 0.5 * model_z_len + real_sphere_info(i,1))) then
                            flag = flag + 1
                        end if
                        if (flag .eq. 3) ghost_sphere_status(j) = 1
                    end do
                end if
                
                ! overlap, new real particle vs old real particles
                do j = 1, i-1
                    distance = sqrt((real_sphere_info(i,2) - real_sphere_info(j,2)) ** 2.0 + &
                                  & (real_sphere_info(i,3) - real_sphere_info(j,3)) ** 2.0 + &
                                  & (real_sphere_info(i,4) - real_sphere_info(j,4)) ** 2.0)
                    if (distance .lt. (real_sphere_info(i,1) + real_sphere_info(j,1))) then
                        is_overlapped = .true.
                        exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new real particle vs old ghost particles
                do j = 1, total_ghost_sphere_num
                    distance = sqrt((real_sphere_info(i,2) - ghost_sphere_info(j,2)) ** 2.0 + &
                                  & (real_sphere_info(i,3) - ghost_sphere_info(j,3)) ** 2.0 + &
                                  & (real_sphere_info(i,4) - ghost_sphere_info(j,4)) ** 2.0)
                    if (distance .lt. (real_sphere_info(i,1) + ghost_sphere_info(j,1))) then
                        is_overlapped = .true.
                        exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old real particles
                do j = 1, 7
                    if (ghost_sphere_status(j) .eq. 1) then
                        do k = 1, i-1
                            distance = sqrt((ghost_sphere_pos(j,1) - real_sphere_info(k,2)) ** 2.0 + &
                                          & (ghost_sphere_pos(j,2) - real_sphere_info(k,3)) ** 2.0 + &
                                          & (ghost_sphere_pos(j,3) - real_sphere_info(k,4)) ** 2.0)
                            if (distance .lt. (real_sphere_info(i,1) + real_sphere_info(k,1))) then
                                is_overlapped = .true.
                                exit
                            end if
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old ghost particles
                do j = 1, 7
                    if (ghost_sphere_status(j) .eq. 1) then
                        do k = 1, total_ghost_sphere_num
                            distance = sqrt((ghost_sphere_pos(j,1) - ghost_sphere_info(k,2)) ** 2.0 + &
                                          & (ghost_sphere_pos(j,2) - ghost_sphere_info(k,3)) ** 2.0 + &
                                          & (ghost_sphere_pos(j,3) - ghost_sphere_info(k,4)) ** 2.0)
                            if (distance .lt. (real_sphere_info(i,1) + ghost_sphere_info(k,1))) then
                                is_overlapped = .true.
                                exit
                            end if
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    ! update the ghost_sphere_info and total_ghost_sphere_num
                    do j = 1, 7
                        if (ghost_sphere_status(j) .eq. 1) then 
                            total_ghost_sphere_num = total_ghost_sphere_num + 1
                            ghost_sphere_info(total_ghost_sphere_num,1) = real_sphere_info(i,1)
                            ghost_sphere_info(total_ghost_sphere_num,2:4) = ghost_sphere_pos(j,1:3)
                        end if
                    end do
                    if (mod(i,nint(0.25 * total_real_sphere_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_real_sphere_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                    i = i + 1
                end if

            end do

            if (allocated(spherical_particle_info)) deallocate(spherical_particle_info)
            allocate(spherical_particle_info(total_real_sphere_num + total_ghost_sphere_num,4))
            do i = 1, total_real_sphere_num
                spherical_particle_info(i,1:4) = real_sphere_info(i,1:4)
            end do
            do i = 1, total_ghost_sphere_num
                spherical_particle_info(total_real_sphere_num+i,1:4) = ghost_sphere_info(i,1:4)
            end do
            write(*,'(5X,a15,i6,a18,i6,a10,i6)') 'Real particles:',total_real_sphere_num,'. Ghost particles:',total_ghost_sphere_num,'. Totally:',total_real_sphere_num+total_ghost_sphere_num
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate shape of ellipsoid particles
        ! ellipsoid_particle_info    size:   [total_particle_num,9]
        ! each row contains:     radius_a   radius_b    radius_c    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 and aspect_ratio2 are used here
        ! radius_a / radius_b = aspect_ratio1
        ! radius_b / radius_c = aspect_ratio2
        ! ONLY the shape of each ellipsoid is determined (x_pos, y_pos and z_pos are undefined)
        subroutine ellipsoid_shape(ellipsoid_particle_info,container_volume)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  ellipsoid_particle_info(:,:)
            real(kind=8),intent(in)     ::  container_volume
            integer(kind=4)             ::  total_particle_num, i, j
            character(len=100)          ::  ellipsoid_radius_file
            real(kind=8),allocatable    ::  particle_needed_vol(:), particle_real_vol(:)
            real(kind=8)                ::  last_particle_vol, r(1), particle_real_vol_in_current_sieve, low_value, high_value, radius_b
            
            if (allocated(particle_needed_vol)) deallocate(particle_needed_vol)
            if (allocated(particle_real_vol)) deallocate(particle_real_vol)
            allocate(particle_needed_vol(total_sieve_num-1))
            allocate(particle_real_vol(total_sieve_num-1))
            
            do i = 1, total_sieve_num-1
                particle_needed_vol(i) = (gradation(2,i+1) - gradation(2,i)) * particle_occupy_percentage * container_volume
            end do
            
            ! randomize
            call random_seed()
 
            ellipsoid_radius_file='ellipsoid_radius_file.tmp'   
            open(unit=100,file=ellipsoid_radius_file,form='formatted',status='replace',action='write')
            particle_real_vol = 0.0
            total_particle_num = 0
            do i = total_sieve_num-1, 1, -1
                particle_real_vol_in_current_sieve = 0.0 
                do while (particle_real_vol(i) .lt. particle_needed_vol(i))
                    ! generate the ellipsoid
                    call random_number(r)
                    low_value = gradation(1,i) / 2.0
                    high_value = gradation(1,i+1) / 2.0
                    radius_b = r(1) * (high_value - low_value) + low_value
                    particle_real_vol_in_current_sieve = particle_real_vol_in_current_sieve + (radius_b * aspect_ratio1) * radius_b * (radius_b / aspect_ratio2) * 3.14159 * 4.0 / 3.0
                    if (particle_real_vol_in_current_sieve .gt. particle_needed_vol(i)) then
                        last_particle_vol = particle_needed_vol(i) - particle_real_vol(i)
                        radius_b = (last_particle_vol * 3.0 / 4.0 / 3.14159 * aspect_ratio2 / aspect_ratio1) ** (1.0 / 3.0)
                    end if
                    write(100,'(e15.5)') radius_b
                    total_particle_num = total_particle_num + 1
                    particle_real_vol(i) = particle_real_vol(i) + (radius_b * aspect_ratio1) * radius_b * (radius_b / aspect_ratio2) * 3.14159 * 4.0 / 3.0
                end do
            end do
            close(unit=100)

            if (allocated(ellipsoid_particle_info)) deallocate(ellipsoid_particle_info)
            allocate(ellipsoid_particle_info(total_particle_num,9))
            ellipsoid_particle_info = 0.0

            open(unit=100,file=ellipsoid_radius_file,form='formatted',status='old',action='read')
            do i = 1, total_particle_num
                read(100,*) ellipsoid_particle_info(i,2)
                ellipsoid_particle_info(i,1) = ellipsoid_particle_info(i,2) * aspect_ratio1
                ellipsoid_particle_info(i,3) = ellipsoid_particle_info(i,2) / aspect_ratio2
            end do
            close(unit=100,status='delete')

            i = 1
            call sort(ellipsoid_particle_info,i,'d')

            ! generate the rotation angle
            do i = 1, total_particle_num

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                ellipsoid_particle_info(i,4) = r(1)

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                ellipsoid_particle_info(i,5) = r(1)

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                ellipsoid_particle_info(i,6) = r(1)
            end do

        end subroutine

        ! generate ellipsoid particles
        ! ellipsoid_particle_info    size:   [total_particle_num,9]
        ! each row contains:     radius_a   radius_b    radius_c    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 and aspect_ratio2 are used here
        ! radius_a / radius_b = aspect_ratio1
        ! radius_b / radius_c = aspect_ratio2
        subroutine ellipsoid_without_periodic_bc(ellipsoid_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  ellipsoid_particle_info(:,:)
            real(kind=8),allocatable    ::  mesh_a_info(:,:), mesh_b_info(:,:)
            real(kind=8),allocatable    ::  ellipsoid_particle_mesh_x_pos(:,:), ellipsoid_particle_mesh_y_pos(:,:), ellipsoid_particle_mesh_z_pos(:,:)
            real(kind=8)                ::  ellipsoid_a_info(9), ellipsoid_b_info(9), r(1), model_x_len, model_y_len, model_z_len
            integer(kind=4)             ::  total_particle_num, longitude_vertex_num, latitude_vertex_num, total_node_num, i, j, k, progress
            logical                     ::  is_out_region, is_overlapped

            model_x_len = container_size1
            model_y_len = container_size2
            model_z_len = container_size3

            call random_seed()

            call ellipsoid_shape(ellipsoid_particle_info,model_x_len*model_y_len*model_z_len)
            total_particle_num = size(ellipsoid_particle_info,1)

            ellipsoid_particle_info(1:total_particle_num,7:9) = 0.0        ! initiate the position of ellipsoid

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            longitude_vertex_num = 10
            latitude_vertex_num = 10
            total_node_num = longitude_vertex_num * latitude_vertex_num

            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(ellipsoid_particle_mesh_x_pos)) deallocate(ellipsoid_particle_mesh_x_pos)
            allocate(ellipsoid_particle_mesh_x_pos(total_particle_num,total_node_num))
            ellipsoid_particle_mesh_x_pos = 0.0

            if (allocated(ellipsoid_particle_mesh_y_pos)) deallocate(ellipsoid_particle_mesh_y_pos)
            allocate(ellipsoid_particle_mesh_y_pos(total_particle_num,total_node_num))
            ellipsoid_particle_mesh_y_pos = 0.0

            if (allocated(ellipsoid_particle_mesh_z_pos)) deallocate(ellipsoid_particle_mesh_z_pos)
            allocate(ellipsoid_particle_mesh_z_pos(total_particle_num,total_node_num))
            ellipsoid_particle_mesh_z_pos = 0.0
            do i = 1, total_particle_num
                do j = 1, 9
                    ellipsoid_a_info(j) = ellipsoid_particle_info(i,j)
                end do
                call generate_mesh_on_ellipsoid_surface(ellipsoid_a_info,longitude_vertex_num,latitude_vertex_num,mesh_a_info)
                do j = 1, total_node_num
                    ellipsoid_particle_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the ellipsoid is located at 0,0,0
                    ellipsoid_particle_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    ellipsoid_particle_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            i = 1
            do while (i .le. total_particle_num)
                is_out_region = .false.
                is_overlapped = .false.

                call random_number(r)
                ellipsoid_particle_info(i,7) = r(1) * model_x_len - 0.5 * model_x_len
                
                call random_number(r)
                ellipsoid_particle_info(i,8) = r(1) * model_y_len - 0.5 * model_y_len

                call random_number(r)
                ellipsoid_particle_info(i,9) = r(1) * model_z_len - 0.5 * model_z_len
                
                do j = 1, total_node_num
                    mesh_a_info(j,1) = ellipsoid_particle_mesh_x_pos(i,j) + ellipsoid_particle_info(i,7)
                    mesh_a_info(j,2) = ellipsoid_particle_mesh_y_pos(i,j) + ellipsoid_particle_info(i,8)
                    mesh_a_info(j,3) = ellipsoid_particle_mesh_z_pos(i,j) + ellipsoid_particle_info(i,9)
                end do
                is_out_region = is_geometry_out_of_box_region(mesh_a_info,-0.5*model_x_len,0.5*model_x_len,-0.5*model_y_len,0.5*model_y_len,-0.5*model_z_len,0.5*model_z_len)
                if (is_out_region .eq. .true.) cycle

                do j = 1, 9
                    ellipsoid_a_info(j) = ellipsoid_particle_info(i,j)
                end do

                do j = 1, i-1
                    do k = 1, 9
                        ellipsoid_b_info(k) = ellipsoid_particle_info(j,k)
                    end do

                    do k = 1, total_node_num
                        mesh_b_info(k,1) = ellipsoid_particle_mesh_x_pos(j,k) + ellipsoid_particle_info(j,7)
                        mesh_b_info(k,2) = ellipsoid_particle_mesh_y_pos(j,k) + ellipsoid_particle_info(j,8)
                        mesh_b_info(k,3) = ellipsoid_particle_mesh_z_pos(j,k) + ellipsoid_particle_info(j,9)
                    end do

                    is_overlapped = is_ellipsoid_overlapped(ellipsoid_a_info,ellipsoid_b_info,mesh_a_info,mesh_b_info)

                    if (is_overlapped .eq. .true.) exit

                end do

                if ((is_out_region .eq. .false.) .and. (is_overlapped .eq. .false.)) then
                    i = i + 1;
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate ellipsoid particles
        ! ellipsoid_particle_info    size:   [total_particle_num,9]
        ! each row contains:     radius_a   radius_b    radius_c    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 and aspect_ratio2 are used here
        ! radius_a / radius_b = aspect_ratio1
        ! radius_b / radius_c = aspect_ratio2
        ! for cylinder container
        subroutine ellipsoid_in_cylinder_container(ellipsoid_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  ellipsoid_particle_info(:,:)
            real(kind=8),allocatable    ::  mesh_a_info(:,:), mesh_b_info(:,:)
            real(kind=8),allocatable    ::  ellipsoid_particle_mesh_x_pos(:,:), ellipsoid_particle_mesh_y_pos(:,:), ellipsoid_particle_mesh_z_pos(:,:)
            real(kind=8)                ::  ellipsoid_a_info(9), ellipsoid_b_info(9), r(1), container_radius, container_height
            integer(kind=4)             ::  total_particle_num, longitude_vertex_num, latitude_vertex_num, total_node_num, i, j, k, progress
            logical                     ::  is_out_region, is_overlapped

            container_radius = container_size1
            container_height = container_size2

            call random_seed()

            call ellipsoid_shape(ellipsoid_particle_info,container_radius**2.0*3.14159*container_height)
            total_particle_num = size(ellipsoid_particle_info,1)

            ellipsoid_particle_info(1:total_particle_num,7:9) = 0.0        ! initiate the position of ellipsoid

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            longitude_vertex_num = 10
            latitude_vertex_num = 10
            total_node_num = longitude_vertex_num * latitude_vertex_num

            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(ellipsoid_particle_mesh_x_pos)) deallocate(ellipsoid_particle_mesh_x_pos)
            allocate(ellipsoid_particle_mesh_x_pos(total_particle_num,total_node_num))
            ellipsoid_particle_mesh_x_pos = 0.0

            if (allocated(ellipsoid_particle_mesh_y_pos)) deallocate(ellipsoid_particle_mesh_y_pos)
            allocate(ellipsoid_particle_mesh_y_pos(total_particle_num,total_node_num))
            ellipsoid_particle_mesh_y_pos = 0.0

            if (allocated(ellipsoid_particle_mesh_z_pos)) deallocate(ellipsoid_particle_mesh_z_pos)
            allocate(ellipsoid_particle_mesh_z_pos(total_particle_num,total_node_num))
            ellipsoid_particle_mesh_z_pos = 0.0
            do i = 1, total_particle_num
                do j = 1, 9
                    ellipsoid_a_info(j) = ellipsoid_particle_info(i,j)
                end do
                call generate_mesh_on_ellipsoid_surface(ellipsoid_a_info,longitude_vertex_num,latitude_vertex_num,mesh_a_info)
                do j = 1, total_node_num
                    ellipsoid_particle_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the ellipsoid is located at 0,0,0
                    ellipsoid_particle_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    ellipsoid_particle_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            i = 1
            do while (i .le. total_particle_num)
                is_out_region = .false.
                is_overlapped = .false.

                call random_number(r)
                ellipsoid_particle_info(i,7) = r(1) * 2 * container_radius - container_radius
                
                call random_number(r)
                ellipsoid_particle_info(i,8) = r(1) * 2 * container_radius - container_radius

                call random_number(r)
                ellipsoid_particle_info(i,9) = r(1) * container_height - 0.5 * container_height
                
                do j = 1, total_node_num
                    mesh_a_info(j,1) = ellipsoid_particle_mesh_x_pos(i,j) + ellipsoid_particle_info(i,7)
                    mesh_a_info(j,2) = ellipsoid_particle_mesh_y_pos(i,j) + ellipsoid_particle_info(i,8)
                    mesh_a_info(j,3) = ellipsoid_particle_mesh_z_pos(i,j) + ellipsoid_particle_info(i,9)
                    if ((mesh_a_info(j,3) .lt. -0.5 * container_height) .or. (mesh_a_info(j,3) .gt. 0.5 * container_height)) then
                        is_out_region = .true.
                        exit
                    end if
                    if (mesh_a_info(j,1) ** 2.0 + mesh_a_info(j,2) ** 2.0 .gt. container_radius ** 2.0) then
                        is_out_region = .true.
                        exit
                    end if
                end do
                if (is_out_region .eq. .true.) cycle

                do j = 1, 9
                    ellipsoid_a_info(j) = ellipsoid_particle_info(i,j)
                end do

                do j = 1, i-1
                    do k = 1, 9
                        ellipsoid_b_info(k) = ellipsoid_particle_info(j,k)
                    end do

                    do k = 1, total_node_num
                        mesh_b_info(k,1) = ellipsoid_particle_mesh_x_pos(j,k) + ellipsoid_particle_info(j,7)
                        mesh_b_info(k,2) = ellipsoid_particle_mesh_y_pos(j,k) + ellipsoid_particle_info(j,8)
                        mesh_b_info(k,3) = ellipsoid_particle_mesh_z_pos(j,k) + ellipsoid_particle_info(j,9)
                    end do

                    is_overlapped = is_ellipsoid_overlapped(ellipsoid_a_info,ellipsoid_b_info,mesh_a_info,mesh_b_info)

                    if (is_overlapped .eq. .true.) exit

                end do

                if ((is_out_region .eq. .false.) .and. (is_overlapped .eq. .false.)) then
                    i = i + 1;
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate ellipsoid particles
        ! ellipsoid_particle_info    size:   [total_particle_num,9]
        ! each row contains:     radius_a   radius_b    radius_c    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 and aspect_ratio2 are used here
        ! radius_a / radius_b = aspect_ratio1
        ! radius_b / radius_c = aspect_ratio2
        subroutine ellipsoid_with_periodic_bc(ellipsoid_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  ellipsoid_particle_info(:,:)
            real(kind=8),allocatable    ::  real_ellipsoid_info(:,:), ghost_ellipsoid_info(:,:)
            real(kind=8),allocatable    ::  real_ellipsoid_mesh_x_pos(:,:), real_ellipsoid_mesh_y_pos(:,:), real_ellipsoid_mesh_z_pos(:,:)
            real(kind=8),allocatable    ::  ghost_ellipsoid_mesh_x_pos(:,:), ghost_ellipsoid_mesh_y_pos(:,:), ghost_ellipsoid_mesh_z_pos(:,:)
            real(kind=8),allocatable    ::  ref_sphere_radius(:), mesh_a_info(:,:), mesh_b_info(:,:)
            real(kind=8)                ::  r(1), distance, ghost_ellipsoid_pos(7,3), x_pos, y_pos, z_pos, ellipsoid_a_info(9), ellipsoid_b_info(9), model_x_len, model_y_len, model_z_len
            integer(kind=4)             ::  i, j, k, ii, flag, progress, ghost_ellipsoid_status(7), total_node_num
            integer(kind=4)             ::  longitude_vertex_num, latitude_vertex_num, total_real_ellipsoid_num, total_ghost_ellipsoid_num
            logical                     ::  is_overlapped, is_out_region, is_part_in_x_range, is_part_in_y_range, is_part_in_z_range

            model_x_len = container_size1
            model_y_len = container_size2
            model_z_len = container_size3

            call random_seed()

            call ellipsoid_shape(real_ellipsoid_info,model_x_len*model_y_len*model_z_len)
            total_real_ellipsoid_num = size(real_ellipsoid_info,1)

            real_ellipsoid_info(1:total_real_ellipsoid_num,7:9) = 0.0        ! initiate the position of ellipsoid

            if (allocated(ghost_ellipsoid_info)) deallocate(ghost_ellipsoid_info)
            allocate(ghost_ellipsoid_info(total_real_ellipsoid_num*7,9))
            ghost_ellipsoid_info = 0.0
            
            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a17)') 'Particle shape done.','(',total_real_ellipsoid_num,'real particles).'

            ! generate mesh of original real ellipsoid
            longitude_vertex_num = 10
            latitude_vertex_num = 10
            total_node_num = longitude_vertex_num * latitude_vertex_num

            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(real_ellipsoid_mesh_x_pos)) deallocate(real_ellipsoid_mesh_x_pos)
            allocate(real_ellipsoid_mesh_x_pos(total_real_ellipsoid_num,total_node_num))
            real_ellipsoid_mesh_x_pos = 0.0

            if (allocated(real_ellipsoid_mesh_y_pos)) deallocate(real_ellipsoid_mesh_y_pos)
            allocate(real_ellipsoid_mesh_y_pos(total_real_ellipsoid_num,total_node_num))
            real_ellipsoid_mesh_y_pos = 0.0

            if (allocated(real_ellipsoid_mesh_z_pos)) deallocate(real_ellipsoid_mesh_z_pos)
            allocate(real_ellipsoid_mesh_z_pos(total_real_ellipsoid_num,total_node_num))
            real_ellipsoid_mesh_z_pos = 0.0
            
            if (allocated(ghost_ellipsoid_mesh_x_pos)) deallocate(ghost_ellipsoid_mesh_x_pos)
            allocate(ghost_ellipsoid_mesh_x_pos(total_real_ellipsoid_num*7,total_node_num))
            ghost_ellipsoid_mesh_x_pos = 0.0

            if (allocated(ghost_ellipsoid_mesh_y_pos)) deallocate(ghost_ellipsoid_mesh_y_pos)
            allocate(ghost_ellipsoid_mesh_y_pos(total_real_ellipsoid_num*7,total_node_num))
            ghost_ellipsoid_mesh_y_pos = 0.0

            if (allocated(ghost_ellipsoid_mesh_z_pos)) deallocate(ghost_ellipsoid_mesh_z_pos)
            allocate(ghost_ellipsoid_mesh_z_pos(total_real_ellipsoid_num*7,total_node_num))
            ghost_ellipsoid_mesh_z_pos = 0.0

            do i = 1, total_real_ellipsoid_num
                do j = 1, 9
                    ellipsoid_a_info(j) = real_ellipsoid_info(i,j)
                end do
                call generate_mesh_on_ellipsoid_surface(ellipsoid_a_info,longitude_vertex_num,latitude_vertex_num,mesh_a_info)
                do j = 1, total_node_num
                    real_ellipsoid_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the ellipsoid is located at 0,0,0
                    real_ellipsoid_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    real_ellipsoid_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            ! generate the position of ellipsoid particles
            total_ghost_ellipsoid_num = 0
            i = 1
            do while (i .le. total_real_ellipsoid_num)

                call random_number(r)
                real_ellipsoid_info(i,7) = r(1) * model_x_len - 0.5 * model_x_len

                call random_number(r)
                real_ellipsoid_info(i,8) = r(1) * model_y_len - 0.5 * model_y_len

                call random_number(r)
                real_ellipsoid_info(i,9) = r(1) * model_z_len - 0.5 * model_z_len

                ! check whether the current real ellipsoid is intersected by the boundary of model
                ! if yes, generate ghost ellipsoids
                is_out_region = .false.
                do j = 1, total_node_num
                    x_pos = real_ellipsoid_mesh_x_pos(i,j) + real_ellipsoid_info(i,7)
                    y_pos = real_ellipsoid_mesh_y_pos(i,j) + real_ellipsoid_info(i,8)
                    z_pos = real_ellipsoid_mesh_z_pos(i,j) + real_ellipsoid_info(i,9)
                    if ((x_pos .lt. -0.5 * model_x_len) .or. (x_pos .gt. 0.5 * model_x_len)) then
                        is_out_region = .true.
                        exit
                    end if 
                    if ((y_pos .lt. -0.5 * model_y_len) .or. (y_pos .gt. 0.5 * model_y_len)) then
                        is_out_region = .true.
                        exit
                    end if
                    if ((z_pos .lt. -0.5 * model_z_len) .or. (z_pos .gt. 0.5 * model_z_len)) then
                        is_out_region = .true.
                        exit
                    end if 
                end do

                ghost_ellipsoid_status = 0
                if (is_out_region .eq. .true.) then
                    ! generate ghost particle based on the current real particle
                    ! and detect whether the ghost particle is intersected with 
                    ! the region
                    ! for ellipsoid, seven possible ghost particles exist
                    if (real_ellipsoid_info(i,7) .gt. 0.0) then
                        x_pos = real_ellipsoid_info(i,7) - model_x_len
                    else
                        x_pos = real_ellipsoid_info(i,7) + model_x_len
                    end if
                    if (real_ellipsoid_info(i,8) .gt. 0.0) then
                        y_pos = real_ellipsoid_info(i,8) - model_y_len
                    else
                        y_pos = real_ellipsoid_info(i,8) + model_y_len
                    end if
                    if (real_ellipsoid_info(i,9) .gt. 0.0) then
                        z_pos = real_ellipsoid_info(i,9) - model_z_len
                    else
                        z_pos = real_ellipsoid_info(i,9) + model_z_len
                    end if

                    ghost_ellipsoid_pos(1,1) = x_pos
                    ghost_ellipsoid_pos(1,2) = real_ellipsoid_info(i,8)
                    ghost_ellipsoid_pos(1,3) = real_ellipsoid_info(i,9)

                    ghost_ellipsoid_pos(2,1) = x_pos
                    ghost_ellipsoid_pos(2,2) = y_pos
                    ghost_ellipsoid_pos(2,3) = real_ellipsoid_info(i,9)

                    ghost_ellipsoid_pos(3,1) = x_pos
                    ghost_ellipsoid_pos(3,2) = y_pos
                    ghost_ellipsoid_pos(3,3) = z_pos

                    ghost_ellipsoid_pos(4,1) = x_pos
                    ghost_ellipsoid_pos(4,2) = real_ellipsoid_info(i,8)
                    ghost_ellipsoid_pos(4,3) = z_pos

                    ghost_ellipsoid_pos(5,1) = real_ellipsoid_info(i,7)
                    ghost_ellipsoid_pos(5,2) = y_pos
                    ghost_ellipsoid_pos(5,3) = real_ellipsoid_info(i,9)

                    ghost_ellipsoid_pos(6,1) = real_ellipsoid_info(i,7)
                    ghost_ellipsoid_pos(6,2) = y_pos
                    ghost_ellipsoid_pos(6,3) = z_pos

                    ghost_ellipsoid_pos(7,1) = real_ellipsoid_info(i,7)
                    ghost_ellipsoid_pos(7,2) = real_ellipsoid_info(i,8)
                    ghost_ellipsoid_pos(7,3) = z_pos

                    do j = 1, 7
                        do k = 1, total_node_num
                            is_part_in_x_range = .false.
                            is_part_in_y_range = .false.
                            is_part_in_z_range = .false.
                            x_pos = real_ellipsoid_mesh_x_pos(i,k) + ghost_ellipsoid_pos(j,1)
                            y_pos = real_ellipsoid_mesh_y_pos(i,k) + ghost_ellipsoid_pos(j,2)
                            z_pos = real_ellipsoid_mesh_z_pos(i,k) + ghost_ellipsoid_pos(j,3)
                            if ((x_pos .gt. -0.5 * model_x_len) .and. (x_pos .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                            if ((y_pos .gt. -0.5 * model_y_len) .and. (y_pos .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                            if ((z_pos .gt. -0.5 * model_z_len) .and. (z_pos .lt. 0.5 * model_z_len)) is_part_in_z_range = .true.
                            if ((is_part_in_x_range .eq. .true.) .and. (is_part_in_y_range .eq. .true.) .and. (is_part_in_z_range .eq. .true.)) then
                                ghost_ellipsoid_status(j) = 1
                                exit
                            end if
                        end do     
                    end do
                end if

                is_overlapped = .false.
                ! overlap, new real particle vs old real particles
                do j = 1, 9
                    ellipsoid_a_info(j) = real_ellipsoid_info(i,j)
                end do
                do j = 1, total_node_num
                    mesh_a_info(j,1) = real_ellipsoid_mesh_x_pos(i,j) + real_ellipsoid_info(i,7)
                    mesh_a_info(j,2) = real_ellipsoid_mesh_y_pos(i,j) + real_ellipsoid_info(i,8)
                    mesh_a_info(j,3) = real_ellipsoid_mesh_z_pos(i,j) + real_ellipsoid_info(i,9)
                end do   
                do j = 1, i-1
                    do k = 1, 9
                        ellipsoid_b_info(k) = real_ellipsoid_info(j,k)
                    end do
                    do k = 1, total_node_num
                        mesh_b_info(k,1) = real_ellipsoid_mesh_x_pos(j,k) + real_ellipsoid_info(j,7)
                        mesh_b_info(k,2) = real_ellipsoid_mesh_y_pos(j,k) + real_ellipsoid_info(j,8)
                        mesh_b_info(k,3) = real_ellipsoid_mesh_z_pos(j,k) + real_ellipsoid_info(j,9)
                    end do
                    is_overlapped = is_ellipsoid_overlapped(ellipsoid_a_info,ellipsoid_b_info,mesh_a_info,mesh_b_info)
                    if (is_overlapped .eq. .true.) exit
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new real particle vs old ghost particles 
                do j = 1, total_ghost_ellipsoid_num
                    do k = 1, 9
                        ellipsoid_b_info(k) = ghost_ellipsoid_info(j,k)
                    end do
                    do k = 1, total_node_num
                        mesh_b_info(k,1) = ghost_ellipsoid_mesh_x_pos(j,k) + ghost_ellipsoid_info(j,7)
                        mesh_b_info(k,2) = ghost_ellipsoid_mesh_y_pos(j,k) + ghost_ellipsoid_info(j,8)
                        mesh_b_info(k,3) = ghost_ellipsoid_mesh_z_pos(j,k) + ghost_ellipsoid_info(j,9)
                    end do
                    is_overlapped = is_ellipsoid_overlapped(ellipsoid_a_info,ellipsoid_b_info,mesh_a_info,mesh_b_info)
                    if (is_overlapped .eq. .true.) exit
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old real particles
                do j = 1, 7
                    if (ghost_ellipsoid_status(j) .eq. 1) then
                        do k = 1, 6
                            ellipsoid_a_info(k) = real_ellipsoid_info(i,k)
                        end do
                        do k = 7, 9
                            ellipsoid_a_info(k) = ghost_ellipsoid_pos(j,k-6)
                        end do
                        do k = 1, total_node_num
                            mesh_a_info(k,1) = real_ellipsoid_mesh_x_pos(i,k) + ghost_ellipsoid_pos(j,1)
                            mesh_a_info(k,2) = real_ellipsoid_mesh_y_pos(i,k) + ghost_ellipsoid_pos(j,2)
                            mesh_a_info(k,3) = real_ellipsoid_mesh_z_pos(i,k) + ghost_ellipsoid_pos(j,3)
                        end do   
                        do k = 1, i-1
                            do ii = 1, 9
                                ellipsoid_b_info(ii) = real_ellipsoid_info(k,ii)
                            end do
                            do ii = 1, total_node_num
                                mesh_b_info(ii,1) = real_ellipsoid_mesh_x_pos(k,ii) + real_ellipsoid_info(k,7)
                                mesh_b_info(ii,2) = real_ellipsoid_mesh_y_pos(k,ii) + real_ellipsoid_info(k,8)
                                mesh_b_info(ii,3) = real_ellipsoid_mesh_z_pos(k,ii) + real_ellipsoid_info(k,9)
                            end do
                            is_overlapped = is_ellipsoid_overlapped(ellipsoid_a_info,ellipsoid_b_info,mesh_a_info,mesh_b_info)
                            if (is_overlapped .eq. .true.) exit
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old ghost particles
                do j = 1, 7
                    if (ghost_ellipsoid_status(j) .eq. 1) then
                        do k = 1, 6
                            ellipsoid_a_info(k) = real_ellipsoid_info(i,k)
                        end do
                        do k = 7, 9
                            ellipsoid_a_info(k) = ghost_ellipsoid_pos(j,k-6)
                        end do
                        do k = 1, total_node_num
                            mesh_a_info(k,1) = real_ellipsoid_mesh_x_pos(i,k) + ghost_ellipsoid_pos(j,1)
                            mesh_a_info(k,2) = real_ellipsoid_mesh_y_pos(i,k) + ghost_ellipsoid_pos(j,2)
                            mesh_a_info(k,3) = real_ellipsoid_mesh_z_pos(i,k) + ghost_ellipsoid_pos(j,3)
                        end do   
                        do k = 1, total_ghost_ellipsoid_num
                            do ii = 1, 9
                                ellipsoid_b_info(ii) = ghost_ellipsoid_info(k,ii)
                            end do
                            do ii = 1, total_node_num
                                mesh_b_info(ii,1) = ghost_ellipsoid_mesh_x_pos(k,ii) + ghost_ellipsoid_info(k,7)
                                mesh_b_info(ii,2) = ghost_ellipsoid_mesh_y_pos(k,ii) + ghost_ellipsoid_info(k,8)
                                mesh_b_info(ii,3) = ghost_ellipsoid_mesh_z_pos(k,ii) + ghost_ellipsoid_info(k,9)
                            end do
                            is_overlapped = is_ellipsoid_overlapped(ellipsoid_a_info,ellipsoid_b_info,mesh_a_info,mesh_b_info)
                            if (is_overlapped .eq. .true.) exit
                        end do                
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    ! update:
                    ! total_ghost_ellipsoid_num
                    ! ghost_ellipsoid_info
                    ! ghost_ellipsoid_mesh_x_pos
                    ! ghost_ellipsoid_mesh_y_pos
                    ! ghost_ellipsoid_mesh_z_pos
                    do j = 1, 7
                        if (ghost_ellipsoid_status(j) .eq. 1) then 
                            total_ghost_ellipsoid_num = total_ghost_ellipsoid_num + 1
                            ghost_ellipsoid_info(total_ghost_ellipsoid_num,1:6) = real_ellipsoid_info(i,1:6)
                            ghost_ellipsoid_info(total_ghost_ellipsoid_num,7:9) = ghost_ellipsoid_pos(j,1:3)
                            ghost_ellipsoid_mesh_x_pos(total_ghost_ellipsoid_num,:) = real_ellipsoid_mesh_x_pos(i,:)
                            ghost_ellipsoid_mesh_y_pos(total_ghost_ellipsoid_num,:) = real_ellipsoid_mesh_y_pos(i,:)
                            ghost_ellipsoid_mesh_z_pos(total_ghost_ellipsoid_num,:) = real_ellipsoid_mesh_z_pos(i,:)
                        end if
                    end do
                    if (mod(i,nint(0.25 * total_real_ellipsoid_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_real_ellipsoid_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                    i = i + 1
                end if

            end do

            if (allocated(ellipsoid_particle_info)) deallocate(ellipsoid_particle_info)
            allocate(ellipsoid_particle_info(total_real_ellipsoid_num + total_ghost_ellipsoid_num,9))
            do i = 1, total_real_ellipsoid_num
                ellipsoid_particle_info(i,1:9) = real_ellipsoid_info(i,1:9)
            end do
            do i = 1, total_ghost_ellipsoid_num
                ellipsoid_particle_info(total_real_ellipsoid_num+i,1:9) = ghost_ellipsoid_info(i,1:9)
            end do
            write(*,'(5X,a15,i6,a18,i6,a10,i6)') 'Real particles:',total_real_ellipsoid_num,'. Ghost particles:',total_ghost_ellipsoid_num,'. Totally:',total_real_ellipsoid_num+total_ghost_ellipsoid_num
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate cylinder particles
        ! cylinder_particle_info    size:   [total_particle_num,8]
        ! each row contains:     radius_a   radius_b   height    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 and aspect_ratio2 are used here
        ! height / (radius_a * 2) = aspect_ratio1
        ! radius_a / radius_b = aspect_ratio2
        subroutine cylinder_without_periodic_bc(cylinder_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  cylinder_particle_info(:,:)
            real(kind=8),allocatable    ::  mesh_a_info(:,:), mesh_b_info(:,:)
            real(kind=8),allocatable    ::  cylinder_particle_mesh_x_pos(:,:), cylinder_particle_mesh_y_pos(:,:), cylinder_particle_mesh_z_pos(:,:)
            real(kind=8)                ::  cylinder_a_info(9), cylinder_b_info(9), r(1), model_x_len, model_y_len, model_z_len
            integer(kind=4)             ::  total_particle_num, total_radius_segment, total_height_segment, total_node_num, i, j, k, progress
            logical                     ::  is_out_region, is_overlapped

            model_x_len = container_size1
            model_y_len = container_size2
            model_z_len = container_size3

            call random_seed()

            call cylinder_shape(cylinder_particle_info,model_x_len*model_y_len*model_z_len)
            total_particle_num = size(cylinder_particle_info,1)

            cylinder_particle_info(1:total_particle_num,7:9) = 0.0        ! initiate the position of cylinder

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            total_radius_segment = 10
            total_height_segment = 10
            total_node_num = total_radius_segment * total_height_segment

            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(cylinder_particle_mesh_x_pos)) deallocate(cylinder_particle_mesh_x_pos)
            allocate(cylinder_particle_mesh_x_pos(total_particle_num,total_node_num))
            cylinder_particle_mesh_x_pos = 0.0

            if (allocated(cylinder_particle_mesh_y_pos)) deallocate(cylinder_particle_mesh_y_pos)
            allocate(cylinder_particle_mesh_y_pos(total_particle_num,total_node_num))
            cylinder_particle_mesh_y_pos = 0.0

            if (allocated(cylinder_particle_mesh_z_pos)) deallocate(cylinder_particle_mesh_z_pos)
            allocate(cylinder_particle_mesh_z_pos(total_particle_num,total_node_num))
            cylinder_particle_mesh_z_pos = 0.0
            do i = 1, total_particle_num
                do j = 1, 9
                    cylinder_a_info(j) = cylinder_particle_info(i,j)
                end do
                call generate_mesh_on_cylinder_surface(cylinder_a_info,total_radius_segment,total_height_segment,mesh_a_info)
                do j = 1, total_node_num
                    cylinder_particle_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the cylinder is located at 0,0,0
                    cylinder_particle_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    cylinder_particle_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            i = 1
            do while (i .le. total_particle_num)

                is_out_region = .false.
                is_overlapped = .false.

                call random_number(r)
                cylinder_particle_info(i,7) = r(1) * model_x_len - 0.5 * model_x_len
                
                call random_number(r)
                cylinder_particle_info(i,8) = r(1) * model_y_len - 0.5 * model_y_len

                call random_number(r)
                cylinder_particle_info(i,9) = r(1) * model_z_len - 0.5 * model_z_len
                
                do j = 1, total_node_num
                    mesh_a_info(j,1) = cylinder_particle_mesh_x_pos(i,j) + cylinder_particle_info(i,7)
                    mesh_a_info(j,2) = cylinder_particle_mesh_y_pos(i,j) + cylinder_particle_info(i,8)
                    mesh_a_info(j,3) = cylinder_particle_mesh_z_pos(i,j) + cylinder_particle_info(i,9)
                end do
                is_out_region = is_geometry_out_of_box_region(mesh_a_info,-0.5*model_x_len,0.5*model_x_len,-0.5*model_y_len,0.5*model_y_len,-0.5*model_z_len,0.5*model_z_len)  
                if (is_out_region .eq. .true.) cycle

                do j = 1, 9
                    cylinder_a_info(j) = cylinder_particle_info(i,j)
                end do

                do j = 1, i-1
                    do k = 1, 9
                        cylinder_b_info(k) = cylinder_particle_info(j,k)
                    end do

                    do k = 1, total_node_num
                        mesh_b_info(k,1) = cylinder_particle_mesh_x_pos(j,k) + cylinder_particle_info(j,7)
                        mesh_b_info(k,2) = cylinder_particle_mesh_y_pos(j,k) + cylinder_particle_info(j,8)
                        mesh_b_info(k,3) = cylinder_particle_mesh_z_pos(j,k) + cylinder_particle_info(j,9)
                    end do

                    is_overlapped = is_cylinder_overlapped(cylinder_a_info,cylinder_b_info,mesh_a_info,mesh_b_info)

                    if (is_overlapped .eq. .true.) exit

                end do

                if ((is_out_region .eq. .false.) .and. (is_overlapped .eq. .false.)) then
                    i = i + 1;
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate cylinder particles
        ! cylinder_particle_info    size:   [total_particle_num,8]
        ! each row contains:     radius_a   radius_b   height    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 and aspect_ratio2 are used here
        ! height / (radius_a * 2) = aspect_ratio1
        ! radius_a / radius_b = aspect_ratio2
        subroutine cylinder_in_cylinder_container(cylinder_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  cylinder_particle_info(:,:)
            real(kind=8),allocatable    ::  mesh_a_info(:,:), mesh_b_info(:,:)
            real(kind=8),allocatable    ::  cylinder_particle_mesh_x_pos(:,:), cylinder_particle_mesh_y_pos(:,:), cylinder_particle_mesh_z_pos(:,:)
            real(kind=8)                ::  cylinder_a_info(9), cylinder_b_info(9), r(1), container_radius, container_height
            integer(kind=4)             ::  total_particle_num, total_radius_segment, total_height_segment, total_node_num, i, j, k, progress
            logical                     ::  is_out_region, is_overlapped

            container_radius = container_size1
            container_height = container_size2

            call random_seed()

            call cylinder_shape(cylinder_particle_info,container_radius**2.0*3.14159*container_height)
            total_particle_num = size(cylinder_particle_info,1)

            cylinder_particle_info(1:total_particle_num,7:9) = 0.0        ! initiate the position of cylinder

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            total_radius_segment = 10
            total_height_segment = 10
            total_node_num = total_radius_segment * total_height_segment

            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(cylinder_particle_mesh_x_pos)) deallocate(cylinder_particle_mesh_x_pos)
            allocate(cylinder_particle_mesh_x_pos(total_particle_num,total_node_num))
            cylinder_particle_mesh_x_pos = 0.0

            if (allocated(cylinder_particle_mesh_y_pos)) deallocate(cylinder_particle_mesh_y_pos)
            allocate(cylinder_particle_mesh_y_pos(total_particle_num,total_node_num))
            cylinder_particle_mesh_y_pos = 0.0

            if (allocated(cylinder_particle_mesh_z_pos)) deallocate(cylinder_particle_mesh_z_pos)
            allocate(cylinder_particle_mesh_z_pos(total_particle_num,total_node_num))
            cylinder_particle_mesh_z_pos = 0.0
            do i = 1, total_particle_num
                do j = 1, 9
                    cylinder_a_info(j) = cylinder_particle_info(i,j)
                end do
                call generate_mesh_on_cylinder_surface(cylinder_a_info,total_radius_segment,total_height_segment,mesh_a_info)
                do j = 1, total_node_num
                    cylinder_particle_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the cylinder is located at 0,0,0
                    cylinder_particle_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    cylinder_particle_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            i = 1
            do while (i .le. total_particle_num)

                is_out_region = .false.
                is_overlapped = .false.

                call random_number(r)
                cylinder_particle_info(i,7) = r(1) * 2 * container_radius - container_radius
                
                call random_number(r)
                cylinder_particle_info(i,8) = r(1) * 2 * container_radius - container_radius

                call random_number(r)
                cylinder_particle_info(i,9) = r(1) * container_height - 0.5 * container_height
                
                do j = 1, total_node_num
                    mesh_a_info(j,1) = cylinder_particle_mesh_x_pos(i,j) + cylinder_particle_info(i,7)
                    mesh_a_info(j,2) = cylinder_particle_mesh_y_pos(i,j) + cylinder_particle_info(i,8)
                    mesh_a_info(j,3) = cylinder_particle_mesh_z_pos(i,j) + cylinder_particle_info(i,9)
                    if ((mesh_a_info(j,3) .lt. -0.5 * container_height) .or. (mesh_a_info(j,3) .gt. 0.5 * container_height)) then
                        is_out_region = .true.
                        exit
                    end if
                    if (mesh_a_info(j,1) ** 2.0 + mesh_a_info(j,2) ** 2.0 .gt. container_radius ** 2.0) then
                        is_out_region = .true.
                        exit
                    end if
                end do
                if (is_out_region .eq. .true.) cycle

                do j = 1, 9
                    cylinder_a_info(j) = cylinder_particle_info(i,j)
                end do

                do j = 1, i-1
                    do k = 1, 9
                        cylinder_b_info(k) = cylinder_particle_info(j,k)
                    end do

                    do k = 1, total_node_num
                        mesh_b_info(k,1) = cylinder_particle_mesh_x_pos(j,k) + cylinder_particle_info(j,7)
                        mesh_b_info(k,2) = cylinder_particle_mesh_y_pos(j,k) + cylinder_particle_info(j,8)
                        mesh_b_info(k,3) = cylinder_particle_mesh_z_pos(j,k) + cylinder_particle_info(j,9)
                    end do

                    is_overlapped = is_cylinder_overlapped(cylinder_a_info,cylinder_b_info,mesh_a_info,mesh_b_info)

                    if (is_overlapped .eq. .true.) exit

                end do

                if ((is_out_region .eq. .false.) .and. (is_overlapped .eq. .false.)) then
                    i = i + 1;
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate cylinder particles
        ! cylinder_particle_info    size:   [total_particle_num,8]
        ! each row contains:     radius_a   radius_b   height    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 and aspect_ratio2 are used here
        ! height / (radius_a * 2) = aspect_ratio1
        ! radius_a / radius_b = aspect_ratio2
        subroutine cylinder_with_periodic_bc(cylinder_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  cylinder_particle_info(:,:)
            real(kind=8),allocatable    ::  mesh_a_info(:,:), mesh_b_info(:,:), real_cylinder_info(:,:), ghost_cylinder_info(:,:)
            real(kind=8),allocatable    ::  real_cylinder_mesh_x_pos(:,:), real_cylinder_mesh_y_pos(:,:), real_cylinder_mesh_z_pos(:,:)
            real(kind=8),allocatable    ::  ghost_cylinder_mesh_x_pos(:,:), ghost_cylinder_mesh_y_pos(:,:), ghost_cylinder_mesh_z_pos(:,:)
            real(kind=8)                ::  cylinder_a_info(9), cylinder_b_info(9), r(1), ghost_cylinder_pos(7,3), x_pos, y_pos, z_pos, model_x_len, model_y_len, model_z_len
            integer(kind=4)             ::  ghost_cylinder_status(7)
            integer(kind=4)             ::  total_real_cylinder_num, total_ghost_cylinder_num, total_radius_segment, total_height_segment
            integer(kind=4)             ::  total_node_num, i, j, k, ii, progress
            logical                     ::  is_out_region, is_overlapped, is_part_in_x_range, is_part_in_y_range, is_part_in_z_range

            model_x_len = container_size1
            model_y_len = container_size2
            model_z_len = container_size3

            call random_seed()

            call cylinder_shape(real_cylinder_info,model_x_len*model_y_len*model_z_len)
            total_real_cylinder_num = size(real_cylinder_info,1)

            real_cylinder_info(1:total_real_cylinder_num,7:9) = 0.0        ! initiate the position of cylinder

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a17)') 'Particle shape done.','(',total_real_cylinder_num,'real particles).'

            total_radius_segment = 10
            total_height_segment = 10
            total_node_num = total_radius_segment * total_height_segment
        
            if (allocated(ghost_cylinder_info)) deallocate(ghost_cylinder_info)
            allocate(ghost_cylinder_info(total_real_cylinder_num*7,9))
            ghost_cylinder_info = 0.0
            
            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(real_cylinder_mesh_x_pos)) deallocate(real_cylinder_mesh_x_pos)
            allocate(real_cylinder_mesh_x_pos(total_real_cylinder_num,total_node_num))
            real_cylinder_mesh_x_pos = 0.0

            if (allocated(real_cylinder_mesh_y_pos)) deallocate(real_cylinder_mesh_y_pos)
            allocate(real_cylinder_mesh_y_pos(total_real_cylinder_num,total_node_num))
            real_cylinder_mesh_y_pos = 0.0

            if (allocated(real_cylinder_mesh_z_pos)) deallocate(real_cylinder_mesh_z_pos)
            allocate(real_cylinder_mesh_z_pos(total_real_cylinder_num,total_node_num))
            real_cylinder_mesh_z_pos = 0.0

            if (allocated(ghost_cylinder_mesh_x_pos)) deallocate(ghost_cylinder_mesh_x_pos)
            allocate(ghost_cylinder_mesh_x_pos(total_real_cylinder_num*7,total_node_num))
            ghost_cylinder_mesh_x_pos = 0.0

            if (allocated(ghost_cylinder_mesh_y_pos)) deallocate(ghost_cylinder_mesh_y_pos)
            allocate(ghost_cylinder_mesh_y_pos(total_real_cylinder_num*7,total_node_num))
            ghost_cylinder_mesh_y_pos = 0.0

            if (allocated(ghost_cylinder_mesh_z_pos)) deallocate(ghost_cylinder_mesh_z_pos)
            allocate(ghost_cylinder_mesh_z_pos(total_real_cylinder_num*7,total_node_num))
            ghost_cylinder_mesh_z_pos = 0.0

            do i = 1, total_real_cylinder_num
                do j = 1, 9
                    cylinder_a_info(j) = real_cylinder_info(i,j)
                end do
                call generate_mesh_on_cylinder_surface(cylinder_a_info,total_radius_segment,total_height_segment,mesh_a_info)
                do j = 1, total_node_num
                    real_cylinder_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the cylinder is located at 0,0,0
                    real_cylinder_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    real_cylinder_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            total_ghost_cylinder_num = 0
            i = 1
            do while (i .le. total_real_cylinder_num)

                call random_number(r)
                real_cylinder_info(i,7) = r(1) * model_x_len - 0.5 * model_x_len
                
                call random_number(r)
                real_cylinder_info(i,8) = r(1) * model_y_len - 0.5 * model_y_len

                call random_number(r)
                real_cylinder_info(i,9) = r(1) * model_z_len - 0.5 * model_z_len
                
                ! check whether the current real cylinder is intersected by the boundary of model
                ! if yes, generate ghost cylinders
                is_out_region = .false.
                do j = 1, total_node_num
                    x_pos = real_cylinder_mesh_x_pos(i,j) + real_cylinder_info(i,7)
                    y_pos = real_cylinder_mesh_y_pos(i,j) + real_cylinder_info(i,8)
                    z_pos = real_cylinder_mesh_z_pos(i,j) + real_cylinder_info(i,9)
                    if ((x_pos .lt. -0.5 * model_x_len) .or. (x_pos .gt. 0.5 * model_x_len)) then
                        is_out_region = .true.
                        exit
                    end if 
                    if ((y_pos .lt. -0.5 * model_y_len) .or. (y_pos .gt. 0.5 * model_y_len)) then
                        is_out_region = .true.
                        exit
                    end if
                    if ((z_pos .lt. -0.5 * model_z_len) .or. (z_pos .gt. 0.5 * model_z_len)) then
                        is_out_region = .true.
                        exit
                    end if 
                end do
                
                ghost_cylinder_status = 0
                if (is_out_region .eq. .true.) then
                    ! generate ghost particle based on the current real particle
                    ! and detect whether the ghost particle is intersected with 
                    ! the region
                    ! for cylinder, seven possible ghost particles exist
                    if (real_cylinder_info(i,7) .gt. 0.0) then
                        x_pos = real_cylinder_info(i,7) - model_x_len
                    else
                        x_pos = real_cylinder_info(i,7) + model_x_len
                    end if
                    if (real_cylinder_info(i,8) .gt. 0.0) then
                        y_pos = real_cylinder_info(i,8) - model_y_len
                    else
                        y_pos = real_cylinder_info(i,8) + model_y_len
                    end if
                    if (real_cylinder_info(i,9) .gt. 0.0) then
                        z_pos = real_cylinder_info(i,9) - model_z_len
                    else
                        z_pos = real_cylinder_info(i,9) + model_z_len
                    end if

                    ghost_cylinder_pos(1,1) = x_pos
                    ghost_cylinder_pos(1,2) = real_cylinder_info(i,8)
                    ghost_cylinder_pos(1,3) = real_cylinder_info(i,9)

                    ghost_cylinder_pos(2,1) = x_pos
                    ghost_cylinder_pos(2,2) = y_pos
                    ghost_cylinder_pos(2,3) = real_cylinder_info(i,9)

                    ghost_cylinder_pos(3,1) = x_pos
                    ghost_cylinder_pos(3,2) = y_pos
                    ghost_cylinder_pos(3,3) = z_pos

                    ghost_cylinder_pos(4,1) = x_pos
                    ghost_cylinder_pos(4,2) = real_cylinder_info(i,8)
                    ghost_cylinder_pos(4,3) = z_pos

                    ghost_cylinder_pos(5,1) = real_cylinder_info(i,7)
                    ghost_cylinder_pos(5,2) = y_pos
                    ghost_cylinder_pos(5,3) = real_cylinder_info(i,9)

                    ghost_cylinder_pos(6,1) = real_cylinder_info(i,7)
                    ghost_cylinder_pos(6,2) = y_pos
                    ghost_cylinder_pos(6,3) = z_pos

                    ghost_cylinder_pos(7,1) = real_cylinder_info(i,7)
                    ghost_cylinder_pos(7,2) = real_cylinder_info(i,8)
                    ghost_cylinder_pos(7,3) = z_pos

                    do j = 1, 7
                        do k = 1, total_node_num
                            is_part_in_x_range = .false.
                            is_part_in_y_range = .false.
                            is_part_in_z_range = .false.
                            x_pos = real_cylinder_mesh_x_pos(i,k) + ghost_cylinder_pos(j,1)
                            y_pos = real_cylinder_mesh_y_pos(i,k) + ghost_cylinder_pos(j,2)
                            z_pos = real_cylinder_mesh_z_pos(i,k) + ghost_cylinder_pos(j,3)
                            if ((x_pos .gt. -0.5 * model_x_len) .and. (x_pos .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                            if ((y_pos .gt. -0.5 * model_y_len) .and. (y_pos .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                            if ((z_pos .gt. -0.5 * model_z_len) .and. (z_pos .lt. 0.5 * model_z_len)) is_part_in_z_range = .true.
                            if ((is_part_in_x_range .eq. .true.) .and. (is_part_in_y_range .eq. .true.) .and. (is_part_in_z_range .eq. .true.)) then
                                ghost_cylinder_status(j) = 1
                                exit
                            end if
                        end do          
                    end do
                end if

                is_overlapped = .false.
                ! overlap, new real particle vs old real particles
                do j = 1, 9
                    cylinder_a_info(j) = real_cylinder_info(i,j)
                end do
                do j = 1, total_node_num
                    mesh_a_info(j,1) = real_cylinder_mesh_x_pos(i,j) + real_cylinder_info(i,7)
                    mesh_a_info(j,2) = real_cylinder_mesh_y_pos(i,j) + real_cylinder_info(i,8)
                    mesh_a_info(j,3) = real_cylinder_mesh_z_pos(i,j) + real_cylinder_info(i,9)
                end do   
                do j = 1, i-1
                    do k = 1, 9
                        cylinder_b_info(k) = real_cylinder_info(j,k)
                    end do
                    do k = 1, total_node_num
                        mesh_b_info(k,1) = real_cylinder_mesh_x_pos(j,k) + real_cylinder_info(j,7)
                        mesh_b_info(k,2) = real_cylinder_mesh_y_pos(j,k) + real_cylinder_info(j,8)
                        mesh_b_info(k,3) = real_cylinder_mesh_z_pos(j,k) + real_cylinder_info(j,9)
                    end do
                    is_overlapped = is_cylinder_overlapped(cylinder_a_info,cylinder_b_info,mesh_a_info,mesh_b_info)
                    if (is_overlapped .eq. .true.) exit
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new real particle vs old ghost particles 
                do j = 1, total_ghost_cylinder_num
                    do k = 1, 9
                        cylinder_b_info(k) = ghost_cylinder_info(j,k)
                    end do
                    do k = 1, total_node_num
                        mesh_b_info(k,1) = ghost_cylinder_mesh_x_pos(j,k) + ghost_cylinder_info(j,7)
                        mesh_b_info(k,2) = ghost_cylinder_mesh_y_pos(j,k) + ghost_cylinder_info(j,8)
                        mesh_b_info(k,3) = ghost_cylinder_mesh_z_pos(j,k) + ghost_cylinder_info(j,9)
                    end do
                    is_overlapped = is_cylinder_overlapped(cylinder_a_info,cylinder_b_info,mesh_a_info,mesh_b_info)
                    if (is_overlapped .eq. .true.) exit
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old real particles
                do j = 1, 7
                    if (ghost_cylinder_status(j) .eq. 1) then
                        do k = 1, 6
                            cylinder_a_info(k) = real_cylinder_info(i,k)
                        end do
                        do k = 7, 9
                            cylinder_a_info(k) = ghost_cylinder_pos(j,k-6)
                        end do
                        do k = 1, total_node_num
                            mesh_a_info(k,1) = real_cylinder_mesh_x_pos(i,k) + ghost_cylinder_pos(j,1)
                            mesh_a_info(k,2) = real_cylinder_mesh_y_pos(i,k) + ghost_cylinder_pos(j,2)
                            mesh_a_info(k,3) = real_cylinder_mesh_z_pos(i,k) + ghost_cylinder_pos(j,3)
                        end do   
                        do k = 1, i-1
                            do ii = 1, 9
                                cylinder_b_info(ii) = real_cylinder_info(k,ii)
                            end do
                            do ii = 1, total_node_num
                                mesh_b_info(ii,1) = real_cylinder_mesh_x_pos(k,ii) + real_cylinder_info(k,7)
                                mesh_b_info(ii,2) = real_cylinder_mesh_y_pos(k,ii) + real_cylinder_info(k,8)
                                mesh_b_info(ii,3) = real_cylinder_mesh_z_pos(k,ii) + real_cylinder_info(k,9)
                            end do
                            is_overlapped = is_cylinder_overlapped(cylinder_a_info,cylinder_b_info,mesh_a_info,mesh_b_info)
                            if (is_overlapped .eq. .true.) exit
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old ghost particles
                do j = 1, 7
                    if (ghost_cylinder_status(j) .eq. 1) then
                        do k = 1, 6
                            cylinder_a_info(k) = real_cylinder_info(i,k)
                        end do
                        do k = 7, 9
                            cylinder_a_info(k) = ghost_cylinder_pos(j,k-6)
                        end do
                        do k = 1, total_node_num
                            mesh_a_info(k,1) = real_cylinder_mesh_x_pos(i,k) + ghost_cylinder_pos(j,1)
                            mesh_a_info(k,2) = real_cylinder_mesh_y_pos(i,k) + ghost_cylinder_pos(j,2)
                            mesh_a_info(k,3) = real_cylinder_mesh_z_pos(i,k) + ghost_cylinder_pos(j,3)
                        end do   
                        do k = 1, total_ghost_cylinder_num
                            do ii = 1, 9
                                cylinder_b_info(ii) = ghost_cylinder_info(k,ii)
                            end do
                            do ii = 1, total_node_num
                                mesh_b_info(ii,1) = ghost_cylinder_mesh_x_pos(k,ii) + ghost_cylinder_info(k,7)
                                mesh_b_info(ii,2) = ghost_cylinder_mesh_y_pos(k,ii) + ghost_cylinder_info(k,8)
                                mesh_b_info(ii,3) = ghost_cylinder_mesh_z_pos(k,ii) + ghost_cylinder_info(k,9)
                            end do
                            is_overlapped = is_cylinder_overlapped(cylinder_a_info,cylinder_b_info,mesh_a_info,mesh_b_info)
                            if (is_overlapped .eq. .true.) exit
                        end do                
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    ! update:
                    ! total_ghost_cylinder_num
                    ! ghost_cylinder_info
                    ! ghost_cylinder_mesh_x_pos
                    ! ghost_cylinder_mesh_y_pos
                    ! ghost_cylinder_mesh_z_pos
                    do j = 1, 7
                        if (ghost_cylinder_status(j) .eq. 1) then 
                            total_ghost_cylinder_num = total_ghost_cylinder_num + 1
                            ghost_cylinder_info(total_ghost_cylinder_num,1:6) = real_cylinder_info(i,1:6)
                            ghost_cylinder_info(total_ghost_cylinder_num,7:9) = ghost_cylinder_pos(j,1:3)
                            ghost_cylinder_mesh_x_pos(total_ghost_cylinder_num,:) = real_cylinder_mesh_x_pos(i,:)
                            ghost_cylinder_mesh_y_pos(total_ghost_cylinder_num,:) = real_cylinder_mesh_y_pos(i,:)
                            ghost_cylinder_mesh_z_pos(total_ghost_cylinder_num,:) = real_cylinder_mesh_z_pos(i,:)
                        end if
                    end do
                    if (mod(i,nint(0.25 * total_real_cylinder_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_real_cylinder_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                    i = i + 1
                end if

            end do

            if (allocated(cylinder_particle_info)) deallocate(cylinder_particle_info)
            allocate(cylinder_particle_info(total_real_cylinder_num + total_ghost_cylinder_num,9))
            do i = 1, total_real_cylinder_num
                cylinder_particle_info(i,1:9) = real_cylinder_info(i,1:9)
            end do
            do i = 1, total_ghost_cylinder_num
                cylinder_particle_info(total_real_cylinder_num+i,1:9) = ghost_cylinder_info(i,1:9)
            end do
            write(*,'(5X,a15,i6,a18,i6,a10,i6)') 'Real particles:',total_real_cylinder_num,'. Ghost particles:',total_ghost_cylinder_num,'. Totally:',total_real_cylinder_num+total_ghost_cylinder_num
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate capsule particles
        ! capsule_particle_info    size:   [total_particle_num,8]
        ! each row contains:     radius   height    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 is used here
        ! height / (radius * 2) = aspect_ratio1     (the height does not include the spherical tips!!!!)
        subroutine capsule_without_periodic_bc(capsule_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  capsule_particle_info(:,:)
            real(kind=8),allocatable    ::  mesh_a_info(:,:), mesh_b_info(:,:)
            real(kind=8),allocatable    ::  capsule_particle_mesh_x_pos(:,:), capsule_particle_mesh_y_pos(:,:), capsule_particle_mesh_z_pos(:,:)
            real(kind=8)                ::  capsule_a_info(8), capsule_b_info(8), r(1), model_x_len, model_y_len, model_z_len
            integer(kind=4)             ::  total_particle_num, total_alpha_segment, total_beta_segment, total_height_segment, total_node_num, i, j, k, progress
            logical                     ::  is_out_region, is_overlapped

            model_x_len = container_size1
            model_y_len = container_size2
            model_z_len = container_size3

            call random_seed()

            call capsule_shape(capsule_particle_info,model_x_len*model_y_len*model_z_len)
            total_particle_num = size(capsule_particle_info,1)

            capsule_particle_info(1:total_particle_num,6:8) = 0.0        ! initiate the position of capsule

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            total_alpha_segment = 10
            total_beta_segment = 5
            total_height_segment = 20
            total_node_num = total_alpha_segment * total_beta_segment * 2.0 + (total_height_segment - 2) * total_alpha_segment

            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(capsule_particle_mesh_x_pos)) deallocate(capsule_particle_mesh_x_pos)
            allocate(capsule_particle_mesh_x_pos(total_particle_num,total_node_num))
            capsule_particle_mesh_x_pos = 0.0

            if (allocated(capsule_particle_mesh_y_pos)) deallocate(capsule_particle_mesh_y_pos)
            allocate(capsule_particle_mesh_y_pos(total_particle_num,total_node_num))
            capsule_particle_mesh_y_pos = 0.0

            if (allocated(capsule_particle_mesh_z_pos)) deallocate(capsule_particle_mesh_z_pos)
            allocate(capsule_particle_mesh_z_pos(total_particle_num,total_node_num))
            capsule_particle_mesh_z_pos = 0.0
            do i = 1, total_particle_num
                do j = 1, 8
                    capsule_a_info(j) = capsule_particle_info(i,j)
                end do
                call generate_mesh_on_capsule_surface(capsule_a_info,total_alpha_segment,total_beta_segment,total_height_segment,mesh_a_info)
                do j = 1, total_node_num
                    capsule_particle_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the capsule is located at 0,0,0
                    capsule_particle_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    capsule_particle_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            i = 1
            do while (i .le. total_particle_num)

                is_out_region = .false.
                is_overlapped = .false.

                call random_number(r)
                capsule_particle_info(i,6) = r(1) * model_x_len - 0.5 * model_x_len
                
                call random_number(r)
                capsule_particle_info(i,7) = r(1) * model_y_len - 0.5 * model_y_len

                call random_number(r)
                capsule_particle_info(i,8) = r(1) * model_z_len - 0.5 * model_z_len
                
                do j = 1, total_node_num
                    mesh_a_info(j,1) = capsule_particle_mesh_x_pos(i,j) + capsule_particle_info(i,6)
                    mesh_a_info(j,2) = capsule_particle_mesh_y_pos(i,j) + capsule_particle_info(i,7)
                    mesh_a_info(j,3) = capsule_particle_mesh_z_pos(i,j) + capsule_particle_info(i,8)
                end do
                is_out_region = is_geometry_out_of_box_region(mesh_a_info,-0.5*model_x_len,0.5*model_x_len,-0.5*model_y_len,0.5*model_y_len,-0.5*model_z_len,0.5*model_z_len)
                if (is_out_region .eq. .true.) cycle

                do j = 1, 8
                    capsule_a_info(j) = capsule_particle_info(i,j)
                end do

                do j = 1, i-1
                    do k = 1, 8
                        capsule_b_info(k) = capsule_particle_info(j,k)
                    end do

                    do k = 1, total_node_num
                        mesh_b_info(k,1) = capsule_particle_mesh_x_pos(j,k) + capsule_particle_info(j,7)
                        mesh_b_info(k,2) = capsule_particle_mesh_y_pos(j,k) + capsule_particle_info(j,8)
                        mesh_b_info(k,3) = capsule_particle_mesh_z_pos(j,k) + capsule_particle_info(j,9)
                    end do

                    is_overlapped = is_capsule_overlapped(capsule_a_info,capsule_b_info,mesh_a_info,mesh_b_info)

                    if (is_overlapped .eq. .true.) exit

                end do

                if ((is_out_region .eq. .false.) .and. (is_overlapped .eq. .false.)) then
                    i = i + 1;
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate capsule particles
        ! capsule_particle_info    size:   [total_particle_num,8]
        ! each row contains:     radius   height    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 is used here
        ! height / (radius * 2) = aspect_ratio1     (the height does not include the spherical tips!!!!)
        ! for cylinder container
        subroutine capsule_in_cylinder_container(capsule_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  capsule_particle_info(:,:)
            real(kind=8),allocatable    ::  mesh_a_info(:,:), mesh_b_info(:,:)
            real(kind=8),allocatable    ::  capsule_particle_mesh_x_pos(:,:), capsule_particle_mesh_y_pos(:,:), capsule_particle_mesh_z_pos(:,:)
            real(kind=8)                ::  capsule_a_info(8), capsule_b_info(8), r(1), container_radius, container_height
            integer(kind=4)             ::  total_particle_num, total_alpha_segment, total_beta_segment, total_height_segment, total_node_num, i, j, k, progress
            logical                     ::  is_out_region, is_overlapped

            container_radius = container_size1
            container_height = container_size2

            call random_seed()

            call capsule_shape(capsule_particle_info,container_radius**2.0*3.14159*container_height)
            total_particle_num = size(capsule_particle_info,1)

            capsule_particle_info(1:total_particle_num,6:8) = 0.0        ! initiate the position of capsule

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            total_alpha_segment = 10
            total_beta_segment = 5
            total_height_segment = 20
            total_node_num = total_alpha_segment * total_beta_segment * 2.0 + (total_height_segment - 2) * total_alpha_segment

            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(capsule_particle_mesh_x_pos)) deallocate(capsule_particle_mesh_x_pos)
            allocate(capsule_particle_mesh_x_pos(total_particle_num,total_node_num))
            capsule_particle_mesh_x_pos = 0.0

            if (allocated(capsule_particle_mesh_y_pos)) deallocate(capsule_particle_mesh_y_pos)
            allocate(capsule_particle_mesh_y_pos(total_particle_num,total_node_num))
            capsule_particle_mesh_y_pos = 0.0

            if (allocated(capsule_particle_mesh_z_pos)) deallocate(capsule_particle_mesh_z_pos)
            allocate(capsule_particle_mesh_z_pos(total_particle_num,total_node_num))
            capsule_particle_mesh_z_pos = 0.0
            do i = 1, total_particle_num
                do j = 1, 8
                    capsule_a_info(j) = capsule_particle_info(i,j)
                end do
                call generate_mesh_on_capsule_surface(capsule_a_info,total_alpha_segment,total_beta_segment,total_height_segment,mesh_a_info)
                do j = 1, total_node_num
                    capsule_particle_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the capsule is located at 0,0,0
                    capsule_particle_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    capsule_particle_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            i = 1
            do while (i .le. total_particle_num)

                is_out_region = .false.
                is_overlapped = .false.

                call random_number(r)
                capsule_particle_info(i,6) = r(1) * 2 * container_radius - container_radius
                
                call random_number(r)
                capsule_particle_info(i,7) = r(1) * 2 * container_radius - container_radius

                call random_number(r)
                capsule_particle_info(i,8) = r(1) * container_height - 0.5 * container_height
                
                do j = 1, total_node_num
                    mesh_a_info(j,1) = capsule_particle_mesh_x_pos(i,j) + capsule_particle_info(i,6)
                    mesh_a_info(j,2) = capsule_particle_mesh_y_pos(i,j) + capsule_particle_info(i,7)
                    mesh_a_info(j,3) = capsule_particle_mesh_z_pos(i,j) + capsule_particle_info(i,8)
                    if ((mesh_a_info(j,3) .lt. -0.5 * container_height) .or. (mesh_a_info(j,3) .gt. 0.5 * container_height)) then
                        is_out_region = .true.
                        exit
                    end if
                    if (mesh_a_info(j,1) ** 2.0 + mesh_a_info(j,2) ** 2.0 .gt. container_radius ** 2.0) then
                        is_out_region = .true.
                        exit
                    end if
                end do
                if (is_out_region .eq. .true.) cycle

                do j = 1, 8
                    capsule_a_info(j) = capsule_particle_info(i,j)
                end do

                do j = 1, i-1
                    do k = 1, 8
                        capsule_b_info(k) = capsule_particle_info(j,k)
                    end do

                    do k = 1, total_node_num
                        mesh_b_info(k,1) = capsule_particle_mesh_x_pos(j,k) + capsule_particle_info(j,7)
                        mesh_b_info(k,2) = capsule_particle_mesh_y_pos(j,k) + capsule_particle_info(j,8)
                        mesh_b_info(k,3) = capsule_particle_mesh_z_pos(j,k) + capsule_particle_info(j,9)
                    end do

                    is_overlapped = is_capsule_overlapped(capsule_a_info,capsule_b_info,mesh_a_info,mesh_b_info)

                    if (is_overlapped .eq. .true.) exit

                end do

                if ((is_out_region .eq. .false.) .and. (is_overlapped .eq. .false.)) then
                    i = i + 1;
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate capsule particles
        ! capsule_particle_info    size:   [total_particle_num,8]
        ! each row contains:     radius   height    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 is used here
        ! height / (radius * 2) = aspect_ratio1     (the height does not include the spherical tips!!!!)
        subroutine capsule_with_periodic_bc(capsule_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  capsule_particle_info(:,:)
            real(kind=8),allocatable    ::  mesh_a_info(:,:), mesh_b_info(:,:), real_capsule_info(:,:), ghost_capsule_info(:,:)
            real(kind=8),allocatable    ::  real_capsule_mesh_x_pos(:,:), real_capsule_mesh_y_pos(:,:), real_capsule_mesh_z_pos(:,:)
            real(kind=8),allocatable    ::  ghost_capsule_mesh_x_pos(:,:), ghost_capsule_mesh_y_pos(:,:), ghost_capsule_mesh_z_pos(:,:)
            real(kind=8)                ::  capsule_a_info(8), capsule_b_info(8), r(1), ghost_capsule_pos(7,3), x_pos, y_pos, z_pos, model_x_len, model_y_len, model_z_len
            integer(kind=4)             ::  total_real_capsule_num, total_ghost_capsule_num, total_alpha_segment, total_beta_segment, total_height_segment
            integer(kind=4)             ::  total_node_num, i, j, k, ii, progress, ghost_capsule_status(7)
            logical                     ::  is_out_region, is_overlapped, is_part_in_x_range, is_part_in_y_range, is_part_in_z_range

            model_x_len = container_size1
            model_y_len = container_size2
            model_z_len = container_size3

            call random_seed()

            call capsule_shape(real_capsule_info,model_x_len*model_y_len*model_z_len)
            total_real_capsule_num = size(real_capsule_info,1)

            real_capsule_info(1:total_real_capsule_num,6:8) = 0.0        ! initiate the position of capsule


            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a17)') 'Particle shape done.','(',total_real_capsule_num,'real particles).'

            total_alpha_segment = 10
            total_beta_segment = 5
            total_height_segment = 20
            total_node_num = total_alpha_segment * total_beta_segment * 2.0 + (total_height_segment - 2) * total_alpha_segment
        
            if (allocated(ghost_capsule_info)) deallocate(ghost_capsule_info)
            allocate(ghost_capsule_info(total_real_capsule_num*7,9))
            ghost_capsule_info = 0.0
            
            if (allocated(mesh_a_info)) deallocate(mesh_a_info)
            allocate(mesh_a_info(total_node_num,3))
            mesh_a_info = 0.0

            if (allocated(mesh_b_info)) deallocate(mesh_b_info)
            allocate(mesh_b_info(total_node_num,3))
            mesh_b_info = 0.0

            if (allocated(real_capsule_mesh_x_pos)) deallocate(real_capsule_mesh_x_pos)
            allocate(real_capsule_mesh_x_pos(total_real_capsule_num,total_node_num))
            real_capsule_mesh_x_pos = 0.0

            if (allocated(real_capsule_mesh_y_pos)) deallocate(real_capsule_mesh_y_pos)
            allocate(real_capsule_mesh_y_pos(total_real_capsule_num,total_node_num))
            real_capsule_mesh_y_pos = 0.0

            if (allocated(real_capsule_mesh_z_pos)) deallocate(real_capsule_mesh_z_pos)
            allocate(real_capsule_mesh_z_pos(total_real_capsule_num,total_node_num))
            real_capsule_mesh_z_pos = 0.0

            if (allocated(ghost_capsule_mesh_x_pos)) deallocate(ghost_capsule_mesh_x_pos)
            allocate(ghost_capsule_mesh_x_pos(total_real_capsule_num*7,total_node_num))
            ghost_capsule_mesh_x_pos = 0.0

            if (allocated(ghost_capsule_mesh_y_pos)) deallocate(ghost_capsule_mesh_y_pos)
            allocate(ghost_capsule_mesh_y_pos(total_real_capsule_num*7,total_node_num))
            ghost_capsule_mesh_y_pos = 0.0

            if (allocated(ghost_capsule_mesh_z_pos)) deallocate(ghost_capsule_mesh_z_pos)
            allocate(ghost_capsule_mesh_z_pos(total_real_capsule_num*7,total_node_num))

            do i = 1, total_real_capsule_num
                do j = 1, 8
                    capsule_a_info(j) = real_capsule_info(i,j)
                end do
                call generate_mesh_on_capsule_surface(capsule_a_info,total_alpha_segment,total_beta_segment,total_height_segment,mesh_a_info)
                do j = 1, total_node_num
                    real_capsule_mesh_x_pos(i,j) = mesh_a_info(j,1)             ! mesh coordinate when the center of the capsule is located at 0,0,0
                    real_capsule_mesh_y_pos(i,j) = mesh_a_info(j,2)
                    real_capsule_mesh_z_pos(i,j) = mesh_a_info(j,3)
                end do
            end do

            total_ghost_capsule_num = 0
            i = 1
            do while (i .le. total_real_capsule_num)

                call random_number(r)
                real_capsule_info(i,6) = r(1) * model_x_len - 0.5 * model_x_len
                
                call random_number(r)
                real_capsule_info(i,7) = r(1) * model_y_len - 0.5 * model_y_len

                call random_number(r)
                real_capsule_info(i,8) = r(1) * model_z_len - 0.5 * model_z_len

                ! check whether the current real capsule is intersected by the boundary of model
                ! if yes, generate ghost capsules
                is_out_region = .false.
                do j = 1, total_node_num
                    x_pos = real_capsule_mesh_x_pos(i,j) + real_capsule_info(i,6)
                    y_pos = real_capsule_mesh_y_pos(i,j) + real_capsule_info(i,7)
                    z_pos = real_capsule_mesh_z_pos(i,j) + real_capsule_info(i,8)
                    if ((x_pos .lt. -0.5 * model_x_len) .or. (x_pos .gt. 0.5 * model_x_len)) then
                        is_out_region = .true.
                        exit
                    end if 
                    if ((y_pos .lt. -0.5 * model_y_len) .or. (y_pos .gt. 0.5 * model_y_len)) then
                        is_out_region = .true.
                        exit
                    end if
                    if ((z_pos .lt. -0.5 * model_z_len) .or. (z_pos .gt. 0.5 * model_z_len)) then
                        is_out_region = .true.
                        exit
                    end if 
                end do
                
                ghost_capsule_status = 0
                if (is_out_region .eq. .true.) then
                    ! generate ghost particle based on the current real particle
                    ! and detect whether the ghost particle is intersected with 
                    ! the region
                    ! for capsule, seven possible ghost particles exist
                    if (real_capsule_info(i,6) .gt. 0.0) then
                        x_pos = real_capsule_info(i,6) - model_x_len
                    else
                        x_pos = real_capsule_info(i,6) + model_x_len
                    end if
                    if (real_capsule_info(i,7) .gt. 0.0) then
                        y_pos = real_capsule_info(i,7) - model_y_len
                    else
                        y_pos = real_capsule_info(i,7) + model_y_len
                    end if
                    if (real_capsule_info(i,8) .gt. 0.0) then
                        z_pos = real_capsule_info(i,8) - model_z_len
                    else
                        z_pos = real_capsule_info(i,8) + model_z_len
                    end if

                    ghost_capsule_pos(1,1) = x_pos
                    ghost_capsule_pos(1,2) = real_capsule_info(i,7)
                    ghost_capsule_pos(1,3) = real_capsule_info(i,8)

                    ghost_capsule_pos(2,1) = x_pos
                    ghost_capsule_pos(2,2) = y_pos
                    ghost_capsule_pos(2,3) = real_capsule_info(i,8)

                    ghost_capsule_pos(3,1) = x_pos
                    ghost_capsule_pos(3,2) = y_pos
                    ghost_capsule_pos(3,3) = z_pos

                    ghost_capsule_pos(4,1) = x_pos
                    ghost_capsule_pos(4,2) = real_capsule_info(i,7)
                    ghost_capsule_pos(4,3) = z_pos

                    ghost_capsule_pos(5,1) = real_capsule_info(i,6)
                    ghost_capsule_pos(5,2) = y_pos
                    ghost_capsule_pos(5,3) = real_capsule_info(i,8)

                    ghost_capsule_pos(6,1) = real_capsule_info(i,6)
                    ghost_capsule_pos(6,2) = y_pos
                    ghost_capsule_pos(6,3) = z_pos

                    ghost_capsule_pos(7,1) = real_capsule_info(i,6)
                    ghost_capsule_pos(7,2) = real_capsule_info(i,7)
                    ghost_capsule_pos(7,3) = z_pos

                    do j = 1, 7
                        do k = 1, total_node_num
                            is_part_in_x_range = .false.
                            is_part_in_y_range = .false.
                            is_part_in_z_range = .false.
                            x_pos = real_capsule_mesh_x_pos(i,k) + ghost_capsule_pos(j,1)
                            y_pos = real_capsule_mesh_y_pos(i,k) + ghost_capsule_pos(j,2)
                            z_pos = real_capsule_mesh_z_pos(i,k) + ghost_capsule_pos(j,3)
                            if ((x_pos .gt. -0.5 * model_x_len) .and. (x_pos .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                            if ((y_pos .gt. -0.5 * model_y_len) .and. (y_pos .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                            if ((z_pos .gt. -0.5 * model_z_len) .and. (z_pos .lt. 0.5 * model_z_len)) is_part_in_z_range = .true.
                            if ((is_part_in_x_range .eq. .true.) .and. (is_part_in_y_range .eq. .true.) .and. (is_part_in_z_range .eq. .true.)) then
                                ghost_capsule_status(j) = 1
                                exit
                            end if
                        end do          
                    end do
                end if

                is_overlapped = .false.
                ! overlap, new real particle vs old real particles
                do j = 1, 8
                    capsule_a_info(j) = real_capsule_info(i,j)
                end do
                do j = 1, total_node_num
                    mesh_a_info(j,1) = real_capsule_mesh_x_pos(i,j) + real_capsule_info(i,6)
                    mesh_a_info(j,2) = real_capsule_mesh_y_pos(i,j) + real_capsule_info(i,7)
                    mesh_a_info(j,3) = real_capsule_mesh_z_pos(i,j) + real_capsule_info(i,8)
                end do   
                do j = 1, i-1
                    do k = 1, 8
                        capsule_b_info(k) = real_capsule_info(j,k)
                    end do
                    do k = 1, total_node_num
                        mesh_b_info(k,1) = real_capsule_mesh_x_pos(j,k) + real_capsule_info(j,6)
                        mesh_b_info(k,2) = real_capsule_mesh_y_pos(j,k) + real_capsule_info(j,7)
                        mesh_b_info(k,3) = real_capsule_mesh_z_pos(j,k) + real_capsule_info(j,8)
                    end do
                    is_overlapped = is_capsule_overlapped(capsule_a_info,capsule_b_info,mesh_a_info,mesh_b_info)
                    if (is_overlapped .eq. .true.) exit
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new real particle vs old ghost particles 
                do j = 1, total_ghost_capsule_num
                    do k = 1, 8
                        capsule_b_info(k) = ghost_capsule_info(j,k)
                    end do
                    do k = 1, total_node_num
                        mesh_b_info(k,1) = ghost_capsule_mesh_x_pos(j,k) + ghost_capsule_info(j,6)
                        mesh_b_info(k,2) = ghost_capsule_mesh_y_pos(j,k) + ghost_capsule_info(j,7)
                        mesh_b_info(k,3) = ghost_capsule_mesh_z_pos(j,k) + ghost_capsule_info(j,8)
                    end do
                    is_overlapped = is_capsule_overlapped(capsule_a_info,capsule_b_info,mesh_a_info,mesh_b_info)
                    if (is_overlapped .eq. .true.) exit
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old real particles
                do j = 1, 7
                    if (ghost_capsule_status(j) .eq. 1) then
                        do k = 1, 5
                            capsule_a_info(k) = real_capsule_info(i,k)
                        end do
                        do k = 6, 8
                            capsule_a_info(k) = ghost_capsule_pos(j,k-5)
                        end do
                        do k = 1, total_node_num
                            mesh_a_info(k,1) = real_capsule_mesh_x_pos(i,k) + ghost_capsule_pos(j,1)
                            mesh_a_info(k,2) = real_capsule_mesh_y_pos(i,k) + ghost_capsule_pos(j,2)
                            mesh_a_info(k,3) = real_capsule_mesh_z_pos(i,k) + ghost_capsule_pos(j,3)
                        end do   
                        do k = 1, i-1
                            do ii = 1, 8
                                capsule_b_info(ii) = real_capsule_info(k,ii)
                            end do
                            do ii = 1, total_node_num
                                mesh_b_info(ii,1) = real_capsule_mesh_x_pos(k,ii) + real_capsule_info(k,6)
                                mesh_b_info(ii,2) = real_capsule_mesh_y_pos(k,ii) + real_capsule_info(k,7)
                                mesh_b_info(ii,3) = real_capsule_mesh_z_pos(k,ii) + real_capsule_info(k,8)
                            end do
                            is_overlapped = is_capsule_overlapped(capsule_a_info,capsule_b_info,mesh_a_info,mesh_b_info)
                            if (is_overlapped .eq. .true.) exit
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old ghost particles
                do j = 1, 7
                    if (ghost_capsule_status(j) .eq. 1) then
                        do k = 1, 5
                            capsule_a_info(k) = real_capsule_info(i,k)
                        end do
                        do k = 6, 8
                            capsule_a_info(k) = ghost_capsule_pos(j,k-5)
                        end do
                        do k = 1, total_node_num
                            mesh_a_info(k,1) = real_capsule_mesh_x_pos(i,k) + ghost_capsule_pos(j,1)
                            mesh_a_info(k,2) = real_capsule_mesh_y_pos(i,k) + ghost_capsule_pos(j,2)
                            mesh_a_info(k,3) = real_capsule_mesh_z_pos(i,k) + ghost_capsule_pos(j,3)
                        end do   
                        do k = 1, total_ghost_capsule_num
                            do ii = 1, 8
                                capsule_b_info(ii) = ghost_capsule_info(k,ii)
                            end do
                            do ii = 1, total_node_num
                                mesh_b_info(ii,1) = ghost_capsule_mesh_x_pos(k,ii) + ghost_capsule_info(k,6)
                                mesh_b_info(ii,2) = ghost_capsule_mesh_y_pos(k,ii) + ghost_capsule_info(k,7)
                                mesh_b_info(ii,3) = ghost_capsule_mesh_z_pos(k,ii) + ghost_capsule_info(k,8)
                            end do
                            is_overlapped = is_capsule_overlapped(capsule_a_info,capsule_b_info,mesh_a_info,mesh_b_info)
                            if (is_overlapped .eq. .true.) exit
                        end do                
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    ! update:
                    ! total_ghost_capsule_num
                    ! ghost_capsule_info
                    ! ghost_capsule_mesh_x_pos
                    ! ghost_capsule_mesh_y_pos
                    ! ghost_capsule_mesh_z_pos
                    do j = 1, 7
                        if (ghost_capsule_status(j) .eq. 1) then 
                            total_ghost_capsule_num = total_ghost_capsule_num + 1
                            ghost_capsule_info(total_ghost_capsule_num,1:5) = real_capsule_info(i,1:5)
                            ghost_capsule_info(total_ghost_capsule_num,6:8) = ghost_capsule_pos(j,1:3)
                            ghost_capsule_mesh_x_pos(total_ghost_capsule_num,:) = real_capsule_mesh_x_pos(i,:)
                            ghost_capsule_mesh_y_pos(total_ghost_capsule_num,:) = real_capsule_mesh_y_pos(i,:)
                            ghost_capsule_mesh_z_pos(total_ghost_capsule_num,:) = real_capsule_mesh_z_pos(i,:)
                        end if
                    end do
                    if (mod(i,nint(0.25 * total_real_capsule_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_real_capsule_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                    i = i + 1
                end if

            end do

            if (allocated(capsule_particle_info)) deallocate(capsule_particle_info)
            allocate(capsule_particle_info(total_real_capsule_num + total_ghost_capsule_num,8))
            do i = 1, total_real_capsule_num
                capsule_particle_info(i,1:8) = real_capsule_info(i,1:8)
            end do
            do i = 1, total_ghost_capsule_num
                capsule_particle_info(total_real_capsule_num+i,1:8) = ghost_capsule_info(i,1:8)
            end do
            write(*,'(5X,a15,i6,a18,i6,a10,i6)') 'Real particles:',total_real_capsule_num,'. Ghost particles:',total_ghost_capsule_num,'. Totally:',total_real_capsule_num+total_ghost_capsule_num
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate shape of cylinder particles
        ! cylinder_particle_info    size:   [total_particle_num,9]
        ! each row contains:     radius_a   radius_b   height    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 and aspect_ratio2 are used here
        ! height / (radius_a * 2) = aspect_ratio1
        ! radius_a / radius_b = aspect_ratio2
        ! ONLY the shape of each cylinder is determined (x_pos, y_pos and z_pos are undefined)
        subroutine cylinder_shape(cylinder_particle_info,container_volume)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  cylinder_particle_info(:,:)
            real(kind=8),intent(in)     ::  container_volume
            real(kind=8)                ::  r(1), low_value, high_value
            integer(kind=4)             ::  total_particle_num, i, j
            character(len=100)          ::  cylinder_radius_file
            real(kind=8),allocatable    ::  particle_needed_vol(:), particle_real_vol(:)
            real(kind=8)                ::  last_particle_vol, particle_real_vol_in_current_sieve, radius_a, radius_b, height
            
            if (allocated(particle_needed_vol)) deallocate(particle_needed_vol)
            if (allocated(particle_real_vol)) deallocate(particle_real_vol)
            allocate(particle_needed_vol(total_sieve_num-1))
            allocate(particle_real_vol(total_sieve_num-1))
            
            do i = 1, total_sieve_num-1
                particle_needed_vol(i) = (gradation(2,i+1) - gradation(2,i)) * particle_occupy_percentage * container_volume
            end do
            
            ! randomize
            call random_seed()
 
            cylinder_radius_file='cylinder_radius_file.tmp'   
            open(unit=100,file=cylinder_radius_file,form='formatted',status='replace',action='write')
            particle_real_vol = 0.0
            total_particle_num = 0
            do i = total_sieve_num-1, 1, -1
                particle_real_vol_in_current_sieve = 0.0 
                do while (particle_real_vol(i) .lt. particle_needed_vol(i))
                    ! generate the cylinder
                    call random_number(r)
                    low_value = gradation(1,i) / 2.0
                    high_value = gradation(1,i+1) / 2.0
                    radius_a = r(1) * (high_value - low_value) + low_value
                    radius_b = radius_a / aspect_ratio2
                    height = radius_a * 2 * aspect_ratio1
                    particle_real_vol_in_current_sieve = particle_real_vol_in_current_sieve + radius_a * radius_b * 3.14159 * height
                    if (particle_real_vol_in_current_sieve .gt. particle_needed_vol(i)) then
                        last_particle_vol = particle_needed_vol(i) - particle_real_vol(i)
                        radius_a = (last_particle_vol * aspect_ratio2 / 2 / 3.14159 / aspect_ratio1) ** (1.0 / 3.0)
                        radius_b = radius_a / aspect_ratio2
                        height = radius_a * 2 * aspect_ratio1
                    end if
                    write(100,'(e15.5)') radius_a
                    total_particle_num = total_particle_num + 1
                    particle_real_vol(i) = particle_real_vol(i) + radius_a * radius_b * 3.14159 * height
                end do
            end do
            close(unit=100)

            if (allocated(cylinder_particle_info)) deallocate(cylinder_particle_info)
            allocate(cylinder_particle_info(total_particle_num,9))
            cylinder_particle_info = 0.0

            open(unit=100,file=cylinder_radius_file,form='formatted',status='old',action='read')
            do i = 1, total_particle_num
                read(100,*) cylinder_particle_info(i,1)
                cylinder_particle_info(i,2) = cylinder_particle_info(i,1) / aspect_ratio2
                cylinder_particle_info(i,3) = cylinder_particle_info(i,1) * 2 * aspect_ratio1
            end do
            close(unit=100,status='delete')

            i = 1
            call sort(cylinder_particle_info,i,'d')

            ! generate the rotation angle
            do i = 1, total_particle_num

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                cylinder_particle_info(i,4) = r(1)

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                cylinder_particle_info(i,5) = r(1)

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                cylinder_particle_info(i,6) = r(1)
            end do

        end subroutine

        ! generate shape of capsule particles
        ! capsule_particle_info    size:   [total_particle_num,8]
        ! each row contains:     radius   height    angle_x     angle_y     angle_z      x_pos       y_pos      z_pos
        ! the aspect_ratio1 is used here
        ! height / (radius * 2) = aspect_ratio1     (the height does not include the spherical tips!!!!)
        ! ONLY the shape of each capsule is determined (x_pos, y_pos and z_pos are undefined)
        subroutine capsule_shape(capsule_particle_info,container_volume)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  capsule_particle_info(:,:)
            real(kind=8),intent(in)     ::  container_volume
            real(kind=8),allocatable    ::  ref_sphere_radius(:)
            real(kind=8)                ::  r(1), low_value, high_value
            character(len=100)          ::  capsule_radius_file
            integer(kind=4)             ::  total_particle_num, i, j
            real(kind=8),allocatable    ::  particle_needed_vol(:), particle_real_vol(:)
            real(kind=8)                ::  last_particle_vol, particle_real_vol_in_current_sieve, radius, height
            
            if (allocated(particle_needed_vol)) deallocate(particle_needed_vol)
            if (allocated(particle_real_vol)) deallocate(particle_real_vol)
            allocate(particle_needed_vol(total_sieve_num-1))
            allocate(particle_real_vol(total_sieve_num-1))
            
            do i = 1, total_sieve_num-1
                particle_needed_vol(i) = (gradation(2,i+1) - gradation(2,i)) * particle_occupy_percentage * container_volume
            end do
            
            ! randomize
            call random_seed()
 
            capsule_radius_file='capsule_radius_file.tmp'   
            open(unit=100,file=capsule_radius_file,form='formatted',status='replace',action='write')
            particle_real_vol = 0.0
            total_particle_num = 0
            do i = total_sieve_num-1, 1, -1
                particle_real_vol_in_current_sieve = 0.0 
                do while (particle_real_vol(i) .lt. particle_needed_vol(i))
                    ! generate the capsule
                    call random_number(r)
                    low_value = gradation(1,i) / 2.0
                    high_value = gradation(1,i+1) / 2.0
                    radius = r(1) * (high_value - low_value) + low_value
                    height = radius * 2 * aspect_ratio1
                    particle_real_vol_in_current_sieve = particle_real_vol_in_current_sieve + radius ** 3.0 * 3.14159 * 4.0 / 3.0 + radius ** 2.0 * 3.14159 * height
                    if (particle_real_vol_in_current_sieve .gt. particle_needed_vol(i)) then
                        last_particle_vol = particle_needed_vol(i) - particle_real_vol(i)
                        radius = (last_particle_vol / (4.0 * 3.14159 / 3.0 + 2.0 * 3.14159 *aspect_ratio1)) ** (1.0 / 3.0)
                        height = radius * 2 * aspect_ratio1
                    end if
                    write(100,'(e15.5)') radius
                    total_particle_num = total_particle_num + 1
                    particle_real_vol(i) = particle_real_vol(i) + radius ** 3.0 * 3.14159 * 4.0 / 3.0 + radius ** 2.0 * 3.14159 * height
                end do
            end do
            close(unit=100)

            if (allocated(capsule_particle_info)) deallocate(capsule_particle_info)
            allocate(capsule_particle_info(total_particle_num,8))
            capsule_particle_info = 0.0

            open(unit=100,file=capsule_radius_file,form='formatted',status='old',action='read')
            do i = 1, total_particle_num
                read(100,*) capsule_particle_info(i,1)
                capsule_particle_info(i,2) = capsule_particle_info(i,1) * 2 * aspect_ratio1
            end do
            close(unit=100,status='delete')

            i = 1
            call sort(capsule_particle_info,i,'d')

            ! generate the rotation angle
            do i = 1, total_particle_num

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                capsule_particle_info(i,3) = r(1)

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                capsule_particle_info(i,4) = r(1)

                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                capsule_particle_info(i,5) = r(1)
            end do

        end subroutine

end module

!dec$ endif