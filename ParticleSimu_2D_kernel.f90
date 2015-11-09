!dec$ if .not. defined(__ParticleSimu_2D_kernel_f90)
!dec$ define __ParticleSimu_2D_kernel_f90

module ParticleSimu_2D_kernel
    implicit none
    contains
        
        ! generate the ref circle for polygon and ellipse, and also the circle itself
        ! ref_circle_radius     :   [total_particle_num]
        ! container_area        :   area of the container
        subroutine get_ref_circle_radius(ref_circle_radius,container_area)
            use ParticleSimu_global_parameter
            implicit none
            real(kind=8),allocatable    ::  ref_circle_radius(:)
            real(kind=8),intent(in)     ::  container_area
            character(len=100)          ::  circle_rad_file
            real(kind=8),allocatable    ::  particle_needed_area(:), particle_real_area(:)
            real(kind=8)                ::  last_particle_area, r(1), particle_real_area_in_current_sieve, low_value, high_value, temp
            integer(kind=4)             ::  total_particle_num, i, j, k
            
            if (allocated(particle_needed_area)) deallocate(particle_needed_area)
            if (allocated(particle_real_area)) deallocate(particle_real_area)
            allocate(particle_needed_area(total_sieve_num-1))
            allocate(particle_real_area(total_sieve_num-1))
            
            do i = 1, total_sieve_num-1
                particle_needed_area(i) = (gradation(2,i+1) - gradation(2,i)) * particle_occupy_percentage * container_area
            end do
            
            ! randomize
            call random_seed()
 
            circle_rad_file='circle_radius.tmp'   
            open(unit=100,file=circle_rad_file,form='formatted',status='replace',action='write')
            particle_real_area = 0.0
            total_particle_num = 0
            do i = total_sieve_num-1, 1, -1
                particle_real_area_in_current_sieve = 0.0 
                do while (particle_real_area(i) .lt. particle_needed_area(i))
                    ! generate the circle
                    call random_number(r)
                    low_value = gradation(1,i) / 2.0
                    high_value = gradation(1,i+1) / 2.0
                    r(1) = r(1) * (high_value - low_value) + low_value
                    particle_real_area_in_current_sieve = particle_real_area_in_current_sieve + r(1) ** 2.0 * 3.14159
                    if (particle_real_area_in_current_sieve .gt. particle_needed_area(i)) then
                        last_particle_area = particle_needed_area(i) - particle_real_area(i)
                        r(1) = (last_particle_area / 3.14159) ** 0.5
                    end if
                    write(100,'(e15.5)') r(1)
                    total_particle_num = total_particle_num + 1
                    particle_real_area(i) = particle_real_area(i) + r(1) ** 2.0 * 3.14159  
                end do
            end do
            close(unit=100)

            ! read radius from file 
            if (allocated(ref_circle_radius)) deallocate(ref_circle_radius)
            allocate(ref_circle_radius(total_particle_num))
            ref_circle_radius = 0.0
            open(unit=100,file=circle_rad_file,form='formatted',status='old',action='read')
            do i = 1, total_particle_num
                read(100,*) ref_circle_radius(i)
            end do
            close(unit=100,status='delete')

            ! sort
            do i = 1, total_particle_num - 1
                k = i
                do j = i+1, total_particle_num
                    if (ref_circle_radius(k) .lt. ref_circle_radius(j)) k = j
                end do
                if (k .ne. i) then
                    temp = ref_circle_radius(k)
                    ref_circle_radius(k) = ref_circle_radius(i)
                    ref_circle_radius(i) = temp
                end if
            end do
        end subroutine

        ! generate circular particles
        ! circular_particle_info    size:   [total_particle_num,3]
        ! each row contains:     radius      x_pos       y_pos
        subroutine circle_without_periodic_bc(circular_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  circular_particle_info(:,:)
            real(kind=8),allocatable    ::  ref_circle_radius(:)
            real(kind=8)                ::  r(1), low_value, high_value, distance, model_x_len, model_y_len
            integer(kind=4)             ::  total_particle_num, i, j, progress
            logical                     ::  is_overlapped

            model_x_len = container_size1
            model_y_len = container_size2

            call random_seed()

            call get_ref_circle_radius(ref_circle_radius,model_x_len*model_y_len)
            total_particle_num = size(ref_circle_radius)

            if (allocated(circular_particle_info)) deallocate(circular_particle_info)
            allocate(circular_particle_info(total_particle_num,3))
            circular_particle_info = 0.0
            do i = 1, total_particle_num
                circular_particle_info(i,1) = ref_circle_radius(i)
            end do
            
            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            ! generate the position of circular particles
            i = 1;
            do while (i .le. total_particle_num)
                is_overlapped = .false.

                call random_number(r)
                circular_particle_info(i,2) = r(1) * (model_x_len - 2 * circular_particle_info(i,1)) + circular_particle_info(i,1)
                circular_particle_info(i,2) = circular_particle_info(i,2) - 0.5 * model_x_len

                call random_number(r)
                circular_particle_info(i,3) = r(1) * (model_y_len - 2 * circular_particle_info(i,1)) + circular_particle_info(i,1)
                circular_particle_info(i,3) = circular_particle_info(i,3) - 0.5 * model_y_len

                do j = 1, i-1
                    distance = sqrt((circular_particle_info(i,2) - circular_particle_info(j,2)) ** 2.0 + (circular_particle_info(i,3) - circular_particle_info(j,3)) ** 2.0)
                    if (distance .lt. (circular_particle_info(i,1) + circular_particle_info(j,1))) then
                        is_overlapped = .true.
                        exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    i = i + 1;
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate circular particles
        ! circular_particle_info    size:   [total_particle_num,3]
        ! each row contains:     radius      x_pos       y_pos
        ! non periodic boundary condition
        ! for circular containder, the periodic
        ! boundary condition is not available
        subroutine circle_in_circle_container(circular_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  circular_particle_info(:,:)
            real(kind=8),allocatable    ::  ref_circle_radius(:)
            real(kind=8)                ::  r(1), low_value, high_value, distance, polar_angle, polar_radius
            integer(kind=4)             ::  total_particle_num, i, j, progress
            logical                     ::  is_overlapped, is_out_region
            real(kind=8)                ::  container_radius

            container_radius = container_size1
            call random_seed()

            call get_ref_circle_radius(ref_circle_radius,container_radius ** 2.0 * 3.14159)
            total_particle_num = size(ref_circle_radius)

            if (allocated(circular_particle_info)) deallocate(circular_particle_info)
            allocate(circular_particle_info(total_particle_num,3))
            circular_particle_info = 0.0
            do i = 1, total_particle_num
                circular_particle_info(i,1) = ref_circle_radius(i)
            end do
            
            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            ! generate the position of circular particles
            i = 1;
            do while (i .le. total_particle_num)
                is_out_region = .false.
                is_overlapped = .false.

                ! fisrt generate position in polar coordinate
                call random_number(r)
                polar_angle = r(1) * 2.0 * 3.14159

                call random_number(r)
                low_value = 0.0
                high_value = container_radius - circular_particle_info(i,1)
                polar_radius = r(1) * (high_value - low_value) + low_value

                circular_particle_info(i,2) = polar_radius * cos(polar_angle)
                circular_particle_info(i,3) = polar_radius * sin(polar_angle)


                do j = 1, i-1
                    distance = sqrt((circular_particle_info(i,2) - circular_particle_info(j,2)) ** 2.0 + (circular_particle_info(i,3) - circular_particle_info(j,3)) ** 2.0)
                    if (distance .lt. (circular_particle_info(i,1) + circular_particle_info(j,1))) then
                        is_overlapped = .true.
                        exit
                    end if
                end do

                if ((is_overlapped .eq. .false.) .and. (is_out_region .eq. .false.)) then
                    i = i + 1;
                    if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate circular particles
        ! circular_particle_info    size:   [total_particle_num,3]
        ! each row contains:     radius      x_pos       y_pos
        subroutine circle_with_periodic_bc(circular_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  circular_particle_info(:,:)
            real(kind=8),allocatable    ::  real_circle_info(:,:), ghost_circle_info(:,:)
            real(kind=8),allocatable    ::  ref_circle_radius(:)
            real(kind=8)                ::  r(1), distance, ghost_circle_pos(3,2), x_pos, y_pos, x1, x2, y1, y2, model_x_len, model_y_len
            integer(kind=4)             ::  total_real_circle_num, total_ghost_circle_num, i, j, k, progress, ghost_circle_status(3)
            logical                     ::  is_overlapped, is_part_in_y_range, is_part_in_x_range

            model_x_len = container_size1
            model_y_len = container_size2

            call random_seed()

            call get_ref_circle_radius(ref_circle_radius,model_x_len*model_y_len)
            total_real_circle_num = size(ref_circle_radius)

            if (allocated(real_circle_info)) deallocate(real_circle_info)
            allocate(real_circle_info(total_real_circle_num,3))
            real_circle_info = 0.0
            do i = 1, total_real_circle_num
                real_circle_info(i,1) = ref_circle_radius(i)
            end do

            if (allocated(ghost_circle_info)) deallocate(ghost_circle_info)
            allocate(ghost_circle_info(total_real_circle_num*3,3))
            ghost_circle_info = 0.0
            
            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a17)') 'Particle shape done.','(',total_real_circle_num,'real particles).'

            ! generate the position of circular particles
            i = 1;
            total_ghost_circle_num = 0
            do while (i .le. total_real_circle_num)
                is_overlapped = .false.

                call random_number(r)
                real_circle_info(i,2) = r(1) * model_x_len
                real_circle_info(i,2) = real_circle_info(i,2) - 0.5 * model_x_len
                
                call random_number(r)
                real_circle_info(i,3) = r(1) * model_y_len
                real_circle_info(i,3) = real_circle_info(i,3) - 0.5 * model_y_len
                
                ! generate ghost particle based on the current real particle
                ! and detect whether the ghost particle is intersected with 
                ! the container
                ! for circle, three possible ghost particles exist
                if (real_circle_info(i,2) .gt. 0.0) then
                    x_pos = real_circle_info(i,2) - model_x_len
                else
                    x_pos = real_circle_info(i,2) + model_x_len
                end if
                if (real_circle_info(i,3) .gt. 0.0) then
                    y_pos = real_circle_info(i,3) - model_y_len
                else
                    y_pos = real_circle_info(i,3) + model_y_len
                end if

                ! all ghost particles have the same shape as the current real particle
                ! so only position is recorded
                ghost_circle_pos(1,1) = x_pos
                ghost_circle_pos(1,2) = real_circle_info(i,3)

                ghost_circle_pos(2,1) = x_pos
                ghost_circle_pos(2,2) = y_pos

                ghost_circle_pos(3,1) = real_circle_info(i,2)
                ghost_circle_pos(3,2) = y_pos

                ghost_circle_status = 0          ! 0 for completely outside the container, 1 for inside or intersected with the container
                do j = 1, 3
                    is_part_in_x_range = .false.
                    is_part_in_y_range = .false.

                    ! x = -0.5 * model_x_len
                    if ((-0.5 * model_x_len - ghost_circle_pos(j,1)) ** 2.0 .le. real_circle_info(i,1) ** 2.0) then
                        y1 = -sqrt(real_circle_info(i,1) ** 2.0 - (-0.5 * model_x_len - ghost_circle_pos(j,1)) ** 2.0) + ghost_circle_pos(j,2)
                        y2 =  sqrt(real_circle_info(i,1) ** 2.0 - (-0.5 * model_x_len - ghost_circle_pos(j,1)) ** 2.0) + ghost_circle_pos(j,2)
                        if ((y1 .gt. -0.5 * model_y_len) .and. (y1 .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                        if ((y2 .gt. -0.5 * model_y_len) .and. (y2 .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                    end if

                    ! x = 0.5 * model_x_len
                    if ((0.5 * model_x_len - ghost_circle_pos(j,1)) ** 2.0 .le. real_circle_info(i,1) ** 2.0) then
                        y1 = -sqrt(real_circle_info(i,1) ** 2.0 - (0.5 * model_x_len - ghost_circle_pos(j,1)) ** 2.0) + ghost_circle_pos(j,2)
                        y2 =  sqrt(real_circle_info(i,1) ** 2.0 - (0.5 * model_x_len - ghost_circle_pos(j,1)) ** 2.0) + ghost_circle_pos(j,2)
                        if ((y1 .gt. -0.5 * model_y_len) .and. (y1 .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                        if ((y2 .gt. -0.5 * model_y_len) .and. (y2 .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                    end if

                    ! y = =0.5 * model_y_len
                    if ((-0.5 * model_y_len - ghost_circle_pos(j,2)) ** 2.0 .le. real_circle_info(i,1) ** 2.0) then
                        x1 = -sqrt(real_circle_info(i,1) ** 2.0 - (-0.5 * model_y_len - ghost_circle_pos(j,2)) ** 2.0) + ghost_circle_pos(j,1)
                        x2 =  sqrt(real_circle_info(i,1) ** 2.0 - (-0.5 * model_y_len - ghost_circle_pos(j,2)) ** 2.0) + ghost_circle_pos(j,1)
                        if ((x1 .gt. -0.5 * model_x_len) .and. (x1 .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                        if ((x2 .gt. -0.5 * model_x_len) .and. (x2 .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                    end if

                    ! y = 0.5 * model_y_len
                    if ((0.5 * model_y_len - ghost_circle_pos(j,2)) ** 2.0 .le. real_circle_info(i,1) ** 2.0) then
                        x1 = -sqrt(real_circle_info(i,1) ** 2.0 - (0.5 * model_y_len - ghost_circle_pos(j,2)) ** 2.0) + ghost_circle_pos(j,1)
                        x2 =  sqrt(real_circle_info(i,1) ** 2.0 - (0.5 * model_y_len - ghost_circle_pos(j,2)) ** 2.0) + ghost_circle_pos(j,1)
                        if ((x1 .gt. -0.5 * model_x_len) .and. (x1 .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                        if ((x2 .gt. -0.5 * model_x_len) .and. (x2 .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                    end if

                    if ((is_part_in_x_range .eq. .true.) .or. (is_part_in_y_range .eq. .true.)) ghost_circle_status(j) = 1
                    
                end do
                
                ! overlap, new real particle vs old real particles
                do j = 1, i-1
                    distance = sqrt((real_circle_info(i,2) - real_circle_info(j,2)) ** 2.0 + (real_circle_info(i,3) - real_circle_info(j,3)) ** 2.0)
                    if (distance .lt. (real_circle_info(i,1) + real_circle_info(j,1))) then
                        is_overlapped = .true.
                        exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new real particle vs old ghost particles
                do j = 1, total_ghost_circle_num
                    distance = sqrt((real_circle_info(i,2) - ghost_circle_info(j,2)) ** 2.0 + (real_circle_info(i,3) - ghost_circle_info(j,3)) ** 2.0)
                    if (distance .lt. (real_circle_info(i,1) + ghost_circle_info(j,1))) then
                        is_overlapped = .true.
                        exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old real particles
                do j = 1, 3
                    if (ghost_circle_status(j) .eq. 1) then
                        do k = 1, i-1
                            distance = sqrt((ghost_circle_pos(j,1) - real_circle_info(k,2)) ** 2.0 + (ghost_circle_pos(j,2) - real_circle_info(k,3)) ** 2.0)
                            if (distance .lt. (real_circle_info(i,1) + real_circle_info(k,1))) then
                                is_overlapped = .true.
                                exit
                            end if
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old ghost particles
                do j = 1, 3
                    if (ghost_circle_status(j) .eq. 1) then
                        do k = 1, total_ghost_circle_num
                            distance = sqrt((ghost_circle_pos(j,1) - ghost_circle_info(k,2)) ** 2.0 + (ghost_circle_pos(j,2) - ghost_circle_info(k,3)) ** 2.0)
                            if (distance .lt. (real_circle_info(i,1) + ghost_circle_info(k,1))) then
                                is_overlapped = .true.
                                exit
                            end if
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    ! update the ghost_circle_info and total_ghost_circle_num
                    do j = 1, 3
                        if (ghost_circle_status(j) .eq. 1) then 
                            total_ghost_circle_num = total_ghost_circle_num + 1
                            ghost_circle_info(total_ghost_circle_num,1) = real_circle_info(i,1)
                            ghost_circle_info(total_ghost_circle_num,2:3) = ghost_circle_pos(j,1:2)
                        end if
                    end do
                    if (mod(i,nint(0.25 * total_real_circle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_real_circle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                    i = i + 1
                end if
            end do

            ! combine real_circle_info and ghost_circle_info into circular_particle_info
            if (allocated(circular_particle_info)) deallocate(circular_particle_info)
            allocate(circular_particle_info(total_real_circle_num + total_ghost_circle_num,3))
            do i = 1, total_real_circle_num
                circular_particle_info(i,1:3) = real_circle_info(i,1:3)
            end do
            do i = 1, total_ghost_circle_num
                circular_particle_info(total_real_circle_num+i,1:3) = ghost_circle_info(i,1:3)
            end do
            write(*,'(5X,a15,i6,a18,i6,a10,i6)') 'Real particles:',total_real_circle_num,'. Ghost particles:',total_ghost_circle_num,'. Totally:',total_real_circle_num+total_ghost_circle_num
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate shape of elliptical particles
        ! elliptical_particle_info    size:   [total_particle_num,5]
        ! each row contains:     radius_a   radius_b    angle      x_pos       y_pos
        ! the aspect_ratio1 is used here
        ! radius_a / radius_b = aspect_ratio1
        ! ONLY the shape of each ellipse is determined (x_pos and y_pos are undefined)
        subroutine ellipse_shape(elliptical_particle_info,container_area)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  elliptical_particle_info(:,:)
            real(kind=8),intent(in)     ::  container_area
            integer(kind=4)             ::  total_particle_num, i, j, k
            character(len=100)          ::  ellipse_major_radius_file
            real(kind=8),allocatable    ::  particle_needed_area(:), particle_real_area(:)
            real(kind=8)                ::  last_particle_area, r(1), particle_real_area_in_current_sieve, low_value, high_value, radius_a

            if (allocated(particle_needed_area)) deallocate(particle_needed_area)
            if (allocated(particle_real_area)) deallocate(particle_real_area)
            allocate(particle_needed_area(total_sieve_num-1))
            allocate(particle_real_area(total_sieve_num-1))
            
            do i = 1, total_sieve_num-1
                particle_needed_area(i) = (gradation(2,i+1) - gradation(2,i)) * particle_occupy_percentage * container_area
            end do
            
            ! randomize
            call random_seed()
 
            ellipse_major_radius_file='ellipse_minor_radius.tmp'   
            open(unit=100,file=ellipse_major_radius_file,form='formatted',status='replace',action='write')
            particle_real_area = 0.0
            total_particle_num = 0
            do i = total_sieve_num-1, 1, -1
                particle_real_area_in_current_sieve = 0.0 
                do while (particle_real_area(i) .lt. particle_needed_area(i))
                    ! generate the ellipse
                    call random_number(r)
                    low_value = gradation(1,i) / 2.0
                    high_value = gradation(1,i+1) / 2.0
                    radius_a = r(1) * (high_value - low_value) + low_value
                    particle_real_area_in_current_sieve = particle_real_area_in_current_sieve + (radius_a / aspect_ratio1) * radius_a * 3.14159
                    if (particle_real_area_in_current_sieve .gt. particle_needed_area(i)) then
                        last_particle_area = particle_needed_area(i) - particle_real_area(i)
                        radius_a = (last_particle_area / 3.14159 * aspect_ratio1) ** 0.5
                    end if
                    write(100,'(e15.5)') radius_a
                    total_particle_num = total_particle_num + 1
                    particle_real_area(i) = particle_real_area(i) + (radius_a / aspect_ratio1) * radius_a * 3.14159
                end do
            end do
            close(unit=100)

            if (allocated(elliptical_particle_info)) deallocate(elliptical_particle_info)
            allocate(elliptical_particle_info(total_particle_num,5))
            elliptical_particle_info = 0.0

            open(unit=100,file=ellipse_major_radius_file,form='formatted',status='old',action='read')
            do i = 1, total_particle_num
                read(100,*) elliptical_particle_info(i,1)
                elliptical_particle_info(i,2) = elliptical_particle_info(i,1) / aspect_ratio1
            end do
            close(unit=100,status='delete')

            k = 1
            call sort(elliptical_particle_info,k,'d')

            ! generate the rotation angle
            do i = 1, total_particle_num
                call random_number(r)
                low_value = 0.0
                high_value = 2.0 * 3.14159
                r(1) = r(1) * (high_value - low_value) + low_value
                elliptical_particle_info(i,3) = r(1)
            end do

        end subroutine

        ! generate elliptical particles
        ! elliptical_particle_info    size:   [total_particle_num,5]
        ! each row contains:     radius_a   radius_b    angle      x_pos       y_pos
        ! the aspect_ratio1 is used here
        ! radius_a / radius_b = aspect_ratio1
        subroutine ellipse_without_periodic_bc(elliptical_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  elliptical_particle_info(:,:)
            real(kind=8)                ::  r(1), rad_a, rad_b, angle, para_m, para_p, para_q, delta, model_x_len, model_y_len
            real(kind=8)                ::  ellipse_a(5), ellipse_b(5)
            integer(kind=4)             ::  flag, total_particle_num, i, j, progress

            model_x_len = container_size1
            model_y_len = container_size2

            call random_seed()

            call ellipse_shape(elliptical_particle_info,model_x_len*model_y_len)
            total_particle_num = size(elliptical_particle_info,1)

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            do i = 1, total_particle_num
                flag = 1
                do while (flag .eq. 1)
                    flag = 0
                    ! generate the random position of ellipse
                    call random_number(r)
                    r(1) = r(1) * model_x_len
                    elliptical_particle_info(i,4) = r(1) - 0.5 * model_x_len
                    call random_number(r)
                    r(1) = r(1) * model_y_len
                    elliptical_particle_info(i,5) = r(1) - 0.5 * model_y_len

                    ! check whether the ellipse is wholly inside the container
                    rad_a = elliptical_particle_info(i,1)
                    rad_b = elliptical_particle_info(i,2)
                    angle = elliptical_particle_info(i,3)
                    ! x = -0.5 * model_x_len
                    para_m = sin(2 * angle) * (rad_b ** 2.0 - rad_a ** 2.0) * (-0.5 * model_x_len - elliptical_particle_info(i,4))
                    para_p = (rad_a ** 2.0) * (cos(angle) ** 2.0) + (rad_b ** 2.0) * (sin(angle) ** 2.0)
                    para_q = (rad_b ** 2.0) * ((-0.5 * model_x_len - elliptical_particle_info(i,4)) ** 2.0) * (cos(angle) ** 2.0) + &
                           & (rad_a ** 2.0) * ((-0.5 * model_x_len - elliptical_particle_info(i,4)) ** 2.0) * (sin(angle) ** 2.0) - &
                           & (rad_a ** 2.0) * (rad_b ** 2.0)
                    delta = para_m ** 2.0 - 4 * para_p * para_q
                    if (delta .ge. 0.0) then
                        flag = 1
                        cycle
                    end if

                    ! x = 0.5 * model_x_len
                    para_m = sin(2 * angle) * (rad_b ** 2.0 - rad_a ** 2.0) * (0.5 * model_x_len - elliptical_particle_info(i,4))
                    para_p = (rad_a ** 2.0) * (cos(angle) ** 2.0) + (rad_b ** 2.0) * (sin(angle) ** 2.0)
                    para_q = (rad_b ** 2.0) * ((0.5 * model_x_len - elliptical_particle_info(i,4)) ** 2.0) * (cos(angle) ** 2.0) + &
                           & (rad_a ** 2.0) * ((0.5 * model_x_len - elliptical_particle_info(i,4)) ** 2.0) * (sin(angle) ** 2.0) - &
                           & (rad_a ** 2.0) * (rad_b ** 2.0)
                    delta = para_m ** 2.0 - 4 * para_p * para_q
                    if (delta .ge. 0.0) then
                        flag = 1
                        cycle
                    end if

                    ! y = -0.5 * model_y_len
                    para_m = 2 * (rad_b ** 2.0 - rad_a ** 2.0) * (-0.5 * model_y_len - elliptical_particle_info(i,5)) * sin(angle) * cos(angle)
                    para_p = (rad_b ** 2.0) * (cos(angle) ** 2.0) + (rad_a ** 2.0) * (sin(angle) ** 2.0)
                    para_q = (rad_b ** 2.0) * ((-0.5 * model_y_len - elliptical_particle_info(i,5)) ** 2.0) * (sin(angle) ** 2.0) + &
                           & (rad_a ** 2.0) * ((-0.5 * model_y_len - elliptical_particle_info(i,5)) ** 2.0) * (cos(angle) ** 2.0) - &
                           & (rad_a ** 2.0) * (rad_b ** 2.0)
                    delta = para_m ** 2.0 - 4 * para_p * para_q
                    if (delta .ge. 0.0) then
                        flag = 1
                        cycle
                    end if

                    ! y = 0.5 * model_y_len
                    para_m = 2 * (rad_b ** 2.0 - rad_a ** 2.0) * (0.5 * model_y_len - elliptical_particle_info(i,5)) * sin(angle) * cos(angle)
                    para_p = (rad_b ** 2.0) * (cos(angle) ** 2.0) + (rad_a ** 2.0) * (sin(angle) ** 2.0)
                    para_q = (rad_b ** 2.0) * ((0.5 * model_y_len - elliptical_particle_info(i,5)) ** 2.0) * (sin(angle) ** 2.0) + &
                           & (rad_a ** 2.0) * ((0.5 * model_y_len - elliptical_particle_info(i,5)) ** 2.0) * (cos(angle) ** 2.0) - &
                           & (rad_a ** 2.0) * (rad_b ** 2.0)
                    delta = para_m ** 2.0 - 4 * para_p * para_q
                    if (delta .ge. 0.0) then
                        flag = 1
                        cycle
                    end if

                    ! seperation check
                    ! prepare ellipse_a
                    ellipse_a(1:5) = elliptical_particle_info(i,1:5)

                    do j = 1, i-1
                        ! prepare ellipse_b
                        ellipse_b(1:5) = elliptical_particle_info(j,1:5)
                        if (is_ellipse_overlapped(ellipse_a,ellipse_b) .eq. .true.) then
                            flag = 1
                            exit
                        end if
                    end do
                end do

                if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                    progress = i / nint(0.25 * total_particle_num) * 25
                    write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate elliptical particles
        ! elliptical_particle_info    size:   [total_particle_num,5]
        ! each row contains:     radius_a   radius_b    angle      x_pos       y_pos
        ! the aspect_ratio1 is used here
        ! radius_a / radius_b = aspect_ratio1
        ! non periodic boundary condition
        ! for circular containder, the periodic
        ! boundary condition is not available
        subroutine ellipse_in_circle_container(elliptical_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  elliptical_particle_info(:,:)
            real(kind=8)                ::  r(1), angle, container_radius
            real(kind=8)                ::  ellipse_a(5), ellipse_b(5), x_pos, y_pos, point(3), trans_rotate(3,3)
            integer(kind=4)             ::  flag, total_particle_num, i, j, progress, angle_segment

            container_radius = container_size1

            call random_seed()

            call ellipse_shape(elliptical_particle_info,container_radius**2.0*3.14159)
            total_particle_num = size(elliptical_particle_info,1)

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            do i = 1, total_particle_num
                flag = 1
                do while (flag .eq. 1)
                    flag = 0
                    ! generate the random position of ellipse
                    call random_number(r)
                    elliptical_particle_info(i,4) = r(1) * 2 * container_radius - container_radius
                    call random_number(r)
                    elliptical_particle_info(i,5) = r(1) * 2 * container_radius - container_radius

                    ! check whether the ellipse is wholly inside the container
                    angle_segment = 30
                    do j = 1, angle_segment
                        angle = 2.0 * 3.14159 / (angle_segment -1) * (j - 1)
                        point(1) = elliptical_particle_info(i,1) * cos(angle)
                        point(2) = elliptical_particle_info(i,2) * sin(angle)
                        point(3) = 1

                        trans_rotate = 0.0
                        trans_rotate(1,1) = cos(elliptical_particle_info(i,3))
                        trans_rotate(1,2) = -sin(elliptical_particle_info(i,3))
                        trans_rotate(2,1) = +sin(elliptical_particle_info(i,3))
                        trans_rotate(2,2) = cos(elliptical_particle_info(i,3))
                        trans_rotate(3,3) = 1.0
                        point = matmul(trans_rotate,point)

                        x_pos = point(1) + elliptical_particle_info(i,4)
                        y_pos = point(2) + elliptical_particle_info(i,5)

                        if ((x_pos) ** 2.0 + (y_pos) ** 2.0 .gt. container_radius ** 2.0) then
                            flag = 1
                            exit
                        end if
                    end do
                    if (flag .eq. 1) cycle

                    ! seperation check
                    ellipse_a(1:5) = elliptical_particle_info(i,1:5)
                    do j = 1, i-1
                        ! prepare ellipse_b
                        ellipse_b(1:5) = elliptical_particle_info(j,1:5)
                        if (is_ellipse_overlapped(ellipse_a,ellipse_b) .eq. .true.) then
                            flag = 1
                            exit
                        end if
                    end do
                end do

                if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                    progress = i / nint(0.25 * total_particle_num) * 25
                    write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                end if

            end do

            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate elliptical particles
        ! elliptical_particle_info    size:   [total_particle_num,5]
        ! each row contains:     radius_a   radius_b    angle      x_pos       y_pos
        ! the aspect_ratio1 is used here
        ! radius_a / radius_b = aspect_ratio1
        subroutine ellipse_with_periodic_bc(elliptical_particle_info)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            real(kind=8),allocatable    ::  elliptical_particle_info(:,:)
            real(kind=8),allocatable    ::  real_ellipse_info(:,:), ghost_ellipse_info(:,:)
            real(kind=8),allocatable    ::  ref_circle_radius(:)
            real(kind=8)                ::  r(1), ghost_ellipse_pos(3,2), x_pos, y_pos, ellipse_a(5), ellipse_b(5), x1, x2, y1, y2, model_x_len, model_y_len
            real(kind=8)                ::  rad_a, rad_b, angle, para_q, para_p, para_m, delta
            integer(kind=4)             ::  total_real_ellipse_num, total_ghost_particle_num, total_particle_num, i, j, k, progress, ghost_ellipse_status(3)
            logical                     ::  is_overlapped, is_part_in_x_range, is_part_in_y_range

            model_x_len = container_size1
            model_y_len = container_size2

            call random_seed()

            call ellipse_shape(real_ellipse_info,model_x_len*model_y_len)
            total_real_ellipse_num = size(real_ellipse_info,1)

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a17)') 'Particle shape done.','(',total_real_ellipse_num,'real particles).'

            if (allocated(ghost_ellipse_info)) deallocate(ghost_ellipse_info)
            allocate(ghost_ellipse_info(total_real_ellipse_num*3,5))
            ghost_ellipse_info = 0.0

            ! generate the position of elliptical particles
            i = 1;
            total_ghost_particle_num = 0
            do while (i .le. total_real_ellipse_num)

                is_overlapped = .false.

                call random_number(r)
                real_ellipse_info(i,4) = r(1) * model_x_len - 0.5 * model_x_len
                
                call random_number(r)
                real_ellipse_info(i,5) = r(1) * model_y_len - 0.5 * model_y_len

                ! generate ghost particle based on the current real particle
                ! and detect whether the ghost particle is intersected with 
                ! the container
                ! for ellipse, three possible ghost particles exist
                if (real_ellipse_info(i,4) .gt. 0.0) then
                    x_pos = real_ellipse_info(i,4) - model_x_len
                else
                    x_pos = real_ellipse_info(i,4) + model_x_len
                end if
                if (real_ellipse_info(i,5) .gt. 0.0) then
                    y_pos = real_ellipse_info(i,5) - model_y_len
                else
                    y_pos = real_ellipse_info(i,5) + model_y_len
                end if

                ghost_ellipse_pos(1,1) = x_pos
                ghost_ellipse_pos(1,2) = real_ellipse_info(i,5)

                ghost_ellipse_pos(2,1) = x_pos
                ghost_ellipse_pos(2,2) = y_pos

                ghost_ellipse_pos(3,1) = real_ellipse_info(i,4)
                ghost_ellipse_pos(3,2) = y_pos

                ghost_ellipse_status = 0          ! 0 for completely outside the container, 1 for inside or intersected with the container
                ! check whether the ghost ellipse is intersected with the container and
                ! update the ghost_ellipse_status
                rad_a = real_ellipse_info(i,1)
                rad_b = real_ellipse_info(i,2)
                angle = real_ellipse_info(i,3)
                do j = 1, 3
                    is_part_in_x_range = .false.
                    is_part_in_y_range = .false.

                    ! x = -0.5 * model_x_len
                    para_m = sin(2 * angle) * (rad_b ** 2.0 - rad_a ** 2.0) * (-0.5 * model_x_len - ghost_ellipse_pos(j,1))
                    para_p = (rad_a ** 2.0) * (cos(angle) ** 2.0) + (rad_b ** 2.0) * (sin(angle) ** 2.0)
                    para_q = (rad_b ** 2.0) * ((-0.5 * model_x_len - ghost_ellipse_pos(j,1)) ** 2.0) * (cos(angle) ** 2.0) + &
                           & (rad_a ** 2.0) * ((-0.5 * model_x_len - ghost_ellipse_pos(j,1)) ** 2.0) * (sin(angle) ** 2.0) - &
                           & (rad_a ** 2.0) * (rad_b ** 2.0)
                    delta = para_m ** 2.0 - 4 * para_p * para_q
                    if (delta .ge. 0.0) then
                        y1 = (-para_m - sqrt(delta)) / 2.0 / para_p
                        y2 = (-para_m + sqrt(delta)) / 2.0 / para_p
                        if ((y1 + ghost_ellipse_pos(j,2) .gt. -0.5 * model_y_len) .and. (y1 + ghost_ellipse_pos(j,2) .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                        if ((y2 + ghost_ellipse_pos(j,2) .gt. -0.5 * model_y_len) .and. (y2 + ghost_ellipse_pos(j,2) .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                    end if

                    ! x = 0.5 * model_x_len
                    para_m = sin(2 * angle) * (rad_b ** 2.0 - rad_a ** 2.0) * (0.5 * model_x_len - ghost_ellipse_pos(j,1))
                    para_p = (rad_a ** 2.0) * (cos(angle) ** 2.0) + (rad_b ** 2.0) * (sin(angle) ** 2.0)
                    para_q = (rad_b ** 2.0) * ((0.5 * model_x_len - ghost_ellipse_pos(j,1)) ** 2.0) * (cos(angle) ** 2.0) + &
                           & (rad_a ** 2.0) * ((0.5 * model_x_len - ghost_ellipse_pos(j,1)) ** 2.0) * (sin(angle) ** 2.0) - &
                           & (rad_a ** 2.0) * (rad_b ** 2.0)
                    delta = para_m ** 2.0 - 4 * para_p * para_q
                    if (delta .ge. 0.0) then
                        y1 = (-para_m - sqrt(delta)) / 2.0 / para_p
                        y2 = (-para_m + sqrt(delta)) / 2.0 / para_p
                        if ((y1 + ghost_ellipse_pos(j,2) .gt. -0.5 * model_y_len) .and. (y1 + ghost_ellipse_pos(j,2) .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                        if ((y2 + ghost_ellipse_pos(j,2) .gt. -0.5 * model_y_len) .and. (y2 + ghost_ellipse_pos(j,2) .lt. 0.5 * model_y_len)) is_part_in_y_range = .true.
                    end if

                    ! y = 0
                    para_m = 2 * (rad_b ** 2.0 - rad_a ** 2.0) * (-0.5 * model_y_len - ghost_ellipse_pos(j,2)) * sin(angle) * cos(angle)
                    para_p = (rad_b ** 2.0) * (cos(angle) ** 2.0) + (rad_a ** 2.0) * (sin(angle) ** 2.0)
                    para_q = (rad_b ** 2.0) * ((-0.5 * model_y_len - ghost_ellipse_pos(j,2)) ** 2.0) * (sin(angle) ** 2.0) + &
                           & (rad_a ** 2.0) * ((-0.5 * model_y_len - ghost_ellipse_pos(j,2)) ** 2.0) * (cos(angle) ** 2.0) - &
                           & (rad_a ** 2.0) * (rad_b ** 2.0)
                    delta = para_m ** 2.0 - 4 * para_p * para_q
                    if (delta .ge. 0.0) then
                        x1 = (-para_m - sqrt(delta)) / 2.0 / para_p
                        x2 = (-para_m + sqrt(delta)) / 2.0 / para_p
                        if ((x1 + ghost_ellipse_pos(j,1) .gt. -0.5 * model_x_len) .and. (x1 + ghost_ellipse_pos(j,1) .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                        if ((x2 + ghost_ellipse_pos(j,1) .gt. -0.5 * model_x_len) .and. (x2 + ghost_ellipse_pos(j,1) .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                    end if

                    ! y = model_y_len
                    para_m = 2 * (rad_b ** 2.0 - rad_a ** 2.0) * (0.5 * model_y_len - ghost_ellipse_pos(j,2)) * sin(angle) * cos(angle)
                    para_p = (rad_b ** 2.0) * (cos(angle) ** 2.0) + (rad_a ** 2.0) * (sin(angle) ** 2.0)
                    para_q = (rad_b ** 2.0) * ((0.5 * model_y_len - ghost_ellipse_pos(j,2)) ** 2.0) * (sin(angle) ** 2.0) + &
                           & (rad_a ** 2.0) * ((0.5 * model_y_len - ghost_ellipse_pos(j,2)) ** 2.0) * (cos(angle) ** 2.0) - &
                           & (rad_a ** 2.0) * (rad_b ** 2.0)
                    delta = para_m ** 2.0 - 4 * para_p * para_q
                    if (delta .ge. 0.0) then
                        x1 = (-para_m - sqrt(delta)) / 2.0 / para_p
                        x2 = (-para_m + sqrt(delta)) / 2.0 / para_p
                        if ((x1 + ghost_ellipse_pos(j,1) .gt. -0.5 * model_x_len) .and. (x1 + ghost_ellipse_pos(j,1) .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                        if ((x2 + ghost_ellipse_pos(j,1) .gt. -0.5 * model_x_len) .and. (x2 + ghost_ellipse_pos(j,1) .lt. 0.5 * model_x_len)) is_part_in_x_range = .true.
                    end if

                    if ((is_part_in_x_range .eq. .true.) .or. (is_part_in_y_range .eq. .true.)) ghost_ellipse_status(j) = 1

                end do
                    
                ! overlap, new real particle vs old real particles
                ellipse_a(1:5) = real_ellipse_info(i,1:5)
                do j = 1, i-1
                    ellipse_b(1:5) = real_ellipse_info(j,1:5)
                    if (is_ellipse_overlapped(ellipse_a,ellipse_b) .eq. .true.) then
                        is_overlapped = .true.
                        exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new real particle vs old ghost particles
                ellipse_a(1:5) = real_ellipse_info(i,1:5)
                do j = 1, total_ghost_particle_num
                    ellipse_b(1:5) = ghost_ellipse_info(j,1:5)
                    if (is_ellipse_overlapped(ellipse_a,ellipse_b) .eq. .true.) then
                        is_overlapped = .true.
                        exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old real particles
                do j = 1, 3
                    if (ghost_ellipse_status(j) .eq. 1) then
                        ellipse_a(1:3) = real_ellipse_info(i,1:3)
                        ellipse_a(4:5) = ghost_ellipse_pos(j,1:2)
                        do k = 1, i-1
                            ellipse_b(1:5) = real_ellipse_info(k,1:5)
                            if (is_ellipse_overlapped(ellipse_a,ellipse_b) .eq. .true.) then
                                is_overlapped = .true.
                                exit
                            end if
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! overlap, new ghost particles vs old ghost particles
                do j = 1, 3
                    if (ghost_ellipse_status(j) .eq. 1) then
                        ellipse_a(1:3) = real_ellipse_info(i,1:3)
                        ellipse_a(4:5) = ghost_ellipse_pos(j,1:2)
                        do k = 1, total_ghost_particle_num
                            ellipse_b(1:5) = ghost_ellipse_info(k,1:5)
                            if (is_ellipse_overlapped(ellipse_a,ellipse_b) .eq. .true.) then
                                is_overlapped = .true.
                                exit
                            end if
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    ! update the ghost_ellipse_info and total_ghost_particle_num
                    do j = 1, 3
                        if (ghost_ellipse_status(j) .eq. 1) then 
                            total_ghost_particle_num = total_ghost_particle_num + 1
                            ghost_ellipse_info(total_ghost_particle_num,1:3) = real_ellipse_info(i,1:3)
                            ghost_ellipse_info(total_ghost_particle_num,4:5) = ghost_ellipse_pos(j,1:2)
                        end if
                    end do
                    if (mod(i,nint(0.25 * total_real_ellipse_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_real_ellipse_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                    i = i + 1
                end if
            end do

            ! combine real_ellipse_info and ghost_ellipse_info into elliptical_particle_info
            if (allocated(elliptical_particle_info)) deallocate(elliptical_particle_info)
            allocate(elliptical_particle_info(total_real_ellipse_num + total_ghost_particle_num,5))
            do i = 1, total_real_ellipse_num
                elliptical_particle_info(i,1:5) = real_ellipse_info(i,1:5)
            end do
            do i = 1, total_ghost_particle_num
                elliptical_particle_info(total_real_ellipse_num+i,1:5) = ghost_ellipse_info(i,1:5)
            end do
            write(*,'(5X,a15,i6,a18,i6,a10,i6)') 'Real particles:',total_real_ellipse_num,'. Ghost particles:',total_ghost_particle_num,'. Totally:',total_real_ellipse_num+total_ghost_particle_num
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

        ! generate polygonal particles
        ! polygon_side_num    size:   [total_particle_num]  contains the number of side for each polygon
        ! polygon_vertex_x_pos  size:   [total_particle_num,max_side_num]  contains the global coordinate for each vertex of each polygon
        ! the aspect_ratio1 is used here
        ! ONLY the shape of each polygon is determined
        subroutine polygon_shape(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos,container_area)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            integer(kind=4),allocatable     ::  polygon_side_num(:)
            real(kind=8),intent(in)         ::  container_area
            real(kind=8),allocatable        ::  polygon_vertex_x_pos(:,:), polygon_vertex_y_pos(:,:), ref_circle_radius(:), polygon_vertex(:,:)
            real(kind=8)                    ::  r(1), low_value, high_value, circle_area
            integer(kind=4)                 ::  i, j, side_num, total_particle_num

            call random_seed()

            call get_ref_circle_radius(ref_circle_radius,container_area)
            total_particle_num = size(ref_circle_radius)

            if (allocated(polygon_side_num)) deallocate(polygon_side_num)
            allocate(polygon_side_num(total_particle_num))
            polygon_side_num = 0
            if (allocated(polygon_vertex_x_pos)) deallocate(polygon_vertex_x_pos)
            allocate(polygon_vertex_x_pos(total_particle_num,max_side_num))
            polygon_vertex_x_pos = 0.0
            if (allocated(polygon_vertex_y_pos)) deallocate(polygon_vertex_y_pos)
            allocate(polygon_vertex_y_pos(total_particle_num,max_side_num))
            polygon_vertex_y_pos = 0.0

            ! cut procedure
            do i = 1, total_particle_num
                call random_number(r)
                low_value = min_side_num * 1.0
                high_value = max_side_num * 1.0
                r(1) = r(1) * (high_value - low_value) + low_value
                side_num = nint(r(1))
                polygon_side_num(i) = side_num
                do j = 1, side_num
                    call random_number(r)
                    low_value = 2 * 3.14159 / side_num * (j - 1)
                    high_value = 2 * 3.14159 / side_num * j
                    r(1) = r(1) * (high_value - low_value) + low_value
                    polygon_vertex_x_pos(i,j) = cos(r(1)) * ref_circle_radius(i)
                    polygon_vertex_y_pos(i,j) = sin(r(1)) * ref_circle_radius(i)
                end do
            end do   
            
            ! extend procedure
            ! Change the shape of polygons according to the given elongation
            do i = 1, total_particle_num
                side_num = polygon_side_num(i)
                if (allocated(polygon_vertex)) deallocate(polygon_vertex)
                allocate(polygon_vertex(2,side_num))
                do j = 1, side_num
                    polygon_vertex(1,j) = polygon_vertex_x_pos(i,j)
                    polygon_vertex(2,j) = polygon_vertex_y_pos(i,j)
                end do
                call extend_polygon(polygon_vertex,aspect_ratio1)      ! will change the value of polygon_vertex
                do j = 1, side_num
                    polygon_vertex_x_pos(i,j) = polygon_vertex(1,j)
                    polygon_vertex_y_pos(i,j) = polygon_vertex(2,j)
                end do
            end do
            
            ! scale procedure
            ! make sure the area of new polygon is the same as that of the refcircle
            ! see document for more detailed information
            do i = 1, total_particle_num
                side_num = polygon_side_num(i)
                if (allocated(polygon_vertex)) deallocate(polygon_vertex)
                allocate(polygon_vertex(2,side_num))
                do j = 1, side_num
                    polygon_vertex(1,j) = polygon_vertex_x_pos(i,j)
                    polygon_vertex(2,j) = polygon_vertex_y_pos(i,j)
                end do
                circle_area = 3.14159 * (ref_circle_radius(i)) ** 2.0
                call expand_polygon(polygon_vertex,circle_area)    ! will change the value of polygon_vertex
                do j = 1, side_num
                    polygon_vertex_x_pos(i,j) = polygon_vertex(1,j)
                    polygon_vertex_y_pos(i,j) = polygon_vertex(2,j)
                end do 
            end do    

            ! generate the position of each polygon

        end subroutine

        ! generate polygonal particles
        ! polygon_side_num    size:   [total_particle_num]  contains the number of side for each polygon
        ! polygon_vertex_x_pos  size:   [total_particle_num,max_side_num]  contains the global coordinate for each vertex of each polygon
        ! the aspect_ratio1 is used here
        subroutine polygon_without_periodic_bc(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            integer(kind=4),allocatable     ::  polygon_side_num(:)
            real(kind=8),allocatable        ::  polygon_vertex_x_pos(:,:), polygon_vertex_y_pos(:,:), polygon_a_vertex(:,:), polygon_b_vertex(:,:), ref_pos(:,:)
            real(kind=8)                    ::  r(1), low_value, high_value, x_pos, y_pos, model_x_len, model_y_len
            integer(kind=4)                 ::  i, j, k, side_num, total_particle_num, flag, progress

            model_x_len = container_size1
            model_y_len = container_size2

            call random_seed()

            call polygon_shape(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos,model_x_len*model_y_len)
            total_particle_num = size(polygon_side_num)

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            !random number generator
            call random_seed()
            
            if (allocated(ref_pos)) deallocate(ref_pos)
            allocate(ref_pos(2,total_particle_num))
            do i = 1, total_particle_num
                flag = 1
                do while (flag .eq. 1)
                    flag = 0
                    ! generate the random position of polygon
                    call random_number(r)
                    r(1) = r(1) * model_x_len
                    ref_pos(1,i) = r(1) - 0.5 * model_x_len
                    call random_number(r)
                    r(1) = r(1) * model_y_len
                    ref_pos(2,i) = r(1) - 0.5 * model_y_len

                    ! detect whether the polygon is outside or intersected by the container
                    side_num = polygon_side_num(i)
                    do j = 1, side_num
                        x_pos = polygon_vertex_x_pos(i,j) + ref_pos(1,i)
                        y_pos = polygon_vertex_y_pos(i,j) + ref_pos(2,i)
                        if ((x_pos .lt. -0.5 * model_x_len) .or. (x_pos .gt. 0.5 * model_x_len)) then
                            flag = 1
                            exit
                        end if
                        if ((y_pos .lt. -0.5 * model_y_len) .or. (y_pos .gt. 0.5 * model_y_len)) then
                            flag = 1
                            exit
                        end if
                    end do

                    if (flag .eq. 1) cycle

                    if (allocated(polygon_a_vertex)) deallocate(polygon_a_vertex)
                    allocate(polygon_a_vertex(2,polygon_side_num(i)))
                    do k = 1, polygon_side_num(i)
                        polygon_a_vertex(1,k) = polygon_vertex_x_pos(i,k) + ref_pos(1,i)
                        polygon_a_vertex(2,k) = polygon_vertex_y_pos(i,k) + ref_pos(2,i)
                    end do

                    ! detect whether the polygons are overlapped with each other
                    do j = 1, i - 1

                        if (allocated(polygon_b_vertex)) deallocate(polygon_b_vertex)
                        allocate(polygon_b_vertex(2,polygon_side_num(j)))
                        do k = 1, polygon_side_num(j)
                            polygon_b_vertex(1,k) = polygon_vertex_x_pos(j,k) + ref_pos(1,j)
                            polygon_b_vertex(2,k) = polygon_vertex_y_pos(j,k) + ref_pos(2,j)
                        end do

                        if (is_polygon_overlapped(polygon_a_vertex,polygon_b_vertex) .eq. .true.) then
                            flag = 1
                            exit
                        end if

                    end do
                end do  

                if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                    progress = i / nint(0.25 * total_particle_num) * 25
                    write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                end if

            end do
            
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

            do i = 1, total_particle_num
                side_num = polygon_side_num(i)
                do j = 1, side_num
                    polygon_vertex_x_pos(i,j) = polygon_vertex_x_pos(i,j) + ref_pos(1,i)
                    polygon_vertex_y_pos(i,j) = polygon_vertex_y_pos(i,j) + ref_pos(2,i)
                end do
            end do

        end subroutine

        ! generate polygonal particles
        ! polygon_side_num    size:   [total_particle_num]  contains the number of side for each polygon
        ! polygon_vertex_x_pos  size:   [total_particle_num,max_side_num]  contains the global coordinate for each vertex of each polygon
        ! the aspect_ratio1 is used here
        subroutine polygon_in_circle_container(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            integer(kind=4),allocatable     ::  polygon_side_num(:)
            real(kind=8),allocatable        ::  polygon_vertex_x_pos(:,:), polygon_vertex_y_pos(:,:), polygon_a_vertex(:,:), polygon_b_vertex(:,:), ref_pos(:,:)
            real(kind=8)                    ::  r(1), low_value, high_value, x_pos, y_pos, container_radius
            integer(kind=4)                 ::  i, j, k, side_num, total_particle_num, flag, progress

            container_radius = container_size1

            call random_seed()

            call polygon_shape(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos,container_radius**2.0*3.14159)
            total_particle_num = size(polygon_side_num)

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a12)') 'Particle shape done.','(',total_particle_num,'particles).'

            !random number generator
            call random_seed()
            
            if (allocated(ref_pos)) deallocate(ref_pos)
            allocate(ref_pos(2,total_particle_num))
            do i = 1, total_particle_num
                flag = 1
                do while (flag .eq. 1)
                    flag = 0
                    ! generate the random position of polygon
                    call random_number(r)
                    ref_pos(1,i) = r(1) * 2 * container_radius - container_radius
                    call random_number(r)
                    ref_pos(2,i) = r(1) * 2 * container_radius - container_radius

                    ! detect whether the polygon is outside or intersected by the container
                    side_num = polygon_side_num(i)
                    do j = 1, side_num
                        x_pos = polygon_vertex_x_pos(i,j) + ref_pos(1,i)
                        y_pos = polygon_vertex_y_pos(i,j) + ref_pos(2,i)
                        if (x_pos ** 2.0 + y_pos ** 2.0 .gt. container_radius ** 2.0) then
                            flag = 1
                            exit
                        end if
                    end do
                    if (flag .eq. 1) cycle

                    if (allocated(polygon_a_vertex)) deallocate(polygon_a_vertex)
                    allocate(polygon_a_vertex(2,polygon_side_num(i)))
                    do k = 1, polygon_side_num(i)
                        polygon_a_vertex(1,k) = polygon_vertex_x_pos(i,k) + ref_pos(1,i)
                        polygon_a_vertex(2,k) = polygon_vertex_y_pos(i,k) + ref_pos(2,i)
                    end do

                    ! detect whether the polygons are overlapped with each other
                    do j = 1, i - 1

                        if (allocated(polygon_b_vertex)) deallocate(polygon_b_vertex)
                        allocate(polygon_b_vertex(2,polygon_side_num(j)))
                        do k = 1, polygon_side_num(j)
                            polygon_b_vertex(1,k) = polygon_vertex_x_pos(j,k) + ref_pos(1,j)
                            polygon_b_vertex(2,k) = polygon_vertex_y_pos(j,k) + ref_pos(2,j)
                        end do

                        if (is_polygon_overlapped(polygon_a_vertex,polygon_b_vertex) .eq. .true.) then
                            flag = 1
                            exit
                        end if

                    end do
                end do  

                if (mod(i,nint(0.25 * total_particle_num)) .eq. 0) then
                    progress = i / nint(0.25 * total_particle_num) * 25
                    write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                end if

            end do
            
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

            do i = 1, total_particle_num
                side_num = polygon_side_num(i)
                do j = 1, side_num
                    polygon_vertex_x_pos(i,j) = polygon_vertex_x_pos(i,j) + ref_pos(1,i)
                    polygon_vertex_y_pos(i,j) = polygon_vertex_y_pos(i,j) + ref_pos(2,i)
                end do
            end do

        end subroutine

        ! generate polygonal particles
        ! polygon_side_num    size:   [total_particle_num]  contains the number of side for each polygon
        ! polygon_vertex_x_pos  size:   [total_particle_num,max_side_num]  contains the global coordinate for each vertex of each polygon
        ! the aspect_ratio1 is used here
        subroutine polygon_with_periodic_bc(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos)
            use ParticleSimu_global_parameter
            use ParticleSimu_general_sub
            implicit none
            integer(kind=4),allocatable     ::  polygon_side_num(:)
            real(kind=8),allocatable        ::  polygon_vertex_x_pos(:,:), polygon_vertex_y_pos(:,:)
            real(kind=8),allocatable        ::  real_polygon_original_vertex_x_pos(:,:), real_polygon_original_vertex_y_pos(:,:)
            real(kind=8),allocatable        ::  ghost_polygon_original_vertex_x_pos(:,:), ghost_polygon_original_vertex_y_pos(:,:)
            integer(kind=4),allocatable     ::  real_polygon_side_num(:), ghost_polygon_side_num(:)
            real(kind=8),allocatable        ::  polygon_a_vertex(:,:), polygon_b_vertex(:,:), real_polygon_center_pos(:,:), ghost_polygon_center_pos(:,:)
            real(kind=8)                    ::  r(1), x_pos, y_pos, local_ghost_polygon_center_pos(3,2), model_x_len, model_y_len
            integer(kind=4)                 ::  i, j, k, ii, total_real_particle_num, total_ghost_particle_num, progress, local_ghost_polygon_status(3)
            logical                         ::  is_overlapped, is_in_x_range, is_in_y_range

            model_x_len = container_size1
            model_y_len = container_size2

            call random_seed()

            call polygon_shape(real_polygon_side_num,real_polygon_original_vertex_x_pos,real_polygon_original_vertex_y_pos,model_x_len*model_y_len)
            total_real_particle_num = size(real_polygon_side_num)

            if (allocated(real_polygon_center_pos)) deallocate(real_polygon_center_pos)
            allocate(real_polygon_center_pos(total_real_particle_num,2))
            real_polygon_center_pos = 0.0

            if (allocated(ghost_polygon_original_vertex_x_pos)) deallocate(ghost_polygon_original_vertex_x_pos)
            allocate(ghost_polygon_original_vertex_x_pos(total_real_particle_num*3,max_side_num))
            ghost_polygon_original_vertex_x_pos = 0.0

            if (allocated(ghost_polygon_original_vertex_y_pos)) deallocate(ghost_polygon_original_vertex_y_pos)
            allocate(ghost_polygon_original_vertex_y_pos(total_real_particle_num*3,max_side_num))
            ghost_polygon_original_vertex_y_pos = 0.0

            if (allocated(ghost_polygon_side_num)) deallocate(ghost_polygon_side_num)
            allocate(ghost_polygon_side_num(total_real_particle_num*3))
            ghost_polygon_side_num = 0

            write(*,*) ''
            write(*,'(5X,a20,a2,i7,a17)') 'Particle shape done.','(',total_real_particle_num,'real particles).'

            !random number generator
            call random_seed()
                   
            total_ghost_particle_num = 0
            if (allocated(real_polygon_center_pos)) deallocate(real_polygon_center_pos)
            allocate(real_polygon_center_pos(total_real_particle_num,2))
            real_polygon_center_pos = 0.0
            if (allocated(ghost_polygon_center_pos)) deallocate(ghost_polygon_center_pos)
            allocate(ghost_polygon_center_pos(total_real_particle_num*3,2))
            ghost_polygon_center_pos = 0.0

            i = 1;
            do while (i .le. total_real_particle_num)
                is_overlapped = .false.
                ! generate the random position of polygon
                call random_number(r)
                r(1) = r(1) * model_x_len
                real_polygon_center_pos(i,1) = r(1) - 0.5 * model_x_len
                call random_number(r)
                r(1) = r(1) * model_y_len
                real_polygon_center_pos(i,2) = r(1) - 0.5 * model_y_len

                ! generate ghost particle based on the current real particle
                ! and detect whether the ghost particle is intersected with 
                ! the container
                ! for polygon, three possible ghost particles exist
                if (real_polygon_center_pos(i,1) .gt. 0.0) then
                    x_pos = real_polygon_center_pos(i,1) - model_x_len
                else
                    x_pos = real_polygon_center_pos(i,1) + model_x_len
                end if
                if (real_polygon_center_pos(i,2) .gt. 0.0) then
                    y_pos = real_polygon_center_pos(i,2) - model_y_len
                else
                    y_pos = real_polygon_center_pos(i,2) + model_y_len
                end if

                local_ghost_polygon_center_pos(1,1) = x_pos
                local_ghost_polygon_center_pos(1,2) = real_polygon_center_pos(i,2)

                local_ghost_polygon_center_pos(2,1) = x_pos
                local_ghost_polygon_center_pos(2,2) = y_pos

                local_ghost_polygon_center_pos(3,1) = real_polygon_center_pos(i,1)
                local_ghost_polygon_center_pos(3,2) = y_pos

                ! check whether the polygon is inside or intersected with the container
                local_ghost_polygon_status = 0
                do j = 1, 3
                    is_in_x_range = .false.
                    is_in_y_range = .false.
                    do k = 1, real_polygon_side_num(i)
                        x_pos = real_polygon_original_vertex_x_pos(i,k) + local_ghost_polygon_center_pos(j,1)
                        y_pos = real_polygon_original_vertex_y_pos(i,k) + local_ghost_polygon_center_pos(j,2)
                        if ((x_pos .gt. -0.5 * model_x_len) .and. (x_pos .lt. 0.5 * model_x_len)) then
                            is_in_x_range = .true.
                        end if
                        if ((y_pos .gt. -0.5 * model_y_len) .and. (y_pos .lt. 0.5 * model_y_len)) then
                            is_in_y_range = .true.
                        end if
                    end do
                    if ((is_in_x_range .eq. .true.) .and. (is_in_y_range .eq. .true.)) local_ghost_polygon_status(j) = 1
                end do

                ! seperation check
                ! new real polygon vs old real polygons
                if (allocated(polygon_a_vertex)) deallocate(polygon_a_vertex)
                allocate(polygon_a_vertex(2,real_polygon_side_num(i)))
                polygon_a_vertex(1,1:real_polygon_side_num(i)) = real_polygon_original_vertex_x_pos(i,1:real_polygon_side_num(i)) + real_polygon_center_pos(i,1)
                polygon_a_vertex(2,1:real_polygon_side_num(i)) = real_polygon_original_vertex_y_pos(i,1:real_polygon_side_num(i)) + real_polygon_center_pos(i,2)

                do j = 1, i - 1
                    if (allocated(polygon_b_vertex)) deallocate(polygon_b_vertex)
                    allocate(polygon_b_vertex(2,real_polygon_side_num(j)))
                    polygon_b_vertex(1,1:real_polygon_side_num(j)) = real_polygon_original_vertex_x_pos(j,1:real_polygon_side_num(j)) + real_polygon_center_pos(j,1)
                    polygon_b_vertex(2,1:real_polygon_side_num(j)) = real_polygon_original_vertex_y_pos(j,1:real_polygon_side_num(j)) + real_polygon_center_pos(j,2)
                    if (is_polygon_overlapped(polygon_a_vertex,polygon_b_vertex) .eq. .true.) then
                        is_overlapped = .true.
                        exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! new real polygon vs old ghost polygons
                if (allocated(polygon_a_vertex)) deallocate(polygon_a_vertex)
                allocate(polygon_a_vertex(2,real_polygon_side_num(i)))
                polygon_a_vertex(1,1:real_polygon_side_num(i)) = real_polygon_original_vertex_x_pos(i,1:real_polygon_side_num(i)) + real_polygon_center_pos(i,1)
                polygon_a_vertex(2,1:real_polygon_side_num(i)) = real_polygon_original_vertex_y_pos(i,1:real_polygon_side_num(i)) + real_polygon_center_pos(i,2)

                do j = 1, total_ghost_particle_num
                    if (allocated(polygon_b_vertex)) deallocate(polygon_b_vertex)
                    allocate(polygon_b_vertex(2,ghost_polygon_side_num(j)))
                    polygon_b_vertex(1,1:ghost_polygon_side_num(j)) = ghost_polygon_original_vertex_x_pos(j,1:ghost_polygon_side_num(j)) + ghost_polygon_center_pos(j,1)
                    polygon_b_vertex(2,1:ghost_polygon_side_num(j)) = ghost_polygon_original_vertex_y_pos(j,1:ghost_polygon_side_num(j)) + ghost_polygon_center_pos(j,2)
                    if (is_polygon_overlapped(polygon_a_vertex,polygon_b_vertex) .eq. .true.) then
                        is_overlapped = .true.
                        exit
                    end if
                end do
                if (is_overlapped .eq. .true.) cycle

                ! new ghost polygon vs old real polygons
                if (allocated(polygon_a_vertex)) deallocate(polygon_a_vertex)
                allocate(polygon_a_vertex(2,real_polygon_side_num(i)))
                do j = 1, 3
                    if (local_ghost_polygon_status(j) .eq. 1) then
                        polygon_a_vertex(1,1:real_polygon_side_num(i)) = real_polygon_original_vertex_x_pos(i,1:real_polygon_side_num(i)) + local_ghost_polygon_center_pos(j,1)
                        polygon_a_vertex(2,1:real_polygon_side_num(i)) = real_polygon_original_vertex_y_pos(i,1:real_polygon_side_num(i)) + local_ghost_polygon_center_pos(j,2)
                        do k = 1, total_real_particle_num
                            if (allocated(polygon_b_vertex)) deallocate(polygon_b_vertex)
                            allocate(polygon_b_vertex(2,real_polygon_side_num(k)))
                            polygon_b_vertex(1,1:real_polygon_side_num(k)) = real_polygon_original_vertex_x_pos(k,1:real_polygon_side_num(k)) + real_polygon_center_pos(k,1)
                            polygon_b_vertex(2,1:real_polygon_side_num(k)) = real_polygon_original_vertex_y_pos(k,1:real_polygon_side_num(k)) + real_polygon_center_pos(k,2)
                            if (is_polygon_overlapped(polygon_a_vertex,polygon_b_vertex) .eq. .true.) then
                                is_overlapped = .true.
                                exit
                            end if
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do

                ! new ghost polygon vs old ghost polygons
                if (allocated(polygon_a_vertex)) deallocate(polygon_a_vertex)
                allocate(polygon_a_vertex(2,real_polygon_side_num(i)))
                do j = 1, 3
                    if (local_ghost_polygon_status(j) .eq. 1) then
                        polygon_a_vertex(1,1:real_polygon_side_num(i)) = real_polygon_original_vertex_x_pos(i,1:real_polygon_side_num(i)) + local_ghost_polygon_center_pos(j,1)
                        polygon_a_vertex(2,1:real_polygon_side_num(i)) = real_polygon_original_vertex_y_pos(i,1:real_polygon_side_num(i)) + local_ghost_polygon_center_pos(j,2)
                        do k = 1, total_ghost_particle_num
                            if (allocated(polygon_b_vertex)) deallocate(polygon_b_vertex)
                            allocate(polygon_b_vertex(2,ghost_polygon_side_num(k)))
                            polygon_b_vertex(1,1:ghost_polygon_side_num(k)) = ghost_polygon_original_vertex_x_pos(k,1:ghost_polygon_side_num(k)) + ghost_polygon_center_pos(k,1)
                            polygon_b_vertex(2,1:ghost_polygon_side_num(k)) = ghost_polygon_original_vertex_y_pos(k,1:ghost_polygon_side_num(k)) + ghost_polygon_center_pos(k,2)
                            if (is_polygon_overlapped(polygon_a_vertex,polygon_b_vertex) .eq. .true.) then
                                is_overlapped = .true.
                                exit
                            end if
                        end do
                        if (is_overlapped .eq. .true.) exit
                    end if
                end do

                if (is_overlapped .eq. .false.) then
                    ! update ghost_polygon_vertex_x_pos, ghost_polygon_vertex_y_pos, ghost_polygon_side_num
                    do j = 1, 3
                        if (local_ghost_polygon_status(j) .eq. 1) then 
                            total_ghost_particle_num = total_ghost_particle_num + 1
                            ghost_polygon_side_num(total_ghost_particle_num) = real_polygon_side_num(i)
                            ghost_polygon_original_vertex_x_pos(total_ghost_particle_num,:) = real_polygon_original_vertex_x_pos(i,:)
                            ghost_polygon_original_vertex_y_pos(total_ghost_particle_num,:) = real_polygon_original_vertex_y_pos(i,:)
                            ghost_polygon_center_pos(total_ghost_particle_num,1:2) = local_ghost_polygon_center_pos(j,1:2)
                        end if
                    end do
                    if (mod(i,nint(0.25 * total_real_particle_num)) .eq. 0) then
                        progress = i / nint(0.25 * total_real_particle_num) * 25
                        write(*,'(5X,a32,i3,a1)') 'Generate particle position....  ',progress,'%'
                    end if
                    i = i + 1;
                end if
            end do

            ! combine real particle info and ghost particle info
            if (allocated(polygon_side_num)) deallocate(polygon_side_num)
            allocate(polygon_side_num(total_real_particle_num + total_ghost_particle_num))
            polygon_side_num = 0

            if (allocated(polygon_vertex_x_pos)) deallocate(polygon_vertex_x_pos)
            allocate(polygon_vertex_x_pos(total_real_particle_num + total_ghost_particle_num,max_side_num))
            polygon_vertex_x_pos = 0.0

            if (allocated(polygon_vertex_y_pos)) deallocate(polygon_vertex_y_pos)
            allocate(polygon_vertex_y_pos(total_real_particle_num + total_ghost_particle_num,max_side_num))
            polygon_vertex_y_pos = 0.0

            do i = 1, total_real_particle_num
                polygon_side_num(i) = real_polygon_side_num(i)
                polygon_vertex_x_pos(i,:) = real_polygon_original_vertex_x_pos(i,:) + real_polygon_center_pos(i,1)
                polygon_vertex_y_pos(i,:) = real_polygon_original_vertex_y_pos(i,:) + real_polygon_center_pos(i,2)
            end do
            do i = 1, total_ghost_particle_num
                polygon_side_num(total_real_particle_num+i) = ghost_polygon_side_num(i)
                polygon_vertex_x_pos(total_real_particle_num+i,:) = ghost_polygon_original_vertex_x_pos(i,:) + ghost_polygon_center_pos(i,1)
                polygon_vertex_y_pos(total_real_particle_num+i,:) = ghost_polygon_original_vertex_y_pos(i,:) + ghost_polygon_center_pos(i,2)
            end do
            write(*,'(5X,a15,i6,a18,i6,a10,i6)') 'Real particles:',total_real_particle_num,'. Ghost particles:',total_ghost_particle_num,'. Totally:',total_real_particle_num+total_ghost_particle_num
            write(*,'(5X,a52)') 'Job is done! Check the output file(s) in the folder.'

        end subroutine

end module

!dec$ endif