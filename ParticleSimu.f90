program main
	use ParticleSimu_global_parameter
	use ParticleSimu_initialization
	use ParticleSimu_2D_kernel
	use ParticleSimu_3D_kernel
	implicit none
	character(len=100)			::	filename
	real(kind=8),allocatable	::	particle_info(:,:), polygon_vertex_x_pos(:,:), polygon_vertex_y_pos(:,:)
	integer(kind=4),allocatable	::	polygon_side_num(:)
	integer(kind=4)				::	total_particle_num, i, j, simu_num
	character(len=1)			::	sign

	filename = 'input.dat'
	call read_input_file(filename)

	! circle, non periodic bc, rectangle container
	if ((particle_type .eq. 1) .and. (container_type .eq. 1) .and. (periodic_bc_on .eq. 0)) then
		call circle_without_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'circle.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(3e20.5)') ((particle_info(i,j),j=1,3),i=1,total_particle_num)
		close(unit=100)
	end if

	! circle, periodic bc, rectangle container
	if ((particle_type .eq. 1) .and. (container_type .eq. 1) .and. (periodic_bc_on .eq. 1)) then
		call circle_with_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'circle.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(3e20.5)') ((particle_info(i,j),j=1,3),i=1,total_particle_num)
		close(unit=100)
	end if

	! circle, non periodic bc, circle container
	if ((particle_type .eq. 1) .and. (container_type .eq. 3)) then
		call circle_in_circle_container(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'circle.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(3e20.5)') ((particle_info(i,j),j=1,3),i=1,total_particle_num)
		close(unit=100)
	end if
	! ///////////////////////////////////////////////////////////////////////////////////////////

	! ellipse, non periodic bc, rectangle container
	if ((particle_type .eq. 2) .and. (container_type .eq. 1) .and. (periodic_bc_on .eq. 0)) then
		call ellipse_without_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'ellipse.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(5e20.5)') ((particle_info(i,j),j=1,5),i=1,total_particle_num)
		close(unit=100)
	end if

	! ellipse, periodic bc, rectangle container
	if ((particle_type .eq. 2) .and. (container_type .eq. 1) .and. (periodic_bc_on .eq. 1)) then
		call ellipse_with_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'ellipse.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(5e20.5)') ((particle_info(i,j),j=1,5),i=1,total_particle_num)
		close(unit=100)
    end if

    ! ellipse, non periodic bc, circle container
	if ((particle_type .eq. 2) .and. (container_type .eq. 3)) then
		call ellipse_in_circle_container(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'ellipse.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(5e20.5)') ((particle_info(i,j),j=1,5),i=1,total_particle_num)
		close(unit=100)
    end if
    ! ///////////////////////////////////////////////////////////////////////////////////////////

    ! polygon, non periodic bc, , rectangle container
	if ((particle_type .eq. 3) .and. (container_type .eq. 1) .and. (periodic_bc_on .eq. 0)) then
		call polygon_without_periodic_bc(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos)
		total_particle_num = size(polygon_side_num)
		filename = 'side_num.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(i10)') (polygon_side_num(i),i=1,total_particle_num)
		close(unit=100)

		filename = 'vertex_x_pos.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		do i = 1, total_particle_num
			do j = 1, max_side_num-1
				write(100,'(e20.5)',advance='no') polygon_vertex_x_pos(i,j)
			end do
			write(100,'(e20.5)',advance='yes') polygon_vertex_x_pos(i,max_side_num)
		end do
		close(unit=100)

		filename = 'vertex_y_pos.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		do i = 1, total_particle_num
			do j = 1, max_side_num-1
				write(100,'(e20.5)',advance='no') polygon_vertex_y_pos(i,j)
			end do
			write(100,'(e20.5)',advance='yes') polygon_vertex_y_pos(i,max_side_num)
		end do
		close(unit=100)
    end if
    
    ! polygon, periodic bc, , rectangle container
	if ((particle_type .eq. 3) .and. (container_type .eq. 1) .and. (periodic_bc_on .eq. 1)) then
		call polygon_with_periodic_bc(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos)
		total_particle_num = size(polygon_side_num)
		filename = 'side_num.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(i10)') (polygon_side_num(i),i=1,total_particle_num)
		close(unit=100)

		filename = 'vertex_x_pos.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		do i = 1, total_particle_num
			do j = 1, max_side_num-1
				write(100,'(e20.5)',advance='no') polygon_vertex_x_pos(i,j)
			end do
			write(100,'(e20.5)',advance='yes') polygon_vertex_x_pos(i,max_side_num)
		end do
		close(unit=100)

		filename = 'vertex_y_pos.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		do i = 1, total_particle_num
			do j = 1, max_side_num-1
				write(100,'(e20.5)',advance='no') polygon_vertex_y_pos(i,j)
			end do
			write(100,'(e20.5)',advance='yes') polygon_vertex_y_pos(i,max_side_num)
		end do
		close(unit=100)
	end if

	! polygon, non periodic bc, circle container
	if ((particle_type .eq. 3) .and. (container_type .eq. 3)) then
		call polygon_in_circle_container(polygon_side_num,polygon_vertex_x_pos,polygon_vertex_y_pos)
		total_particle_num = size(polygon_side_num)
		filename = 'side_num.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(i10)') (polygon_side_num(i),i=1,total_particle_num)
		close(unit=100)

		filename = 'vertex_x_pos.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		do i = 1, total_particle_num
			do j = 1, max_side_num-1
				write(100,'(e20.5)',advance='no') polygon_vertex_x_pos(i,j)
			end do
			write(100,'(e20.5)',advance='yes') polygon_vertex_x_pos(i,max_side_num)
		end do
		close(unit=100)

		filename = 'vertex_y_pos.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		do i = 1, total_particle_num
			do j = 1, max_side_num-1
				write(100,'(e20.5)',advance='no') polygon_vertex_y_pos(i,j)
			end do
			write(100,'(e20.5)',advance='yes') polygon_vertex_y_pos(i,max_side_num)
		end do
		close(unit=100)
    end if
    ! ///////////////////////////////////////////////////////////////////////////////////////////

	! sphere, non periodic bc, box container
	if ((particle_type .eq. 4) .and. (container_type .eq. 2) .and. (periodic_bc_on .eq. 0)) then
		call sphere_without_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'sphere.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(4e20.5)') ((particle_info(i,j),j=1,4),i=1,total_particle_num)
		close(unit=100)
	end if

	! sphere, periodic bc, box container
	if ((particle_type .eq. 4) .and. (container_type .eq. 2) .and. (periodic_bc_on .eq. 1)) then
		call sphere_with_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'sphere.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(4e20.5)') ((particle_info(i,j),j=1,4),i=1,total_particle_num)
		close(unit=100)
    end if

    ! sphere, cylinder container
	if ((particle_type .eq. 4) .and. (container_type .eq. 4)) then
		call sphere_in_cylinder_container(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'sphere.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(4e20.5)') ((particle_info(i,j),j=1,4),i=1,total_particle_num)
		close(unit=100)
	end if
    ! ///////////////////////////////////////////////////////////////////////////////////////////

    ! ellipsoid, non periodic bc, box container
	if ((particle_type .eq. 5) .and. (container_type .eq. 2) .and. (periodic_bc_on .eq. 0)) then
		call ellipsoid_without_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'ellipsoid.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(9e20.5)') ((particle_info(i,j),j=1,9),i=1,total_particle_num)
		close(unit=100)
    end if
    
    ! ellipsoid, periodic bc, box container
	if ((particle_type .eq. 5) .and. (container_type .eq. 2) .and. (periodic_bc_on .eq. 1)) then
		call ellipsoid_with_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'ellipsoid.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(9e20.5)') ((particle_info(i,j),j=1,9),i=1,total_particle_num)
		close(unit=100)
	end if

	! ellipsoid, non periodic bc, cylinder container
	if ((particle_type .eq. 5) .and. (container_type .eq. 4)) then
		call ellipsoid_in_cylinder_container(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'ellipsoid.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(9e20.5)') ((particle_info(i,j),j=1,9),i=1,total_particle_num)
		close(unit=100)
    end if
	! ///////////////////////////////////////////////////////////////////////////////////////////

	! cylinder, non periodic bc, box container
	if ((particle_type .eq. 6) .and. (container_type .eq. 2) .and. (periodic_bc_on .eq. 0)) then
		call cylinder_without_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'cylinder.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(9e20.5)') ((particle_info(i,j),j=1,9),i=1,total_particle_num)
		close(unit=100)
	end if

	! cylinder, periodic bc, box container
	if ((particle_type .eq. 6) .and. (container_type .eq. 2) .and. (periodic_bc_on .eq. 1)) then
		call cylinder_with_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'cylinder.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(9e20.5)') ((particle_info(i,j),j=1,9),i=1,total_particle_num)
		close(unit=100)
    end if

    ! cylinder, non periodic bc, cylinder container
	if ((particle_type .eq. 6) .and. (container_type .eq. 4)) then
		call cylinder_in_cylinder_container(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'cylinder.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(9e20.5)') ((particle_info(i,j),j=1,9),i=1,total_particle_num)
		close(unit=100)
	end if
    ! ///////////////////////////////////////////////////////////////////////////////////////////

    ! capsule, non periodic bc, box container
	if ((particle_type .eq. 7) .and. (container_type .eq. 2) .and. (periodic_bc_on .eq. 0)) then
		call capsule_without_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'capsule.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(8e20.5)') ((particle_info(i,j),j=1,8),i=1,total_particle_num)
		close(unit=100)
    end if
    
    ! capsule, periodic bc, box container
	if ((particle_type .eq. 7) .and. (container_type .eq. 2) .and. (periodic_bc_on .eq. 1)) then
		call capsule_with_periodic_bc(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'capsule.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(8e20.5)') ((particle_info(i,j),j=1,8),i=1,total_particle_num)
		close(unit=100)
	end if

	! capsule, non periodic bc, cylinder container
	if ((particle_type .eq. 7) .and. (container_type .eq. 4)) then
		call capsule_in_cylinder_container(particle_info)
		total_particle_num = size(particle_info,1)
		filename = 'capsule.dat'
		open(unit=100,file=filename,form='formatted',status='replace',action='write')
		write(100,'(8e20.5)') ((particle_info(i,j),j=1,8),i=1,total_particle_num)
		close(unit=100)
    end if
	! ///////////////////////////////////////////////////////////////////////////////////////////

	write(*,'(5X,a22)') 'Press any key to quit.'
	pause

end program