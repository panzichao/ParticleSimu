!dec$ if .not. defined(__ParticleSimu_initialization_f90)
!dec$ define __ParticleSimu_initialization_f90

module ParticleSimu_initialization
    implicit none
    contains
    
        !****************************************************************************
        !  Initiate global parameters by reading input file
        !****************************************************************************
        subroutine read_input_file(filename)
            use ParticleSimu_global_parameter
            implicit none
            character(len=*),intent(in) ::  filename
            integer(kind=4)             ::  i
            
            open(unit=100,file=filename,form='formatted',status='old',action='read')
            read(100,*) spatial_dimension

            ! container_type = 1     for rectangle
            !                  2     for box
            !                  3     for circle
            !                  4     for cylinder
            ! when container_type = 1
            ! container_size1 = x_len, container_size2 = y_len
            ! when container_type = 2
            ! container_size1 = x_len, container_size2 = y_len, container_size3 = z_len
            ! when container_type = 3
            ! container_size1 = radius
            ! when container_type = 4, container_size2 = height
            ! container_size1 = radius
            ! CAUTION: when container_type = 3 or 4, the periodic boundary condition is not available
            read(100,*) container_type
            read(100,*) container_size1, container_size2, container_size3
            read(100,*) particle_occupy_percentage
            read(100,*) total_sieve_num
            if (allocated(gradation)) deallocate(gradation)
            allocate(gradation(2,total_sieve_num))
            gradation = 0.0
            read(100,*) (gradation(1,i),i=1,total_sieve_num)
            read(100,*) (gradation(2,i),i=1,total_sieve_num)
            read(100,*) particle_type
            read(100,*) aspect_ratio1, aspect_ratio2
            read(100,*) min_side_num, max_side_num
            read(100,*) periodic_bc_on
            close(unit=100)
        end subroutine
        
end module

!dec$ endif