!dec$ if .not. defined(__ParticleSimu_global_parameter_f90)
!dec$ define __ParticleSimu_global_parameter_f90

module ParticleSimu_global_parameter
    implicit none
    integer(kind=4),save            ::  spatial_dimension
    integer(kind=4),save            ::  container_type
    real(kind=8),save               ::  container_size1
    real(kind=8),save               ::  container_size2
    real(kind=8),save               ::  container_size3
    real(kind=8),save               ::  particle_occupy_percentage  !   area ratio and volumn ratio
    integer(kind=4),save            ::  total_sieve_num
    real(kind=8),save,allocatable   ::  gradation(:,:)
    integer(kind=4),save            ::  particle_type
    real(kind=8),save               ::  aspect_ratio1
    real(kind=8),save               ::  aspect_ratio2
    integer(kind=4),save            ::  min_side_num
    integer(kind=4),save            ::  max_side_num
    integer(kind=4),save            ::  periodic_bc_on
                
end module

!dec$ endif