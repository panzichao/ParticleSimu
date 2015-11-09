# make file for ParticleSimu
# Start of the makefile
ParticleSimu: ParticleSimu.obj ParticleSimu_initialization.obj ParticleSimu_global_parameter.obj ParticleSimu_2D_kernel.obj ParticleSimu_3D_kernel.obj ParticleSimu_general_sub.obj
	ifort ParticleSimu.obj ParticleSimu_initialization.obj ParticleSimu_global_parameter.obj ParticleSimu_2D_kernel.obj ParticleSimu_3D_kernel.obj ParticleSimu_general_sub.obj

ParticleSimu.obj: ParticleSimu_initialization.mod ParticleSimu_global_parameter.mod ParticleSimu_general_sub.mod ParticleSimu_2D_kernel.mod ParticleSimu_3D_kernel.mod ParticleSimu.f90
	ifort /c ParticleSimu.f90

ParticleSimu_initialization.mod: ParticleSimu_initialization.obj ParticleSimu_initialization.f90
	ifort /c ParticleSimu_initialization.f90
ParticleSimu_initialization.obj: ParticleSimu_global_parameter.mod ParticleSimu_initialization.f90
	ifort /c ParticleSimu_initialization.f90

ParticleSimu_global_parameter.mod: ParticleSimu_global_parameter.obj ParticleSimu_global_parameter.f90
	ifort /c ParticleSimu_global_parameter.f90
ParticleSimu_global_parameter.obj: ParticleSimu_global_parameter.f90
	ifort /c ParticleSimu_global_parameter.f90

ParticleSimu_2D_kernel.mod: ParticleSimu_2D_kernel.obj ParticleSimu_2D_kernel.f90
	ifort /c ParticleSimu_2D_kernel.f90
ParticleSimu_2D_kernel.obj: ParticleSimu_global_parameter.mod ParticleSimu_general_sub.mod ParticleSimu_2D_kernel.f90
	ifort /c ParticleSimu_2D_kernel.f90

ParticleSimu_3D_kernel.mod: ParticleSimu_3D_kernel.obj ParticleSimu_3D_kernel.f90
	ifort /c ParticleSimu_3D_kernel.f90
ParticleSimu_3D_kernel.obj: ParticleSimu_global_parameter.mod ParticleSimu_general_sub.mod ParticleSimu_3D_kernel.f90
	ifort /c ParticleSimu_3D_kernel.f90

ParticleSimu_general_sub.mod: ParticleSimu_general_sub.obj ParticleSimu_general_sub.f90
	ifort /c ParticleSimu_general_sub.f90
ParticleSimu_general_sub.obj: ParticleSimu_general_sub.f90
	ifort /c ParticleSimu_general_sub.f90

clean:
	main.obj

 # End of the makefile