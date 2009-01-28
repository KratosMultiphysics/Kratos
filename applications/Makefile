#####################################################################
#
#       MAKEFILE FOR Kratos ApplicationS
#		Authors: Carlos Labra and Riccardo Rossi
#
#####################################################################

ifeq "$(target)" ""
	target  = release;
endif

installapplications:  installale installstructure installfluid installconvdiff installpfem installexternalsolvers

cleanapplications: cleanale cleanstructure cleanfluid cleanconvdiff cleanpfem cleanexternalsolvers

all: ale structure fluid convdiff pfem

installale:
	cd ALEapplication/Linux; \
	make new target=$(target)

ale:
	cd ALEapplication/Linux; \
	make all target=$(target)

cleanale:
	cd ALEapplication/Linux; \
	make clean target=$(target) \
	make cleanlibs


	
installstructure:
	cd structural_application/Linux; \
	make new target=$(target)

structure:
	cd structural_application/Linux; \
	make all target=$(target)

cleanstructure:
	cd structural_application/Linux; \
	make clean target=$(target) \
	make cleanlibs;


	
installfluid:
	cd incompressible_fluid_application/Linux; \
	make new target=$(target)

fluid:
	cd incompressible_fluid_application/Linux; \
	make all target=$(target)

cleanfluid:
	cd incompressible_fluid_application/Linux; \
	make clean target=$(target) \
	make cleanlibs;



installconvdiff:
	cd convection_diffusion_application/Linux; \
	make new target=$(target)

convdiff:
	cd convection_diffusion_application/Linux; \
	make all target=$(target)

cleanconvdiff:
	cd convection_diffusion_application/Linux; \
	make clean target=$(target) \
	make cleanlibs



installpfem:
	cd PFEMapplication/Linux; \
	make new target=$(target)

pfem:
	cd PFEMapplication/Linux; \
	make all target=$(target)

cleanpfem:
	cd PFEMapplication/Linux; \
	make clean target=$(target) \
	make cleanlibs

installexternalsolvers:
	cd ExternalSolversApplication/Linux; \
	make new target=$(target)

externalsolvers:
	cd ExternalSolversApplication/Linux; \
	make all target=$(target)

cleanexternalsolvers:
	cd ExternalSolversApplication/Linux; \
	make clean target=$(target) \
	make cleanlibs
