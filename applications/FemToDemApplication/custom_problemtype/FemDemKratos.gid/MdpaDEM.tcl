
proc WriteMdpaDEM { basename dir problemtypedir } {

    ## Start MDPA file
	set basenameDEM $basename
	append basenameDEM "DEM"
    set filename [file join $dir ${basenameDEM}.mdpa]
    set FileVar [open $filename w]

    puts $FileVar ""
    puts $FileVar "Begin Properties 1"
	puts $FileVar "  PARTICLE_DENSITY                        [GiD_AccessValue get gendata Density]"
	puts $FileVar "  YOUNG_MODULUS                           [GiD_AccessValue get gendata Young_Modulus]"
	puts $FileVar "  POISSON_RATIO                           [GiD_AccessValue get gendata Poisson_Ratio]"
	puts $FileVar "  FRICTION                       [GiD_AccessValue get gendata Particle_Friction]"
	puts $FileVar "  PARTICLE_COHESION                       [GiD_AccessValue get gendata Cohesion]"
	puts $FileVar "  COEFFICIENT_OF_RESTITUTION              [GiD_AccessValue get gendata Coefficion_of_Restitution]"
	puts $FileVar "  PARTICLE_MATERIAL                       [GiD_AccessValue get gendata Color]"
	puts $FileVar "  ROLLING_FRICTION                        [GiD_AccessValue get gendata Rolling_Friction]"
	puts $FileVar "  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME  DEM_D_Linear_viscous_Coulomb"
	puts $FileVar "  DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME     DEMContinuumConstitutiveLaw"
    puts $FileVar "End Properties"
    puts $FileVar ""
    close $FileVar

	# create empty mdpa 
	set basenameCluster $basename
	append basenameCluster "DEM_Clusters"
	set filename [file join $dir ${basenameCluster}.mdpa]
    set FileVar [open $filename w]
	close $FileVar
	
	set basenameDEMFEM $basename
	append basenameDEMFEM "DEM_FEM_boundary"
	set filename [file join $dir ${basenameDEMFEM}.mdpa]
    set FileVar [open $filename w]
	close $FileVar
	
	set basenameDEMinlet $basename
	append basenameDEMinlet "DEM_Inlet"
	set filename [file join $dir ${basenameDEMinlet}.mdpa]
    set FileVar [open $filename w]
	close $FileVar

}
