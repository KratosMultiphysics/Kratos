proc WriteMdpa { basename dir problemtypedir } {

    ## Source auxiliar procedures
    source [file join $problemtypedir MdpaAuxProcs.tcl]

    ## Start MDPA file
    set filename [file join $dir ${basename}.mdpa]
    set FileVar [open $filename w]

    ## ModelPart Data
    #puts $FileVar "Begin ModelPartData"
    #puts $FileVar "  // VARIABLE_NAME value"
    #puts $FileVar "End ModelPartData"
    #puts $FileVar ""
    #puts $FileVar ""

    ## Tables
    set TableId 0
    set TableDict [dict create]
    # Solid_Displacement
    ConstraintVectorTable FileVar TableId TableDict Solid_Displacement DISPLACEMENT
    # Fluid_Pressure
    PressureTable FileVar TableId TableDict Fluid_Pressure WATER_PRESSURE
    # Force
    VectorTable FileVar TableId TableDict Force FORCE
    # Face_Load
    VectorTable FileVar TableId TableDict Face_Load FACE_LOAD
    # Normal_Load
    NormalTangentialTable FileVar TableId TableDict Normal_Load NORMAL_CONTACT_STRESS TANGENTIAL_CONTACT_STRESS
    # Normal_Fluid_Flux
    ScalarTable FileVar TableId TableDict Normal_Fluid_Flux NORMAL_FLUID_FLUX
    # Interface_Face_Load
    VectorTable FileVar TableId TableDict Interface_Face_Load FACE_LOAD
    # Interface_Normal_Fluid_Flux
    ScalarTable FileVar TableId TableDict Interface_Normal_Fluid_Flux NORMAL_FLUID_FLUX
    # Body_Acceleration
    VectorTable FileVar TableId TableDict Body_Acceleration VOLUME_ACCELERATION
    puts $FileVar ""

    ## Properties
    set PropertyId 0
    set PropertyDict [dict create]
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 3] eq "LinearElastic3DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "  CONSTITUTIVE_LAW_NAME LinearElastic3DLaw"
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  PERMEABILITY_ZZ [lindex [lindex $Groups $i] 13]"
            puts $FileVar "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $FileVar "  PERMEABILITY_YZ [lindex [lindex $Groups $i] 15]"
            puts $FileVar "  PERMEABILITY_ZX [lindex [lindex $Groups $i] 16]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif { ([lindex [lindex $Groups $i] 3] eq "LinearElasticPlaneStrain2DLaw") || ([lindex [lindex $Groups $i] 3] eq "LinearElasticPlaneStress2DLaw")} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "  CONSTITUTIVE_LAW_NAME [lindex [lindex $Groups $i] 3]"
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $FileVar "  THICKNESS [lindex [lindex $Groups $i] 18]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJuDamage3DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            if {[GiD_AccessValue get gendata Non-local_Damage] eq true} {
                puts $FileVar "  CONSTITUTIVE_LAW_NAME SimoJuNonlocalDamage3DLaw"
            } else {
                puts $FileVar "  CONSTITUTIVE_LAW_NAME SimoJuLocalDamage3DLaw"
            }
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  PERMEABILITY_ZZ [lindex [lindex $Groups $i] 13]"
            puts $FileVar "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $FileVar "  PERMEABILITY_YZ [lindex [lindex $Groups $i] 15]"
            puts $FileVar "  PERMEABILITY_ZX [lindex [lindex $Groups $i] 16]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $FileVar "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 19]"
            puts $FileVar "  STRENGTH_RATIO [lindex [lindex $Groups $i] 20]"
            puts $FileVar "  FRACTURE_ENERGY [lindex [lindex $Groups $i] 21]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif {([lindex [lindex $Groups $i] 3] eq "SimoJuDamagePlaneStrain2DLaw") || ([lindex [lindex $Groups $i] 3] eq "SimoJuDamagePlaneStress2DLaw")} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            if {[GiD_AccessValue get gendata Non-local_Damage] eq true} {
                if {[lindex [lindex $Groups $i] 3] eq "SimoJuDamagePlaneStrain2DLaw"} {
                    puts $FileVar "  CONSTITUTIVE_LAW_NAME SimoJuNonlocalDamagePlaneStrain2DLaw"
                } else {
                    puts $FileVar "  CONSTITUTIVE_LAW_NAME SimoJuNonlocalDamagePlaneStress2DLaw"
                }
            } else {
                if {[lindex [lindex $Groups $i] 3] eq "SimoJuDamagePlaneStrain2DLaw"} {
                    puts $FileVar "  CONSTITUTIVE_LAW_NAME SimoJuLocalDamagePlaneStrain2DLaw"
                } else {
                    puts $FileVar "  CONSTITUTIVE_LAW_NAME SimoJuLocalDamagePlaneStress2DLaw"
                }
            }
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $FileVar "  THICKNESS [lindex [lindex $Groups $i] 18]"
            puts $FileVar "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 19]"
            puts $FileVar "  STRENGTH_RATIO [lindex [lindex $Groups $i] 20]"
            puts $FileVar "  FRACTURE_ENERGY [lindex [lindex $Groups $i] 21]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif {[lindex [lindex $Groups $i] 3] eq "ModifiedMisesDamage3DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "  CONSTITUTIVE_LAW_NAME ModifiedMisesNonlocalDamage3DLaw"
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  PERMEABILITY_ZZ [lindex [lindex $Groups $i] 13]"
            puts $FileVar "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $FileVar "  PERMEABILITY_YZ [lindex [lindex $Groups $i] 15]"
            puts $FileVar "  PERMEABILITY_ZX [lindex [lindex $Groups $i] 16]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $FileVar "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 19]"
            puts $FileVar "  STRENGTH_RATIO [lindex [lindex $Groups $i] 20]"
            puts $FileVar "  RESIDUAL_STRENGTH [lindex [lindex $Groups $i] 22]"
            puts $FileVar "  SOFTENING_SLOPE [lindex [lindex $Groups $i] 23]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif {[lindex [lindex $Groups $i] 3] eq "ModifiedMisesDamagePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 3] eq "ModifiedMisesDamagePlaneStress2DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            if {[lindex [lindex $Groups $i] 3] eq "ModifiedMisesDamagePlaneStrain2DLaw"} {
                puts $FileVar "  CONSTITUTIVE_LAW_NAME ModifiedMisesNonlocalDamagePlaneStrain2DLaw"
            } else {
                puts $FileVar "  CONSTITUTIVE_LAW_NAME ModifiedMisesNonlocalDamagePlaneStress2DLaw"
            }
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  PERMEABILITY_XX [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  PERMEABILITY_YY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  PERMEABILITY_XY [lindex [lindex $Groups $i] 14]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 17]"
            puts $FileVar "  THICKNESS [lindex [lindex $Groups $i] 18]"
            puts $FileVar "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 19]"
            puts $FileVar "  STRENGTH_RATIO [lindex [lindex $Groups $i] 20]"
            puts $FileVar "  RESIDUAL_STRENGTH [lindex [lindex $Groups $i] 22]"
            puts $FileVar "  SOFTENING_SLOPE [lindex [lindex $Groups $i] 23]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        }
    }
    # Interface_Part
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 4] eq "BilinearCohesive3DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "  CONSTITUTIVE_LAW_NAME BilinearCohesive3DLaw"
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  TRANSVERSAL_PERMEABILITY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 13]"
            puts $FileVar "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 15]"
            puts $FileVar "  MINIMUM_JOINT_WIDTH [lindex [lindex $Groups $i] 16]"
            puts $FileVar "  CRITICAL_DISPLACEMENT [lindex [lindex $Groups $i] 17]"
            puts $FileVar "  YIELD_STRESS [lindex [lindex $Groups $i] 18]"
            puts $FileVar "  FRICTION_COEFFICIENT [lindex [lindex $Groups $i] 19]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif {[lindex [lindex $Groups $i] 4] eq "BilinearCohesivePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 4] eq "BilinearCohesivePlaneStress2DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "  CONSTITUTIVE_LAW_NAME BilinearCohesive2DLaw"
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  TRANSVERSAL_PERMEABILITY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 13]"
            puts $FileVar "  THICKNESS [lindex [lindex $Groups $i] 14]"
            puts $FileVar "  DAMAGE_THRESHOLD [lindex [lindex $Groups $i] 15]"
            puts $FileVar "  MINIMUM_JOINT_WIDTH [lindex [lindex $Groups $i] 16]"
            puts $FileVar "  CRITICAL_DISPLACEMENT [lindex [lindex $Groups $i] 17]"
            puts $FileVar "  YIELD_STRESS [lindex [lindex $Groups $i] 18]"
            puts $FileVar "  FRICTION_COEFFICIENT [lindex [lindex $Groups $i] 19]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif {[lindex [lindex $Groups $i] 4] eq "ExponentialCohesive3DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "  CONSTITUTIVE_LAW_NAME ExponentialCohesive3DLaw"
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  TRANSVERSAL_PERMEABILITY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 13]"
            puts $FileVar "  MINIMUM_JOINT_WIDTH [lindex [lindex $Groups $i] 16]"
            puts $FileVar "  YIELD_STRESS [lindex [lindex $Groups $i] 18]"
            puts $FileVar "  FRACTURE_ENERGY [lindex [lindex $Groups $i] 22]"
            puts $FileVar "  SHEAR_FRACTURE_ENERGY [lindex [lindex $Groups $i] 23]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif {[lindex [lindex $Groups $i] 4] eq "ExponentialCohesivePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 4] eq "ExponentialCohesivePlaneStress2DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            puts $FileVar "  CONSTITUTIVE_LAW_NAME ExponentialCohesive2DLaw"
            puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 5]"
            puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 6]"
            puts $FileVar "  DENSITY_SOLID [lindex [lindex $Groups $i] 7]"
            puts $FileVar "  DENSITY_WATER [lindex [lindex $Groups $i] 8]"
            puts $FileVar "  POROSITY [lindex [lindex $Groups $i] 9]"
            puts $FileVar "  BULK_MODULUS_SOLID [lindex [lindex $Groups $i] 10]"
            puts $FileVar "  BULK_MODULUS_FLUID [lindex [lindex $Groups $i] 11]"
            puts $FileVar "  TRANSVERSAL_PERMEABILITY [lindex [lindex $Groups $i] 12]"
            puts $FileVar "  DYNAMIC_VISCOSITY [lindex [lindex $Groups $i] 13]"
            puts $FileVar "  THICKNESS [lindex [lindex $Groups $i] 14]"
            puts $FileVar "  MINIMUM_JOINT_WIDTH [lindex [lindex $Groups $i] 16]"
            puts $FileVar "  YIELD_STRESS [lindex [lindex $Groups $i] 18]"
            puts $FileVar "  FRACTURE_ENERGY [lindex [lindex $Groups $i] 22]"
            puts $FileVar "  SHEAR_FRACTURE_ENERGY [lindex [lindex $Groups $i] 23]"
            puts $FileVar "End Properties"
            puts $FileVar ""
        }
    }
    puts $FileVar ""

    ## Nodes
    set Nodes [GiD_Info Mesh Nodes]
    puts $FileVar "Begin Nodes"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
        # puts $FileVar "  [lindex $Nodes $i]  [lindex $Nodes [expr { $i+1 }]] [lindex $Nodes [expr { $i+2 }]] [lindex $Nodes [expr { $i+3 }]]"
        puts -nonewline $FileVar "  [lindex $Nodes $i]  "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+1 }]]]
        puts -nonewline $FileVar " "
        puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+2 }]]]
        puts -nonewline $FileVar " "
        puts $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+3 }]]]
    }
    puts $FileVar "End Nodes"
    puts $FileVar ""
    puts $FileVar ""

    ## Elements
    set FIC [GiD_AccessValue get gendata FIC_Stabilization]
    set IsQuadratic [GiD_Info Project Quadratic]
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    if {$IsQuadratic eq 0} {
        if {$FIC eq false} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

                # UPwSmallStrainElement2D3N
                WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainElement2D3N $BodyElemsProp Triangle2D3Connectivities
                # UPwSmallStrainElement2D4N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # UPwSmallStrainElement3D4N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # UPwSmallStrainElement3D8N
                WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                # Elements Property
                set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

                # UPwSmallStrainFICElement2D3N
                WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainFICElement2D3N $BodyElemsProp Triangle2D3Connectivities
                # UPwSmallStrainFICElement2D4N
                WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainFICElement2D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # UPwSmallStrainFICElement3D4N
                WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainFICElement3D4N $BodyElemsProp Quadrilateral2D4Connectivities
                # UPwSmallStrainFICElement3D8N
                WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainFICElement3D8N $BodyElemsProp Hexahedron3D8Connectivities
            }
        }
    } elseif {$IsQuadratic eq 1} {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

            # SmallStrainUPwDiffOrderElement2D6N
            WriteElements FileVar [lindex $Groups $i] triangle SmallStrainUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
            # SmallStrainUPwDiffOrderElement2D8N
            WriteElements FileVar [lindex $Groups $i] quadrilateral SmallStrainUPwDiffOrderElement2D8N $BodyElemsProp Hexahedron3D8Connectivities
            # SmallStrainUPwDiffOrderElement3D10N
            WriteElements FileVar [lindex $Groups $i] tetrahedra SmallStrainUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
            # SmallStrainUPwDiffOrderElement3D20N
            WriteElements FileVar [lindex $Groups $i] hexahedra SmallStrainUPwDiffOrderElement3D20N $BodyElemsProp Hexahedron3D20Connectivities
        }
    } else {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            # Elements Property
            set BodyElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]

            # SmallStrainUPwDiffOrderElement2D6N
            WriteElements FileVar [lindex $Groups $i] triangle SmallStrainUPwDiffOrderElement2D6N $BodyElemsProp Triangle2D6Connectivities
            # SmallStrainUPwDiffOrderElement2D9N
            WriteElements FileVar [lindex $Groups $i] quadrilateral SmallStrainUPwDiffOrderElement2D9N $BodyElemsProp Quadrilateral2D9Connectivities
            # SmallStrainUPwDiffOrderElement3D10N
            WriteElements FileVar [lindex $Groups $i] tetrahedra SmallStrainUPwDiffOrderElement3D10N $BodyElemsProp Tetrahedron3D10Connectivities
            # SmallStrainUPwDiffOrderElement3D27N
            WriteElements FileVar [lindex $Groups $i] hexahedra SmallStrainUPwDiffOrderElement3D27N $BodyElemsProp Hexahedron3D27Connectivities
        }
    }
    # Interface_Part
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 3] eq false} {
            # Elements Property
            set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            # UPwSmallStrainInterfaceElement2D4N
            WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainInterfaceElement2D4N $InterfaceElemsProp Quadrilateral2D4Connectivities
            # UPwSmallStrainInterfaceElement3D6N
            WriteElements FileVar [lindex $Groups $i] prism UPwSmallStrainInterfaceElement3D6N $InterfaceElemsProp PrismInterface3D6Connectivities
            # UPwSmallStrainInterfaceElement3D8N
            WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainInterfaceElement3D8N $InterfaceElemsProp HexahedronInterface3D8Connectivities
        } else {
            # Elements Property
            set LinkInterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
            # UPwSmallStrainLinkInterfaceElement2D4N
            WriteElements FileVar [lindex $Groups $i] quadrilateral UPwSmallStrainLinkInterfaceElement2D4N $LinkInterfaceElemsProp QuadrilateralInterface2D4Connectivities
            WriteElements FileVar [lindex $Groups $i] triangle UPwSmallStrainLinkInterfaceElement2D4N $LinkInterfaceElemsProp TriangleInterface2D4Connectivities
            # UPwSmallStrainLinkInterfaceElement3D6N
            WriteElements FileVar [lindex $Groups $i] prism UPwSmallStrainLinkInterfaceElement3D6N $LinkInterfaceElemsProp Triangle2D6Connectivities
            WriteElements FileVar [lindex $Groups $i] tetrahedra UPwSmallStrainLinkInterfaceElement3D6N $LinkInterfaceElemsProp TetrahedronInterface3D6Connectivities
            # UPwSmallStrainLinkInterfaceElement3D8N
            WriteElements FileVar [lindex $Groups $i] hexahedra UPwSmallStrainLinkInterfaceElement3D8N $LinkInterfaceElemsProp Hexahedron3D8Connectivities
        }
    }
    # PropagationUnion (InterfaceElement)
    if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
        # UPwSmallStrainInterfaceElement3D6N
        set PropUnionElementList [WritePropUnionElements FileVar $InterfaceElemsProp]
    }
    puts $FileVar ""

    ## Conditions
    set ConditionId 0
    set ConditionDict [dict create]
    set Dim [GiD_AccessValue get gendata Domain_Size]
    # Force
    set Groups [GiD_Info conditions Force groups]
    if {$Dim eq 2} {
        # UPwForceCondition2D1N
        WriteNodalConditions FileVar ConditionId ConditionDict $Groups UPwForceCondition2D1N $BodyElemsProp
    } else {
        # UPwForceCondition3D1N
        WriteNodalConditions FileVar ConditionId ConditionDict $Groups UPwForceCondition3D1N $BodyElemsProp
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    if {$Dim eq 2} {
        if {$IsQuadratic eq 0} {
            # UPwFaceLoadCondition2D2N
            WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwFaceLoadCondition2D2N $PropertyDict
        } else {
            # LineLoadDiffOrderCondition2D3N
            WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineLoadDiffOrderCondition2D3N $PropertyDict
        }
    } else {
        if {$IsQuadratic eq 0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # UPwFaceLoadCondition3D3N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwFaceLoadCondition3D3N $PropertyDict
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwFaceLoadCondition3D3N $PropertyDict
                # UPwFaceLoadCondition3D4N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwFaceLoadCondition3D4N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        } elseif {$IsQuadratic eq 1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceLoadDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceLoadDiffOrderCondition3D6N $PropertyDict
                # SurfaceLoadDiffOrderCondition3D8N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceLoadDiffOrderCondition3D8N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceLoadDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceLoadDiffOrderCondition3D6N $PropertyDict
                # SurfaceLoadDiffOrderCondition3D9N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceLoadDiffOrderCondition3D9N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        }
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    if {$Dim eq 2} {
        if {$IsQuadratic eq 0} {
            # UPwNormalFaceLoadCondition2D2N
            WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFaceLoadCondition2D2N $PropertyDict
        } else {
            # LineNormalLoadDiffOrderCondition2D3N
            WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineNormalLoadDiffOrderCondition2D3N $PropertyDict
        }
    } else {
        if {$IsQuadratic eq 0} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # UPwNormalFaceLoadCondition3D3N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFaceLoadCondition3D3N $PropertyDict
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFaceLoadCondition3D3N $PropertyDict
                # UpwNormalFaceLoadCondition3D4N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UpwNormalFaceLoadCondition3D4N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        } elseif {$IsQuadratic eq 1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceNormalLoadDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalLoadDiffOrderCondition3D6N $PropertyDict
                # SurfaceNormalLoadDiffOrderCondition3D8N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalLoadDiffOrderCondition3D8N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceNormalLoadDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalLoadDiffOrderCondition3D6N $PropertyDict
                # SurfaceNormalLoadDiffOrderCondition3D9N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalLoadDiffOrderCondition3D9N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        }
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    if {$Dim eq 2} {
        if {$IsQuadratic eq 0} {
            if {$FIC eq false} {
                # UPwNormalFluxCondition2D2N
                WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFluxCondition2D2N $PropertyDict
            } else {
                # UPwNormalFluxFICCondition2D2N
                WriteFaceConditions FileVar ConditionId ConditionDict $Groups UPwNormalFluxFICCondition2D2N $PropertyDict
            }
        } else {
            # LineNormalFluidFluxDiffOrderCondition2D3N
            WriteFaceConditions FileVar ConditionId ConditionDict $Groups LineNormalFluidFluxDiffOrderCondition2D3N $PropertyDict
        }
    } else {
        if {$IsQuadratic eq 0} {
            if {$FIC eq false} {
                for {set i 0} {$i < [llength $Groups]} {incr i} {
                    set MyConditionList [list]
                    # UPwNormalFluxCondition3D3N
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFluxCondition3D3N $PropertyDict
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFluxCondition3D3N $PropertyDict
                    # UPwNormalFluxCondition3D4N
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwNormalFluxCondition3D4N $PropertyDict
                    dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
                }
            } else {
                for {set i 0} {$i < [llength $Groups]} {incr i} {
                    set MyConditionList [list]
                    # UPwNormalFluxFICCondition3D3N
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra UPwNormalFluxFICCondition3D3N $PropertyDict
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] prism UPwNormalFluxFICCondition3D3N $PropertyDict
                    # UPwNormalFluxFICCondition3D4N
                    WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra UPwNormalFluxFICCondition3D4N $PropertyDict
                    dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
                }
            }
        } elseif {$IsQuadratic eq 1} {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceNormalFluidFluxDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalFluidFluxDiffOrderCondition3D6N $PropertyDict
                # SurfaceNormalFluidFluxDiffOrderCondition3D8N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalFluidFluxDiffOrderCondition3D8N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        } else {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                set MyConditionList [list]
                # SurfaceNormalFluidFluxDiffOrderCondition3D6N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] tetrahedra SurfaceNormalFluidFluxDiffOrderCondition3D6N $PropertyDict
                # SurfaceNormalFluidFluxDiffOrderCondition3D9N
                WriteTypeFaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] hexahedra SurfaceNormalFluidFluxDiffOrderCondition3D9N $PropertyDict
                dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
            }
        }
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set MyConditionList [list]
        # UPwFaceLoadInterfaceCondition2D2N
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] linear UPwFaceLoadInterfaceCondition2D2N $InterfaceElemsProp Line2D2Connectivities
        # UPwFaceLoadInterfaceCondition3D4N
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] triangle UPwFaceLoadInterfaceCondition3D4N $InterfaceElemsProp TriangleInterface3D4Connectivities
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] quadrilateral UPwFaceLoadInterfaceCondition3D4N $InterfaceElemsProp QuadrilateralInterface3D4Connectivities
        dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set MyConditionList [list]
        # UPwNormalFluxInterfaceCondition2D2N
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] linear UPwNormalFluxInterfaceCondition2D2N $InterfaceElemsProp Line2D2Connectivities
        # UPwNormalFluxInterfaceCondition3D4N
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] triangle UPwNormalFluxInterfaceCondition3D4N $InterfaceElemsProp TriangleInterface3D4Connectivities
        WriteInterfaceConditions FileVar ConditionId MyConditionList [lindex $Groups $i] quadrilateral UPwNormalFluxInterfaceCondition3D4N $InterfaceElemsProp QuadrilateralInterface3D4Connectivities
        dict set ConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
    }

    # Periodic_Bars
    set IsPeriodic [GiD_AccessValue get gendata Periodic_Interface_Conditions]
    if {$IsPeriodic eq true} {
        set PeriodicBarsDict [dict create]
        set Groups [GiD_Info conditions Interface_Part groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            if {[lindex [lindex $Groups $i] 20] eq true} {
                # Elements Property
                set InterfaceElemsProp [dict get $PropertyDict [lindex [lindex $Groups $i] 1]]
                set ConditionList [list]
                # InterfaceElement2D4N
                SavePeriodicBarsFromIE2D4N PeriodicBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp
                # InterfaceElement3D6N
                SavePeriodicBarsFromIE3D6N PeriodicBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp
                # InterfaceElement3D8N
                SavePeriodicBarsFromIE3D8N PeriodicBarsDict ConditionId ConditionList [lindex $Groups $i] $InterfaceElemsProp

                dict set ConditionDict Periodic_Bars_[lindex [lindex $Groups $i] 1] $ConditionList
            }
        }

        if {[dict size $PeriodicBarsDict] > 0} {
            puts $FileVar "Begin Conditions PeriodicCondition"
            dict for {Name PeriodicBar} $PeriodicBarsDict {
                puts $FileVar "  [dict get $PeriodicBar Id]  [dict get $PeriodicBar PropertyId]  [dict get $PeriodicBar Connectivities]"
            }
            puts $FileVar "End Conditions"
            puts $FileVar ""
        }
    }

    puts $FileVar ""

    ## SubModelParts
    # Body_Part
    WriteElementSubmodelPart FileVar Body_Part
    # Interface_Part
    WriteElementSubmodelPart FileVar Interface_Part
    # PropagationUnion (InterfaceElement)
    if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
        WritePropUnionElementSubmodelPart FileVar $PropUnionElementList
    }
    # Solid_Displacement
    WriteConstraintSubmodelPart FileVar Solid_Displacement $TableDict
    # Fluid_Pressure
    WriteConstraintSubmodelPart FileVar Fluid_Pressure $TableDict
    # Force
    WriteLoadSubmodelPart FileVar Force $TableDict $ConditionDict
    # Face_Load
    WriteLoadSubmodelPart FileVar Face_Load $TableDict $ConditionDict
    # Normal_Load
    WriteLoadSubmodelPart FileVar Normal_Load $TableDict $ConditionDict
    # Normal_Fluid_Flux
    WriteLoadSubmodelPart FileVar Normal_Fluid_Flux $TableDict $ConditionDict
    # Interface_Face_Load
    WriteLoadSubmodelPart FileVar Interface_Face_Load $TableDict $ConditionDict
    # Interface_Normal_Fluid_Flux
    WriteLoadSubmodelPart FileVar Interface_Normal_Fluid_Flux $TableDict $ConditionDict
    # Body_Acceleration
    WriteConstraintSubmodelPart FileVar Body_Acceleration $TableDict

    # Periodic_Bars
    if {$IsPeriodic eq true} {
        WritePeriodicBarsSubmodelPart FileVar Interface_Part $ConditionDict
    }

    close $FileVar

    return $TableDict
}