proc WritePoroMaterials { basename dir problemtypedir } {

    ## Source auxiliar procedures
    # source [file join $problemtypedir PoroMaterialsAuxProcs.tcl]

    ## Start PoroMaterials.json file
    set filename [file join $dir PoroMaterials.json]
    set FileVar [open $filename w]

    puts $FileVar "\{"

    puts $FileVar "    \"properties\": \[\{"

    # Total number of groups with properties
    set Groups [GiD_Info conditions Body_Part groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Interface_Part groups]
    incr NumGroups [llength $Groups]
    set iGroup 0

    set PropertyId 0
    # set PropertyDict [dict create]
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 3] eq "LinearElasticSolid3DLaw"} {
            incr PropertyId
            puts $FileVar "        \"model_part_name\": \"PorousModelPart.[lindex [lindex $Groups $i] 1]\","
            puts $FileVar "        \"properties_id\": $PropertyId,"
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            if { ([GiD_AccessValue get gendata Initial_Stresses] eq false) || (([GiD_AccessValue get gendata Initial_Stresses] eq true) && ([GiD_AccessValue get gendata Mode] eq "save")) } {
                puts $FileVar "  CONSTITUTIVE_LAW_NAME LinearElasticSolid3DLaw"
            } else {
                puts $FileVar "  CONSTITUTIVE_LAW_NAME HistoryLinearElastic3DLaw"
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
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif { ([lindex [lindex $Groups $i] 3] eq "LinearElasticPlaneStrainSolid2DLaw") || ([lindex [lindex $Groups $i] 3] eq "LinearElasticPlaneStressSolid2DLaw")} {
            incr PropertyId
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            puts $FileVar "Begin Properties $PropertyId"
            if { ([GiD_AccessValue get gendata Initial_Stresses] eq false) || (([GiD_AccessValue get gendata Initial_Stresses] eq true) && ([GiD_AccessValue get gendata Mode] eq "save")) } {
                puts $FileVar "  CONSTITUTIVE_LAW_NAME [lindex [lindex $Groups $i] 3]"
            } else {
                if {[lindex [lindex $Groups $i] 3] eq "LinearElasticPlaneStrainSolid2DLaw"} {
                    puts $FileVar "  CONSTITUTIVE_LAW_NAME HistoryLinearElasticPlaneStrain2DLaw"
                } else {
                    puts $FileVar "  CONSTITUTIVE_LAW_NAME HistoryLinearElasticPlaneStress2DLaw"
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
            puts $FileVar "End Properties"
            puts $FileVar ""
        } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJuDamage3DLaw"} {
            incr PropertyId
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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
            # dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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




    ## Processes
    puts $FileVar "    \"processes\": \{"
    ## constraints_process_list
    set Groups [GiD_Info conditions Solid_Displacement groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    incr NumGroups [llength $Groups]
    set iGroup 0
    puts $FileVar "        \"constraints_process_list\": \[\{"
    # Solid_Displacement
    set Groups [GiD_Info conditions Solid_Displacement groups]
    WriteConstraintVectorProcess FileVar iGroup $Groups volumes DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups surfaces DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups lines DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups points DISPLACEMENT $TableDict $NumGroups
    # Note: it is important to write processes in the following order to account for intersections between conditions
    # Fluid_Pressure
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    WritePressureConstraintProcess FileVar iGroup $Groups volumes WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups surfaces WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups lines WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups points WATER_PRESSURE $TableDict $NumGroups
    ## loads_process_list
    set Groups [GiD_Info conditions Force groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Face_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Normal_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Body_Acceleration groups]
    incr NumGroups [llength $Groups]
    if {$NumGroups > 0} {
        set iGroup 0
        puts $FileVar "        \"loads_process_list\": \[\{"
        # Force
        set Groups [GiD_Info conditions Force groups]
        WriteLoadVectorProcess FileVar iGroup $Groups FORCE $TableDict $NumGroups
        # Face_Load
        set Groups [GiD_Info conditions Face_Load groups]
        WriteLoadVectorProcess FileVar iGroup $Groups FACE_LOAD $TableDict $NumGroups
        # Normal_Load
        set Groups [GiD_Info conditions Normal_Load groups]
        WriteNormalLoadProcess FileVar iGroup $Groups NORMAL_CONTACT_STRESS $TableDict $NumGroups
        # Normal_Fluid_Flux
        set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
        WriteLoadScalarProcess FileVar iGroup $Groups NORMAL_FLUID_FLUX $TableDict $NumGroups
        # Interface_Face_Load
        set Groups [GiD_Info conditions Interface_Face_Load groups]
        WriteLoadVectorProcess FileVar iGroup $Groups FACE_LOAD $TableDict $NumGroups
        # Interface_Normal_Fluid_Flux
        set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
        WriteLoadScalarProcess FileVar iGroup $Groups NORMAL_FLUID_FLUX $TableDict $NumGroups
        # Body_Acceleration
        set Groups [GiD_Info conditions Body_Acceleration groups]
        WriteLoadVectorProcess FileVar iGroup $Groups VOLUME_ACCELERATION $TableDict $NumGroups
    } else {
        puts $FileVar "        \"loads_process_list\":       \[\],"
    }
    ## auxiliar_process_list
    set NumGroups 0
    if {$IsPeriodic eq true} {
        set Groups [GiD_Info conditions Interface_Part groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            if {[lindex [lindex $Groups $i] 20] eq true} {
                incr NumGroups
            }
        }
    }
    if {$NumGroups > 0} {
        set iGroup 0
        puts $FileVar "        \"auxiliar_process_list\": \[\{"
        # Periodic_Bars
        if {$IsPeriodic eq true} {
            set Groups [GiD_Info conditions Interface_Part groups]
            WritePeriodicInterfaceProcess FileVar iGroup $Groups $NumGroups
        }
    } else {
        puts $FileVar "        \"auxiliar_process_list\": \[\]"
    }

    puts $FileVar "    \}"
    puts $FileVar "\}"

    close $FileVar
}
