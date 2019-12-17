proc WritePoroMaterials { basename dir problemtypedir } {

    ## Source auxiliar procedures
    # source [file join $problemtypedir PoroMaterialsAuxProcs.tcl]

    ## Start PoroMaterials.json file
    set filename [file join $dir PoroMaterials.json]
    set FileVar [open $filename w]

    puts $FileVar "\{"



    set PropertyId 0
    set PropertyDict [dict create]
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 3] eq "LinearElasticSolid3DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
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


    ## problem_domain_sub_model_part_list
    set PutStrings \[
    # Body_Part
    AppendGroupNames PutStrings Body_Part
    # Interface_Part
    AppendGroupNames PutStrings Interface_Part
    if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
        append PutStrings \" PropagationUnion_3d_6 \" \]
    } else {
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
    }
    puts $FileVar "        \"problem_domain_sub_model_part_list\": $PutStrings,"
    ## processes_sub_model_part_list
    set PutStrings \[
    # Solid_Displacement
    AppendGroupNames PutStrings Solid_Displacement
    # Fluid_Pressure
    AppendGroupNames PutStrings Fluid_Pressure
    # Force
    AppendGroupNames PutStrings Force
    # Face_Load
    AppendGroupNames PutStrings Face_Load
    # Normal_Load
    AppendGroupNames PutStrings Normal_Load
    # Normal_Fluid_Flux
    AppendGroupNames PutStrings Normal_Fluid_Flux
    # Interface_Face_Load
    AppendGroupNames PutStrings Interface_Face_Load
    # Interface_Normal_Fluid_Flux
    AppendGroupNames PutStrings Interface_Normal_Fluid_Flux
    # Body_Acceleration
    AppendGroupNames PutStrings Body_Acceleration
    # Periodic_Bars
    if {$IsPeriodic eq true} {
        set Groups [GiD_Info conditions Interface_Part groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            if {[lindex [lindex $Groups $i] 20] eq true} {
                append PutStrings \" Periodic_Bars_[lindex [lindex $Groups $i] 1] \" ,
            }
        }
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar "        \"processes_sub_model_part_list\":      $PutStrings,"
    ## body_domain_sub_model_part_list
    set PutStrings \[
    AppendGroupNames PutStrings Body_Part
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    if {[GiD_AccessValue get gendata Strategy_Type] eq "arc_length"} {
        puts $FileVar "        \"body_domain_sub_model_part_list\":    $PutStrings,"
        ## loads_sub_model_part_list
        set PutStrings \[
        set iGroup 0
        # Force
        AppendGroupNamesWithNum PutStrings iGroup Force
        # Face_Load
        AppendGroupNamesWithNum PutStrings iGroup Face_Load
        # Normal_Load
        AppendGroupNamesWithNum PutStrings iGroup Normal_Load
        # Normal_Fluid_Flux
        AppendGroupNamesWithNum PutStrings iGroup Normal_Fluid_Flux
        # Interface_Face_Load
        AppendGroupNamesWithNum PutStrings iGroup Interface_Face_Load
        # Interface_Normal_Fluid_Flux
        AppendGroupNamesWithNum PutStrings iGroup Interface_Normal_Fluid_Flux
        # Body_Acceleration
        AppendGroupNamesWithNum PutStrings iGroup Body_Acceleration
        if {$iGroup > 0} {
            set PutStrings [string trimright $PutStrings ,]
        }
        append PutStrings \]
        puts $FileVar "        \"loads_sub_model_part_list\":          $PutStrings,"
        ## loads_variable_list
        set PutStrings \[
        # Force
        AppendGroupVariables PutStrings Force FORCE
        # Face_Load
        AppendGroupVariables PutStrings Face_Load FACE_LOAD
        # Normal_Load
        AppendGroupVariables PutStrings Normal_Load NORMAL_CONTACT_STRESS
        # Normal_Fluid_Flux
        AppendGroupVariables PutStrings Normal_Fluid_Flux NORMAL_FLUID_FLUX
        # Interface_Face_Load
        AppendGroupVariables PutStrings Interface_Face_Load FACE_LOAD
        # Interface_Normal_Fluid_Flux
        AppendGroupVariables PutStrings Interface_Normal_Fluid_Flux NORMAL_FLUID_FLUX
        # Body_Acceleration
        AppendGroupVariables PutStrings Body_Acceleration VOLUME_ACCELERATION
        if {$iGroup > 0} {
            set PutStrings [string trimright $PutStrings ,]
        }
        append PutStrings \]
        puts $FileVar "        \"loads_variable_list\":                $PutStrings"
        puts $FileVar "    \},"
    } else {
        puts $FileVar "        \"body_domain_sub_model_part_list\":    $PutStrings"
        puts $FileVar "    \},"
    }

    ## Output processes
    puts $FileVar "    \"output_processes\": \{"
    puts $FileVar "        \"gid_output\": \[\{"
    puts $FileVar "            \"python_module\": \"gid_output_process\","
    puts $FileVar "            \"kratos_module\": \"KratosMultiphysics\","
    puts $FileVar "            \"process_name\": \"GiDOutputProcess\","
    puts $FileVar "            \"Parameters\":    \{"
    puts $FileVar "                \"model_part_name\": \"PorousModelPart.porous_computational_model_part\","
    puts $FileVar "                \"output_name\": \"$basename\","
    puts $FileVar "                \"postprocess_parameters\": \{"
    puts $FileVar "                    \"result_file_configuration\": \{"
    puts $FileVar "                        \"gidpost_flags\":       \{"
    puts $FileVar "                            \"WriteDeformedMeshFlag\": \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $FileVar "                            \"WriteConditionsFlag\":   \"[GiD_AccessValue get gendata Write_conditions]\","
    if { ([GiD_AccessValue get gendata Fracture_Propagation] eq true) || ($IsPeriodic eq true) } {
        puts $FileVar "                            \"GiDPostMode\":           \"GiD_PostAscii\","
        puts $FileVar "                            \"MultiFileFlag\":         \"MultipleFiles\""
        puts $FileVar "                        \},"
        puts $FileVar "                        \"file_label\":          \"time\","
    } else {
        puts $FileVar "                            \"GiDPostMode\":           \"[GiD_AccessValue get gendata GiD_post_mode]\","
        puts $FileVar "                            \"MultiFileFlag\":         \"[GiD_AccessValue get gendata Multi_file_flag]\""
        puts $FileVar "                        \},"
        puts $FileVar "                        \"file_label\":          \"[GiD_AccessValue get gendata File_label]\","
    }
    puts $FileVar "                        \"output_control_type\": \"[GiD_AccessValue get gendata Output_control_type]\","
    puts $FileVar "                        \"output_frequency\":    [GiD_AccessValue get gendata Output_frequency],"
    puts $FileVar "                        \"body_output\":         [GiD_AccessValue get gendata Body_output],"
    puts $FileVar "                        \"node_output\":         [GiD_AccessValue get gendata Node_output],"
    puts $FileVar "                        \"skin_output\":         [GiD_AccessValue get gendata Skin_output],"
    puts $FileVar "                        \"plane_output\":        \[\],"
    # nodal_results
    set PutStrings \[
    set iGroup 0
    AppendOutputVariables PutStrings iGroup Write_Solid_Displacement DISPLACEMENT
    AppendOutputVariables PutStrings iGroup Write_Fluid_Pressure WATER_PRESSURE
    if {[GiD_AccessValue get gendata Write_Reactions] eq true} {
        incr iGroup
        append PutStrings \" REACTION \" , \" REACTION_WATER_PRESSURE \" ,
    }
    AppendOutputVariables PutStrings iGroup Write_Force FORCE
    AppendOutputVariables PutStrings iGroup Write_Face_Load FACE_LOAD
    AppendOutputVariables PutStrings iGroup Write_Normal_Load NORMAL_CONTACT_STRESS
    AppendOutputVariables PutStrings iGroup Write_Tangential_Load TANGENTIAL_CONTACT_STRESS
    AppendOutputVariables PutStrings iGroup Write_Normal_Fluid_Flux NORMAL_FLUID_FLUX
    AppendOutputVariables PutStrings iGroup Write_Body_Acceleration VOLUME_ACCELERATION
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        incr iGroup
        append PutStrings \" PARTITION_INDEX \" ,
    }
    # Nodal smoothed variables
    if {[GiD_AccessValue get gendata Nodal_Smoothing] eq true} {
        AppendOutputVariables PutStrings iGroup Write_Effective_Stress NODAL_CAUCHY_STRESS_TENSOR
        AppendOutputVariables PutStrings iGroup Write_Damage NODAL_DAMAGE_VARIABLE
        AppendOutputVariables PutStrings iGroup Write_Joint_Width NODAL_JOINT_WIDTH
        AppendOutputVariables PutStrings iGroup Write_Damage NODAL_JOINT_DAMAGE
    }
    AppendOutputVariables PutStrings iGroup Write_Initial_Stress INITIAL_STRESS_TENSOR
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $FileVar "                        \"nodal_results\":       $PutStrings,"
    # gauss_point_results
    set PutStrings \[
    set iGroup 0
    AppendOutputVariables PutStrings iGroup Write_Strain GREEN_LAGRANGE_STRAIN_TENSOR
    AppendOutputVariables PutStrings iGroup Write_Effective_Stress CAUCHY_STRESS_TENSOR
    AppendOutputVariables PutStrings iGroup Write_Total_Stress TOTAL_STRESS_TENSOR
    AppendOutputVariables PutStrings iGroup Write_Von_Mises_Stress VON_MISES_STRESS
    AppendOutputVariables PutStrings iGroup Write_Fluid_Flux FLUID_FLUX_VECTOR
    AppendOutputVariables PutStrings iGroup Write_Permeability PERMEABILITY_MATRIX
    AppendOutputVariables PutStrings iGroup Write_Damage DAMAGE_VARIABLE
    AppendOutputVariables PutStrings iGroup Write_Joint_Width JOINT_WIDTH
    AppendOutputVariables PutStrings iGroup Write_Local_Stress_Vector LOCAL_STRESS_VECTOR
    AppendOutputVariables PutStrings iGroup Write_Local_Relative_Displacement LOCAL_RELATIVE_DISPLACEMENT_VECTOR
    AppendOutputVariables PutStrings iGroup Write_Local_Fluid_Flux LOCAL_FLUID_FLUX_VECTOR
    AppendOutputVariables PutStrings iGroup Write_Local_Permeability LOCAL_PERMEABILITY_MATRIX
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $FileVar "                        \"gauss_point_results\": $PutStrings"
    puts $FileVar "                    \},"
    puts $FileVar "                    \"point_data_configuration\":  \[\]"
    puts $FileVar "                \}"
    puts $FileVar "            \}"
    puts $FileVar "        \}\]"
    puts $FileVar "    \},"

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
