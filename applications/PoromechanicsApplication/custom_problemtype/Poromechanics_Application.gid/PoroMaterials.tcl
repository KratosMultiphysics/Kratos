proc WritePoroMaterials { basename dir problemtypedir PropertyId } {

    upvar $PropertyId MyPropertyId

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
    incr NumGroups $MyPropertyId
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr MyPropertyId
        puts $FileVar "        \"model_part_name\": \"PorousModelPart.[lindex [lindex $Groups $i] 1]\","
        puts $FileVar "        \"properties_id\": $MyPropertyId,"
        puts $FileVar "        \"Material\": \{"
        puts $FileVar "            \"constitutive_law\": \{"
        if {[lindex [lindex $Groups $i] 3] eq "LinearElasticSolid3DLaw"} {
            if { ([GiD_AccessValue get gendata Initial_Stresses] eq false) || (([GiD_AccessValue get gendata Initial_Stresses] eq true) && ([GiD_AccessValue get gendata Mode] eq "save")) } {
                puts $FileVar "                \"name\": \"[lindex [lindex $Groups $i] 3]\""
            } else {
                puts $FileVar "                \"name\": \"HistoryLinearElastic3DLaw\""
            }
        } elseif {[lindex [lindex $Groups $i] 3] eq "LinearElasticPlaneStrainSolid2DLaw"} {
            if { ([GiD_AccessValue get gendata Initial_Stresses] eq false) || (([GiD_AccessValue get gendata Initial_Stresses] eq true) && ([GiD_AccessValue get gendata Mode] eq "save")) } {
                puts $FileVar "                \"name\": \"[lindex [lindex $Groups $i] 3]\""
            } else {
                puts $FileVar "                \"name\": \"HistoryLinearElasticPlaneStrain2DLaw\""
            }
        } elseif {[lindex [lindex $Groups $i] 3] eq "LinearElasticPlaneStressSolid2DLaw"} {
            if { ([GiD_AccessValue get gendata Initial_Stresses] eq false) || (([GiD_AccessValue get gendata Initial_Stresses] eq true) && ([GiD_AccessValue get gendata Mode] eq "save")) } {
                puts $FileVar "                \"name\": \"[lindex [lindex $Groups $i] 3]\""
            } else {
                puts $FileVar "                \"name\": \"HistoryLinearElasticPlaneStress2DLaw\""
            }
        } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJuDamage3DLaw"} {
            if {[GiD_AccessValue get gendata Non-local_Damage] eq true} {
                puts $FileVar "                \"name\": \"SimoJuNonlocalDamage3DLaw\""
            } else {
                puts $FileVar "                \"name\": \"SimoJuLocalDamage3DLaw\""
            }
        } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJuDamagePlaneStrain2DLaw"} {
            if {[GiD_AccessValue get gendata Non-local_Damage] eq true} {
                puts $FileVar "                \"name\": \"SimoJuNonlocalDamagePlaneStrain2DLaw\""
            } else {
                puts $FileVar "                \"name\": \"SimoJuLocalDamagePlaneStrain2DLaw\""
            }
        } elseif {[lindex [lindex $Groups $i] 3] eq "SimoJuDamagePlaneStress2DLaw"} {
            if {[GiD_AccessValue get gendata Non-local_Damage] eq true} {
                puts $FileVar "                \"name\": \"SimoJuNonlocalDamagePlaneStress2DLaw\""
            } else {
                puts $FileVar "                \"name\": \"SimoJuLocalDamagePlaneStress2DLaw\""
            }
        } elseif {[lindex [lindex $Groups $i] 3] eq "ModifiedMisesDamage3DLaw"} {
            puts $FileVar "                \"name\": \"ModifiedMisesNonlocalDamage3DLaw\""
        } elseif {[lindex [lindex $Groups $i] 3] eq "ModifiedMisesDamagePlaneStrain2DLaw"} {
            puts $FileVar "                \"name\": \"ModifiedMisesNonlocalDamagePlaneStrain2DLaw\""
        } elseif {[lindex [lindex $Groups $i] 3] eq "ModifiedMisesDamagePlaneStress2DLaw"} {
            puts $FileVar "                \"name\": \"ModifiedMisesNonlocalDamagePlaneStress2DLaw\""
        }
        puts $FileVar "            \},"
        puts $FileVar "            \"Variables\": \{"
        puts $FileVar "                \"YOUNG_MODULUS\": [lindex [lindex $Groups $i] 4],"
        puts $FileVar "                \"POISSON_RATIO\": [lindex [lindex $Groups $i] 5],"
        puts $FileVar "                \"DENSITY_SOLID\": [lindex [lindex $Groups $i] 6],"
        puts $FileVar "                \"DENSITY_WATER\": [lindex [lindex $Groups $i] 7],"
        puts $FileVar "                \"POROSITY\": [lindex [lindex $Groups $i] 8],"
        puts $FileVar "                \"BULK_MODULUS_SOLID\": [lindex [lindex $Groups $i] 9],"
        puts $FileVar "                \"BULK_MODULUS_FLUID\": [lindex [lindex $Groups $i] 10],"
        puts $FileVar "                \"PERMEABILITY_XX\": [lindex [lindex $Groups $i] 11],"
        puts $FileVar "                \"PERMEABILITY_YY\": [lindex [lindex $Groups $i] 12],"
        puts $FileVar "                \"PERMEABILITY_ZZ\": [lindex [lindex $Groups $i] 13],"
        puts $FileVar "                \"PERMEABILITY_XY\": [lindex [lindex $Groups $i] 14],"
        puts $FileVar "                \"PERMEABILITY_YZ\": [lindex [lindex $Groups $i] 15],"
        puts $FileVar "                \"PERMEABILITY_ZX\": [lindex [lindex $Groups $i] 16],"
        puts $FileVar "                \"DYNAMIC_VISCOSITY\": [lindex [lindex $Groups $i] 17],"
        puts $FileVar "                \"THICKNESS\": [lindex [lindex $Groups $i] 18],"
        puts $FileVar "                \"DAMAGE_THRESHOLD\": [lindex [lindex $Groups $i] 19],"
        puts $FileVar "                \"STRENGTH_RATIO\": [lindex [lindex $Groups $i] 20],"
        puts $FileVar "                \"FRACTURE_ENERGY\": [lindex [lindex $Groups $i] 21],"
        puts $FileVar "                \"RESIDUAL_STRENGTH\": [lindex [lindex $Groups $i] 22],"
        puts $FileVar "                \"SOFTENING_SLOPE\": [lindex [lindex $Groups $i] 23]"
        puts $FileVar "            \},"
        puts $FileVar "            \"Tables\": \{\}"
        puts $FileVar "        \}"
        if {$MyPropertyId < $NumGroups} {
            puts $FileVar "    \},\{"
        } else {
            puts $FileVar "    \}\]"
        }
    }
    # Interface_Part
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr MyPropertyId
        puts $FileVar "        \"model_part_name\": \"PorousModelPart.[lindex [lindex $Groups $i] 1]\","
        puts $FileVar "        \"properties_id\": $MyPropertyId,"
        puts $FileVar "        \"Material\": \{"
        puts $FileVar "            \"constitutive_law\": \{"
        if {[lindex [lindex $Groups $i] 4] eq "BilinearCohesive3DLaw" || [lindex [lindex $Groups $i] 4] eq "ExponentialCohesive3DLaw"} {
            puts $FileVar "                \"name\": \"[lindex [lindex $Groups $i] 4]\""
        } elseif {[lindex [lindex $Groups $i] 4] eq "BilinearCohesivePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 4] eq "BilinearCohesivePlaneStress2DLaw"} {
            puts $FileVar "                \"name\": \"BilinearCohesive2DLaw\""
        } elseif {[lindex [lindex $Groups $i] 4] eq "ExponentialCohesivePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 4] eq "ExponentialCohesivePlaneStress2DLaw"} {
            puts $FileVar "                \"name\": \"ExponentialCohesive2DLaw\""
        }
        puts $FileVar "            \},"
        puts $FileVar "            \"Variables\": \{"
        puts $FileVar "                \"YOUNG_MODULUS\": [lindex [lindex $Groups $i] 5],"
        puts $FileVar "                \"POISSON_RATIO\": [lindex [lindex $Groups $i] 6],"
        puts $FileVar "                \"DENSITY_SOLID\": [lindex [lindex $Groups $i] 7],"
        puts $FileVar "                \"DENSITY_WATER\": [lindex [lindex $Groups $i] 8],"
        puts $FileVar "                \"POROSITY\": [lindex [lindex $Groups $i] 9],"
        puts $FileVar "                \"BULK_MODULUS_SOLID\": [lindex [lindex $Groups $i] 10],"
        puts $FileVar "                \"BULK_MODULUS_FLUID\": [lindex [lindex $Groups $i] 11],"
        puts $FileVar "                \"TRANSVERSAL_PERMEABILITY\": [lindex [lindex $Groups $i] 12],"
        puts $FileVar "                \"DYNAMIC_VISCOSITY\": [lindex [lindex $Groups $i] 13],"
        puts $FileVar "                \"THICKNESS\": [lindex [lindex $Groups $i] 14],"
        puts $FileVar "                \"DAMAGE_THRESHOLD\": [lindex [lindex $Groups $i] 15],"
        puts $FileVar "                \"MINIMUM_JOINT_WIDTH\": [lindex [lindex $Groups $i] 16],"
        puts $FileVar "                \"CRITICAL_DISPLACEMENT\": [lindex [lindex $Groups $i] 17],"
        puts $FileVar "                \"YIELD_STRESS\": [lindex [lindex $Groups $i] 18],"
        puts $FileVar "                \"FRICTION_COEFFICIENT\": [lindex [lindex $Groups $i] 19],"
        puts $FileVar "                \"FRACTURE_ENERGY\": [lindex [lindex $Groups $i] 22],"
        puts $FileVar "                \"SHEAR_FRACTURE_ENERGY\": [lindex [lindex $Groups $i] 23],"
        puts $FileVar "                \"STATE_VARIABLE\": [lindex [lindex $Groups $i] 24]"
        puts $FileVar "            \},"
        puts $FileVar "            \"Tables\": \{\}"
        puts $FileVar "        \}"
        if {$MyPropertyId < $NumGroups} {
            puts $FileVar "    \},\{"
        } else {
            puts $FileVar "    \}\]"
        }
    }

    puts $FileVar "\}"

    close $FileVar
}
