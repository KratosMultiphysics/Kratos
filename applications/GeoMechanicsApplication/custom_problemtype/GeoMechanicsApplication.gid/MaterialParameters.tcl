proc WriteMaterialParameters {basename dir problemtypedir TableDict} {

    ## Start ProjectParameters.json file
    set filename [file join $dir MaterialParameters.json]
    set FileVar [open $filename w]
    puts $FileVar "\{"

    ## problem_data
    set PropertyId 0
    set PropertyDict [dict create]

    puts $FileVar "   \"properties\": \[\{"

    set SolutionType [GiD_AccessValue get gendata Solution_Type]
    set Dim [GiD_AccessValue get gendata Domain_Size]

    # Soil_two_phase part
    set Groups [GiD_Info conditions Soil_two_phase groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {$SolutionType eq "K0-Procedure"} {
            if {$Dim eq 3} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"LinearElasticK03DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  0.0,"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  0.0,"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  0.0,"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  1.0e-3,"
                if {[lindex [lindex $Groups $i] 28] eq "Y"} {
                  set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 28] eq "Z"} {
                  set PutStrings 2
                } else {
                  set PutStrings 0
                }
                puts $FileVar "              \"K0_MAIN_DIRECTION\"        :  $PutStrings,"
                puts $FileVar "              \"K0_VALUE_XX\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"K0_VALUE_YY\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31],"
            } else {
                # 2D soil elements
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"LinearElasticPlaneStrainK02DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  0.0,"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  1.0e-3,"
                puts $FileVar "              \"THICKNESS\"                :  1.0,"
                if {[lindex [lindex $Groups $i] 28] eq "Y"} {
                  set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 28] eq "Z"} {
                  set PutStrings 2
                } else {
                  set PutStrings 0
                }
                puts $FileVar "              \"K0_MAIN_DIRECTION\"        :  $PutStrings,"
                puts $FileVar "              \"K0_VALUE_XX\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"K0_VALUE_YY\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31],"
            }
        } else {
            if {[lindex [lindex $Groups $i] 3] eq "SmallStrainIsotropicPlasticity3D"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"SmallStrainIsotropicPlasticity3D[lindex [lindex $Groups $i] 4][lindex [lindex $Groups $i] 5]\", "
                puts $FileVar "              \"yield_surface\"    :  \"[lindex [lindex $Groups $i] 4]\", "
                puts $FileVar "              \"plastic_potential\":  \"[lindex [lindex $Groups $i] 5]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
                puts $FileVar "              \"FRICTION_ANGLE\"           :  [lindex [lindex $Groups $i] 20],"
                puts $FileVar "              \"DILATANCY_ANGLE\"          :  [lindex [lindex $Groups $i] 21],"
                puts $FileVar "              \"HARDENING_CURVE\"          :  [lindex [lindex $Groups $i] 22],"
                puts $FileVar "              \"FRACTURE_ENERGY\"          :  [lindex [lindex $Groups $i] 23],"
                puts $FileVar "              \"YIELD_STRESS_COMPRESSION\" :  [lindex [lindex $Groups $i] 24],"
                puts $FileVar "              \"YIELD_STRESS_TENSION\"     :  [lindex [lindex $Groups $i] 25],"
                puts $FileVar "              \"MAXIMUM_STRESS_POSITION\"  :  [lindex [lindex $Groups $i] 26],"
                puts $FileVar "              \"MAXIMUM_STRESS\"           :  [lindex [lindex $Groups $i] 27],"
            } elseif {[lindex [lindex $Groups $i] 3] eq "LinearElastic3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
            } elseif {[lindex [lindex $Groups $i] 3] eq "GeoLinearElasticPlaneStrain2DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
                puts $FileVar "              \"THICKNESS\"                :  1.0,"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUDSM3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 32]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 34],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUDSM2DPlaneStrainLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 32]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 34],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUMAT3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 35]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 37],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 38],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 89]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+89}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUMAT2DPlaneStrainLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 35]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 37],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 38],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 89]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+89}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            }
        }
        if {[lindex [lindex $Groups $i] 147] eq true} {
            puts $FileVar "              \"BIOT_COEFFICIENT\"         :  [lindex [lindex $Groups $i] 148],"
        }
        puts $FileVar "              \"RETENTION_LAW\"                    : \"[lindex [lindex $Groups $i] 140]\","
        puts $FileVar "              \"SATURATED_SATURATION\"             :  [lindex [lindex $Groups $i] 141],"
        puts $FileVar "              \"RESIDUAL_SATURATION\"              :  [lindex [lindex $Groups $i] 142],"
        puts $FileVar "              \"VAN_GENUCHTEN_AIR_ENTRY_PRESSURE\" :  [lindex [lindex $Groups $i] 143],"
        puts $FileVar "              \"VAN_GENUCHTEN_GN\"                 :  [lindex [lindex $Groups $i] 144],"
        puts $FileVar "              \"VAN_GENUCHTEN_GL\"                 :  [lindex [lindex $Groups $i] 145],"
        puts $FileVar "              \"MINIMUM_RELATIVE_PERMEABILITY\"    :  [lindex [lindex $Groups $i] 146]"

        puts $FileVar "         \},"
        puts $FileVar "         \"Tables\": \{\}"
        puts $FileVar "      \}"
    }

    # Soil_drained part
    set Groups [GiD_Info conditions Soil_drained groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {$SolutionType eq "K0-Procedure"} {
            if {$Dim eq 3} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"LinearElasticK03DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  0.0,"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  0.0,"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  0.0,"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  1.0e-3,"
                if {[lindex [lindex $Groups $i] 28] eq "Y"} {
                  set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 28] eq "Z"} {
                  set PutStrings 2
                } else {
                  set PutStrings 0
                }
                puts $FileVar "              \"K0_MAIN_DIRECTION\"        :  $PutStrings,"
                puts $FileVar "              \"K0_VALUE_XX\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"K0_VALUE_YY\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31],"

            } else {
                # 2D elements
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"LinearElasticPlaneStrainK02DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  1.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  0.0,"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  1.0e-3,"
                puts $FileVar "              \"THICKNESS\"                :  1.0,"
                if {[lindex [lindex $Groups $i] 28] eq "Y"} {
                  set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 28] eq "Z"} {
                  set PutStrings 2
                } else {
                  set PutStrings 0
                }
                puts $FileVar "              \"K0_MAIN_DIRECTION\"        :  $PutStrings,"
                puts $FileVar "              \"K0_VALUE_XX\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"K0_VALUE_YY\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31],"
            }
        } else {
            if {[lindex [lindex $Groups $i] 3] eq "SmallStrainIsotropicPlasticity3D"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"SmallStrainIsotropicPlasticity3D[lindex [lindex $Groups $i] 4][lindex [lindex $Groups $i] 5]\", "
                puts $FileVar "              \"yield_surface\"    :  \"[lindex [lindex $Groups $i] 4]\", "
                puts $FileVar "              \"plastic_potential\":  \"[lindex [lindex $Groups $i] 5]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
                puts $FileVar "              \"FRICTION_ANGLE\"           :  [lindex [lindex $Groups $i] 20],"
                puts $FileVar "              \"DILATANCY_ANGLE\"          :  [lindex [lindex $Groups $i] 21],"
                puts $FileVar "              \"HARDENING_CURVE\"          :  [lindex [lindex $Groups $i] 22],"
                puts $FileVar "              \"FRACTURE_ENERGY\"          :  [lindex [lindex $Groups $i] 23],"
                puts $FileVar "              \"YIELD_STRESS_COMPRESSION\" :  [lindex [lindex $Groups $i] 24],"
                puts $FileVar "              \"YIELD_STRESS_TENSION\"     :  [lindex [lindex $Groups $i] 25],"
                puts $FileVar "              \"MAXIMUM_STRESS_POSITION\"  :  [lindex [lindex $Groups $i] 26],"
                puts $FileVar "              \"MAXIMUM_STRESS\"           :  [lindex [lindex $Groups $i] 27],"
            } elseif {[lindex [lindex $Groups $i] 3] eq "LinearElastic3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
            } elseif {[lindex [lindex $Groups $i] 3] eq "GeoLinearElasticPlaneStrain2DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
                puts $FileVar "              \"THICKNESS\"                :  1.0,"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUDSM3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 32]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 34],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUDSM2DPlaneStrainLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 32]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 34],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUMAT3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 35]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 37],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 38],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 89]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+89}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUMAT2DPlaneStrainLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 35]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 37],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 38],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 89]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+89}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            }
        }
        if {[lindex [lindex $Groups $i] 147] eq true} {
            puts $FileVar "              \"BIOT_COEFFICIENT\"         :  [lindex [lindex $Groups $i] 148],"
        }
        puts $FileVar "              \"RETENTION_LAW\"                    : \"[lindex [lindex $Groups $i] 140]\","
        puts $FileVar "              \"SATURATED_SATURATION\"             :  [lindex [lindex $Groups $i] 141],"
        puts $FileVar "              \"RESIDUAL_SATURATION\"              :  [lindex [lindex $Groups $i] 142],"
        puts $FileVar "              \"VAN_GENUCHTEN_AIR_ENTRY_PRESSURE\" :  [lindex [lindex $Groups $i] 143],"
        puts $FileVar "              \"VAN_GENUCHTEN_GN\"                 :  [lindex [lindex $Groups $i] 144],"
        puts $FileVar "              \"VAN_GENUCHTEN_GL\"                 :  [lindex [lindex $Groups $i] 145],"
        puts $FileVar "              \"MINIMUM_RELATIVE_PERMEABILITY\"    :  [lindex [lindex $Groups $i] 146]"

        puts $FileVar "         \},"
        puts $FileVar "         \"Tables\": \{\}"
        puts $FileVar "      \}"
    }

    # Soil_undrained part
    set Groups [GiD_Info conditions Soil_undrained groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {$SolutionType eq "K0-Procedure"} {
            if {$Dim eq 3} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"LinearElasticK03DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  2.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  4.5e-30,"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  4.5e-30,"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  4.5e-30,"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  0.0,"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  0.0,"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  0.0,"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  1.0e-3,"
                if {[lindex [lindex $Groups $i] 28] eq "Y"} {
                  set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 28] eq "Z"} {
                  set PutStrings 2
                } else {
                  set PutStrings 0
                }
                puts $FileVar "              \"K0_MAIN_DIRECTION\"        :  $PutStrings,"
                puts $FileVar "              \"K0_VALUE_XX\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"K0_VALUE_YY\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31],"
            } else {
                #2D elements
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"LinearElasticPlaneStrainK02DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  2.0e-30,"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  4.5e-30,"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  4.5e-30,"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  0.0,"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  1.0e-3,"
                puts $FileVar "              \"THICKNESS\"                :  1.0,"
                if {[lindex [lindex $Groups $i] 28] eq "Y"} {
                  set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 28] eq "Z"} {
                  set PutStrings 2
                } else {
                  set PutStrings 0
                }
                puts $FileVar "              \"K0_MAIN_DIRECTION\"        :  $PutStrings,"
                puts $FileVar "              \"K0_VALUE_XX\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"K0_VALUE_YY\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31],"
            }
        } else {    
            if {[lindex [lindex $Groups $i] 3] eq "SmallStrainIsotropicPlasticity3D"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"SmallStrainIsotropicPlasticity3D[lindex [lindex $Groups $i] 4][lindex [lindex $Groups $i] 5]\", "
                puts $FileVar "              \"yield_surface\"    :  \"[lindex [lindex $Groups $i] 4]\", "
                puts $FileVar "              \"plastic_potential\":  \"[lindex [lindex $Groups $i] 5]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
                puts $FileVar "              \"FRICTION_ANGLE\"           :  [lindex [lindex $Groups $i] 20],"
                puts $FileVar "              \"DILATANCY_ANGLE\"          :  [lindex [lindex $Groups $i] 21],"
                puts $FileVar "              \"HARDENING_CURVE\"          :  [lindex [lindex $Groups $i] 22],"
                puts $FileVar "              \"FRACTURE_ENERGY\"          :  [lindex [lindex $Groups $i] 23],"
                puts $FileVar "              \"YIELD_STRESS_COMPRESSION\" :  [lindex [lindex $Groups $i] 24],"
                puts $FileVar "              \"YIELD_STRESS_TENSION\"     :  [lindex [lindex $Groups $i] 25],"
                puts $FileVar "              \"MAXIMUM_STRESS_POSITION\"  :  [lindex [lindex $Groups $i] 26],"
                puts $FileVar "              \"MAXIMUM_STRESS\"           :  [lindex [lindex $Groups $i] 27],"
            } elseif {[lindex [lindex $Groups $i] 3] eq "LinearElastic3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
            } elseif {[lindex [lindex $Groups $i] 3] eq "GeoLinearElasticPlaneStrain2DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"
                puts $FileVar "              \"THICKNESS\"                :  1.0,"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUDSM3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 32]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 34],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUDSM2DPlaneStrainLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 32]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 34],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUMAT3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 35]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 37],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 38],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 89]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+89}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUMAT2DPlaneStrainLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 35]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 37],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 38],"

                set nParameters [expr {[lindex [lindex $Groups $i] 38]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+38}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 89]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+89}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

            }
        }
        if {[lindex [lindex $Groups $i] 147] eq true} {
            puts $FileVar "              \"BIOT_COEFFICIENT\"         :  [lindex [lindex $Groups $i] 148],"
        }
        puts $FileVar "              \"RETENTION_LAW\"                    : \"[lindex [lindex $Groups $i] 140]\","
        puts $FileVar "              \"SATURATED_SATURATION\"             :  [lindex [lindex $Groups $i] 141],"
        puts $FileVar "              \"RESIDUAL_SATURATION\"              :  [lindex [lindex $Groups $i] 142],"
        puts $FileVar "              \"VAN_GENUCHTEN_AIR_ENTRY_PRESSURE\" :  [lindex [lindex $Groups $i] 143],"
        puts $FileVar "              \"VAN_GENUCHTEN_GN\"                 :  [lindex [lindex $Groups $i] 144],"
        puts $FileVar "              \"VAN_GENUCHTEN_GL\"                 :  [lindex [lindex $Groups $i] 145],"
        puts $FileVar "              \"MINIMUM_RELATIVE_PERMEABILITY\"    :  [lindex [lindex $Groups $i] 146]"

        puts $FileVar "         \},"
        puts $FileVar "         \"Tables\": \{\}"
        puts $FileVar "      \}"
    }

    # Non_porous part
    set Groups [GiD_Info conditions Non_porous groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {$SolutionType eq "K0-Procedure"} {
            if {$Dim eq 3} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"LinearElasticK03DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"    :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"    :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY\"          :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"         :  [lindex [lindex $Groups $i] 9],"
                if {[lindex [lindex $Groups $i] 18] eq "Y"} {
                  set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 18] eq "Z"} {
                  set PutStrings 2
                } else {
                  set PutStrings 0
                }
                puts $FileVar "              \"K0_MAIN_DIRECTION\":  $PutStrings,"
                puts $FileVar "              \"K0_VALUE_XX\"      :  [lindex [lindex $Groups $i] 19],"
                puts $FileVar "              \"K0_VALUE_YY\"      :  [lindex [lindex $Groups $i] 20],"
                puts $FileVar "              \"K0_VALUE_ZZ\"      :  [lindex [lindex $Groups $i] 21]"
                
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            } else {
                # 2D soil elements
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"LinearElasticPlaneStrainK02DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"    :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"    :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY\"          :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"         :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"THICKNESS\"        :  1.0,"
                if {[lindex [lindex $Groups $i] 18] eq "Y"} {
                  set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 18] eq "Z"} {
                  set PutStrings 2
                } else {
                  set PutStrings 0
                }
                puts $FileVar "              \"K0_MAIN_DIRECTION\":  $PutStrings,"
                puts $FileVar "              \"K0_VALUE_XX\"      :  [lindex [lindex $Groups $i] 19],"
                puts $FileVar "              \"K0_VALUE_YY\"      :  [lindex [lindex $Groups $i] 20],"
                puts $FileVar "              \"K0_VALUE_ZZ\"      :  [lindex [lindex $Groups $i] 21]"
                
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            }
        } else {
    
            if {[lindex [lindex $Groups $i] 3] eq "SmallStrainIsotropicPlasticity3D"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"SmallStrainIsotropicPlasticity3D[lindex [lindex $Groups $i] 4][lindex [lindex $Groups $i] 5]\", "
                puts $FileVar "              \"yield_surface\"    :  \"[lindex [lindex $Groups $i] 4]\", "
                puts $FileVar "              \"plastic_potential\":  \"[lindex [lindex $Groups $i] 5]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY\"                  :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"FRICTION_ANGLE\"           :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"DILATANCY_ANGLE\"          :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"HARDENING_CURVE\"          :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"FRACTURE_ENERGY\"          :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"YIELD_STRESS_COMPRESSION\" :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"YIELD_STRESS_TENSION\"     :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"MAXIMUM_STRESS_POSITION\"  :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"MAXIMUM_STRESS\"           :  [lindex [lindex $Groups $i] 17]"
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            } elseif {[lindex [lindex $Groups $i] 3] eq "LinearElastic3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"    :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"    :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY\"          :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"         :  [lindex [lindex $Groups $i] 9]"
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"

            } elseif {[lindex [lindex $Groups $i] 3] eq "GeoLinearElasticPlaneStrain2DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"    :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"POISSON_RATIO\"    :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY\"          :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"         :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"THICKNESS\"        :  1.0"
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUDSM3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 22]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 23],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 24],"

                set nParameters [expr {[lindex [lindex $Groups $i] 28]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+28}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                puts $FileVar "              \"DENSITY\"                  :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9]"

                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            
            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUDSM2DPlaneStrainLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 22]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 23],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 24],"

                set nParameters [expr {[lindex [lindex $Groups $i] 28]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+28}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                puts $FileVar "              \"DENSITY\"                  :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9]"

                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUMAT3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 25]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 27],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 28],"

                set nParameters [expr {[lindex [lindex $Groups $i] 28]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+28}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 79]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+79}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                puts $FileVar "              \"DENSITY\"                  :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9]"

                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            
            } elseif {[lindex [lindex $Groups $i] 3] eq "SmallStrainUMAT2DPlaneStrainLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"DENSITY\"                  :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 25]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 27],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 28],"

                set nParameters [expr {[lindex [lindex $Groups $i] 28]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+28}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 79]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+79}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            }
        }
    }


    # Soil_Groundwater_Flow part
    set Groups [GiD_Info conditions Soil_Groundwater_Flow groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 3] eq "LinearElastic3DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            if {$PropertyId > 1} {
                puts $FileVar "   \},\{"
            }
            puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $FileVar "      \"properties_id\":           $PropertyId,"
            puts $FileVar "      \"Material\": \{"
            puts $FileVar "          \"constitutive_law\": \{"
            puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
            puts $FileVar "          \},"
            puts $FileVar "          \"Variables\": \{"
            puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 4],"
            puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 5],"
            puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 6],"
            puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 7],"
            puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 8],"
            puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 9],"
            puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 10],"
            puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 11],"
            puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 12],"
            puts $FileVar "              \"PERMEABILITY_ZZ\"          :  [lindex [lindex $Groups $i] 13],"
            puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 14],"
            puts $FileVar "              \"PERMEABILITY_YZ\"          :  [lindex [lindex $Groups $i] 15],"
            puts $FileVar "              \"PERMEABILITY_ZX\"          :  [lindex [lindex $Groups $i] 16],"
            puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 17],"
        } elseif {[lindex [lindex $Groups $i] 3] eq "GeoLinearElasticPlaneStrain2DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            if {$PropertyId > 1} {
                puts $FileVar "   \},\{"
            }
            puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $FileVar "      \"properties_id\":           $PropertyId,"
            puts $FileVar "      \"Material\": \{"
            puts $FileVar "          \"constitutive_law\": \{"
            puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 3]\" "
            puts $FileVar "          \},"
            puts $FileVar "          \"Variables\": \{"
            puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 4],"
            puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 5],"
            puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 6],"
            puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 7],"
            puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 8],"
            puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 9],"
            puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 10],"
            puts $FileVar "              \"PERMEABILITY_XX\"          :  [lindex [lindex $Groups $i] 11],"
            puts $FileVar "              \"PERMEABILITY_YY\"          :  [lindex [lindex $Groups $i] 12],"
            puts $FileVar "              \"PERMEABILITY_XY\"          :  [lindex [lindex $Groups $i] 14],"
            puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 17],"
            puts $FileVar "              \"THICKNESS\"                :  1.0,"
        }
        if {[lindex [lindex $Groups $i] 25] eq true} {
            puts $FileVar "              \"BIOT_COEFFICIENT\"         :  [lindex [lindex $Groups $i] 26],"
        }
        puts $FileVar "              \"RETENTION_LAW\"                    : \"[lindex [lindex $Groups $i] 18]\","
        puts $FileVar "              \"SATURATED_SATURATION\"             :  [lindex [lindex $Groups $i] 19],"
        puts $FileVar "              \"RESIDUAL_SATURATION\"              :  [lindex [lindex $Groups $i] 20],"
        puts $FileVar "              \"VAN_GENUCHTEN_AIR_ENTRY_PRESSURE\" :  [lindex [lindex $Groups $i] 21],"
        puts $FileVar "              \"VAN_GENUCHTEN_GN\"                 :  [lindex [lindex $Groups $i] 22],"
        puts $FileVar "              \"VAN_GENUCHTEN_GL\"                 :  [lindex [lindex $Groups $i] 23],"
        puts $FileVar "              \"MINIMUM_RELATIVE_PERMEABILITY\"    :  [lindex [lindex $Groups $i] 24]"

        puts $FileVar "         \},"
        puts $FileVar "         \"Tables\": \{\}"
        puts $FileVar "      \}"
    }


    # Beam part
    set Groups [GiD_Info conditions Beam groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
       incr PropertyId
       dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
       if {$PropertyId > 1} {
           puts $FileVar "   \},\{"
       }
       puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
       puts $FileVar "      \"properties_id\":           $PropertyId,"
       puts $FileVar "      \"Material\": \{"
       puts $FileVar "          \"constitutive_law\": \{"
       puts $FileVar "              \"name\"             :  \"KratosMultiphysics.StructuralMechanicsApplication.BeamConstitutiveLaw\" "
       puts $FileVar "          \},"
       puts $FileVar "          \"Variables\": \{"
       puts $FileVar "              \"YOUNG_MODULUS\"     :  [lindex [lindex $Groups $i] 5],"
       puts $FileVar "              \"POISSON_RATIO\"     :  [lindex [lindex $Groups $i] 6],"
       puts $FileVar "              \"DENSITY\"           :  [lindex [lindex $Groups $i] 7],"
       puts $FileVar "              \"CROSS_AREA\"        :  [lindex [lindex $Groups $i] 8],"
       puts $FileVar "              \"I22\"               :  [lindex [lindex $Groups $i] 9],"
       puts $FileVar "              \"I33\"               :  [lindex [lindex $Groups $i] 10],"
       puts $FileVar "              \"TORSIONAL_INERTIA\" :  [lindex [lindex $Groups $i] 11]"
       puts $FileVar "         \},"
       puts $FileVar "         \"Tables\": \{\}"
       puts $FileVar "      \}"
    }

    # Shell_thin_corotational part
    set Groups [GiD_Info conditions Shell_thin_corotational groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
       incr PropertyId
       dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
       if {$PropertyId > 1} {
           puts $FileVar "   \},\{"
       }
       puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
       puts $FileVar "      \"properties_id\":           $PropertyId,"
       puts $FileVar "      \"Material\": \{"
       puts $FileVar "          \"constitutive_law\": \{"
       puts $FileVar "              \"name\"             :  \"KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw\" "
       puts $FileVar "          \},"
       puts $FileVar "          \"Variables\": \{"
       puts $FileVar "              \"YOUNG_MODULUS\"    :  [lindex [lindex $Groups $i] 4],"
       puts $FileVar "              \"POISSON_RATIO\"    :  [lindex [lindex $Groups $i] 5],"
       puts $FileVar "              \"DENSITY\"          :  [lindex [lindex $Groups $i] 6],"
       puts $FileVar "              \"THICKNESS\"        :  [lindex [lindex $Groups $i] 7]"
       puts $FileVar "         \},"
       puts $FileVar "         \"Tables\": \{\}"
       puts $FileVar "      \}"
    }

    # Shell_thick_corotational part
    set Groups [GiD_Info conditions Shell_thick_corotational groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
       incr PropertyId
       dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
       if {$PropertyId > 1} {
           puts $FileVar "   \},\{"
       }
       puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
       puts $FileVar "      \"properties_id\":           $PropertyId,"
       puts $FileVar "      \"Material\": \{"
       puts $FileVar "          \"constitutive_law\": \{"
       puts $FileVar "              \"name\"             :  \"KratosMultiphysics.StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw\" "
       puts $FileVar "          \},"
       puts $FileVar "          \"Variables\": \{"
       puts $FileVar "              \"YOUNG_MODULUS\"    :  [lindex [lindex $Groups $i] 4],"
       puts $FileVar "              \"POISSON_RATIO\"    :  [lindex [lindex $Groups $i] 5],"
       puts $FileVar "              \"DENSITY\"          :  [lindex [lindex $Groups $i] 6],"
       puts $FileVar "              \"THICKNESS\"        :  [lindex [lindex $Groups $i] 7]"
       puts $FileVar "         \},"
       puts $FileVar "         \"Tables\": \{\}"
       puts $FileVar "      \}"
    }

    # Truss part
    set Groups [GiD_Info conditions Truss groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
       incr PropertyId
       dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
       if {$PropertyId > 1} {
           puts $FileVar "   \},\{"
       }
       puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
       puts $FileVar "      \"properties_id\":           $PropertyId,"
       puts $FileVar "      \"Material\": \{"
       puts $FileVar "          \"constitutive_law\": \{"
       puts $FileVar "              \"name\"             :  \"KratosMultiphysics.StructuralMechanicsApplication.TrussConstitutiveLaw\" "
       puts $FileVar "          \},"
       puts $FileVar "          \"Variables\": \{"
       puts $FileVar "              \"YOUNG_MODULUS\"       :  [lindex [lindex $Groups $i] 4],"
       puts $FileVar "              \"POISSON_RATIO\"       :  [lindex [lindex $Groups $i] 5],"
       puts $FileVar "              \"DENSITY\"             :  [lindex [lindex $Groups $i] 6],"
       puts $FileVar "              \"CROSS_AREA\"          :  [lindex [lindex $Groups $i] 7],"
       puts $FileVar "              \"TRUSS_PRESTRESS_PK2\" :  [lindex [lindex $Groups $i] 8]"
       puts $FileVar "         \},"
       puts $FileVar "         \"Tables\": \{\}"
       puts $FileVar "      \}"
    }

    # Anchor part
    set Groups [GiD_Info conditions Anchor groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
       incr PropertyId
       dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
       if {$PropertyId > 1} {
           puts $FileVar "   \},\{"
       }
       puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
       puts $FileVar "      \"properties_id\":           $PropertyId,"
       puts $FileVar "      \"Material\": \{"
       puts $FileVar "          \"constitutive_law\": \{"
       puts $FileVar "              \"name\"             :  \"KratosMultiphysics.StructuralMechanicsApplication.TrussConstitutiveLaw\" "
       puts $FileVar "          \},"
       puts $FileVar "          \"Variables\": \{"
       puts $FileVar "              \"YOUNG_MODULUS\"       :  [lindex [lindex $Groups $i] 4],"
       puts $FileVar "              \"POISSON_RATIO\"       :  [lindex [lindex $Groups $i] 5],"
       puts $FileVar "              \"DENSITY\"             :  [lindex [lindex $Groups $i] 6],"
       puts $FileVar "              \"CROSS_AREA\"          :  [lindex [lindex $Groups $i] 7],"
       puts $FileVar "              \"TRUSS_PRESTRESS_PK2\" :  [lindex [lindex $Groups $i] 8]"
       puts $FileVar "         \},"
       puts $FileVar "         \"Tables\": \{\}"
       puts $FileVar "      \}"
    }

    # Interface drained part
    set Groups [GiD_Info conditions Interface_drained groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 4] eq "BilinearCohesive3DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            if {$PropertyId > 1} {
                puts $FileVar "   \},\{"
            }
            puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $FileVar "      \"properties_id\":           $PropertyId,"
            puts $FileVar "      \"Material\": \{"
            puts $FileVar "          \"constitutive_law\": \{"
            puts $FileVar "              \"name\"             :  \"BilinearCohesive3DLaw\" "
            puts $FileVar "          \},"
            puts $FileVar "          \"Variables\": \{"
            puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
            puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 5],"
            puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 6],"
            puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
            puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
            puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
            puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
            puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
            puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
            puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
            puts $FileVar "              \"DAMAGE_THRESHOLD\"         :  [lindex [lindex $Groups $i] 15],"
            puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16],"
            puts $FileVar "              \"CRITICAL_DISPLACEMENT\"    :  [lindex [lindex $Groups $i] 17],"
            puts $FileVar "              \"CONSIDER_GAP_CLOSURE\"     :  [lindex [lindex $Groups $i] 135],"
            puts $FileVar "              \"YIELD_STRESS\"             :  [lindex [lindex $Groups $i] 18],"
            puts $FileVar "              \"FRICTION_COEFFICIENT\"     :  [lindex [lindex $Groups $i] 19]"
            puts $FileVar "         \},"
            puts $FileVar "         \"Tables\": \{\}"
            puts $FileVar "      \}"
        } elseif {[lindex [lindex $Groups $i] 4] eq "BilinearCohesivePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 4] eq "BilinearCohesivePlaneStress2DLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            if {$PropertyId > 1} {
                puts $FileVar "   \},\{"
            }
            puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $FileVar "      \"properties_id\":           $PropertyId,"
            puts $FileVar "      \"Material\": \{"
            puts $FileVar "          \"constitutive_law\": \{"
            puts $FileVar "              \"name\"             :  \"BilinearCohesive2DLaw\" "
            puts $FileVar "          \},"
            puts $FileVar "          \"Variables\": \{"
            puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
            puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 5],"
            puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 6],"
            puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
            puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
            puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
            puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
            puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
            puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
            puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
            puts $FileVar "              \"THICKNESS\"                :  [lindex [lindex $Groups $i] 14],"
            puts $FileVar "              \"DAMAGE_THRESHOLD\"         :  [lindex [lindex $Groups $i] 15],"
            puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16],"
            puts $FileVar "              \"CRITICAL_DISPLACEMENT\"    :  [lindex [lindex $Groups $i] 17],"
            puts $FileVar "              \"CONSIDER_GAP_CLOSURE\"     :  [lindex [lindex $Groups $i] 135],"
            puts $FileVar "              \"YIELD_STRESS\"             :  [lindex [lindex $Groups $i] 18],"
            puts $FileVar "              \"FRICTION_COEFFICIENT\"     :  [lindex [lindex $Groups $i] 19]"
            puts $FileVar "         \},"
            puts $FileVar "         \"Tables\": \{\}"
            puts $FileVar "      \}"
        } elseif {[lindex [lindex $Groups $i] 4] eq "SmallStrainUDSM2DInterfaceLaw" || [lindex [lindex $Groups $i] 4] eq "SmallStrainUDSM3DInterfaceLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            if {$PropertyId > 1} {
                puts $FileVar "   \},\{"
            }
            puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $FileVar "      \"properties_id\":           $PropertyId,"
            puts $FileVar "      \"Material\": \{"
            puts $FileVar "          \"constitutive_law\": \{"
            puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 4]\" "
            puts $FileVar "          \},"
            puts $FileVar "          \"Variables\": \{"
            puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
            puts $FileVar "              \"CONSIDER_GAP_CLOSURE\"     :  [lindex [lindex $Groups $i] 135],"

            puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 20]\","
            puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 21],"
            puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 22],"

            set nParameters [expr {[lindex [lindex $Groups $i] 26]}]
            puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

            for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                set j [expr {$iParam+1}]
                set k [expr {$j+26}]
                if {$j eq $nParameters} {
                    puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                } else {
                    puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                }
            }

            puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
            puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
            puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
            puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
            puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
            puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
            puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
            puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16]"

            puts $FileVar "         \},"
            puts $FileVar "         \"Tables\": \{\}"
            puts $FileVar "      \}"
        } elseif {[lindex [lindex $Groups $i] 4] eq "SmallStrainUMAT2DInterfaceLaw" || [lindex [lindex $Groups $i] 4] eq "SmallStrainUMAT3DInterfaceLaw"} {
            incr PropertyId
            dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
            if {$PropertyId > 1} {
                puts $FileVar "   \},\{"
            }
            puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $FileVar "      \"properties_id\":           $PropertyId,"
            puts $FileVar "      \"Material\": \{"
            puts $FileVar "          \"constitutive_law\": \{"
            puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 4]\" "
            puts $FileVar "          \},"
            puts $FileVar "          \"Variables\": \{"
            puts $FileVar "              \"IGNORE_UNDRAINED\"         :  true,"
            puts $FileVar "              \"CONSIDER_GAP_CLOSURE\"     :  [lindex [lindex $Groups $i] 135],"

            puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 20]\","
            puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 21],"
            puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 22],"

            set nParameters [expr {[lindex [lindex $Groups $i] 26]}]
            puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

            for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                set j [expr {$iParam+1}]
                set k [expr {$j+26}]
                if {$j eq $nParameters} {
                    puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                } else {
                    puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                }
            }

            set nStateVariables [expr {[lindex [lindex $Groups $i] 77]}]
            puts $FileVar "              \"STATE_VARIABLES\"          :  \["

            for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                set j [expr {$iStateVar+1}]
                set k [expr {$j+77}]
                if {$j eq $nStateVariables} {
                    puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                } else {
                    puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                }
            }


            puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
            puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
            puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
            puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
            puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
            puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
            puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
            puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16]"

            puts $FileVar "         \},"
            puts $FileVar "         \"Tables\": \{\}"
            puts $FileVar "      \}"
        }
    }


    # Interface two-phase or undrained part
    set interface_Groups [list [GiD_Info conditions Interface_two_phase groups] [GiD_Info conditions Interface_undrained groups]]
    foreach Groups $interface_Groups {
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            if {[lindex [lindex $Groups $i] 4] eq "BilinearCohesive3DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"BilinearCohesive3DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 5],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"DAMAGE_THRESHOLD\"         :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"CRITICAL_DISPLACEMENT\"    :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"CONSIDER_GAP_CLOSURE\"     :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"YIELD_STRESS\"             :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"FRICTION_COEFFICIENT\"     :  [lindex [lindex $Groups $i] 19]"
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            } elseif {[lindex [lindex $Groups $i] 4] eq "BilinearCohesivePlaneStrain2DLaw" || [lindex [lindex $Groups $i] 4] eq "BilinearCohesivePlaneStress2DLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"BilinearCohesive2DLaw\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"YOUNG_MODULUS\"            :  [lindex [lindex $Groups $i] 5],"
                puts $FileVar "              \"POISSON_RATIO\"            :  [lindex [lindex $Groups $i] 6],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"THICKNESS\"                :  [lindex [lindex $Groups $i] 14],"
                puts $FileVar "              \"DAMAGE_THRESHOLD\"         :  [lindex [lindex $Groups $i] 15],"
                puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16],"
                puts $FileVar "              \"CRITICAL_DISPLACEMENT\"    :  [lindex [lindex $Groups $i] 17],"
                puts $FileVar "              \"CONSIDER_GAP_CLOSURE\"     :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"YIELD_STRESS\"             :  [lindex [lindex $Groups $i] 18],"
                puts $FileVar "              \"FRICTION_COEFFICIENT\"     :  [lindex [lindex $Groups $i] 19]"
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            } elseif {[lindex [lindex $Groups $i] 4] eq "SmallStrainUDSM2DInterfaceLaw" || [lindex [lindex $Groups $i] 4] eq "SmallStrainUDSM3DInterfaceLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 4]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 20]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 21],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 22],"

                set nParameters [expr {[lindex [lindex $Groups $i] 26]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+26}]
                    if {$j eq $nParameters} {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                       puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                puts $FileVar "              \"CONSIDER_GAP_CLOSURE\"     :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16]"

                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            } elseif {[lindex [lindex $Groups $i] 4] eq "SmallStrainUMAT2DInterfaceLaw" || [lindex [lindex $Groups $i] 4] eq "SmallStrainUMAT3DInterfaceLaw"} {
                incr PropertyId
                dict set PropertyDict [lindex [lindex $Groups $i] 1] $PropertyId
                if {$PropertyId > 1} {
                    puts $FileVar "   \},\{"
                }
                puts $FileVar "      \"model_part_name\":         \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $FileVar "      \"properties_id\":           $PropertyId,"
                puts $FileVar "      \"Material\": \{"
                puts $FileVar "          \"constitutive_law\": \{"
                puts $FileVar "              \"name\"             :  \"[lindex [lindex $Groups $i] 4]\" "
                puts $FileVar "          \},"
                puts $FileVar "          \"Variables\": \{"
                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 20]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 21],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 22],"

                set nParameters [expr {[lindex [lindex $Groups $i] 26]}]
                puts $FileVar "              \"UMAT_PARAMETERS\"          :  \["

                for {set iParam 0} {$iParam < $nParameters} {incr iParam} {
                    set j [expr {$iParam+1}]
                    set k [expr {$j+26}]
                    if {$j eq $nParameters} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                set nStateVariables [expr {[lindex [lindex $Groups $i] 77]}]
                puts $FileVar "              \"STATE_VARIABLES\"          :  \["

                for {set iStateVar 0} {$iStateVar < $nStateVariables} {incr iStateVar} {
                    set j [expr {$iStateVar+1}]
                    set k [expr {$j+77}]
                    if {$j eq $nStateVariables} {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k]\],"
                    } else {
                        puts $FileVar "                                              [lindex [lindex $Groups $i] $k],"
                    }
                }

                puts $FileVar "              \"CONSIDER_GAP_CLOSURE\"     :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
                puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
                puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
                puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
                puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
                puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16]"

                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            }
        }
    }

    puts $FileVar "   \}\]"
    puts $FileVar "\}"

    close $FileVar
}

