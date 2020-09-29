proc WriteStructuralMaterials {basename dir problemtypedir TableDict} {

    ## Start ProjectParameters.json file
    set filename [file join $dir MaterialParameters.json]
    set FileVar [open $filename w]
    puts $FileVar "\{"

    ## problem_data
    set PropertyId 0
    set PropertyDict [dict create]

    puts $FileVar "   \"properties\": \[\{"

    set IsK0 [GiD_AccessValue get gendata Solution_Type]
    set Dim [GiD_AccessValue get gendata Domain_Size]
    
    # Soil_two_phase part
    set Groups [GiD_Info conditions Soil_two_phase groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {$IsK0 eq "K0-Procedure"} {
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
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31]"
                
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
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31]"
                
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
                puts $FileVar "              \"MAXIMUM_STRESS\"           :  [lindex [lindex $Groups $i] 27]"
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
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19]"
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
                puts $FileVar "              \"THICKNESS\"                :  1.0"
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88]"
                
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88]"
                        
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88],"

                puts $FileVar "              \"NUMBER_OF_UMAT_STATE_VARIABLES\":  [lindex [lindex $Groups $i] 89],"

                puts $FileVar "              \"STATE_VARIABLE_1\"         :  [lindex [lindex $Groups $i] 90],"
                puts $FileVar "              \"STATE_VARIABLE_2\"         :  [lindex [lindex $Groups $i] 91],"
                puts $FileVar "              \"STATE_VARIABLE_3\"         :  [lindex [lindex $Groups $i] 92],"
                puts $FileVar "              \"STATE_VARIABLE_4\"         :  [lindex [lindex $Groups $i] 93],"
                puts $FileVar "              \"STATE_VARIABLE_5\"         :  [lindex [lindex $Groups $i] 94],"
                puts $FileVar "              \"STATE_VARIABLE_6\"         :  [lindex [lindex $Groups $i] 95],"
                puts $FileVar "              \"STATE_VARIABLE_7\"         :  [lindex [lindex $Groups $i] 96],"
                puts $FileVar "              \"STATE_VARIABLE_8\"         :  [lindex [lindex $Groups $i] 97],"
                puts $FileVar "              \"STATE_VARIABLE_9\"         :  [lindex [lindex $Groups $i] 98],"
                puts $FileVar "              \"STATE_VARIABLE_10\"         :  [lindex [lindex $Groups $i] 99],"
                puts $FileVar "              \"STATE_VARIABLE_11\"         :  [lindex [lindex $Groups $i] 100],"
                
                puts $FileVar "              \"STATE_VARIABLE_12\"         :  [lindex [lindex $Groups $i] 101],"
                puts $FileVar "              \"STATE_VARIABLE_13\"         :  [lindex [lindex $Groups $i] 102],"
                puts $FileVar "              \"STATE_VARIABLE_14\"         :  [lindex [lindex $Groups $i] 103],"
                puts $FileVar "              \"STATE_VARIABLE_15\"         :  [lindex [lindex $Groups $i] 104],"
                puts $FileVar "              \"STATE_VARIABLE_16\"         :  [lindex [lindex $Groups $i] 105],"
                puts $FileVar "              \"STATE_VARIABLE_17\"         :  [lindex [lindex $Groups $i] 106],"
                puts $FileVar "              \"STATE_VARIABLE_18\"         :  [lindex [lindex $Groups $i] 107],"
                puts $FileVar "              \"STATE_VARIABLE_19\"         :  [lindex [lindex $Groups $i] 108],"
                puts $FileVar "              \"STATE_VARIABLE_20\"         :  [lindex [lindex $Groups $i] 109],"
                puts $FileVar "              \"STATE_VARIABLE_21\"         :  [lindex [lindex $Groups $i] 110],"
                
                puts $FileVar "              \"STATE_VARIABLE_22\"         :  [lindex [lindex $Groups $i] 111],"
                puts $FileVar "              \"STATE_VARIABLE_23\"         :  [lindex [lindex $Groups $i] 112],"
                puts $FileVar "              \"STATE_VARIABLE_24\"         :  [lindex [lindex $Groups $i] 113],"
                puts $FileVar "              \"STATE_VARIABLE_25\"         :  [lindex [lindex $Groups $i] 114],"
                puts $FileVar "              \"STATE_VARIABLE_26\"         :  [lindex [lindex $Groups $i] 115],"
                puts $FileVar "              \"STATE_VARIABLE_27\"         :  [lindex [lindex $Groups $i] 116],"
                puts $FileVar "              \"STATE_VARIABLE_28\"         :  [lindex [lindex $Groups $i] 117],"
                puts $FileVar "              \"STATE_VARIABLE_29\"         :  [lindex [lindex $Groups $i] 118],"
                puts $FileVar "              \"STATE_VARIABLE_30\"         :  [lindex [lindex $Groups $i] 119],"
                puts $FileVar "              \"STATE_VARIABLE_31\"         :  [lindex [lindex $Groups $i] 120],"
                
                puts $FileVar "              \"STATE_VARIABLE_32\"         :  [lindex [lindex $Groups $i] 121],"
                puts $FileVar "              \"STATE_VARIABLE_33\"         :  [lindex [lindex $Groups $i] 122],"
                puts $FileVar "              \"STATE_VARIABLE_34\"         :  [lindex [lindex $Groups $i] 123],"
                puts $FileVar "              \"STATE_VARIABLE_35\"         :  [lindex [lindex $Groups $i] 124],"
                puts $FileVar "              \"STATE_VARIABLE_36\"         :  [lindex [lindex $Groups $i] 125],"
                puts $FileVar "              \"STATE_VARIABLE_37\"         :  [lindex [lindex $Groups $i] 126],"
                puts $FileVar "              \"STATE_VARIABLE_38\"         :  [lindex [lindex $Groups $i] 127],"
                puts $FileVar "              \"STATE_VARIABLE_39\"         :  [lindex [lindex $Groups $i] 128],"
                puts $FileVar "              \"STATE_VARIABLE_40\"         :  [lindex [lindex $Groups $i] 129],"
                puts $FileVar "              \"STATE_VARIABLE_41\"         :  [lindex [lindex $Groups $i] 130],"
                
                puts $FileVar "              \"STATE_VARIABLE_42\"         :  [lindex [lindex $Groups $i] 131],"
                puts $FileVar "              \"STATE_VARIABLE_43\"         :  [lindex [lindex $Groups $i] 132],"
                puts $FileVar "              \"STATE_VARIABLE_44\"         :  [lindex [lindex $Groups $i] 133],"
                puts $FileVar "              \"STATE_VARIABLE_45\"         :  [lindex [lindex $Groups $i] 134],"
                puts $FileVar "              \"STATE_VARIABLE_46\"         :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"STATE_VARIABLE_47\"         :  [lindex [lindex $Groups $i] 136],"
                puts $FileVar "              \"STATE_VARIABLE_48\"         :  [lindex [lindex $Groups $i] 137],"
                puts $FileVar "              \"STATE_VARIABLE_49\"         :  [lindex [lindex $Groups $i] 138],"
                puts $FileVar "              \"STATE_VARIABLE_50\"         :  [lindex [lindex $Groups $i] 139]"

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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88]"
                
                puts $FileVar "              \"NUMBER_OF_UMAT_STATE_VARIABLES\":  [lindex [lindex $Groups $i] 89],"

                puts $FileVar "              \"STATE_VARIABLE_1\"         :  [lindex [lindex $Groups $i] 90],"
                puts $FileVar "              \"STATE_VARIABLE_2\"         :  [lindex [lindex $Groups $i] 91],"
                puts $FileVar "              \"STATE_VARIABLE_3\"         :  [lindex [lindex $Groups $i] 92],"
                puts $FileVar "              \"STATE_VARIABLE_4\"         :  [lindex [lindex $Groups $i] 93],"
                puts $FileVar "              \"STATE_VARIABLE_5\"         :  [lindex [lindex $Groups $i] 94],"
                puts $FileVar "              \"STATE_VARIABLE_6\"         :  [lindex [lindex $Groups $i] 95],"
                puts $FileVar "              \"STATE_VARIABLE_7\"         :  [lindex [lindex $Groups $i] 96],"
                puts $FileVar "              \"STATE_VARIABLE_8\"         :  [lindex [lindex $Groups $i] 97],"
                puts $FileVar "              \"STATE_VARIABLE_9\"         :  [lindex [lindex $Groups $i] 98],"
                puts $FileVar "              \"STATE_VARIABLE_10\"         :  [lindex [lindex $Groups $i] 99],"
                puts $FileVar "              \"STATE_VARIABLE_11\"         :  [lindex [lindex $Groups $i] 100],"
                
                puts $FileVar "              \"STATE_VARIABLE_12\"         :  [lindex [lindex $Groups $i] 101],"
                puts $FileVar "              \"STATE_VARIABLE_13\"         :  [lindex [lindex $Groups $i] 102],"
                puts $FileVar "              \"STATE_VARIABLE_14\"         :  [lindex [lindex $Groups $i] 103],"
                puts $FileVar "              \"STATE_VARIABLE_15\"         :  [lindex [lindex $Groups $i] 104],"
                puts $FileVar "              \"STATE_VARIABLE_16\"         :  [lindex [lindex $Groups $i] 105],"
                puts $FileVar "              \"STATE_VARIABLE_17\"         :  [lindex [lindex $Groups $i] 106],"
                puts $FileVar "              \"STATE_VARIABLE_18\"         :  [lindex [lindex $Groups $i] 107],"
                puts $FileVar "              \"STATE_VARIABLE_19\"         :  [lindex [lindex $Groups $i] 108],"
                puts $FileVar "              \"STATE_VARIABLE_20\"         :  [lindex [lindex $Groups $i] 109],"
                puts $FileVar "              \"STATE_VARIABLE_21\"         :  [lindex [lindex $Groups $i] 110],"
                
                puts $FileVar "              \"STATE_VARIABLE_22\"         :  [lindex [lindex $Groups $i] 111],"
                puts $FileVar "              \"STATE_VARIABLE_23\"         :  [lindex [lindex $Groups $i] 112],"
                puts $FileVar "              \"STATE_VARIABLE_24\"         :  [lindex [lindex $Groups $i] 113],"
                puts $FileVar "              \"STATE_VARIABLE_25\"         :  [lindex [lindex $Groups $i] 114],"
                puts $FileVar "              \"STATE_VARIABLE_26\"         :  [lindex [lindex $Groups $i] 115],"
                puts $FileVar "              \"STATE_VARIABLE_27\"         :  [lindex [lindex $Groups $i] 116],"
                puts $FileVar "              \"STATE_VARIABLE_28\"         :  [lindex [lindex $Groups $i] 117],"
                puts $FileVar "              \"STATE_VARIABLE_29\"         :  [lindex [lindex $Groups $i] 118],"
                puts $FileVar "              \"STATE_VARIABLE_30\"         :  [lindex [lindex $Groups $i] 119],"
                puts $FileVar "              \"STATE_VARIABLE_31\"         :  [lindex [lindex $Groups $i] 120],"
                
                puts $FileVar "              \"STATE_VARIABLE_32\"         :  [lindex [lindex $Groups $i] 121],"
                puts $FileVar "              \"STATE_VARIABLE_33\"         :  [lindex [lindex $Groups $i] 122],"
                puts $FileVar "              \"STATE_VARIABLE_34\"         :  [lindex [lindex $Groups $i] 123],"
                puts $FileVar "              \"STATE_VARIABLE_35\"         :  [lindex [lindex $Groups $i] 124],"
                puts $FileVar "              \"STATE_VARIABLE_36\"         :  [lindex [lindex $Groups $i] 125],"
                puts $FileVar "              \"STATE_VARIABLE_37\"         :  [lindex [lindex $Groups $i] 126],"
                puts $FileVar "              \"STATE_VARIABLE_38\"         :  [lindex [lindex $Groups $i] 127],"
                puts $FileVar "              \"STATE_VARIABLE_39\"         :  [lindex [lindex $Groups $i] 128],"
                puts $FileVar "              \"STATE_VARIABLE_40\"         :  [lindex [lindex $Groups $i] 129],"
                puts $FileVar "              \"STATE_VARIABLE_41\"         :  [lindex [lindex $Groups $i] 130],"
                
                puts $FileVar "              \"STATE_VARIABLE_42\"         :  [lindex [lindex $Groups $i] 131],"
                puts $FileVar "              \"STATE_VARIABLE_43\"         :  [lindex [lindex $Groups $i] 132],"
                puts $FileVar "              \"STATE_VARIABLE_44\"         :  [lindex [lindex $Groups $i] 133],"
                puts $FileVar "              \"STATE_VARIABLE_45\"         :  [lindex [lindex $Groups $i] 134],"
                puts $FileVar "              \"STATE_VARIABLE_46\"         :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"STATE_VARIABLE_47\"         :  [lindex [lindex $Groups $i] 136],"
                puts $FileVar "              \"STATE_VARIABLE_48\"         :  [lindex [lindex $Groups $i] 137],"
                puts $FileVar "              \"STATE_VARIABLE_49\"         :  [lindex [lindex $Groups $i] 138],"
                puts $FileVar "              \"STATE_VARIABLE_50\"         :  [lindex [lindex $Groups $i] 139]"

                
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            }
        }
    }

    # Soil_drained part
    set Groups [GiD_Info conditions Soil_drained groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {$IsK0 eq "K0-Procedure"} {
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
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31]"
                
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"

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
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31]"
                
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
                puts $FileVar "              \"MAXIMUM_STRESS\"           :  [lindex [lindex $Groups $i] 27]"
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
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19]"
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
                puts $FileVar "              \"THICKNESS\"                :  1.0"
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88]"
                
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88]"
                
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88],"

                puts $FileVar "              \"NUMBER_OF_UMAT_STATE_VARIABLES\":  [lindex [lindex $Groups $i] 89],"

                puts $FileVar "              \"STATE_VARIABLE_1\"         :  [lindex [lindex $Groups $i] 90],"
                puts $FileVar "              \"STATE_VARIABLE_2\"         :  [lindex [lindex $Groups $i] 91],"
                puts $FileVar "              \"STATE_VARIABLE_3\"         :  [lindex [lindex $Groups $i] 92],"
                puts $FileVar "              \"STATE_VARIABLE_4\"         :  [lindex [lindex $Groups $i] 93],"
                puts $FileVar "              \"STATE_VARIABLE_5\"         :  [lindex [lindex $Groups $i] 94],"
                puts $FileVar "              \"STATE_VARIABLE_6\"         :  [lindex [lindex $Groups $i] 95],"
                puts $FileVar "              \"STATE_VARIABLE_7\"         :  [lindex [lindex $Groups $i] 96],"
                puts $FileVar "              \"STATE_VARIABLE_8\"         :  [lindex [lindex $Groups $i] 97],"
                puts $FileVar "              \"STATE_VARIABLE_9\"         :  [lindex [lindex $Groups $i] 98],"
                puts $FileVar "              \"STATE_VARIABLE_10\"        :  [lindex [lindex $Groups $i] 99],"
                
                puts $FileVar "              \"STATE_VARIABLE_11\"        :  [lindex [lindex $Groups $i] 100],"
                puts $FileVar "              \"STATE_VARIABLE_12\"        :  [lindex [lindex $Groups $i] 101],"
                puts $FileVar "              \"STATE_VARIABLE_13\"        :  [lindex [lindex $Groups $i] 102],"
                puts $FileVar "              \"STATE_VARIABLE_14\"        :  [lindex [lindex $Groups $i] 103],"
                puts $FileVar "              \"STATE_VARIABLE_15\"        :  [lindex [lindex $Groups $i] 104],"
                puts $FileVar "              \"STATE_VARIABLE_16\"        :  [lindex [lindex $Groups $i] 105],"
                puts $FileVar "              \"STATE_VARIABLE_17\"        :  [lindex [lindex $Groups $i] 106],"
                puts $FileVar "              \"STATE_VARIABLE_18\"        :  [lindex [lindex $Groups $i] 107],"
                puts $FileVar "              \"STATE_VARIABLE_19\"        :  [lindex [lindex $Groups $i] 108],"
                puts $FileVar "              \"STATE_VARIABLE_20\"        :  [lindex [lindex $Groups $i] 109],"
                
                puts $FileVar "              \"STATE_VARIABLE_21\"        :  [lindex [lindex $Groups $i] 110],"
                puts $FileVar "              \"STATE_VARIABLE_22\"        :  [lindex [lindex $Groups $i] 111],"
                puts $FileVar "              \"STATE_VARIABLE_23\"        :  [lindex [lindex $Groups $i] 112],"
                puts $FileVar "              \"STATE_VARIABLE_24\"        :  [lindex [lindex $Groups $i] 113],"
                puts $FileVar "              \"STATE_VARIABLE_25\"        :  [lindex [lindex $Groups $i] 114],"
                puts $FileVar "              \"STATE_VARIABLE_26\"        :  [lindex [lindex $Groups $i] 115],"
                puts $FileVar "              \"STATE_VARIABLE_27\"        :  [lindex [lindex $Groups $i] 116],"
                puts $FileVar "              \"STATE_VARIABLE_28\"        :  [lindex [lindex $Groups $i] 117],"
                puts $FileVar "              \"STATE_VARIABLE_29\"        :  [lindex [lindex $Groups $i] 118],"
                puts $FileVar "              \"STATE_VARIABLE_30\"        :  [lindex [lindex $Groups $i] 119],"
                
                puts $FileVar "              \"STATE_VARIABLE_31\"        :  [lindex [lindex $Groups $i] 120],"
                puts $FileVar "              \"STATE_VARIABLE_32\"        :  [lindex [lindex $Groups $i] 121],"
                puts $FileVar "              \"STATE_VARIABLE_33\"        :  [lindex [lindex $Groups $i] 122],"
                puts $FileVar "              \"STATE_VARIABLE_34\"        :  [lindex [lindex $Groups $i] 123],"
                puts $FileVar "              \"STATE_VARIABLE_35\"        :  [lindex [lindex $Groups $i] 124],"
                puts $FileVar "              \"STATE_VARIABLE_36\"        :  [lindex [lindex $Groups $i] 125],"
                puts $FileVar "              \"STATE_VARIABLE_37\"        :  [lindex [lindex $Groups $i] 126],"
                puts $FileVar "              \"STATE_VARIABLE_38\"        :  [lindex [lindex $Groups $i] 127],"
                puts $FileVar "              \"STATE_VARIABLE_39\"        :  [lindex [lindex $Groups $i] 128],"
                puts $FileVar "              \"STATE_VARIABLE_40\"        :  [lindex [lindex $Groups $i] 129],"
                
                puts $FileVar "              \"STATE_VARIABLE_41\"        :  [lindex [lindex $Groups $i] 130],"
                puts $FileVar "              \"STATE_VARIABLE_42\"        :  [lindex [lindex $Groups $i] 131],"
                puts $FileVar "              \"STATE_VARIABLE_43\"        :  [lindex [lindex $Groups $i] 132],"
                puts $FileVar "              \"STATE_VARIABLE_44\"        :  [lindex [lindex $Groups $i] 133],"
                puts $FileVar "              \"STATE_VARIABLE_45\"        :  [lindex [lindex $Groups $i] 134],"
                puts $FileVar "              \"STATE_VARIABLE_46\"        :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"STATE_VARIABLE_47\"        :  [lindex [lindex $Groups $i] 136],"
                puts $FileVar "              \"STATE_VARIABLE_48\"        :  [lindex [lindex $Groups $i] 137],"
                puts $FileVar "              \"STATE_VARIABLE_49\"        :  [lindex [lindex $Groups $i] 138],"
                puts $FileVar "              \"STATE_VARIABLE_50\"        :  [lindex [lindex $Groups $i] 139]"

                
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88],"

                puts $FileVar "              \"NUMBER_OF_UMAT_STATE_VARIABLES\":  [lindex [lindex $Groups $i] 89],"

                puts $FileVar "              \"STATE_VARIABLE_1\"         :  [lindex [lindex $Groups $i] 90],"
                puts $FileVar "              \"STATE_VARIABLE_2\"         :  [lindex [lindex $Groups $i] 91],"
                puts $FileVar "              \"STATE_VARIABLE_3\"         :  [lindex [lindex $Groups $i] 92],"
                puts $FileVar "              \"STATE_VARIABLE_4\"         :  [lindex [lindex $Groups $i] 93],"
                puts $FileVar "              \"STATE_VARIABLE_5\"         :  [lindex [lindex $Groups $i] 94],"
                puts $FileVar "              \"STATE_VARIABLE_6\"         :  [lindex [lindex $Groups $i] 95],"
                puts $FileVar "              \"STATE_VARIABLE_7\"         :  [lindex [lindex $Groups $i] 96],"
                puts $FileVar "              \"STATE_VARIABLE_8\"         :  [lindex [lindex $Groups $i] 97],"
                puts $FileVar "              \"STATE_VARIABLE_9\"         :  [lindex [lindex $Groups $i] 98],"
                puts $FileVar "              \"STATE_VARIABLE_10\"        :  [lindex [lindex $Groups $i] 99],"
                
                puts $FileVar "              \"STATE_VARIABLE_11\"        :  [lindex [lindex $Groups $i] 100],"
                puts $FileVar "              \"STATE_VARIABLE_12\"        :  [lindex [lindex $Groups $i] 101],"
                puts $FileVar "              \"STATE_VARIABLE_13\"        :  [lindex [lindex $Groups $i] 102],"
                puts $FileVar "              \"STATE_VARIABLE_14\"        :  [lindex [lindex $Groups $i] 103],"
                puts $FileVar "              \"STATE_VARIABLE_15\"        :  [lindex [lindex $Groups $i] 104],"
                puts $FileVar "              \"STATE_VARIABLE_16\"        :  [lindex [lindex $Groups $i] 105],"
                puts $FileVar "              \"STATE_VARIABLE_17\"        :  [lindex [lindex $Groups $i] 106],"
                puts $FileVar "              \"STATE_VARIABLE_18\"        :  [lindex [lindex $Groups $i] 107],"
                puts $FileVar "              \"STATE_VARIABLE_19\"        :  [lindex [lindex $Groups $i] 108],"
                puts $FileVar "              \"STATE_VARIABLE_20\"        :  [lindex [lindex $Groups $i] 109],"
                
                puts $FileVar "              \"STATE_VARIABLE_21\"        :  [lindex [lindex $Groups $i] 110],"
                puts $FileVar "              \"STATE_VARIABLE_22\"        :  [lindex [lindex $Groups $i] 111],"
                puts $FileVar "              \"STATE_VARIABLE_23\"        :  [lindex [lindex $Groups $i] 112],"
                puts $FileVar "              \"STATE_VARIABLE_24\"        :  [lindex [lindex $Groups $i] 113],"
                puts $FileVar "              \"STATE_VARIABLE_25\"        :  [lindex [lindex $Groups $i] 114],"
                puts $FileVar "              \"STATE_VARIABLE_26\"        :  [lindex [lindex $Groups $i] 115],"
                puts $FileVar "              \"STATE_VARIABLE_27\"        :  [lindex [lindex $Groups $i] 116],"
                puts $FileVar "              \"STATE_VARIABLE_28\"        :  [lindex [lindex $Groups $i] 117],"
                puts $FileVar "              \"STATE_VARIABLE_29\"        :  [lindex [lindex $Groups $i] 118],"
                puts $FileVar "              \"STATE_VARIABLE_30\"        :  [lindex [lindex $Groups $i] 119],"
                
                puts $FileVar "              \"STATE_VARIABLE_31\"        :  [lindex [lindex $Groups $i] 120],"
                puts $FileVar "              \"STATE_VARIABLE_32\"        :  [lindex [lindex $Groups $i] 121],"
                puts $FileVar "              \"STATE_VARIABLE_33\"        :  [lindex [lindex $Groups $i] 122],"
                puts $FileVar "              \"STATE_VARIABLE_34\"        :  [lindex [lindex $Groups $i] 123],"
                puts $FileVar "              \"STATE_VARIABLE_35\"        :  [lindex [lindex $Groups $i] 124],"
                puts $FileVar "              \"STATE_VARIABLE_36\"        :  [lindex [lindex $Groups $i] 125],"
                puts $FileVar "              \"STATE_VARIABLE_37\"        :  [lindex [lindex $Groups $i] 126],"
                puts $FileVar "              \"STATE_VARIABLE_38\"        :  [lindex [lindex $Groups $i] 127],"
                puts $FileVar "              \"STATE_VARIABLE_39\"        :  [lindex [lindex $Groups $i] 128],"
                puts $FileVar "              \"STATE_VARIABLE_40\"        :  [lindex [lindex $Groups $i] 129],"
                
                puts $FileVar "              \"STATE_VARIABLE_41\"        :  [lindex [lindex $Groups $i] 130],"
                puts $FileVar "              \"STATE_VARIABLE_42\"        :  [lindex [lindex $Groups $i] 131],"
                puts $FileVar "              \"STATE_VARIABLE_43\"        :  [lindex [lindex $Groups $i] 132],"
                puts $FileVar "              \"STATE_VARIABLE_44\"        :  [lindex [lindex $Groups $i] 133],"
                puts $FileVar "              \"STATE_VARIABLE_45\"        :  [lindex [lindex $Groups $i] 134],"
                puts $FileVar "              \"STATE_VARIABLE_46\"        :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"STATE_VARIABLE_47\"        :  [lindex [lindex $Groups $i] 136],"
                puts $FileVar "              \"STATE_VARIABLE_48\"        :  [lindex [lindex $Groups $i] 137],"
                puts $FileVar "              \"STATE_VARIABLE_49\"        :  [lindex [lindex $Groups $i] 138],"
                puts $FileVar "              \"STATE_VARIABLE_50\"        :  [lindex [lindex $Groups $i] 139]"

                
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            }        
        }
    }

    # Soil_undrained part
    set Groups [GiD_Info conditions Soil_undrained groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {$IsK0 eq "K0-Procedure"} {
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
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31]"
                
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
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
                puts $FileVar "              \"K0_VALUE_ZZ\"              :  [lindex [lindex $Groups $i] 31]"
                
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
                puts $FileVar "              \"MAXIMUM_STRESS\"           :  [lindex [lindex $Groups $i] 27]"
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
                puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 19]"
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
                puts $FileVar "              \"THICKNESS\"                :  1.0"
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],
                "
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88]"
                
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88]"
                
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],
                "
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88]"

                puts $FileVar "              \"NUMBER_OF_UMAT_STATE_VARIABLES\":  [lindex [lindex $Groups $i] 89],"

                puts $FileVar "              \"STATE_VARIABLE_1\"         :  [lindex [lindex $Groups $i] 90],"
                puts $FileVar "              \"STATE_VARIABLE_2\"         :  [lindex [lindex $Groups $i] 91],"
                puts $FileVar "              \"STATE_VARIABLE_3\"         :  [lindex [lindex $Groups $i] 92],"
                puts $FileVar "              \"STATE_VARIABLE_4\"         :  [lindex [lindex $Groups $i] 93],"
                puts $FileVar "              \"STATE_VARIABLE_5\"         :  [lindex [lindex $Groups $i] 94],"
                puts $FileVar "              \"STATE_VARIABLE_6\"         :  [lindex [lindex $Groups $i] 95],"
                puts $FileVar "              \"STATE_VARIABLE_7\"         :  [lindex [lindex $Groups $i] 96],"
                puts $FileVar "              \"STATE_VARIABLE_8\"         :  [lindex [lindex $Groups $i] 97],"
                puts $FileVar "              \"STATE_VARIABLE_9\"         :  [lindex [lindex $Groups $i] 98],"
                puts $FileVar "              \"STATE_VARIABLE_10\"         :  [lindex [lindex $Groups $i] 99],"
                
                puts $FileVar "              \"STATE_VARIABLE_11\"         :  [lindex [lindex $Groups $i] 100],"
                puts $FileVar "              \"STATE_VARIABLE_12\"         :  [lindex [lindex $Groups $i] 101],"
                puts $FileVar "              \"STATE_VARIABLE_13\"         :  [lindex [lindex $Groups $i] 102],"
                puts $FileVar "              \"STATE_VARIABLE_14\"         :  [lindex [lindex $Groups $i] 103],"
                puts $FileVar "              \"STATE_VARIABLE_15\"         :  [lindex [lindex $Groups $i] 104],"
                puts $FileVar "              \"STATE_VARIABLE_16\"         :  [lindex [lindex $Groups $i] 105],"
                puts $FileVar "              \"STATE_VARIABLE_17\"         :  [lindex [lindex $Groups $i] 106],"
                puts $FileVar "              \"STATE_VARIABLE_18\"         :  [lindex [lindex $Groups $i] 107],"
                puts $FileVar "              \"STATE_VARIABLE_19\"         :  [lindex [lindex $Groups $i] 108],"
                puts $FileVar "              \"STATE_VARIABLE_20\"         :  [lindex [lindex $Groups $i] 109],"
                
                puts $FileVar "              \"STATE_VARIABLE_21\"         :  [lindex [lindex $Groups $i] 110],"
                puts $FileVar "              \"STATE_VARIABLE_22\"         :  [lindex [lindex $Groups $i] 111],"
                puts $FileVar "              \"STATE_VARIABLE_23\"         :  [lindex [lindex $Groups $i] 112],"
                puts $FileVar "              \"STATE_VARIABLE_24\"         :  [lindex [lindex $Groups $i] 113],"
                puts $FileVar "              \"STATE_VARIABLE_25\"         :  [lindex [lindex $Groups $i] 114],"
                puts $FileVar "              \"STATE_VARIABLE_26\"         :  [lindex [lindex $Groups $i] 115],"
                puts $FileVar "              \"STATE_VARIABLE_27\"         :  [lindex [lindex $Groups $i] 116],"
                puts $FileVar "              \"STATE_VARIABLE_28\"         :  [lindex [lindex $Groups $i] 117],"
                puts $FileVar "              \"STATE_VARIABLE_29\"         :  [lindex [lindex $Groups $i] 118],"
                puts $FileVar "              \"STATE_VARIABLE_30\"         :  [lindex [lindex $Groups $i] 119],"
                
                puts $FileVar "              \"STATE_VARIABLE_31\"         :  [lindex [lindex $Groups $i] 120],"
                puts $FileVar "              \"STATE_VARIABLE_32\"         :  [lindex [lindex $Groups $i] 121],"
                puts $FileVar "              \"STATE_VARIABLE_33\"         :  [lindex [lindex $Groups $i] 122],"
                puts $FileVar "              \"STATE_VARIABLE_34\"         :  [lindex [lindex $Groups $i] 123],"
                puts $FileVar "              \"STATE_VARIABLE_35\"         :  [lindex [lindex $Groups $i] 124],"
                puts $FileVar "              \"STATE_VARIABLE_36\"         :  [lindex [lindex $Groups $i] 125],"
                puts $FileVar "              \"STATE_VARIABLE_37\"         :  [lindex [lindex $Groups $i] 126],"
                puts $FileVar "              \"STATE_VARIABLE_38\"         :  [lindex [lindex $Groups $i] 127],"
                puts $FileVar "              \"STATE_VARIABLE_39\"         :  [lindex [lindex $Groups $i] 128],"
                puts $FileVar "              \"STATE_VARIABLE_40\"         :  [lindex [lindex $Groups $i] 129],"
                
                puts $FileVar "              \"STATE_VARIABLE_41\"         :  [lindex [lindex $Groups $i] 130],"
                puts $FileVar "              \"STATE_VARIABLE_42\"         :  [lindex [lindex $Groups $i] 131],"
                puts $FileVar "              \"STATE_VARIABLE_43\"         :  [lindex [lindex $Groups $i] 132],"
                puts $FileVar "              \"STATE_VARIABLE_44\"         :  [lindex [lindex $Groups $i] 133],"
                puts $FileVar "              \"STATE_VARIABLE_45\"         :  [lindex [lindex $Groups $i] 134],"
                puts $FileVar "              \"STATE_VARIABLE_46\"         :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"STATE_VARIABLE_47\"         :  [lindex [lindex $Groups $i] 136],"
                puts $FileVar "              \"STATE_VARIABLE_48\"         :  [lindex [lindex $Groups $i] 137],"
                puts $FileVar "              \"STATE_VARIABLE_49\"         :  [lindex [lindex $Groups $i] 138],"
                puts $FileVar "              \"STATE_VARIABLE_50\"         :  [lindex [lindex $Groups $i] 139]"

                
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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 78],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 79],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 88],"

                puts $FileVar "              \"NUMBER_OF_UMAT_STATE_VARIABLES\":  [lindex [lindex $Groups $i] 89],"

                puts $FileVar "              \"STATE_VARIABLE_1\"         :  [lindex [lindex $Groups $i] 90],"
                puts $FileVar "              \"STATE_VARIABLE_2\"         :  [lindex [lindex $Groups $i] 91],"
                puts $FileVar "              \"STATE_VARIABLE_3\"         :  [lindex [lindex $Groups $i] 92],"
                puts $FileVar "              \"STATE_VARIABLE_4\"         :  [lindex [lindex $Groups $i] 93],"
                puts $FileVar "              \"STATE_VARIABLE_5\"         :  [lindex [lindex $Groups $i] 94],"
                puts $FileVar "              \"STATE_VARIABLE_6\"         :  [lindex [lindex $Groups $i] 95],"
                puts $FileVar "              \"STATE_VARIABLE_7\"         :  [lindex [lindex $Groups $i] 96],"
                puts $FileVar "              \"STATE_VARIABLE_8\"         :  [lindex [lindex $Groups $i] 97],"
                puts $FileVar "              \"STATE_VARIABLE_9\"         :  [lindex [lindex $Groups $i] 98],"
                puts $FileVar "              \"STATE_VARIABLE_10\"         :  [lindex [lindex $Groups $i] 99],"
                
                puts $FileVar "              \"STATE_VARIABLE_11\"         :  [lindex [lindex $Groups $i] 100],"
                puts $FileVar "              \"STATE_VARIABLE_12\"         :  [lindex [lindex $Groups $i] 101],"
                puts $FileVar "              \"STATE_VARIABLE_13\"         :  [lindex [lindex $Groups $i] 102],"
                puts $FileVar "              \"STATE_VARIABLE_14\"         :  [lindex [lindex $Groups $i] 103],"
                puts $FileVar "              \"STATE_VARIABLE_15\"         :  [lindex [lindex $Groups $i] 104],"
                puts $FileVar "              \"STATE_VARIABLE_16\"         :  [lindex [lindex $Groups $i] 105],"
                puts $FileVar "              \"STATE_VARIABLE_17\"         :  [lindex [lindex $Groups $i] 106],"
                puts $FileVar "              \"STATE_VARIABLE_18\"         :  [lindex [lindex $Groups $i] 107],"
                puts $FileVar "              \"STATE_VARIABLE_19\"         :  [lindex [lindex $Groups $i] 108],"
                puts $FileVar "              \"STATE_VARIABLE_20\"         :  [lindex [lindex $Groups $i] 109],"
                
                puts $FileVar "              \"STATE_VARIABLE_21\"         :  [lindex [lindex $Groups $i] 110],"
                puts $FileVar "              \"STATE_VARIABLE_22\"         :  [lindex [lindex $Groups $i] 111],"
                puts $FileVar "              \"STATE_VARIABLE_23\"         :  [lindex [lindex $Groups $i] 112],"
                puts $FileVar "              \"STATE_VARIABLE_24\"         :  [lindex [lindex $Groups $i] 113],"
                puts $FileVar "              \"STATE_VARIABLE_25\"         :  [lindex [lindex $Groups $i] 114],"
                puts $FileVar "              \"STATE_VARIABLE_26\"         :  [lindex [lindex $Groups $i] 115],"
                puts $FileVar "              \"STATE_VARIABLE_27\"         :  [lindex [lindex $Groups $i] 116],"
                puts $FileVar "              \"STATE_VARIABLE_28\"         :  [lindex [lindex $Groups $i] 117],"
                puts $FileVar "              \"STATE_VARIABLE_29\"         :  [lindex [lindex $Groups $i] 118],"
                puts $FileVar "              \"STATE_VARIABLE_30\"         :  [lindex [lindex $Groups $i] 119],"
                
                puts $FileVar "              \"STATE_VARIABLE_31\"         :  [lindex [lindex $Groups $i] 120],"
                puts $FileVar "              \"STATE_VARIABLE_32\"         :  [lindex [lindex $Groups $i] 121],"
                puts $FileVar "              \"STATE_VARIABLE_33\"         :  [lindex [lindex $Groups $i] 122],"
                puts $FileVar "              \"STATE_VARIABLE_34\"         :  [lindex [lindex $Groups $i] 123],"
                puts $FileVar "              \"STATE_VARIABLE_35\"         :  [lindex [lindex $Groups $i] 124],"
                puts $FileVar "              \"STATE_VARIABLE_36\"         :  [lindex [lindex $Groups $i] 125],"
                puts $FileVar "              \"STATE_VARIABLE_37\"         :  [lindex [lindex $Groups $i] 126],"
                puts $FileVar "              \"STATE_VARIABLE_38\"         :  [lindex [lindex $Groups $i] 127],"
                puts $FileVar "              \"STATE_VARIABLE_39\"         :  [lindex [lindex $Groups $i] 128],"
                puts $FileVar "              \"STATE_VARIABLE_40\"         :  [lindex [lindex $Groups $i] 129],"
                
                puts $FileVar "              \"STATE_VARIABLE_41\"         :  [lindex [lindex $Groups $i] 130],"
                puts $FileVar "              \"STATE_VARIABLE_42\"         :  [lindex [lindex $Groups $i] 131],"
                puts $FileVar "              \"STATE_VARIABLE_43\"         :  [lindex [lindex $Groups $i] 132],"
                puts $FileVar "              \"STATE_VARIABLE_44\"         :  [lindex [lindex $Groups $i] 133],"
                puts $FileVar "              \"STATE_VARIABLE_45\"         :  [lindex [lindex $Groups $i] 134],"
                puts $FileVar "              \"STATE_VARIABLE_46\"         :  [lindex [lindex $Groups $i] 135],"
                puts $FileVar "              \"STATE_VARIABLE_47\"         :  [lindex [lindex $Groups $i] 136],"
                puts $FileVar "              \"STATE_VARIABLE_48\"         :  [lindex [lindex $Groups $i] 137],"
                puts $FileVar "              \"STATE_VARIABLE_49\"         :  [lindex [lindex $Groups $i] 138],"
                puts $FileVar "              \"STATE_VARIABLE_50\"         :  [lindex [lindex $Groups $i] 139]"

                
                
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            }
        }
    }

    # Non_porous part
    set Groups [GiD_Info conditions Non_porous groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {$IsK0 eq "K0-Procedure"} {
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
                puts $FileVar "              \"DENSITY\"                  :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 22]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 23],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 24],"

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 31],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 32],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 34],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 35],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 36],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 37],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 38],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 78]"

                
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
                puts $FileVar "              \"DENSITY\"                  :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 22]\","
                puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 23],"
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 24],"

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 31],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 32],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 34],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 35],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 36],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 37],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 38],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 78]"
            
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
                puts $FileVar "              \"DENSITY\"                  :  [lindex [lindex $Groups $i] 8],"
                puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"

                puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 25]\","
                puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 27],"

                puts $FileVar "              \"NUMBER_OF_UMAT_PARAMETERS\":  [lindex [lindex $Groups $i] 28],"

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 31],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 32],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 34],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 35],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 36],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 37],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 38],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 78],"

                puts $FileVar "              \"NUMBER_OF_UMAT_STATE_VARIABLES\":  [lindex [lindex $Groups $i] 79],"

                puts $FileVar "              \"STATE_VARIABLE_1\"         :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"STATE_VARIABLE_2\"         :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"STATE_VARIABLE_3\"         :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"STATE_VARIABLE_4\"         :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"STATE_VARIABLE_5\"         :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"STATE_VARIABLE_6\"         :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"STATE_VARIABLE_7\"         :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"STATE_VARIABLE_8\"         :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"STATE_VARIABLE_9\"         :  [lindex [lindex $Groups $i] 88],"
                puts $FileVar "              \"STATE_VARIABLE_10\"         :  [lindex [lindex $Groups $i] 89],"
                
                puts $FileVar "              \"STATE_VARIABLE_11\"         :  [lindex [lindex $Groups $i] 90],"
                puts $FileVar "              \"STATE_VARIABLE_12\"         :  [lindex [lindex $Groups $i] 91],"
                puts $FileVar "              \"STATE_VARIABLE_13\"         :  [lindex [lindex $Groups $i] 92],"
                puts $FileVar "              \"STATE_VARIABLE_14\"         :  [lindex [lindex $Groups $i] 93],"
                puts $FileVar "              \"STATE_VARIABLE_15\"         :  [lindex [lindex $Groups $i] 94],"
                puts $FileVar "              \"STATE_VARIABLE_16\"         :  [lindex [lindex $Groups $i] 95],"
                puts $FileVar "              \"STATE_VARIABLE_17\"         :  [lindex [lindex $Groups $i] 96],"
                puts $FileVar "              \"STATE_VARIABLE_18\"         :  [lindex [lindex $Groups $i] 97],"
                puts $FileVar "              \"STATE_VARIABLE_19\"         :  [lindex [lindex $Groups $i] 98],"
                puts $FileVar "              \"STATE_VARIABLE_20\"         :  [lindex [lindex $Groups $i] 99],"
                
                puts $FileVar "              \"STATE_VARIABLE_21\"         :  [lindex [lindex $Groups $i] 100],"
                puts $FileVar "              \"STATE_VARIABLE_22\"         :  [lindex [lindex $Groups $i] 101],"
                puts $FileVar "              \"STATE_VARIABLE_23\"         :  [lindex [lindex $Groups $i] 102],"
                puts $FileVar "              \"STATE_VARIABLE_24\"         :  [lindex [lindex $Groups $i] 103],"
                puts $FileVar "              \"STATE_VARIABLE_25\"         :  [lindex [lindex $Groups $i] 104],"
                puts $FileVar "              \"STATE_VARIABLE_26\"         :  [lindex [lindex $Groups $i] 105],"
                puts $FileVar "              \"STATE_VARIABLE_27\"         :  [lindex [lindex $Groups $i] 106],"
                puts $FileVar "              \"STATE_VARIABLE_28\"         :  [lindex [lindex $Groups $i] 107],"
                puts $FileVar "              \"STATE_VARIABLE_29\"         :  [lindex [lindex $Groups $i] 108],"
                puts $FileVar "              \"STATE_VARIABLE_30\"         :  [lindex [lindex $Groups $i] 109],"
                
                puts $FileVar "              \"STATE_VARIABLE_31\"         :  [lindex [lindex $Groups $i] 110],"
                puts $FileVar "              \"STATE_VARIABLE_32\"         :  [lindex [lindex $Groups $i] 111],"
                puts $FileVar "              \"STATE_VARIABLE_33\"         :  [lindex [lindex $Groups $i] 112],"
                puts $FileVar "              \"STATE_VARIABLE_34\"         :  [lindex [lindex $Groups $i] 113],"
                puts $FileVar "              \"STATE_VARIABLE_35\"         :  [lindex [lindex $Groups $i] 114],"
                puts $FileVar "              \"STATE_VARIABLE_36\"         :  [lindex [lindex $Groups $i] 115],"
                puts $FileVar "              \"STATE_VARIABLE_37\"         :  [lindex [lindex $Groups $i] 116],"
                puts $FileVar "              \"STATE_VARIABLE_38\"         :  [lindex [lindex $Groups $i] 117],"
                puts $FileVar "              \"STATE_VARIABLE_39\"         :  [lindex [lindex $Groups $i] 118],"
                puts $FileVar "              \"STATE_VARIABLE_40\"         :  [lindex [lindex $Groups $i] 119],"
                
                puts $FileVar "              \"STATE_VARIABLE_41\"         :  [lindex [lindex $Groups $i] 120],"
                puts $FileVar "              \"STATE_VARIABLE_42\"         :  [lindex [lindex $Groups $i] 121],"
                puts $FileVar "              \"STATE_VARIABLE_43\"         :  [lindex [lindex $Groups $i] 122],"
                puts $FileVar "              \"STATE_VARIABLE_44\"         :  [lindex [lindex $Groups $i] 123],"
                puts $FileVar "              \"STATE_VARIABLE_45\"         :  [lindex [lindex $Groups $i] 124],"
                puts $FileVar "              \"STATE_VARIABLE_46\"         :  [lindex [lindex $Groups $i] 125],"
                puts $FileVar "              \"STATE_VARIABLE_47\"         :  [lindex [lindex $Groups $i] 126],"
                puts $FileVar "              \"STATE_VARIABLE_48\"         :  [lindex [lindex $Groups $i] 127],"
                puts $FileVar "              \"STATE_VARIABLE_49\"         :  [lindex [lindex $Groups $i] 128],"
                puts $FileVar "              \"STATE_VARIABLE_50\"         :  [lindex [lindex $Groups $i] 129],"

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

                puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 29],"
                puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 30],"
                puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 31],"
                puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 32],"
                puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 33],"
                puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 34],"
                puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 35],"
                puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 36],"
                puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 37],"
                puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 38],"
                
                puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 39],"
                puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 40],"
                puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 41],"
                puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 42],"
                puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 43],"
                puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 44],"
                puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 45],"
                puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 46],"
                puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 47],"
                puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 48],"
                
                puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 49],"
                puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 50],"
                puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 51],"
                puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 52],"
                puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 53],"
                puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 54],"
                puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 55],"
                puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 56],"
                puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 57],"
                puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 58],"
                
                puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 59],"
                puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 60],"
                puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 61],"
                puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 62],"
                puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 63],"
                puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 64],"
                puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 65],"
                puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 66],"
                puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 67],"
                puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 68],"
                
                puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 69],"
                puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 70],"
                puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 71],"
                puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 72],"
                puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 73],"
                puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 74],"
                puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 75],"
                puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 76],"
                puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 77],"
                puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 78],"

                puts $FileVar "              \"NUMBER_OF_UMAT_STATE_VARIABLES\":  [lindex [lindex $Groups $i] 79],"

                puts $FileVar "              \"STATE_VARIABLE_1\"         :  [lindex [lindex $Groups $i] 80],"
                puts $FileVar "              \"STATE_VARIABLE_2\"         :  [lindex [lindex $Groups $i] 81],"
                puts $FileVar "              \"STATE_VARIABLE_3\"         :  [lindex [lindex $Groups $i] 82],"
                puts $FileVar "              \"STATE_VARIABLE_4\"         :  [lindex [lindex $Groups $i] 83],"
                puts $FileVar "              \"STATE_VARIABLE_5\"         :  [lindex [lindex $Groups $i] 84],"
                puts $FileVar "              \"STATE_VARIABLE_6\"         :  [lindex [lindex $Groups $i] 85],"
                puts $FileVar "              \"STATE_VARIABLE_7\"         :  [lindex [lindex $Groups $i] 86],"
                puts $FileVar "              \"STATE_VARIABLE_8\"         :  [lindex [lindex $Groups $i] 87],"
                puts $FileVar "              \"STATE_VARIABLE_9\"         :  [lindex [lindex $Groups $i] 88],"
                puts $FileVar "              \"STATE_VARIABLE_10\"         :  [lindex [lindex $Groups $i] 89],"
                
                puts $FileVar "              \"STATE_VARIABLE_11\"         :  [lindex [lindex $Groups $i] 90],"
                puts $FileVar "              \"STATE_VARIABLE_12\"         :  [lindex [lindex $Groups $i] 91],"
                puts $FileVar "              \"STATE_VARIABLE_13\"         :  [lindex [lindex $Groups $i] 92],"
                puts $FileVar "              \"STATE_VARIABLE_14\"         :  [lindex [lindex $Groups $i] 93],"
                puts $FileVar "              \"STATE_VARIABLE_15\"         :  [lindex [lindex $Groups $i] 94],"
                puts $FileVar "              \"STATE_VARIABLE_16\"         :  [lindex [lindex $Groups $i] 95],"
                puts $FileVar "              \"STATE_VARIABLE_17\"         :  [lindex [lindex $Groups $i] 96],"
                puts $FileVar "              \"STATE_VARIABLE_18\"         :  [lindex [lindex $Groups $i] 97],"
                puts $FileVar "              \"STATE_VARIABLE_19\"         :  [lindex [lindex $Groups $i] 98],"
                puts $FileVar "              \"STATE_VARIABLE_20\"         :  [lindex [lindex $Groups $i] 99],"
                
                puts $FileVar "              \"STATE_VARIABLE_21\"         :  [lindex [lindex $Groups $i] 100],"
                puts $FileVar "              \"STATE_VARIABLE_22\"         :  [lindex [lindex $Groups $i] 101],"
                puts $FileVar "              \"STATE_VARIABLE_23\"         :  [lindex [lindex $Groups $i] 102],"
                puts $FileVar "              \"STATE_VARIABLE_24\"         :  [lindex [lindex $Groups $i] 103],"
                puts $FileVar "              \"STATE_VARIABLE_25\"         :  [lindex [lindex $Groups $i] 104],"
                puts $FileVar "              \"STATE_VARIABLE_26\"         :  [lindex [lindex $Groups $i] 105],"
                puts $FileVar "              \"STATE_VARIABLE_27\"         :  [lindex [lindex $Groups $i] 106],"
                puts $FileVar "              \"STATE_VARIABLE_28\"         :  [lindex [lindex $Groups $i] 107],"
                puts $FileVar "              \"STATE_VARIABLE_29\"         :  [lindex [lindex $Groups $i] 108],"
                puts $FileVar "              \"STATE_VARIABLE_30\"         :  [lindex [lindex $Groups $i] 109],"
                
                puts $FileVar "              \"STATE_VARIABLE_31\"         :  [lindex [lindex $Groups $i] 110],"
                puts $FileVar "              \"STATE_VARIABLE_32\"         :  [lindex [lindex $Groups $i] 111],"
                puts $FileVar "              \"STATE_VARIABLE_33\"         :  [lindex [lindex $Groups $i] 112],"
                puts $FileVar "              \"STATE_VARIABLE_34\"         :  [lindex [lindex $Groups $i] 113],"
                puts $FileVar "              \"STATE_VARIABLE_35\"         :  [lindex [lindex $Groups $i] 114],"
                puts $FileVar "              \"STATE_VARIABLE_36\"         :  [lindex [lindex $Groups $i] 115],"
                puts $FileVar "              \"STATE_VARIABLE_37\"         :  [lindex [lindex $Groups $i] 116],"
                puts $FileVar "              \"STATE_VARIABLE_38\"         :  [lindex [lindex $Groups $i] 117],"
                puts $FileVar "              \"STATE_VARIABLE_39\"         :  [lindex [lindex $Groups $i] 118],"
                puts $FileVar "              \"STATE_VARIABLE_40\"         :  [lindex [lindex $Groups $i] 119],"
                
                puts $FileVar "              \"STATE_VARIABLE_41\"         :  [lindex [lindex $Groups $i] 120],"
                puts $FileVar "              \"STATE_VARIABLE_42\"         :  [lindex [lindex $Groups $i] 121],"
                puts $FileVar "              \"STATE_VARIABLE_43\"         :  [lindex [lindex $Groups $i] 122],"
                puts $FileVar "              \"STATE_VARIABLE_44\"         :  [lindex [lindex $Groups $i] 123],"
                puts $FileVar "              \"STATE_VARIABLE_45\"         :  [lindex [lindex $Groups $i] 124],"
                puts $FileVar "              \"STATE_VARIABLE_46\"         :  [lindex [lindex $Groups $i] 125],"
                puts $FileVar "              \"STATE_VARIABLE_47\"         :  [lindex [lindex $Groups $i] 126],"
                puts $FileVar "              \"STATE_VARIABLE_48\"         :  [lindex [lindex $Groups $i] 127],"
                puts $FileVar "              \"STATE_VARIABLE_49\"         :  [lindex [lindex $Groups $i] 128],"
                puts $FileVar "              \"STATE_VARIABLE_50\"         :  [lindex [lindex $Groups $i] 129]"
            
                puts $FileVar "         \},"
                puts $FileVar "         \"Tables\": \{\}"
                puts $FileVar "      \}"
            }
        }
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
       puts $FileVar "              \"YOUNG_MODULUS\"    :  [lindex [lindex $Groups $i] 4],"
       puts $FileVar "              \"POISSON_RATIO\"    :  [lindex [lindex $Groups $i] 5],"
       puts $FileVar "              \"DENSITY\"          :  [lindex [lindex $Groups $i] 6],"
       puts $FileVar "              \"CROSS_AREA\"       :  [lindex [lindex $Groups $i] 7],"
       puts $FileVar "              \"I33\"              :  [lindex [lindex $Groups $i] 8]"
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

    # Cable part
    set Groups [GiD_Info conditions Cable groups]
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

    # Interface drained and interface undrained part
	set interface_Groups [list   [GiD_Info conditions Interface_drained groups] [GiD_Info conditions Interface_undrained groups]]
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
				puts $FileVar "              \"DENSITY_SOLID\"            :  [lindex [lindex $Groups $i] 7],"
				puts $FileVar "              \"DENSITY_WATER\"            :  [lindex [lindex $Groups $i] 8],"
				puts $FileVar "              \"POROSITY\"                 :  [lindex [lindex $Groups $i] 9],"
				puts $FileVar "              \"BULK_MODULUS_SOLID\"       :  [lindex [lindex $Groups $i] 10],"
				puts $FileVar "              \"BULK_MODULUS_FLUID\"       :  [lindex [lindex $Groups $i] 11],"
				puts $FileVar "              \"TRANSVERSAL_PERMEABILITY\" :  [lindex [lindex $Groups $i] 12],"
				puts $FileVar "              \"DYNAMIC_VISCOSITY\"        :  [lindex [lindex $Groups $i] 13],"
				puts $FileVar "              \"MINIMUM_JOINT_WIDTH\"      :  [lindex [lindex $Groups $i] 16],"

				puts $FileVar "              \"UDSM_NAME\"                :  \"[lindex [lindex $Groups $i] 20]\","
				puts $FileVar "              \"UDSM_NUMBER\"              :  [lindex [lindex $Groups $i] 21],"
				puts $FileVar "              \"IS_FORTRAN_UDSM\"          :  [lindex [lindex $Groups $i] 22],"

				puts $FileVar "              \"PARAMETER_1\"              :  [lindex [lindex $Groups $i] 27],"
				puts $FileVar "              \"PARAMETER_2\"              :  [lindex [lindex $Groups $i] 28],"
				puts $FileVar "              \"PARAMETER_3\"              :  [lindex [lindex $Groups $i] 29],"
				puts $FileVar "              \"PARAMETER_4\"              :  [lindex [lindex $Groups $i] 30],"
				puts $FileVar "              \"PARAMETER_5\"              :  [lindex [lindex $Groups $i] 31],"
				puts $FileVar "              \"PARAMETER_6\"              :  [lindex [lindex $Groups $i] 32],"
				puts $FileVar "              \"PARAMETER_7\"              :  [lindex [lindex $Groups $i] 33],"
				puts $FileVar "              \"PARAMETER_8\"              :  [lindex [lindex $Groups $i] 34],"
				puts $FileVar "              \"PARAMETER_9\"              :  [lindex [lindex $Groups $i] 35],"
				puts $FileVar "              \"PARAMETER_10\"              :  [lindex [lindex $Groups $i] 36],"
				
				puts $FileVar "              \"PARAMETER_11\"              :  [lindex [lindex $Groups $i] 37],"
				puts $FileVar "              \"PARAMETER_12\"              :  [lindex [lindex $Groups $i] 38],"
				puts $FileVar "              \"PARAMETER_13\"              :  [lindex [lindex $Groups $i] 39],"
				puts $FileVar "              \"PARAMETER_14\"              :  [lindex [lindex $Groups $i] 40],"
				puts $FileVar "              \"PARAMETER_15\"              :  [lindex [lindex $Groups $i] 41],"
				puts $FileVar "              \"PARAMETER_16\"              :  [lindex [lindex $Groups $i] 42],"
				puts $FileVar "              \"PARAMETER_17\"              :  [lindex [lindex $Groups $i] 43],"
				puts $FileVar "              \"PARAMETER_18\"              :  [lindex [lindex $Groups $i] 44],"
				puts $FileVar "              \"PARAMETER_19\"              :  [lindex [lindex $Groups $i] 45],"
				puts $FileVar "              \"PARAMETER_20\"              :  [lindex [lindex $Groups $i] 46],"
				
				puts $FileVar "              \"PARAMETER_21\"              :  [lindex [lindex $Groups $i] 47],"
				puts $FileVar "              \"PARAMETER_22\"              :  [lindex [lindex $Groups $i] 48],"
				puts $FileVar "              \"PARAMETER_23\"              :  [lindex [lindex $Groups $i] 49],"
				puts $FileVar "              \"PARAMETER_24\"              :  [lindex [lindex $Groups $i] 50],"
				puts $FileVar "              \"PARAMETER_25\"              :  [lindex [lindex $Groups $i] 51],"
				puts $FileVar "              \"PARAMETER_26\"              :  [lindex [lindex $Groups $i] 52],"
				puts $FileVar "              \"PARAMETER_27\"              :  [lindex [lindex $Groups $i] 53],"
				puts $FileVar "              \"PARAMETER_28\"              :  [lindex [lindex $Groups $i] 54],"
				puts $FileVar "              \"PARAMETER_29\"              :  [lindex [lindex $Groups $i] 55],"
				puts $FileVar "              \"PARAMETER_30\"              :  [lindex [lindex $Groups $i] 56],"
				
				puts $FileVar "              \"PARAMETER_31\"              :  [lindex [lindex $Groups $i] 57],"
				puts $FileVar "              \"PARAMETER_32\"              :  [lindex [lindex $Groups $i] 58],"
				puts $FileVar "              \"PARAMETER_33\"              :  [lindex [lindex $Groups $i] 59],"
				puts $FileVar "              \"PARAMETER_34\"              :  [lindex [lindex $Groups $i] 60],"
				puts $FileVar "              \"PARAMETER_35\"              :  [lindex [lindex $Groups $i] 61],"
				puts $FileVar "              \"PARAMETER_36\"              :  [lindex [lindex $Groups $i] 62],"
				puts $FileVar "              \"PARAMETER_37\"              :  [lindex [lindex $Groups $i] 63],"
				puts $FileVar "              \"PARAMETER_38\"              :  [lindex [lindex $Groups $i] 64],"
				puts $FileVar "              \"PARAMETER_39\"              :  [lindex [lindex $Groups $i] 65],"
				puts $FileVar "              \"PARAMETER_40\"              :  [lindex [lindex $Groups $i] 66],"
				
				puts $FileVar "              \"PARAMETER_41\"              :  [lindex [lindex $Groups $i] 67],"
				puts $FileVar "              \"PARAMETER_42\"              :  [lindex [lindex $Groups $i] 68],"
				puts $FileVar "              \"PARAMETER_43\"              :  [lindex [lindex $Groups $i] 69],"
				puts $FileVar "              \"PARAMETER_44\"              :  [lindex [lindex $Groups $i] 70],"
				puts $FileVar "              \"PARAMETER_45\"              :  [lindex [lindex $Groups $i] 71],"
				puts $FileVar "              \"PARAMETER_46\"              :  [lindex [lindex $Groups $i] 72],"
				puts $FileVar "              \"PARAMETER_47\"              :  [lindex [lindex $Groups $i] 73],"
				puts $FileVar "              \"PARAMETER_48\"              :  [lindex [lindex $Groups $i] 74],"
				puts $FileVar "              \"PARAMETER_49\"              :  [lindex [lindex $Groups $i] 75],"
				puts $FileVar "              \"PARAMETER_50\"              :  [lindex [lindex $Groups $i] 76]"
				
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

