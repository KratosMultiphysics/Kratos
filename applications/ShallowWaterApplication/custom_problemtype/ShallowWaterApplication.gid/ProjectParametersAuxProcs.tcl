
proc WriteGiDOutputProcess {FileVar ModelPartName FileName OutputVars} {
    upvar $FileVar MyFileVar
    set OutputVars [string trimright $OutputVars ,]

    if {[GiD_AccessValue get gendata Multi_file_flag] eq "SingleFile"} {
        set NodalNonhistoricalVars \"DISPLACEMENT\"
    } elseif {[GiD_AccessValue get gendata Multi_file_flag] eq "MultipleFiles"} {
        set NodalNonhistoricalVars ""
    } else {
        set NodalNonhistoricalVars "something_went_wrong_in_the_problemtype"
    }

    puts $MyFileVar "            \"kratos_module\"        : \"KratosMultiphysics\","
    puts $MyFileVar "            \"python_module\"        : \"gid_output_process\","
    puts $MyFileVar "            \"Parameters\"           : \{"
    puts $MyFileVar "                \"model_part_name\"        : \"$ModelPartName\","
    puts $MyFileVar "                \"output_name\"            : \"$FileName\","
    puts $MyFileVar "                \"postprocess_parameters\" : \{"
    puts $MyFileVar "                    \"result_file_configuration\" : \{"
    puts $MyFileVar "                        \"gidpost_flags\"         : \{"
    puts $MyFileVar "                            \"GiDPostMode\"           : \"[GiD_AccessValue get gendata GiD_post_mode]\","
    puts $MyFileVar "                            \"WriteDeformedMeshFlag\" : \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $MyFileVar "                            \"WriteConditionsFlag\"   : \"[GiD_AccessValue get gendata Write_conditions]\","
    puts $MyFileVar "                            \"MultiFileFlag\"         : \"[GiD_AccessValue get gendata Multi_file_flag]\""
    puts $MyFileVar "                        \},"
    puts $MyFileVar "                        \"output_control_type\"   : \"[GiD_AccessValue get gendata Output_control_type]\","
    puts $MyFileVar "                        \"output_interval\"       : [GiD_AccessValue get gendata Output_interval],"
    puts $MyFileVar "                        \"body_output\"           : true,"
    puts $MyFileVar "                        \"node_output\"           : false,"
    puts $MyFileVar "                        \"nodal_results\"         : \[$OutputVars\],"
    puts $MyFileVar "                        \"gauss_point_results\"   : \[\],"
    puts $MyFileVar "                        \"nodal_nonhistorical_results\" : \[$NodalNonhistoricalVars\]"
    puts $MyFileVar "                    \}"
    puts $MyFileVar "                \}"
    puts $MyFileVar "            \}"
}

proc WriteTopographyProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\"   : \"set_topography_process\","
            puts $MyFileVar "            \"kratos_module\"   : \"KratosMultiphysics.ShallowWaterApplication\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"model_part_name\" : \"model_part.[lindex [lindex $Groups $i] 1]\","
            if {[lindex [lindex $Groups $i] 3] eq "From_digital_model"} {
                puts $MyFileVar "                \"value\"           : \"z\""
            } elseif {[lindex [lindex $Groups $i] 3] eq "By_function"} {
                puts $MyFileVar "                \"value\"           : \"[lindex [lindex $Groups $i] 4]\""
            } else {
                puts $MyFileVar "                \"value\"           : \"0.0\""
            }
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            }
        }
    }
}

proc WriteBottomFrictionProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\"   : \"process_factory\","
            puts $MyFileVar "            \"kratos_module\"   : \"KratosMultiphysics\","
            puts $MyFileVar "            \"process_name\"    : \"ApplyConstantScalarValueProcess\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"model_part_name\" : \"model_part.[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"variable_name\"   : \"MANNING\","
            puts $MyFileVar "                \"value\"           : [lindex [lindex $Groups $i] 3],"
            puts $MyFileVar "                \"is_fixed\"        : \"false\""
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            }
        }
    }
}

proc WriteInitialWaterLevelProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\"   : \"set_initial_water_level_process\","
            puts $MyFileVar "            \"kratos_module\"   : \"KratosMultiphysics.ShallowWaterApplication\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"model_part_name\" : \"model_part.[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"variable_name\"   : \"[lindex [lindex $Groups $i] 3]\","
            puts $MyFileVar "                \"value\"           : \"[lindex [lindex $Groups $i] 4]\""
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            }
        }
    }
}

proc WriteSlipConditionProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\"   : \"apply_slip_process\","
            puts $MyFileVar "            \"kratos_module\"   : \"KratosMultiphysics.ShallowWaterApplication\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"model_part_name\" : \"model_part.[lindex [lindex $Groups $i] 1]\""
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            }
        }
    }
}

proc WriteConstantScalarConditionProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\"   : \"process_factory\","
            puts $MyFileVar "            \"kratos_module\"   : \"KratosMultiphysics\","
            puts $MyFileVar "            \"process_name\"    : \"ApplyConstantScalarValueProcess\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"model_part_name\" : \"model_part.[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"variable_name\"   : \"HEIGHT\","
            puts $MyFileVar "                \"value\"           : [lindex [lindex $Groups $i] 3],"
            puts $MyFileVar "                \"is_fixed\"        : [lindex [lindex $Groups $i] 4]"
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            }
        }
    }
}

proc WriteConstantVectorConditionProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\"   : \"process_factory\","
            puts $MyFileVar "            \"kratos_module\"   : \"KratosMultiphysics\","
            puts $MyFileVar "            \"process_name\"    : \"ApplyConstantVectorValueProcess\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"model_part_name\" : \"model_part.[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"variable_name\"   : \"MOMENTUM\","
            puts $MyFileVar "                \"modulus\"         : [lindex [lindex $Groups $i] 3],"
            set pi 3.1415926535897931
            set alpha [expr {[lindex [lindex $Groups $i] 4] * $pi / 180 }]
            set ux [expr {cos($alpha)}]
            set uy [expr {sin($alpha)}]
            puts $MyFileVar "                \"direction\"       : \[$ux, $uy, 0.0\],"
            puts $MyFileVar "                \"is_fixed_x\"      : [lindex [lindex $Groups $i] 5],"
            puts $MyFileVar "                \"is_fixed_y\"      : [lindex [lindex $Groups $i] 6]"
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            }
        }
    }
}

proc WriteVisualizationMeshProcess {FileVar ModelPartName1 ModelPartName2} {
    upvar $FileVar MyFileVar

    if {[GiD_AccessValue get gendata Multi_file_flag] eq "SingleFile"} {
        set TheMeshDeformationMode "use_nodal_displacement"
    } elseif {[GiD_AccessValue get gendata Multi_file_flag] eq "MultipleFiles"} {
        set TheMeshDeformationMode "use_z_coordinate"
    } else {
        set TheMeshDeformationMode "something_went_wrong_in_the_problemtype"
    }

    puts $MyFileVar "            \"kratos_module\"        : \"KratosMultiphysics.ShallowWaterApplication\","
    puts $MyFileVar "            \"python_module\"        : \"visualization_mesh_process\","
    puts $MyFileVar "            \"Parameters\"           : \{"
    puts $MyFileVar "                \"model_part_name\"                : \"$ModelPartName1\","
    puts $MyFileVar "                \"topographic_model_part_name\"    : \"$ModelPartName2\","
    puts $MyFileVar "                \"create_topographic_model_part\"  : [GiD_AccessValue get gendata Print_topography_as_separate_output],"
    puts $MyFileVar "                \"use_properties_as_dry_wet_flag\" : false,"
    puts $MyFileVar "                \"mesh_deformation_mode\"          : \"$TheMeshDeformationMode\","
    puts $MyFileVar "                \"topography_variable\"            : \"TOPOGRAPHY\","
    puts $MyFileVar "                \"free_surface_variable\"          : \"FREE_SURFACE_ELEVATION\","
    puts $MyFileVar "                \"nodal_variables_to_transfer\"    : \[\"TOPOGRAPHY\"\]"
    puts $MyFileVar "            \}"
}
