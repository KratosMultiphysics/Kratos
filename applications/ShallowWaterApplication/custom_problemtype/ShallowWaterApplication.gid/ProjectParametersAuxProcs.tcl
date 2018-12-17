
proc WriteInitialWaterLevelProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\"   : \"set_initial_water_level_process\","
            puts $MyFileVar "            \"kratos_module\"   : \"KratosMultiphysics.ShallowWaterApplication\","
            puts $MyFileVar "            \"process_name\"    : \"SetInitialWaterLevelProcess\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"mesh_id\"         : 0,"
            puts $MyFileVar "                \"model_part_name\" : \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"variable_name\"   : \"[lindex [lindex $Groups $i] 3]\","
            puts $MyFileVar "                \"value\"           : \"[lindex [lindex $Groups $i] 4]\""
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            } else {
                puts $MyFileVar "        \}\],"
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
            puts $MyFileVar "            \"process_name\"    : \"ApplySlipProcess\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"model_part_name\" : \"[lindex [lindex $Groups $i] 1]\""
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            } else {
                puts $MyFileVar "        \}\],"
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
            puts $MyFileVar "                \"model_part_name\" : \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"variable_name\"   : \"HEIGHT\","
            puts $MyFileVar "                \"value\"           : [lindex [lindex $Groups $i] 3],"
            puts $MyFileVar "                \"is_fixed\"        : [lindex [lindex $Groups $i] 4]"
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            } else {
                puts $MyFileVar "        \}\],"
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
            puts $MyFileVar "                \"model_part_name\" : \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"variable_name\"   : \"[lindex [lindex $Groups $i] 3]\","
            puts $MyFileVar "                \"modulus\"         : [lindex [lindex $Groups $i] 4],"
            set pi 3.1415926535897931
            set alpha [expr {[lindex [lindex $Groups $i] 5] * $pi / 180 }]
            set ux [expr {cos($alpha)}]
            set uy [expr {sin($alpha)}]
            puts $MyFileVar "                \"direction\"       : \[$ux, $uy, 0.0\],"
            puts $MyFileVar "                \"is_fixed_x\"      : [lindex [lindex $Groups $i] 6],"
            puts $MyFileVar "                \"is_fixed_y\"      : [lindex [lindex $Groups $i] 7]"
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            } else {
                puts $MyFileVar "        \}\],"
            }
        }
    }
}

proc WriteBathymetryProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\"   : \"set_bathymetry_process\","
            puts $MyFileVar "            \"kratos_module\"   : \"KratosMultiphysics.ShallowWaterApplication\","
            puts $MyFileVar "            \"process_name\"    : \"SetBathymetryProcess\","
            puts $MyFileVar "            \"Parameters\"      : \{"
            puts $MyFileVar "                \"model_part_name\" : \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"variable_name\"   : \"BATHYMETRY\","
            if {[lindex [lindex $Groups $i] 3] eq "From_digital_model"} {
                puts $MyFileVar "                \"value\"           : \"z\""
            } elseif {[lindex [lindex $Groups $i] 3] eq "Expression"} {
                puts $MyFileVar "                \"value\"           : \"[lindex [lindex $Groups $i] 4]\""
            } else {
                puts $MyFileVar "                \"value\"           : \"0.0\""
            }
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            } else {
                puts $MyFileVar "        \}\]"
            }
        }
    }
}
