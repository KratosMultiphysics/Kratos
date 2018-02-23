
proc WriteInitialWaterLevelProcess {FileVar GroupNum Groups EntityType NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "        \"python_module\"   : \"set_initial_water_level_process\","
            puts $MyFileVar "        \"kratos_module\"   : \"KratosMultiphysics.ShallowWaterApplication\","
            puts $MyFileVar "        \"process_name\"    : \"SetInitialWaterLevelProcess\","
            puts $MyFileVar "        \"Parameters\"      : \{"
            puts $MyFileVar "            \"mesh_id\"         : 0,"
            puts $MyFileVar "            \"model_part_name\" : \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"variable_name\"   : \"[lindex [lindex $Groups $i] 3]\","
            puts $MyFileVar "            \"value\"           : [lindex [lindex $Groups $i] 4]"
            puts $MyFileVar "        \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "    \},\{"
            } else {
                puts $MyFileVar "    \}\],"
            }
        }
    }
}