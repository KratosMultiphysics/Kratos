proc AppendGroupNames {String CondName} {
    upvar $String MyString

    set Groups [GiD_Info conditions $CondName groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append MyString \" [lindex [lindex $Groups $i] 1] \" ,
    }
}

#-------------------------------------------------------------------------------

proc AppendOutputVariables {String GroupNum QuestionName VarName} {
    upvar $String MyString
    upvar $GroupNum MyGroupNum

    if {[GiD_AccessValue get gendata $QuestionName] eq true} {
        incr MyGroupNum
        append MyString \" $VarName \" ,
    }
}

#-------------------------------------------------------------------------------

proc WriteConstraintVectorProcess {FileVar GroupNum Groups EntityType VarName TableDict NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]

        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "        \"python_module\": \"apply_vector_constraint_table_process\","
            puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.FluidTransportApplication\","
            puts $MyFileVar "        \"process_name\":  \"ApplyVectorConstraintTableProcess\","
            puts $MyFileVar "        \"Parameters\":    \{"
            puts $MyFileVar "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"variable_name\":   \"$VarName\","
            puts $MyFileVar "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 13]\],"
            puts $MyFileVar "            \"is_fixed\":        \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 10],[lindex [lindex $Groups $i] 15]\],"
            puts $MyFileVar "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
            puts $MyFileVar "            \"table\":           \[[dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table1],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table2]\]"
            puts $MyFileVar "        \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "    \},\{"
            } else {
                puts $MyFileVar "    \}\],"
            }
        }
    }
}

#-------------------------------------------------------------------------------

proc WritePressureConstraintProcess {FileVar GroupNum Groups EntityType VarName TableDict NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum


    for {set i 0} {$i < [llength $Groups]} {incr i} {

        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "        \"python_module\": \"apply_scalar_constraint_table_process\","
            puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.FluidTransportApplication\","
            puts $MyFileVar "        \"process_name\":  \"ApplyScalarConstraintTableProcess\","
            puts $MyFileVar "        \"Parameters\":    \{"
            puts $MyFileVar "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"variable_name\":        \"$VarName\","
            puts $MyFileVar "            \"is_fixed\":             [lindex [lindex $Groups $i] 4],"
            puts $MyFileVar "            \"value\":                [lindex [lindex $Groups $i] 3],"
            puts $MyFileVar "            \"table\":                [dict get $TableDict [lindex [lindex $Groups $i] 1] Table0]"
            puts $MyFileVar "        \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "    \},\{"
            } else {
                puts $MyFileVar "    \}\],"
            }
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteTempConstraintProcess {FileVar GroupNum Groups VarName TableDict NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum


    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr MyGroupNum
        puts $MyFileVar "        \"python_module\": \"apply_scalar_constraint_table_process\","
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.FluidTransportApplication\","
        puts $MyFileVar "        \"process_name\":  \"ApplyScalarConstraintTableProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        puts $MyFileVar "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
        puts $MyFileVar "            \"variable_name\":        \"$VarName\","
        puts $MyFileVar "            \"value\":                [lindex [lindex $Groups $i] 3],"
        puts $MyFileVar "            \"table\":                [dict get $TableDict [lindex [lindex $Groups $i] 1] Table0]"
        puts $MyFileVar "        \}"
        if {$MyGroupNum < $NumGroups} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\],"
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteLoadScalarProcess {FileVar GroupNum Groups VarName TableDict NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr MyGroupNum
        puts $MyFileVar "        \"python_module\": \"apply_scalar_constraint_table_process\","
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.FluidTransportApplication\","
        puts $MyFileVar "        \"process_name\":  \"ApplyScalarConstraintTableProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        puts $MyFileVar "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
        puts $MyFileVar "            \"variable_name\":   \"$VarName\","
        puts $MyFileVar "            \"value\":           [lindex [lindex $Groups $i] 3],"
        puts $MyFileVar "            \"table\":           [dict get $TableDict [lindex [lindex $Groups $i] 1] Table0]"
        puts $MyFileVar "        \}"
        if {$MyGroupNum < $NumGroups} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\]"
        }
    }
}