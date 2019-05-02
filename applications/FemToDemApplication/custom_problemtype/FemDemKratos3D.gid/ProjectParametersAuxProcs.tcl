proc AppendGroupName {String CondName} {
    upvar $String MyString
    
    set Groups [GiD_Info conditions $CondName groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append MyString \" [lindex [lindex $Groups $i] 1] \" 
    }
}

proc AppendGroupNames {String CondName} {
    upvar $String MyString
    
    set Groups [GiD_Info conditions $CondName groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append MyString \" [lindex [lindex $Groups $i] 1] \" ,
    }
}

#-------------------------------------------------------------------------------

proc AppendGroupNamesWithNum {String GroupNum CondName} {
    upvar $String MyString
    upvar $GroupNum MyGroupNum
    
    set Groups [GiD_Info conditions $CondName groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr MyGroupNum
        append MyString \" [lindex [lindex $Groups $i] 1] \" ,
    }
}

#-------------------------------------------------------------------------------

proc AppendGroupVariables {String CondName VarName} {
    upvar $String MyString
    
    set Groups [GiD_Info conditions $CondName groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append MyString \" $VarName \" ,
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
            puts $MyFileVar "        \"python_module\": \"assign_vector_components_to_nodes_process\","
            puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.SolidMechanicsApplication\","
            puts $MyFileVar "        \"process_name\":  \"AssignVectorComponentsToNodesProcess\","
            puts $MyFileVar "        \"Parameters\":    \{"
            #puts $MyFileVar "            \"mesh_id\":         0,"
            puts $MyFileVar "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"variable_name\":   \"$VarName\","
            #puts $MyFileVar "            \"constrained\":     true,"
            puts $MyFileVar "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
            puts $MyFileVar "            \"interval\":        \[[lindex [lindex $Groups $i] 18],[lindex [lindex $Groups $i] 19]\]"
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
            puts $MyFileVar "        \"python_module\": \"apply_pore_pressure_table_process\","
            puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $MyFileVar "        \"process_name\":  \"ApplyPorePressureTableProcess\","
            puts $MyFileVar "        \"Parameters\":    \{"
            puts $MyFileVar "            \"mesh_id\":              0,"
            puts $MyFileVar "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"variable_name\":        \"$VarName\","
            puts $MyFileVar "            \"is_fixed\":             [lindex [lindex $Groups $i] 8],"
            puts $MyFileVar "            \"value\":                [lindex [lindex $Groups $i] 4],"
            puts $MyFileVar "            \"table\":                [dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],"
            if {[lindex [lindex $Groups $i] 3] eq "Hydrostatic"} {
                set PutStrings true
            } else {
                set PutStrings false
            }
            puts $MyFileVar "            \"hydrostatic\":          $PutStrings,"
            if {[lindex [lindex $Groups $i] 5] eq "Y"} {
                set PutStrings 2
            } elseif {[lindex [lindex $Groups $i] 5] eq "Z"} {
                set PutStrings 3
            } else {
                set PutStrings 1
            }
            puts $MyFileVar "            \"gravity_direction\":    $PutStrings,"
            puts $MyFileVar "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
            puts $MyFileVar "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
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

proc WriteLoadVectorProcess {FileVar GroupNum Groups VarName TableDict NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr MyGroupNum
        puts $MyFileVar "        \"python_module\": \"assign_modulus_and_direction_to_conditions_process\","
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.SolidMechanicsApplication\","
        puts $MyFileVar "        \"process_name\":  \"AssignModulusAndDirectionToConditionsProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        #puts $MyFileVar "            \"mesh_id\":         0,"
        puts $MyFileVar "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
        puts $MyFileVar "            \"variable_name\":   \"$VarName\","
        puts $MyFileVar "            \"modulus\":             [lindex [lindex $Groups $i] 3],"
        puts $MyFileVar "            \"direction\":           \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 13]\],"
        puts $MyFileVar "            \"interval\":            \[[lindex [lindex $Groups $i] 16],[lindex [lindex $Groups $i] 17]\]"
        puts $MyFileVar "        \}"
        if {$MyGroupNum < $NumGroups} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\],"
        }
    }
}
#-------------------------------------------------------------------------------

proc WriteGLoadVectorProcess {FileVar GroupNum Groups VarName TableDict NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr MyGroupNum
        puts $MyFileVar "        \"python_module\": \"assign_modulus_and_direction_to_nodes_process\","
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.SolidMechanicsApplication\","
        puts $MyFileVar "        \"process_name\":  \"AssignModulusAndDirectionToNodesProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        #puts $MyFileVar "            \"mesh_id\":         0,"
        puts $MyFileVar "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
        puts $MyFileVar "            \"variable_name\":   \"$VarName\","
        puts $MyFileVar "            \"modulus\":             [lindex [lindex $Groups $i] 3],"
        puts $MyFileVar "            \"direction\":           \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 13]\],"
        puts $MyFileVar "            \"interval\":            \[[lindex [lindex $Groups $i] 16],[lindex [lindex $Groups $i] 17]\]"
        puts $MyFileVar "        \}"

        if {$MyGroupNum < $NumGroups} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\],"
        }
    }
}
#-------------------------------------------------------------------------------

proc WriteNormalLoadProcess {FileVar GroupNum Groups VarName TableDict NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr MyGroupNum
        puts $MyFileVar "        \"python_module\": \"assign_scalar_to_conditions_process\","
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.SolidMechanicsApplication\","
        puts $MyFileVar "        \"process_name\":  \"AssignScalarToConditionsProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        #puts $MyFileVar "            \"mesh_id\":              0,"
        puts $MyFileVar "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
        puts $MyFileVar "            \"variable_name\":        \"$VarName\","
        puts $MyFileVar "            \"value\":                 [lindex [lindex $Groups $i] 5],"
        puts $MyFileVar "            \"interval\":             \[[lindex [lindex $Groups $i] 15],[lindex [lindex $Groups $i] 16]\]"
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
        puts $MyFileVar "        \"python_module\": \"apply_load_scalar_table_process\","
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
        puts $MyFileVar "        \"process_name\":  \"ApplyLoadScalarTableProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        puts $MyFileVar "            \"mesh_id\":         0,"
        puts $MyFileVar "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
        puts $MyFileVar "            \"variable_name\":   \"$VarName\","
        puts $MyFileVar "            \"value\":           [lindex [lindex $Groups $i] 3],"
        if {[GiD_AccessValue get gendata Strategy_Type] eq "Arc-Length"} {
            puts $MyFileVar "            \"table\":           0"
        } else {
            puts $MyFileVar "            \"table\":           [dict get $TableDict [lindex [lindex $Groups $i] 1] Table0]"
        }
        puts $MyFileVar "        \}"
        if {$MyGroupNum < $NumGroups} {
            puts $MyFileVar "    \},\{"
        } else {
            puts $MyFileVar "    \}\]"
        }
    }
}

#-------------------------------------------------------------------------------

proc WritePeriodicInterfaceProcess {FileVar GroupNum Groups NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 20] eq true} {
            incr MyGroupNum
            puts $MyFileVar "        \"python_module\": \"periodic_interface_activation_process\","
            puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $MyFileVar "        \"process_name\":  \"PeriodicInterfaceActivationProcess\","
            puts $MyFileVar "        \"Parameters\":    \{"
            puts $MyFileVar "            \"mesh_id\":         0,"
            puts $MyFileVar "            \"model_part_name\": \"Periodic_Bars_[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"dimension\":       [GiD_AccessValue get gendata Domain_Size],"
            puts $MyFileVar "            \"von_mises_limit\": [lindex [lindex $Groups $i] 21]"
            puts $MyFileVar "        \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "    \},\{"
            } else {
                puts $MyFileVar "    \}\]"
            }
        }
    }
}