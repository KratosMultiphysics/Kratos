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
            puts $MyFileVar "        \"python_module\": \"apply_vector_constraint_table_process\","
            puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
            puts $MyFileVar "        \"process_name\":  \"ApplyVectorConstraintTableProcess\","
            puts $MyFileVar "        \"Parameters\":    \{"
            puts $MyFileVar "            \"model_part_name\": \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"variable_name\":   \"$VarName\","
            puts $MyFileVar "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 13]\],"
            puts $MyFileVar "            \"is_fixed\":        \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 10],[lindex [lindex $Groups $i] 15]\],"
            puts $MyFileVar "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
            puts $MyFileVar "            \"table\":           \[[dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table1],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table2]\]"
            puts $MyFileVar "        \}"
            
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "    \},\{"
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
            if {[lindex [lindex $Groups $i] 3] ne "Interpolate_Line"} {
                incr MyGroupNum
                puts $MyFileVar "        \"python_module\": \"apply_scalar_constraint_table_process\","
                puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
                puts $MyFileVar "        \"process_name\":  \"ApplyScalarConstraintTableProcess\","
                puts $MyFileVar "        \"Parameters\":    \{"
                puts $MyFileVar "            \"model_part_name\":      \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $MyFileVar "            \"variable_name\":        \"$VarName\","
                puts $MyFileVar "            \"is_fixed\":             [lindex [lindex $Groups $i] 22],"
                puts $MyFileVar "            \"fluid_pressure_type\": \"[lindex [lindex $Groups $i] 3]\","

                if {[lindex [lindex $Groups $i] 5] eq "Y"} {
                    set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 5] eq "Z"} {
                    set PutStrings 2
                } else {
                    set PutStrings 0
                }
                
                if {[lindex [lindex $Groups $i] 3] eq "Hydrostatic"} {
                    puts $MyFileVar "            \"gravity_direction\":    $PutStrings,"
                    puts $MyFileVar "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                    puts $MyFileVar "            \"table\":                [dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],"
                    puts $MyFileVar "            \"pressure_tension_cut_off\":      [lindex [lindex $Groups $i] 8],"
                    puts $MyFileVar "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                } elseif {[lindex [lindex $Groups $i] 3] eq "Uniform"} {
                    puts $MyFileVar "            \"value\":                [lindex [lindex $Groups $i] 4],"
                    puts $MyFileVar "            \"table\":                [dict get $TableDict [lindex [lindex $Groups $i] 1] Table0]"
                } elseif {[lindex [lindex $Groups $i] 3] eq "Phreatic_Line"} {
                    puts $MyFileVar "            \"gravity_direction\":    $PutStrings,"

                    if {[lindex [lindex $Groups $i] 9] eq "Y"} {
                        set PutStrings 1
                    } elseif {[lindex [lindex $Groups $i] 9] eq "Z"} {
                        set PutStrings 2
                    } else {
                        set PutStrings 0
                    }
                    puts $MyFileVar "            \"out_of_plane_direction\":    $PutStrings,"
                    puts $MyFileVar "            \"first_reference_coordinate\" :    \[[lindex [lindex $Groups $i] 11],[lindex [lindex $Groups $i] 12],[lindex [lindex $Groups $i] 13]\],"
                    puts $MyFileVar "            \"second_reference_coordinate\":    \[[lindex [lindex $Groups $i] 15],[lindex [lindex $Groups $i] 16],[lindex [lindex $Groups $i] 17]\],"
                    puts $MyFileVar "            \"table\":                \[[dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table1]\],"
                    puts $MyFileVar "            \"pressure_tension_cut_off\":      [lindex [lindex $Groups $i] 8],"
                    puts $MyFileVar "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                } elseif {[lindex [lindex $Groups $i] 3] eq "Phreatic_Surface"} {
                    puts $MyFileVar "            \"first_reference_coordinate\" :    \[[lindex [lindex $Groups $i] 11],[lindex [lindex $Groups $i] 12],[lindex [lindex $Groups $i] 13]\],"
                    puts $MyFileVar "            \"second_reference_coordinate\":    \[[lindex [lindex $Groups $i] 15],[lindex [lindex $Groups $i] 16],[lindex [lindex $Groups $i] 17]\],"
                    puts $MyFileVar "            \"third_reference_coordinate\" :    \[[lindex [lindex $Groups $i] 19],[lindex [lindex $Groups $i] 20],[lindex [lindex $Groups $i] 21]\],"
                    puts $MyFileVar "            \"table\":                \[[dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table1],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table2]\],"
                    puts $MyFileVar "            \"pressure_tension_cut_off\":      [lindex [lindex $Groups $i] 8],"
                    puts $MyFileVar "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                }

                puts $MyFileVar "        \}"
                if {$MyGroupNum < $NumGroups} {
                    puts $MyFileVar "    \},\{"
                }
        
            }
        }
    }

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            if {[lindex [lindex $Groups $i] 3] eq "Interpolate_Line"} {
                incr MyGroupNum
                puts $MyFileVar "        \"python_module\": \"apply_scalar_constraint_table_process\","
                puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
                puts $MyFileVar "        \"process_name\":  \"ApplyScalarConstraintTableProcess\","
                puts $MyFileVar "        \"Parameters\":    \{"
                puts $MyFileVar "            \"model_part_name\":      \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
                puts $MyFileVar "            \"variable_name\":        \"$VarName\","
                puts $MyFileVar "            \"is_fixed\":             [lindex [lindex $Groups $i] 22],"
                puts $MyFileVar "            \"table\":                [dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],"
                puts $MyFileVar "            \"fluid_pressure_type\": \"[lindex [lindex $Groups $i] 3]\","

                if {[lindex [lindex $Groups $i] 5] eq "Y"} {
                    set PutStrings 1
                } elseif {[lindex [lindex $Groups $i] 5] eq "Z"} {
                    set PutStrings 2
                } else {
                    set PutStrings 0
                }

                puts $MyFileVar "            \"gravity_direction\":    $PutStrings,"
                
                if {[lindex [lindex $Groups $i] 3] eq "Interpolate_Line"} {
                    if {[lindex [lindex $Groups $i] 9] eq "Y"} {
                        set PutStrings 1
                    } elseif {[lindex [lindex $Groups $i] 9] eq "Z"} {
                        set PutStrings 2
                    } else {
                        set PutStrings 0
                    }
                    puts $MyFileVar "            \"out_of_plane_direction\":    $PutStrings"
                }

                puts $MyFileVar "        \}"
                if {$MyGroupNum < $NumGroups} {
                    puts $MyFileVar "    \},\{"
                }
            }
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteExcavationConstraintProcess {FileVar GroupNum Groups EntityType VarName NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "        \"python_module\": \"apply_excavation_process\","
            puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
            puts $MyFileVar "        \"process_name\":  \"ApplyExcavationProcess\","
            puts $MyFileVar "        \"Parameters\":    \{"
            puts $MyFileVar "            \"model_part_name\":      \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"variable_name\":        \"$VarName\","
            puts $MyFileVar "            \"deactivate_soil_part\":             [lindex [lindex $Groups $i] 3]"
            puts $MyFileVar "        \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "    \},\{"
            } 
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteResultVectorProcess {FileVar GroupNum Groups EntityType VarName NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum
    
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] $EntityType]
        if {[llength $Entities] > 0} {
            incr MyGroupNum
            puts $MyFileVar "        \"python_module\": \"apply_write_result_vector_process\","
            puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
            puts $MyFileVar "        \"process_name\":  \"ApplyWriteVectorProcess\","
            puts $MyFileVar "        \"Parameters\":    \{"
            puts $MyFileVar "            \"model_part_name\": \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "            \"variable_name\":   \"$VarName\","
            puts $MyFileVar "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 5]\],"
            puts $MyFileVar "            \"append_file\":     [lindex [lindex $Groups $i] 6]"
            puts $MyFileVar "        \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "    \},\{"
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
        puts $MyFileVar "        \"python_module\": \"apply_vector_constraint_table_process\","
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
        puts $MyFileVar "        \"process_name\":  \"ApplyVectorConstraintTableProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        puts $MyFileVar "            \"model_part_name\": \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
        puts $MyFileVar "            \"variable_name\":   \"$VarName\","
        puts $MyFileVar "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 7],[lindex [lindex $Groups $i] 11]\],"
        puts $MyFileVar "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
        if {[GiD_AccessValue get gendata Strategy_Type] eq "Arc-Length"} {
            puts $MyFileVar "            \"table\":           \[0,0,0\]"
        } else {
            puts $MyFileVar "            \"table\":           \[[dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table1],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table2]\]"
        }
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
        puts $MyFileVar "        \"python_module\": \"apply_normal_load_table_process\","
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
        puts $MyFileVar "        \"process_name\":  \"ApplyNormalLoadTableProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        puts $MyFileVar "            \"model_part_name\":      \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
        puts $MyFileVar "            \"variable_name\":        \"$VarName\","
        puts $MyFileVar "            \"active\":               \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 11]\],"
        puts $MyFileVar "            \"value\":                \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 12]\],"
        if {[GiD_AccessValue get gendata Strategy_Type] eq "Arc-Length"} {
            puts $MyFileVar "            \"table\":                \[0,0\],"
        } else {
            puts $MyFileVar "            \"table\":                \[[dict get $TableDict [lindex [lindex $Groups $i] 1] Table0],[dict get $TableDict [lindex [lindex $Groups $i] 1] Table1]\],"
        }
        puts $MyFileVar "            \"fluid_pressure_type\": \"[lindex [lindex $Groups $i] 4]\","

        if {[lindex [lindex $Groups $i] 6] eq "Y"} {
            set PutStrings 1
        } elseif {[lindex [lindex $Groups $i] 6] eq "Z"} {
            set PutStrings 2
        } else {
            set PutStrings 0
        }
        puts $MyFileVar "            \"gravity_direction\":    $PutStrings,"
        puts $MyFileVar "            \"reference_coordinate\": [lindex [lindex $Groups $i] 7],"
        puts $MyFileVar "            \"specific_weight\":      [lindex [lindex $Groups $i] 8]"
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
        puts $MyFileVar "        \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
        puts $MyFileVar "        \"process_name\":  \"ApplyScalarConstraintTableProcess\","
        puts $MyFileVar "        \"Parameters\":    \{"
        puts $MyFileVar "            \"model_part_name\": \"PorousDomain.[lindex [lindex $Groups $i] 1]\","
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
            puts $MyFileVar "    \}\],"
        }
    }
}

#-------------------------------------------------------------------------------
proc WriteGapClosureInterfaceProcess {FileVar GroupNum Groups NumGroups} {
    upvar $FileVar MyFileVar
    upvar $GroupNum MyGroupNum

    for {set i 0} {$i < [llength $Groups]} {incr i} {
        if {[lindex [lindex $Groups $i] 135] eq true} {
            incr MyGroupNum
            puts $MyFileVar "            \"python_module\": \"gap_closure_interface_activation_process\","
            puts $MyFileVar "            \"kratos_module\": \"KratosMultiphysics.GeoMechanicsApplication\","
            puts $MyFileVar "            \"Parameters\":    \{"
            puts $MyFileVar "                \"model_part_name\": \"PorousDomain.Gap_Closure_Bars_[lindex [lindex $Groups $i] 1]\","
            puts $MyFileVar "                \"consider_gap_closure\": [lindex [lindex $Groups $i] 135],"
            puts $MyFileVar "                \"gap_width_threshold\": [lindex [lindex $Groups $i] 136]"
            puts $MyFileVar "            \}"
            if {$MyGroupNum < $NumGroups} {
                puts $MyFileVar "        \},\{"
            } else {
                puts $MyFileVar "        \}\]"
            }
        }
    }
}