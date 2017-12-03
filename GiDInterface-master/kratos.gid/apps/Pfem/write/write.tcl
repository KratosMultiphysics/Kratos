namespace eval Pfem::write {
    variable remesh_domains_dict
    variable bodies_list
}

proc Pfem::write::Init { } {
    variable remesh_domains_dict
    set remesh_domains [dict create ]
    variable bodies_list
    set bodies_list [list ]
    Solid::write::AddValidApps "Pfem"
}

proc Pfem::write::writeParametersEvent { } {
    write::WriteJSON [getParametersDict]
}

# Model Part Blocks
proc Pfem::write::writeModelPartEvent { } {
    set parts_un_list [GetPartsUN]
    foreach part_un $parts_un_list {
        write::initWriteData $part_un "PFEM_Materials"
    }
    
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    write::writeMaterials "Pfem"
    
    write::writeNodalCoordinates
    foreach part_un $parts_un_list {
        write::initWriteData $part_un "PFEM_Materials"
        write::writeElementConnectivities
    }
    Solid::write::writeConditions
    Pfem::write::writeMeshes
}

proc Pfem::write::writeMeshes { } {
    
    foreach part_un [GetPartsUN] {
        write::initWriteData $part_un "PFEM_Materials"
        write::writePartMeshes
    }
    # Solo Malla , no en conditions
    writeNodalConditions "PFEM_NodalConditions"
    
    # A Condition y a meshes-> salvo lo que no tenga topologia
    Solid::write::writeLoads
}


proc Pfem::write::writeNodalConditions { keyword } {
    write::writeNodalConditions $keyword
    return ""
    
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute $keyword]/container/blockdata"
    set groups [$root selectNodes $xp1]
    foreach group $groups {
        set cid [[$group parent] @n]
        set groupid [$group @name]
        set groupid [write::GetWriteGroupName $groupid]
        # Aqui hay que gestionar la escritura de los bodies
        # Una opcion es crear un megagrupo temporal con esa informacion, mandar a pintar, y luego borrar el grupo.
        # Otra opcion es no escribir el submodelpart. Ya tienen las parts y el project parameters tiene el conformado de los bodies
        ::write::writeGroupMesh $cid $groupid "nodal"
    }
}

proc Pfem::write::GetPartsUN { } {
    customlib::UpdateDocument
    set lista [list ]
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute "PFEM_Bodies"]/blockdata/condition"
    set i 0
    foreach part_node [$root selectNodes $xp1] {
        if {![$part_node hasAttribute "un"]} {
            set un "PFEM_Part$i"
            while {[spdAux::getRoute $un] ne ""} {
                incr i
                set un "PFEM_Part$i"
            }
            $part_node setAttribute un $un
            spdAux::setRoute $un [$part_node toXPath]
        }
        lappend lista [get_domnode_attribute $part_node un]
    }
    customlib::UpdateDocument
    return $lista
}

# Custom files (Copy python scripts, write materials file...)
proc Pfem::write::writeCustomFilesEvent { } {
    Solid::write::WriteMaterialsFile
    
    write::CopyFileIntoModel "python/RunPFEM.py"
    write::RenameFileInModel "RunPFEM.py" "MainKratos.py"
    
    #write::RenameFileInModel "ProjectParameters.json" "ProjectParameters.py"
}


proc Pfem::write::getConditionsParametersDict {un {condition_type "Condition"}} {
    
    set root [customlib::GetBaseRoot]
    
    set bcCondsDict [list ]
    
    set xp1 "[spdAux::getRoute $un]/condition/group"
    set groups [$root selectNodes $xp1]
    if {$groups eq ""} {
        set xp1 "[spdAux::getRoute $un]/group"
        set groups [$root selectNodes $xp1]
    }
    foreach group $groups {
        set groupName [$group @n]
        set cid [[$group parent] @n]
        set groupName [write::GetWriteGroupName $groupName]
        set groupId [::write::getMeshId $cid $groupName]
        set condId [[$group parent] @n]
        if {$condition_type eq "Condition"} {
            set condition [::Model::getCondition $condId]
        } {
            set condition [::Model::getNodalConditionbyId $condId]
        }
        set processName [$condition getProcessName]
        set process [::Model::GetProcess $processName]
        set processDict [dict create]
        set paramDict [dict create]
        dict set paramDict model_part_name $groupId
        
        set process_attributes [$process getAttributes]
        set process_parameters [$process getInputs]
        
        dict set process_attributes process_name [dict get $process_attributes n]
        dict unset process_attributes n
        dict unset process_attributes pn
        
        set processDict [dict merge $processDict $process_attributes]
        if {[$condition hasAttribute VariableName]} {
            set variable_name [$condition getAttribute VariableName]
            # "lindex" is a rough solution. Look for a better one.
            if {$variable_name ne ""} {dict set paramDict variable_name [lindex $variable_name 0]}
        }
        foreach {inputName in_obj} $process_parameters {
            set in_type [$in_obj getType]
            if {$in_type eq "vector"} {
                set vector_type [$in_obj getAttribute "vectorType"]
                if {$vector_type eq "bool"} {
                    set ValX [expr [get_domnode_attribute [$group find n ${inputName}X] v] ? True : False]
                    set ValY [expr [get_domnode_attribute [$group find n ${inputName}Y] v] ? True : False]
                    set ValZ [expr False]
                    if {[$group find n ${inputName}Z] ne ""} {set ValZ [expr [get_domnode_attribute [$group find n ${inputName}Z] v] ? True : False]}
                    dict set paramDict $inputName [list $ValX $ValY $ValZ]
                } {
                    if {[$in_obj getAttribute "enabled"] in [list "1" "0"]} {
                        foreach i [list "X" "Y" "Z"] {
                            if {[expr [get_domnode_attribute [$group find n Enabled_$i] v] ] ne "Yes"} {
                                set Val$i null
                            } else {
                                set printed 0
                                if {[$in_obj getAttribute "function"] eq "1"} {
                                    if {[get_domnode_attribute [$group find n "ByFunction$i"] v]  eq "Yes"} {
                                        set funcinputName "${i}function_$inputName"
                                        set value [get_domnode_attribute [$group find n $funcinputName] v]
                                        set Val$i $value
                                        set printed 1
                                    }
                                }
                                if {!$printed} {
                                    set value [expr [gid_groups_conds::convert_value_to_default [$group find n ${inputName}$i] ] ]
                                    set Val$i $value
                                }
                            }
                        }
                    } elseif {$vector_type eq "tablefile" || $vector_type eq "file"} {
                        set ValX "[get_domnode_attribute [$group find n ${inputName}X] v]"
                        set ValY "[get_domnode_attribute [$group find n ${inputName}Y] v]"
                        set ValZ "0"
                        if {[$group find n ${inputName}Z] ne ""} {set ValZ "[get_domnode_attribute [$group find n ${inputName}Z] v]"}
                    } else {
                        set ValX [expr [gid_groups_conds::convert_value_to_default [$group find n ${inputName}X] ] ]
                        set ValY [expr [gid_groups_conds::convert_value_to_default [$group find n ${inputName}Y] ] ]
                        set ValZ [expr 0.0]
                        if {[$group find n ${inputName}Z] ne ""} {set ValZ [expr [gid_groups_conds::convert_value_to_default [$group find n ${inputName}Z] ]]}
                    }
                    dict set paramDict $inputName [list $ValX $ValY $ValZ]
                }
            } elseif {$in_type eq "double" || $in_type eq "integer"} {
                set printed 0
                if {[$in_obj getAttribute "function"] eq "1"} {
                    if {[get_domnode_attribute [$group find n "ByFunction"] v]  eq "Yes"} {
                        set funcinputName "function_$inputName"
                        set value [get_domnode_attribute [$group find n $funcinputName] v]
                        dict set paramDict $inputName $value
                        set printed 1
                    }
                }
                if {!$printed} {
                    set value [gid_groups_conds::convert_value_to_default [$group find n $inputName]]
                    #set value [get_domnode_attribute [$group find n $inputName] v]
                    dict set paramDict $inputName [expr $value]
                }
            } elseif {$in_type eq "bool"} {
                set value [get_domnode_attribute [$group find n $inputName] v]
                set value [expr $value ? True : False]
                dict set paramDict $inputName [expr $value]
            } elseif {$in_type eq "tablefile"} {
                set value [get_domnode_attribute [$group find n $inputName] v]
                dict set paramDict $inputName $value
            } else {
                if {[get_domnode_attribute [$group find n $inputName] state] ne "hidden" } {
                    set value [get_domnode_attribute [$group find n $inputName] v]
                    dict set paramDict $inputName $value
                }
            }
        }
        if {[$group find n Interval] ne ""} {dict set paramDict interval [write::getInterval  [get_domnode_attribute [$group find n Interval] v]] }
        dict set processDict Parameters $paramDict
        lappend bcCondsDict $processDict
    }
    return $bcCondsDict
}

Pfem::write::Init
