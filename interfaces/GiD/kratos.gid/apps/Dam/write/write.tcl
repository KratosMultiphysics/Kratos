namespace eval Dam::write {
    variable nodalwrite
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    variable TableDict
}

proc Dam::write::Init { } {
    # Namespace variables inicialization
    variable nodalwrite
    set nodalwrite 0
    
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    set ConditionsDictGroupIterators [dict create]
    set NodalConditionsGroup [list ]
    
    # key = file path
    # value = id table
    variable TableDict
    set TableDict [dict create]
}


proc Dam::write::writeCustomFilesEvent { } {
    
    write::CopyFileIntoModel "python/dam_thermo_mechanic_script.py"
    write::RenameFileInModel "dam_thermo_mechanic_script.py" "MainKratos.py"
    
    #write::RenameFileInModel "ProjectParameters.json" "ProjectParameters.py"
}

# MDPA Blocks

proc Dam::write::writeModelPartEvent { } {
    write::initWriteData "DamParts" "DamMaterials"
    
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    
    UpdateMaterials
    write::writeMaterials
    Dam::write:::writeTables
    write::writeNodalCoordinates
    write::writeElementConnectivities
    writeConditions
    writeMeshes
    #writeCustomBlock
}

proc Dam::write::UpdateMaterials { } {
    set matdict [write::getMatDict]
    foreach {mat props} $matdict {
        set constlaw [dict get $props ConstitutiveLaw]
        # Modificar la ley constitutiva
        set newconstlaw $constlaw
        if {$constlaw eq "BilinearCohesive2DPlaneStress"} {set newconstlaw "BilinearCohesive2DLaw"}
        if {$constlaw eq "BilinearCohesive2DPlaneStrain"} {
            dict set matdict $mat THICKNESS  1.0000E+00
            set newconstlaw "BilinearCohesive2DLaw"
            }
        dict set matdict $mat CONSTITUTIVE_LAW_NAME $newconstlaw
    }
    write::setMatDict $matdict
}

proc Dam::write::writeConditions { } {
    variable ConditionsDictGroupIterators
    set ConditionsDictGroupIterators [write::writeConditions "DamLoads"]
}

proc Dam::write::writeMeshes { } {
    
    write::writePartMeshes
    
    # Solo Malla , no en conditions
    writeNodalConditions "DamNodalConditions"
    
    # A Condition y a meshes-> salvo lo que no tenga topologia
    writeLoads
}

proc Dam::write::writeNodalConditions { keyword } {
    variable TableDict
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute $keyword]/condition/group"
    set groups [$root selectNodes $xp1]
    if {$groups eq ""} {
        set xp1 "[spdAux::getRoute $keyword]/group"
        set groups [$root selectNodes $xp1]
    }
    foreach group $groups {
        set condid [[$group parent] @n]
        set groupid [$group @n]
        set tableid ""
        if {[dict exists $TableDict $condid $groupid]} {
            set tableid [dict get $TableDict $condid $groupid tableid]
        }
        ::write::writeGroupMesh $condid $groupid "nodal" "" $tableid
    }
}

proc Dam::write::writeLoads { } {
    variable TableDict
    variable nodalwrite
    variable ConditionsDictGroupIterators
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute "DamLoads"]/condition/group"
    foreach group [$root selectNodes $xp1] {
        set condid [get_domnode_attribute [$group parent] n]
        set groupid [get_domnode_attribute $group n]
        #W "Writing mesh of Load $condid $groupid"
        set tableid ""
        if {[dict exists $TableDict $condid $groupid]} {
            set tableid [dict get $TableDict $condid $groupid tableid]
        }
        #W "table $tableid"
        if {$groupid in [dict keys $ConditionsDictGroupIterators]} {
            ::write::writeGroupMesh [[$group parent] @n] $groupid "Conditions" [dict get $ConditionsDictGroupIterators $groupid] $tableid
        } else {
            ::write::writeGroupMesh [[$group parent] @n] $groupid "nodal" "" $tableid
        }
    }
}

proc Dam::write::getVariableParametersDict {un {condition_type "Condition"}} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    
    set xp1 "[spdAux::getRoute $un]/condition/group"
    set groups [$root selectNodes $xp1]

    foreach group $groups {
        set groupName [$group @n]
        #W "GROUP $groupName"
        set cid [[$group parent] @n]
        set groupId [::write::getMeshId $cid $groupName]
        set condId [[$group parent] @n]
        if {$condition_type eq "Condition"} {
            set condition [::Model::getCondition $condId]
        } {
            set condition [::Model::getNodalConditionbyId $condId]
        }
        #W "Condition = $condition"
        catch {
            set variable_name [$condition getAttribute VariableName]
            #W $variable_name
            # "lindex" is a rough solution. Look for a better one.
            if {$variable_name ne ""} {dict set paramDict variable_name [lindex $variable_name 0]}
        }        
        
        return $variable_name
    }
}

proc Dam::write::GetTableidFromFileid { filename } {
    variable TableDict
    foreach condid [dict keys $TableDict] {
        foreach groupid [dict keys [dict get $TableDict $condid]] {
            if {[dict get $TableDict $condid $groupid fileid] eq $filename} {
                return [dict get $TableDict $condid $groupid tableid]
            }
        }
    }
    return 0
}

proc Dam::write::writeTables { } {
    variable TableDict
    foreach table [GetPrinTables] {
        lassign $table tableid fileid condid groupid
        dict set TableDict $condid $groupid tableid $tableid
        dict set TableDict $condid $groupid fileid $fileid
        write::WriteString "Begin Table $tableid"
        set data [GidUtils::ReadFile $fileid]
        write::WriteString [string map {; { }} $data]
        write::WriteString "End Table"
        write::WriteString ""
    }
}

proc Dam::write::GetPrinTables {} {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set listaTablas [list ]
    set listaFiles [list ]
    set num 1
    set origins [list "DamLoads" "DamNodalConditions"]
    foreach unique_name $origins {
        set xpathCond "[spdAux::getRoute $unique_name]/condition/group/value\[@type='tablefile'\]"
        foreach node [$root selectNodes $xpathCond] {
            set fileid [get_domnode_attribute $node v]
            set groupid [get_domnode_attribute [$node parent] n]
            set condid [get_domnode_attribute [[$node parent] parent] n]
            if {$fileid ne "" && $fileid ni $listaFiles} {
                lappend listaTablas [list $num $fileid $condid $groupid]
                lappend listaFiles $fileid
                incr num
            }
        }
    }
    return $listaTablas
}



Dam::write::Init
