namespace eval Dam::write {
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    variable TableDict
    
    variable ThermalSubModelPartDict
}

proc Dam::write::Init { } {
    # Namespace variables inicialization
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    set ConditionsDictGroupIterators [dict create]
    set NodalConditionsGroup [list ]
    
    variable TableDict
    catch {unset TableDict}
    set TableDict [dict create]
}

proc Dam::write::writeCustomFilesEvent { } {
    
    set damTypeofProblem [write::getValue DamTypeofProblem]
    set damSelfweight [write::getValue DamSelfweight ConsiderSelf]
    set damConstructionProcess [write::getValue DamConstructionProcess Activate_construction]
    
    if {$damTypeofProblem eq "Acoustic"} {
        write::CopyFileIntoModel "python/dam_acoustic_script.py"
        write::RenameFileInModel "dam_acoustic_script.py" "MainKratos.py"
    } elseif {$damTypeofProblem eq "Modal-Analysis" } {
        write::CopyFileIntoModel "python/dam_eigen_script.py"
        write::RenameFileInModel "dam_eigen_script.py" "MainKratos.py"
    } elseif {$damSelfweight eq "Yes" } {
        write::CopyFileIntoModel "python/dam_main_selfweight.py"
        write::RenameFileInModel "dam_main_selfweight.py" "MainKratos.py"
    } elseif {$damConstructionProcess eq "True" } { 
        write::CopyFileIntoModel "python/dam_main_construction.py"
        write::RenameFileInModel "dam_main_construction.py" "MainKratos.py"
    } else {
        write::CopyFileIntoModel "python/dam_main.py"
        write::RenameFileInModel "dam_main.py" "MainKratos.py"
    }
    
}

# MDPA Blocks
proc Dam::write::writeModelPartEvent { } {
    write::initWriteData "DamParts" "DamMaterials"
    
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    
    Dam::write::UpdateMaterials
    write::writeMaterials
    Dam::write::writeTables
    write::writeNodalCoordinates
    write::writeElementConnectivities
    
    set damTypeofProblem [write::getValue DamTypeofProblem]
    if {$damTypeofProblem eq "Thermo-Mechanical" || $damTypeofProblem eq "UP_Thermo-Mechanical"} {
        Dam::write::writeThermalElements
    }
    
    Dam::write::writeConditions
    Dam::write::writeMeshes

    # Creation of special mdpa for computing an extra problem just considering selfweight
    set damSelfweight [write::getValue DamSelfweight ConsiderSelf]
    if {$damSelfweight eq "Yes" } {
   
        Dam::write:writeExtraMdpaSelfWeight
    }
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
    set pairs [lsort -increasing -index end [dict values $ConditionsDictGroupIterators] ]
    set index [lindex [lindex [lsort -integer -index 0 $pairs] end] end]
    if {$index eq ""} {
        set index 0
    }
    
    set ThermalConditionGroups [write::writeConditions "DamThermalLoads" $index]
    set ConditionsDictGroupIterators [dict merge $ConditionsDictGroupIterators $ThermalConditionGroups]
}

proc Dam::write::writeMeshes { } {
    
    write::writePartMeshes
    
    set damTypeofProblem [write::getValue DamTypeofProblem]
    if {$damTypeofProblem eq "Thermo-Mechanical" || $damTypeofProblem eq "UP_Thermo-Mechanical"} {
        Dam::write::ThermalSubModelPart
    }
    
    # Solo Malla , no en conditions
    writeNodalConditions "DamNodalConditions"
    
    # A Condition y a meshes-> salvo lo que no tenga topologia
    writeLoads "DamLoads"
    writeLoads "DamThermalLoads"
}



proc Dam::write::writeNodalConditions { keyword } {
    variable TableDict
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute $keyword]/condition/group"
    set groups [$root selectNodes $xp1]
    if {$groups eq ""} {
        set xp1 "[spdAux::getRoute $keyword]/group"
        set groups [$root selectNodes $xp1]
    }
    foreach group $groups {
        set condid [[$group parent] @n]
        set groupid [$group @n]
        set groupid [write::GetWriteGroupName $groupid]
        set tableid [list ]
        if {[dict exists $TableDict $condid $groupid]} {
            set groupdict [dict get $TableDict $condid $groupid]
            foreach valueid [dict keys $groupdict] {
                lappend tableid [dict get $groupdict $valueid tableid]
            }
        }
        ::write::writeGroupMesh $condid $groupid "nodal" "" $tableid
    }
}

proc Dam::write::writeLoads { baseUN } {
    variable TableDict
    variable ConditionsDictGroupIterators
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute $baseUN]/condition/group"
    foreach group [$root selectNodes $xp1] {
        set condid [get_domnode_attribute [$group parent] n]
        set groupid [get_domnode_attribute $group n]
        set groupid [write::GetWriteGroupName $groupid]
        set tableid [list ]
        if {[dict exists $TableDict $condid $groupid]} {
            set groupdict [dict get $TableDict $condid $groupid]
            foreach valueid [dict keys $groupdict] {
                lappend tableid [dict get $groupdict $valueid tableid]
            }
        }
        if {$groupid in [dict keys $ConditionsDictGroupIterators]} {
            ::write::writeGroupMesh [[$group parent] @n] $groupid "Conditions" [dict get $ConditionsDictGroupIterators $groupid] $tableid
        } else {
            ::write::writeGroupMesh [[$group parent] @n] $groupid "nodal" "" $tableid
        }
    }
}

proc Dam::write::getVariableNameList {un {condition_type "Condition"}} {
    set xp1 "[spdAux::getRoute $un]/condition/group"
    set groups [[customlib::GetBaseRoot] selectNodes $xp1]
    
    set variable_list [list ]
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
        set variable_name [$condition getAttribute VariableName]
        if {$variable_name ne ""} {lappend variable_list [lindex $variable_name 0]}  
    }
    return $variable_list
}

proc Dam::write::GetTableidFromFileid { filename } {
    variable TableDict
    foreach condid [dict keys $TableDict] {
        foreach groupid [dict keys [dict get $TableDict $condid]] {
            foreach valueid [dict keys [dict get $TableDict $condid $groupid]] {
                if {[dict get $TableDict $condid $groupid $valueid fileid] eq $filename} {
                    return [dict get $TableDict $condid $groupid $valueid tableid]
                }
            }
        }
    }
    return 0
}

proc Dam::write::writeTables { } {
    variable TableDict
    set printed_tables [list ]
    foreach table [GetPrinTables] {
        lassign $table tableid fileid condid groupid valueid
        dict set TableDict $condid $groupid $valueid tableid $tableid
        dict set TableDict $condid $groupid $valueid fileid $fileid
        if {$tableid ni $printed_tables} {
            lappend printed_tables $tableid
            write::WriteString "Begin Table $tableid TIME VALUE"
            if {[string index $fileid 0] eq "."} {
                set modelname [GiD_Info project ModelName]
                set filename [string range $fileid 2 end]
                set fileid [file join "$modelname.gid" $filename]
            }
            set data [GidUtils::ReadFile $fileid]
            write::WriteString [string map {; { }} $data]
            write::WriteString "End Table"
            write::WriteString ""
        }
    }
}

proc Dam::write::GetPrinTables {} {
    
    set root [customlib::GetBaseRoot]
    FileSelector::CopyFilesIntoModel [file join [GiD_Info project ModelName] ".gid"]
    set listaTablas [list ]
    set listaFiles [list ]
    set num 0
    set origins [list "DamLoads" "DamThermalLoads" "DamNodalConditions"]
    foreach unique_name $origins {
        set xpathCond "[spdAux::getRoute $unique_name]/condition/group/value\[@type='tablefile'\]"
        foreach node [$root selectNodes $xpathCond] {
            set fileid [get_domnode_attribute $node v]
            set valueid [get_domnode_attribute $node n]
            set groupid [get_domnode_attribute [$node parent] n]
            set condid [get_domnode_attribute [[$node parent] parent] n]
            #W $condid
            if {$fileid ni [list "" "- No file"]} {
                if {$fileid ni $listaFiles} {
                    lappend listaFiles $fileid
                    incr num
                    set tableid $num
                } else {
                    set tableid 0
                    foreach table $listaTablas {
                        lassign $table tableid2 fileid2 condid2 groupid2 valueid2
                        if {$fileid2 eq $fileid} {set tableid $tableid2; break}
                    }
                }
                #W "$tableid $fileid $condid $groupid $valueid"
                lappend listaTablas [list $tableid $fileid $condid $groupid $valueid]
            }
        }
    }
    return $listaTablas
}

#-------------------------------------------------------------------------------

proc Dam::write::writeThermalElements {} {
    
    set ThermalGroups [list]
    
    set mat_dict [write::getMatDict]
    foreach part_name [dict keys $mat_dict] {
        if {[[Model::getElement [dict get $mat_dict $part_name Element]] getAttribute "ElementType"] eq "Solid"} {
            lappend ThermalGroups $part_name
        }
    }
    
    set ElementId [GiD_Info Mesh MaxNumElements]
    variable ThermalSubModelPartDict
    set ThermalSubModelPartDict [dict create]
    
    for {set i 0} {$i < [llength $ThermalGroups]} {incr i} {
        
        set ElementList [list]
        
        # EulerianConvDiff2D
        Dam::write::writeThermalConnectivities [lindex $ThermalGroups $i] triangle EulerianConvDiff2D "Dam::write::Triangle2D3Connectivities" ElementId ElementList
        # EulerianConvDiff2D4N
        Dam::write::writeThermalConnectivities [lindex $ThermalGroups $i] quadrilateral EulerianConvDiff2D4N "Dam::write::Quadrilateral2D4Connectivities" ElementId ElementList
        # EulerianConvDiff3D
        Dam::write::writeThermalConnectivities [lindex $ThermalGroups $i] tetrahedra EulerianConvDiff3D "Dam::write::Quadrilateral2D4Connectivities" ElementId ElementList
        # EulerianConvDiff3D8N
        Dam::write::writeThermalConnectivities [lindex $ThermalGroups $i] hexahedra EulerianConvDiff3D8N "Dam::write::Hexahedron3D8Connectivities" ElementId ElementList
        
        dict set ThermalSubModelPartDict [lindex $ThermalGroups $i] Elements $ElementList
        dict set ThermalSubModelPartDict [lindex $ThermalGroups $i] SubModelPartName "Thermal_Part_Auto_[expr {$i+1}]"
    }
    
    
}

proc Dam::write::writeThermalConnectivities {Group ElemType ElemName ConnectivityType ElementId ElementList} {
    set Entities [GiD_EntitiesGroups get $Group elements -element_type $ElemType]
    if {[llength $Entities] > 0} {
        upvar $ElementId MyElementId
        upvar $ElementList MyElementList
        
        write::WriteString "Begin Elements $ElemName // GUI group identifier: $Group"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            incr MyElementId
            lappend MyElementList $MyElementId
            write::WriteString "  $MyElementId  0  [$ConnectivityType [lindex $Entities $j]]"
        }
        write::WriteString "End Elements"
        write::WriteString ""
    }
}

proc Dam::write::Triangle2D3Connectivities { ElemId } {
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"
}


proc Dam::write::Quadrilateral2D4Connectivities { ElemId } {
    
    #Note: It is the same for the Tethrahedron3D4
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
        [lindex $ElementInfo 6]"
}

proc Dam::write::Hexahedron3D8Connectivities { ElemId } {
    
    #It is the same for Quadrilateral2D8
    
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
        [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
        [lindex $ElementInfo 9] [lindex $ElementInfo 10]"
}

#-------------------------------------------------------------------------------


proc Dam::write::ThermalSubModelPart { } {
    
    variable ThermalSubModelPartDict
    
    dict for {Group ThermalPart} $ThermalSubModelPartDict {
        
        write::WriteString "Begin SubModelPart [dict get $ThermalPart SubModelPartName] // Group $Group // Subtree Parts"
        # Nodes
        set ThermalNodes [GiD_EntitiesGroups get $Group nodes]
        write::WriteString "  Begin SubModelPartNodes"
        for {set i 0} {$i < [llength $ThermalNodes]} {incr i} {
            write::WriteString "    [lindex $ThermalNodes $i]"
        }
        write::WriteString "  End SubModelPartNodes"
        # Elements
        set ThermalElements [dict get $ThermalPart Elements]
        write::WriteString "  Begin SubModelPartElements"
        for {set i 0} {$i < [llength $ThermalElements]} {incr i} {
            write::WriteString "    [lindex $ThermalElements $i]"
        }
        write::WriteString "  End SubModelPartElements"
        # Conditions
        write::WriteString "  Begin SubModelPartConditions"
        write::WriteString "  End SubModelPartConditions"
        write::WriteString "End SubModelPart"
        write::WriteString ""
    }
}

#-------------------------------------------------------------------------------

proc Dam::write::getSubModelPartThermalNames { } {
    
    set submodelThermalPartsNames [list]
    
    variable ThermalSubModelPartDict
    dict for {Group ThermalPart} $ThermalSubModelPartDict {
        lappend submodelThermalPartsNames [dict get $ThermalPart SubModelPartName]  
    }
    
    return $submodelThermalPartsNames
}

#-------------------------------------------------------------------------------

# Processes for extra mdpa for selfweight calculations
proc Dam::write:writeExtraMdpaSelfWeight { } {

    write::OpenFile "[file tail selfweight].mdpa"
        
    write::initWriteData "DamParts" "DamMaterials"
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    Dam::write::UpdateMaterialsSelfweight
    write::writeMaterials
    write::writeNodalCoordinates
    write::writeElementConnectivities
    Dam::write::writePartMeshes
    Dam::write::writeNodalConditionsSelfWeight "DamNodalConditions"
    write::CloseFile      
}

proc Dam::write::writePartMeshes { } {
    foreach group [write::getPartsGroupsId] {
        Dam::write::writeGroupMesh Parts $group "Elements"
    }
}

proc Dam::write::UpdateMaterialsSelfweight { } {
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
        if {([string first 3D $constlaw] != -1) && ([string first Bilinear $constlaw] == -1)} {
            set newconstlaw "LinearElastic3DLaw"
        }
        if {[string first Bilinear $constlaw] == -1 } {
            if {[string first Strain $constlaw] != -1} {
                set newconstlaw "LinearElasticPlaneStrain2DLaw"
            }
            if {[string first Stress $constlaw] != -1} {
                set newconstlaw "LinearElasticPlaneStress2DLaw"
            }
        }

        dict set matdict $mat CONSTITUTIVE_LAW_NAME $newconstlaw
    }
    write::setMatDict $matdict
}

proc Dam::write::writeNodalConditionsSelfWeight { keyword } {
    variable TableDict
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute $keyword]/condition/group"
    set groups [$root selectNodes $xp1]
    if {$groups eq ""} {
        set xp1 "[spdAux::getRoute $keyword]/group"
        set groups [$root selectNodes $xp1]
    }
    foreach group $groups {
        set condid [[$group parent] @n]
        if {$condid eq "DISPLACEMENT"} {
            set groupid [$group @n]
            set groupid [write::GetWriteGroupName $groupid]
            set tableid [list ]
            if {[dict exists $TableDict $condid $groupid]} {
                set groupdict [dict get $TableDict $condid $groupid]
                foreach valueid [dict keys $groupdict] {
                    lappend tableid [dict get $groupdict $valueid tableid]
                }
            }
            Dam::write::writeGroupMesh $condid $groupid "nodal" "" $tableid
        }
    }
}

# what can be: nodal, Elements, Conditions or Elements&Conditions
proc Dam::write::writeGroupMesh { cid group {what "Elements"} {iniend ""} {tableid_list ""} } {
    variable meshes
    set meshes [dict create]

    set what [split $what "&"]
    set gtn [write::GetConfigurationAttribute groups_type_name]
    set group [write::GetWriteGroupName $group]
    if {![dict exists $meshes [list $cid ${group}]]} {
        set mid [expr [llength [dict keys $meshes]] +1]
        if {$gtn ne "Mesh"} {
            set good_name [write::transformGroupName $group]
            set mid "${cid}_${good_name}"
        }
        dict set meshes [list $cid ${group}] $mid
        set gdict [dict create]
        set f "%10i\n"
        set f [subst $f]
        dict set gdict $group $f
        write::WriteString "Begin $gtn $mid // Group $group // Subtree $cid"
        if {$tableid_list ne ""} {
            write::WriteString "    Begin SubModelPartTables"
            foreach tableid $tableid_list {
                write::WriteString "    $tableid"
            }
            write::WriteString "    End SubModelPartTables"
        }
        write::WriteString "    Begin ${gtn}Nodes"
        GiD_WriteCalculationFile nodes -sorted $gdict
        write::WriteString "    End ${gtn}Nodes"
        write::WriteString "    Begin ${gtn}Elements"
        if {"Elements" in $what} {
            GiD_WriteCalculationFile elements -sorted $gdict
        }
        write::WriteString "    End ${gtn}Elements"
        write::WriteString "    Begin ${gtn}Conditions"
        if {"Conditions" in $what} {
            #GiD_WriteCalculationFile elements -sorted $gdict
            if {$iniend ne ""} {
                #W $iniend
                foreach {ini end} $iniend {
                    for {set i $ini} {$i<=$end} {incr i} {
                        write::WriteString [format %10d $i]
                    }
                }
            }
        }
        write::WriteString "    End ${gtn}Conditions"
        write::WriteString "End $gtn"
    }
}

Dam::write::Init
