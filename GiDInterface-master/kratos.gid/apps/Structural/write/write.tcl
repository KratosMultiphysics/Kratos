namespace eval Structural::write {
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    variable writeAttributes
}

proc Structural::write::Init { } {
    variable ConditionsDictGroupIterators
    variable NodalConditionsGroup
    set ConditionsDictGroupIterators [dict create]
    set NodalConditionsGroup [list ]
    
    variable writeAttributes
    set writeAttributes [dict create]
    SetAttribute validApps [list "Structural"]
    SetAttribute writeCoordinatesByGroups 0
    SetAttribute properties_location json 
    SetAttribute parts_un STParts
    SetAttribute materials_un STMaterials
    SetAttribute conditions_un STLoads
    SetAttribute nodal_conditions_un STNodalConditions
    SetAttribute materials_file "StructuralMaterials.json"
    SetAttribute main_script_file "KratosStructural.py"
}

proc Structural::write::GetAttribute {att} {
    variable writeAttributes
    return [dict get $writeAttributes $att]
}

proc Structural::write::GetAttributes {} {
    variable writeAttributes
    return $writeAttributes
}

proc Structural::write::SetAttribute {att val} {
    variable writeAttributes
    dict set writeAttributes $att $val
}

proc Structural::write::AddAttribute {att val} {
    variable writeAttributes
    dict append writeAttributes $att $val]
}

proc Structural::write::AddAttributes {configuration} {
    variable writeAttributes
    set writeAttributes [dict merge $writeAttributes $configuration]
}

proc Structural::write::AddValidApps {appList} {
    AddAttribute validApps $appList
}

proc Structural::write::writeCustomFilesEvent { } {
    WriteMaterialsFile
    
    write::SetParallelismConfiguration
    
    set orig_name [GetAttribute main_script_file]
    write::CopyFileIntoModel [file join "python" $orig_name ]
    write::RenameFileInModel $orig_name "MainKratos.py"
}

proc Structural::write::SetCoordinatesByGroups {value} {
    SetAttribute writeCoordinatesByGroups $value
}

proc Structural::write::ApplyConfiguration { } {
    variable writeAttributes
    write::SetConfigurationAttributes $writeAttributes
}

# MDPA Blocks
proc Structural::write::writeModelPartEvent { } {
    variable ConditionsDictGroupIterators
    initLocalWriteConfiguration
    write::initWriteConfiguration [GetAttributes]
    
    # Headers
    write::writeModelPartData
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"

    # Materials
    # write::writeMaterials [GetAttribute validApps 

    # Nodal coordinates (1: Print only Structural nodes <inefficient> | 0: the whole mesh <efficient>)
    if {[GetAttribute writeCoordinatesByGroups]} {write::writeNodalCoordinatesOnParts} {write::writeNodalCoordinates}
    
    # Element connectivities (Groups on STParts)
    write::writeElementConnectivities

    # Local Axes
    Structural::write::writeLocalAxes

    # Nodal conditions and conditions
    writeConditions

    # SubmodelParts
    writeMeshes

    # Custom SubmodelParts
    set basicConds [write::writeBasicSubmodelParts [getLastConditionId]]
    set ConditionsDictGroupIterators [dict merge $ConditionsDictGroupIterators $basicConds]

}


proc Structural::write::writeConditions { } {
    variable ConditionsDictGroupIterators
    set ConditionsDictGroupIterators [write::writeConditions [GetAttribute conditions_un] ]
}

proc Structural::write::writeMeshes { } {
    
    write::writePartMeshes
    
    # Solo Malla , no en conditions
    write::writeNodalConditions [GetAttribute nodal_conditions_un]
    
    # A Condition y a meshes-> salvo lo que no tenga topologia
    writeLoads
}

proc Structural::write::writeLoads { } {
    variable ConditionsDictGroupIterators
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute [GetAttribute conditions_un]]/condition/group"
    foreach group [$root selectNodes $xp1] {
        set groupid [$group @n]
        set groupid [write::GetWriteGroupName $groupid]
        #W "Writing mesh of Load $groupid"
        if {$groupid in [dict keys $ConditionsDictGroupIterators]} {
            ::write::writeGroupMesh [[$group parent] @n] $groupid "Conditions" [dict get $ConditionsDictGroupIterators $groupid]
        } else {
            ::write::writeGroupMesh [[$group parent] @n] $groupid "nodal"
        }
    }
}

proc Structural::write::writeCustomBlock { } {
    write::WriteString "Begin Custom"
    write::WriteString "Custom write for Structural, any app can call me, so be careful!"
    write::WriteString "End Custom"
    write::WriteString ""
}

proc Structural::write::getLastConditionId { } { 
    variable ConditionsDictGroupIterators
    set top 1
    if {$ConditionsDictGroupIterators ne ""} {
        foreach {group iters} $ConditionsDictGroupIterators {
            set top [expr max($top,[lindex $iters 1])]
        }
    }
    return $top
}

# Custom files
proc Structural::write::WriteMaterialsFile { } {
    write::writePropertiesJsonFile [GetAttribute parts_un] [GetAttribute materials_file]
}

proc Structural::write::GetUsedElements { {get "Objects"} } {
    set xp1 "[spdAux::getRoute [GetAttribute parts_un]]/group"
    set lista [list ]
    foreach gNode [[customlib::GetBaseRoot] selectNodes $xp1] {
        set elem_name [get_domnode_attribute [$gNode selectNodes ".//value\[@n='Element']"] v]
        set e [Model::getElement $elem_name]
        if {$get eq "Name"} { set e [$e getName] }
        lappend lista $e
    }
    return $lista
}

proc Structural::write::writeLocalAxes { } {
    set xp1 "[spdAux::getRoute [GetAttribute parts_un]]/group"
    foreach gNode [[customlib::GetBaseRoot] selectNodes $xp1] {
        set elem_name [get_domnode_attribute [$gNode selectNodes ".//value\[@n='Element']"] v]
        set e [Model::getElement $elem_name]
        if {[write::isBooleanTrue [$e getAttribute "RequiresLocalAxes"]]} { 
            set group [$gNode @n]
            if {[GiD_EntitiesGroups get $group elements -count -element_type linear]} {
                write::WriteString "Begin ElementalData LOCAL_AXIS_2 // Element: $elem_name // Groups: $group"
                foreach line [GiD_EntitiesGroups get $group elements -element_type linear] {
                    set raw [lindex [lindex [GiD_Info conditions -localaxesmat line_Local_axes mesh $line] 0] 3]
                    set y0 [lindex $raw 1]
                    set y1 [lindex $raw 4]
                    set y2 [lindex $raw 7]
                    write::WriteString [format "%5d \[3\](%14.10f, %14.10f, %14.10f)" $line $y0 $y1 $y2]
                }
                write::WriteString "End ElementalData"
                write::WriteString ""
            }
        }
    }
}

proc Structural::write::initLocalWriteConfiguration { } {
    
    if {[usesContact]} {
         SetAttribute main_script_file "KratosContactStructural.py"
    }
}

proc Structural::write::usesContact { } {
    set result_node [[customlib::GetBaseRoot] selectNodes "[spdAux::getRoute STNodalConditions]/condition\[@n = 'CONTACT'\]/group"]
    
    if {$result_node ne ""} {
        return 1
    } {
        return 0
    }
}

Structural::write::Init
