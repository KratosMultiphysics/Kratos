proc DEM::write::WriteMDPAWalls { } {
    # Headers
    write::writeModelPartData

    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"

    # Nodal coordinates (only for Walls <inefficient> )
    write::writeNodalCoordinatesOnGroups [GetWallsGroups]
    
    # Nodal conditions and conditions
    writeConditions

    # SubmodelParts
    writeConditionMeshes
}

proc DEM::write::writeConditions { } {
    variable ConditionsDictGroupIterators
    set ConditionsDictGroupIterators [write::writeConditions [GetAttribute conditions_un] ]
}

proc DEM::write::GetWallsGroups { } {
    set groups [list ]
    set xp1 "[spdAux::getRoute [GetAttribute conditions_un]]/condition\[@n = 'DEM-FEM-Wall'\]/group"
    foreach group [[customlib::GetBaseRoot] selectNodes $xp1] {
        set groupid [$group @n]
        lappend groups [write::GetWriteGroupName $groupid]
    }
    return $groups
}
proc DEM::write::GetConditionsGroups { } {
    set groups [list ]
    set xp1 "[spdAux::getRoute [GetAttribute conditions_un]]/condition/group"
    foreach group [[customlib::GetBaseRoot] selectNodes $xp1] {
        set groupid [$group @n]
        lappend groups [write::GetWriteGroupName $groupid]
    }
    return $groups
}

proc DEM::write::writeConditionMeshes { } {
    variable ConditionsDictGroupIterators
    foreach groupid [GetWallsGroups] {
        if {$groupid in [dict keys $ConditionsDictGroupIterators]} {
            ::write::writeGroupMesh "DEM-FEM-Wall" $groupid "Conditions" [dict get $ConditionsDictGroupIterators $groupid]
        } 
    }
}

