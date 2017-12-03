proc DEM::write::WriteMDPAInlet { } {
    # Headers
    write::writeModelPartData

    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"

    # Nodal coordinates (only for DEM Parts <inefficient> )
    write::writeNodalCoordinatesOnGroups [GetInletGroups]
    
    # SubmodelParts
    writeInletMeshes
}

proc DEM::write::GetInletGroups { } {
    set groups [list ]
    set xp1 "[spdAux::getRoute [GetAttribute conditions_un]]/condition\[@n = 'Inlet'\]/group"
    foreach group [[customlib::GetBaseRoot] selectNodes $xp1] {
        set groupid [$group @n]
        lappend groups [write::GetWriteGroupName $groupid]
    }
    return $groups
}

proc DEM::write::writeInletMeshes { } {
    foreach groupid [DEM::write::GetInletGroups] {
        ::write::writeGroupMesh Inlet $groupid "nodal"
    }
}