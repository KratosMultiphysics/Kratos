proc DEM::write::WriteMDPAParts { } {
    # Headers
    write::writeModelPartData

    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
    write::writeMaterials [GetAttribute validApps]

    # Nodal coordinates (only for DEM Parts <inefficient> )
    write::writeNodalCoordinatesOnParts
    
    # Element connectivities (Groups on STParts)
    write::writeElementConnectivities

    # Element radius
    writeSphereRadius

    # SubmodelParts
    write::writePartMeshes
    writeVelocityMeshes
}

proc DEM::write::writeSphereRadius { } {
    set root [customlib::GetBaseRoot]
    set xp1 "[spdAux::getRoute [GetAttribute parts_un]]/group"
    foreach group [$root selectNodes $xp1] {
        set groupid [$group @n]
        set grouppid [write::GetWriteGroupName $groupid]
        write::WriteString "Begin NodalData RADIUS // GUI group identifier: $grouppid"
        GiD_WriteCalculationFile connectivities [dict create $groupid "%.0s %10d 0 %10g\n"]
        write::WriteString "End NodalData"
        write::WriteString ""
    }
}

proc DEM::write::GetNodalConditionsGroups { {include_cond 0} } {
    set groups [list ]
    set xp1 "[spdAux::getRoute [GetAttribute nodal_conditions_un]]/condition/group"
    foreach group [[customlib::GetBaseRoot] selectNodes $xp1] {
        set groupid [$group @n]
        if {$include_cond} {lappend groups [[$group parent] @n]}
        lappend groups [write::GetWriteGroupName $groupid]
    }
    return $groups
}

proc DEM::write::writeVelocityMeshes { } {
    foreach {cid groupid} [DEM::write::GetNodalConditionsGroups 1] {
        ::write::writeGroupMesh $cid $groupid "nodal"
    }
}
