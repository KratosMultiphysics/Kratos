namespace eval DEM::write {
    variable writeAttributes
}

proc DEM::write::Init { } {    
    variable writeAttributes
    set writeAttributes [dict create]
    SetAttribute validApps [list "DEM"]
    SetAttribute writeCoordinatesByGroups 1
    SetAttribute properties_location py 
    SetAttribute parts_un DEMParts
    SetAttribute materials_un DEMMaterials
    SetAttribute conditions_un DEMConditions
    SetAttribute nodal_conditions_un DEMNodalConditions
    SetAttribute materials_file "DEMMaterials.json"
    SetAttribute main_script_file "KratosDEM.py"
}

# Attributes block
proc DEM::write::GetAttribute {att} {
    variable writeAttributes
    return [dict get $writeAttributes $att]
}

proc DEM::write::SetAttribute {att val} {
    variable writeAttributes
    dict set writeAttributes $att $val
}

proc DEM::write::AddAttribute {att val} {
    variable writeAttributes
    dict append writeAttributes $att $val]
}

proc DEM::write::AddAttributes {configuration} {
    variable writeAttributes
    set writeAttributes [dict merge $writeAttributes $configuration]
}

# MultiApp events
proc DEM::write::AddValidApps {appList} {
    AddAttribute validApps $appList
}

proc DEM::write::SetCoordinatesByGroups {value} {
    SetAttribute writeCoordinatesByGroups $value
}

proc DEM::write::ApplyConfiguration { } {
    variable writeAttributes
    write::SetConfigurationAttributes $writeAttributes
}

# MDPA Blocks
proc DEM::write::writeModelPartEvent { } {
    variable writeAttributes
    write::initWriteConfiguration $writeAttributes
    
    # MDPA Parts
    WriteMDPAParts
    write::CloseFile

    # MDPA Walls
    write::OpenFile "[file tail [GiD_Info project ModelName]]_FEM_boundary.mdpa"
    WriteMDPAWalls
    write::CloseFile

    # MDPA Inlet
    write::OpenFile "[file tail [GiD_Info project ModelName]]_Inlet.mdpa"
    WriteMDPAInlet
    write::CloseFile

    # MDPA Walls
    write::OpenFile "[file tail [GiD_Info project ModelName]]_Clusters.mdpa"
    WriteMDPAClusters
    write::CloseFile
}

proc DEM::write::writeCustomFilesEvent { } {
    set orig_name [GetAttribute main_script_file]
    write::CopyFileIntoModel [file join "python" $orig_name ]
    
    write::RenameFileInModel $orig_name "MainKratos.py"
}

DEM::write::Init