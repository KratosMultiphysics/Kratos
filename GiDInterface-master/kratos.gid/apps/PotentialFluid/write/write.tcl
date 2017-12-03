namespace eval PotentialFluid::write {
    variable writeAttributes
}

proc PotentialFluid::write::Init { } {
    # Namespace variables inicialization
    SetAttribute parts_un FLParts
    SetAttribute nodal_conditions_un FLNodalConditions
    SetAttribute conditions_un FLBC
    SetAttribute materials_un PTFLMaterials
    SetAttribute drag_un FLDrags
    SetAttribute writeCoordinatesByGroups 0
    SetAttribute validApps [list "Fluid" "PotentialFluid"]
    SetAttribute main_script_file "KratosFluid.py"
    SetAttribute materials_file "FluidMaterials.json"
}

# Events
proc PotentialFluid::write::writeModelPartEvent { } {
    Fluid::write::AddValidApps "PotentialFluid"
    set err [Fluid::write::Validate]
    if {$err ne ""} {error $err}
    write::initWriteConfiguration [GetAttributes]
    write::writeModelPartData
    Fluid::write::writeProperties
    write::writeMaterials [::Fluid::GetAttribute validApps]
    write::writeNodalCoordinatesOnParts
    write::writeElementConnectivities
    Fluid::write::writeConditions
    Fluid::write::writeMeshes
}
proc PotentialFluid::write::writeCustomFilesEvent { } {
    write::CopyFileIntoModel "python/KratosPotentialFlow.py"
    write::RenameFileInModel "KratosPotentialFlow.py" "MainKratos.py"
}


proc PotentialFluid::write::GetAttribute {att} {
    variable writeAttributes
    return [dict get $writeAttributes $att]
}

proc PotentialFluid::write::GetAttributes {} {
    variable writeAttributes
    return $writeAttributes
}

proc PotentialFluid::write::SetAttribute {att val} {
    variable writeAttributes
    dict set writeAttributes $att $val
}

proc PotentialFluid::write::AddAttribute {att val} {
    variable writeAttributes
    dict append writeAttributes $att $val]
}

proc PotentialFluid::write::AddAttributes {configuration} {
    variable writeAttributes
    set writeAttributes [dict merge $writeAttributes $configuration]
}


PotentialFluid::write::Init
