namespace eval EmbeddedFluid::write {
    variable writeAttributes
}

proc EmbeddedFluid::write::Init { } {
    # Namespace variables inicialization
    SetAttribute parts_un FLParts
    SetAttribute nodal_conditions_un FLNodalConditions
    SetAttribute conditions_un FLBC
    SetAttribute materials_un EMBFLMaterials
    SetAttribute writeCoordinatesByGroups 0
    SetAttribute validApps [list "Fluid" "EmbeddedFluid"]
    SetAttribute main_script_file "KratosFluid.py"
    SetAttribute materials_file "FluidMaterials.json"
}

# Events
proc EmbeddedFluid::write::writeModelPartEvent { } {
    # Fluid::write::AddValidApps "EmbeddedFluid"
    set err [Fluid::write::Validate]
    if {$err ne ""} {error $err}
    write::initWriteConfiguration [GetAttributes]
    write::writeModelPartData
    Fluid::write::writeProperties
    write::writeMaterials [GetAttribute validApps]
    write::writeNodalCoordinatesOnParts
    write::writeElementConnectivities
    Fluid::write::writeConditions
    Fluid::write::writeMeshes
    writeDistances
}
proc EmbeddedFluid::write::writeCustomFilesEvent { } {
    write::CopyFileIntoModel "python/KratosFluid.py"
    write::RenameFileInModel "KratosFluid.py" "MainKratos.py"
}

proc EmbeddedFluid::write::writeDistances { } {
    set must_write [write::getValue EMBFLDistanceSettings ReadingMode]
    if {$must_write eq "from_mdpa"} {
        set go 0
        set distfilepath [file join [write::GetConfigurationAttribute dir] "[file tail [GiD_Info project modelname] ].post.res"]
        if {![file exists $distfilepath]} {error [= "Distances file does not exist. Please check meshing parameters and mesh again."]}
        set a [open $distfilepath r]
        write::WriteString "Begin NodalData DISTANCE"
        while {[gets $a line]>=0} {
            if {$line eq "End Values"} {set go 0}
            if {$go == 2} {
                lassign [split $line " "] node dist
                write::WriteString "$node 0 $dist"
            }
            if {$line eq {Result "Distance" "Signed distance" 1 Scalar OnNodes}} {incr go 1}
            if {$line eq "Values"} {incr go 1}

        }
        write::WriteString "End NodalData"
        close $a
    }
}

proc EmbeddedFluid::write::GetAttribute {att} {
    variable writeAttributes
    return [dict get $writeAttributes $att]
}

proc EmbeddedFluid::write::GetAttributes {} {
    variable writeAttributes
    return $writeAttributes
}

proc EmbeddedFluid::write::SetAttribute {att val} {
    variable writeAttributes
    dict set writeAttributes $att $val
}

EmbeddedFluid::write::Init
