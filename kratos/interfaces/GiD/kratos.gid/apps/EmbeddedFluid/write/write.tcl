namespace eval EmbeddedFluid::write {

}

proc EmbeddedFluid::write::Init { } {
    # Namespace variables inicialization
}

# Events
proc EmbeddedFluid::write::writeModelPartEvent { } {
    Fluid::write::AddValidApps "EmbeddedFluid"
    set err [Fluid::write::Validate]
    if {$err ne ""} {error $err}
    write::initWriteData $Fluid::write::PartsUN "EMBFLMaterials"
    write::writeModelPartData
    Fluid::write::writeProperties
    write::writeMaterials $Fluid::write::validApps
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
    set must_write [write::getValue EMBFLDistanceReading ReadingMode]
    if {$must_write eq "from_mdpa"} {
        set go 0
        set distfilepath [file join $::write::dir "[file tail [GiD_Info project modelname] ].post.res"]
        set a [open $distfilepath r]
        write::WriteString "Begin NodalData DISTANCE"
        while {[gets $a line]>=0} {
            if {$line eq "End Values"} {set go 0}
            if {$go == 2} {
                write::WriteString $line
            }
            if {$line eq {Result "Distance" "Signed distance" 1 Scalar OnNodes}} {incr go 1}
            if {$line eq "Values"} {incr go 1}
            
        }
        write::WriteString "End NodalData"
        close $a
    }
}

EmbeddedFluid::write::Init
