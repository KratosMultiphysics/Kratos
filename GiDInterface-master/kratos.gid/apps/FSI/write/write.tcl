namespace eval FSI::write {

}

proc FSI::write::Init { } {

}

# Events
proc FSI::write::writeModelPartEvent { } {
    set filename "[file tail [GiD_Info project ModelName]]"
    
    ::Fluid::write::SetCoordinatesByGroups 1
    write::writeAppMDPA Fluid
    write::RenameFileInModel "$filename.mdpa" "${filename}_Fluid.mdpa"
    
    ::Structural::write::SetCoordinatesByGroups 1
    write::writeAppMDPA Structural
    write::RenameFileInModel "$filename.mdpa" "${filename}_Structural.mdpa"
}

proc FSI::write::writeCustomFilesEvent { } {
    Structural::write::WriteMaterialsFile
    
    write::CopyFileIntoModel "python/KratosFSI.py"
    write::RenameFileInModel "KratosFSI.py" "MainKratos.py"
}


FSI::write::Init
