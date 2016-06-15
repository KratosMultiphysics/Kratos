namespace eval FSI::write {

}

proc FSI::write::Init { } {

}

# Events
proc FSI::write::writeModelPartEvent { } {
    set filename "[file tail [GiD_Info project ModelName]]"
    write::writeAppMDPA Structural
    write::RenameFileInModel "$filename.mdpa" "${filename}_Structural.mdpa"
    write::writeAppMDPA Fluid
    write::RenameFileInModel "$filename.mdpa" "${filename}_Fluid.mdpa"
}

proc FSI::write::writeCustomFilesEvent { } {
    write::CopyFileIntoModel "python/MAIN_FILE_FSI.py"
    write::RenameFileInModel "MAIN_FILE_FSI.py" "MainKratos.py"
}


FSI::write::Init
