namespace eval Pfem::write {
}

proc Pfem::write::Init { } {
    # Namespace variables inicialization
}

# Project Parameters
proc Pfem::write::writeParametersEvent { } {
    write::WriteString "Project parameters file"

}

# Model Part Blocks
proc Pfem::write::writeModelPartEvent { } {
    write::initWriteData "SMParts" "SMMaterials"
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
}

# Custom files (Copy python scripts, write materials file...)
proc Pfem::write::writeCustomFilesEvent { } {

}


Pfem::write::Init
