namespace eval BaseApp::write {
}

proc BaseApp::write::Init { } {
    # Namespace variables inicialization
}

# Project Parameters
proc BaseApp::write::writeParametersEvent { } {
    write::WriteString "Project parameters file"

}

# Model Part Blocks
proc BaseApp::write::writeModelPartEvent { } {
    write::initWriteData "SMParts" "SMMaterials"
    write::WriteString "Begin Properties 0"
    write::WriteString "End Properties"
}

# Custom files (Copy python scripts, write materials file...)
proc BaseApp::write::writeCustomFilesEvent { } {

}


BaseApp::write::Init
