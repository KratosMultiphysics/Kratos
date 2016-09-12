namespace eval Pfem::write {
}

proc Pfem::write::Init { } {
    Solid::write::AddValidApps "Pfem"
}

proc Pfem::write::writeParametersEvent { } {
    write::WriteJSON [getParametersDict]
}

# Model Part Blocks
proc Pfem::write::writeModelPartEvent { } {
    Solid::write::writeModelPartEvent
}


# Custom files (Copy python scripts, write materials file...)
proc Pfem::write::writeCustomFilesEvent { } {

}


Pfem::write::Init
