namespace eval Structural::write {

}

proc Structural::write::Init { } {
    
}


proc Structural::write::writeCustomFilesEvent { } {
    return [Solid::write::writeCustomFilesEvent]
}

# MDPA Blocks

proc Structural::write::writeModelPartEvent { } {
    Solid::write::writeModelPartEvent
}

# Project Parameters
proc Structural::write::writeParametersEvent { } {
    Solid::write::writeParametersEvent
}

Structural::write::Init
