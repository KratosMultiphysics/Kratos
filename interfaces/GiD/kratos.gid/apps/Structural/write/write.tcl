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
    set project_parameters_dict [::Solid::write::getParametersDict]
    dict set project_parameters_dict solver_settings rotation_dofs [UsingRotationDofElements]
    write::WriteJSON $project_parameters_dict
}

proc Structural::write::UsingRotationDofElements { } {
    set doc $gid_groups_conds::doc
    set root [$doc documentElement]
    set xp1 "[spdAux::getRoute STParts]/group/value\[@n='Element'\]"
    set elements [$root selectNodes $xp1]
    set bool false
    foreach element_node $elements {
        set elemid [$element_node @v]
        set elem [Model::getElement $elemid]
        if {[$elem getAttribute "RotationDofs"]} {set bool true; break}
    }
    
    return $bool
}



Structural::write::Init
