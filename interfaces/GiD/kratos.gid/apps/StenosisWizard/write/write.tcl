namespace eval StenosisWizard::write {

}

proc StenosisWizard::write::Init { } {
    
}


proc StenosisWizard::write::writeCustomFilesEvent { } {
    return [Fluid::write::writeCustomFilesEvent]
}

# MDPA Blocks

proc StenosisWizard::write::writeModelPartEvent { } {
    Fluid::write::writeModelPartEvent
}

# Project Parameters
proc StenosisWizard::write::writeParametersEvent { } {
    set project_parameters_dict [::Fluid::write::getParametersDict]
    write::WriteJSON $project_parameters_dict
}

StenosisWizard::write::Init
