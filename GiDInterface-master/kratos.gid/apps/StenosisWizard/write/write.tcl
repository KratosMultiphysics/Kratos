namespace eval StenosisWizard::write {

}

proc StenosisWizard::write::Init { } {
    
}


proc StenosisWizard::write::writeCustomFilesEvent { } {
    write::CopyFileIntoModel "../Fluid/python/KratosFluid.py"
    write::RenameFileInModel "KratosFluid.py" "MainKratos.py"
}

# MDPA Blocks

proc StenosisWizard::write::writeModelPartEvent { } {
    Fluid::write::AddValidApps StenosisWizard
    write::writeAppMDPA Fluid
}

# Project Parameters
proc StenosisWizard::write::writeParametersEvent { } {
    set project_parameters_dict [::Fluid::write::getParametersDict]
    write::WriteJSON $project_parameters_dict
}

StenosisWizard::write::Init
