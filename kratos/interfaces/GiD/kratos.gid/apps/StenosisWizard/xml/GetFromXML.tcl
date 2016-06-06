namespace eval StenosisWizard::xml {
     variable dir
}

proc StenosisWizard::xml::Init { } {
     variable dir
     Model::InitVariables dir $StenosisWizard::dir
    
     Model::getSolutionStrategies Strategies.xml
     Model::getProcesses Processes.xml
     Model::getElements Elements.xml
     Model::getNodalConditions NodalConditions.xml
     Model::getConstitutiveLaws ConstitutiveLaws.xml
     Model::getConditions Conditions.xml
     Model::getSolvers Solvers.xml
     
}

proc StenosisWizard::xml::getUniqueName {name} {
    return StenWiz$name
}

proc ::StenosisWizard::xml::MultiAppEvent {args} {
   if {$args eq "init"} {
     spdAux::parseRoutes
     spdAux::ConvertAllUniqueNames FL StenWiz
   }
}

StenosisWizard::xml::Init
