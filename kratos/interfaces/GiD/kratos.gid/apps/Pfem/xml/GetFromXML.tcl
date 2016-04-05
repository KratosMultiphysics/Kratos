namespace eval Pfem::xml {
     variable dir
}

proc Pfem::xml::Init { } {
    variable dir
    Model::InitVariables dir $Pfem::dir
    
    Model::getSolutionStrategies strategydefinition.xml
    Model::getElements Elements.xml
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    Model::getProcesses Processes.xml
    Model::getConditions Conditions.xml
    Model::getSolvers Solvers.xml
}

proc Pfem::xml::getUniqueName {name} {
    return BA$name
}

Pfem::xml::Init
