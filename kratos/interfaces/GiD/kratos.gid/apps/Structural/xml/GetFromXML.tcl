namespace eval Structural::xml {
     variable dir
}

proc Structural::xml::Init { } {
     variable dir
     Model::InitVariables dir $Structural::dir
    
     Model::getSolutionStrategies strategydefinition.xml
     Model::getElements Elements.xml
     Model::getConstitutiveLaws ConstitutiveLaws.xml
     Model::getProcesses Processes.xml
     Model::getConditions Conditions.xml
     Model::getSolvers Solvers.xml
}

proc Structural::xml::getUniqueName {name} {
    return ST$name
}

Structural::xml::Init
