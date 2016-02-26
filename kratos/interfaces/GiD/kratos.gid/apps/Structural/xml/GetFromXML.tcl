namespace eval Structural::xml {
     variable dir
}

proc Structural::xml::Init { } {
    variable dir
    Model::InitVariables dir $Structural::dir
    
    getSolutionStrategies
    getElements
    getConditions
    getConstitutiveLaws
    getSolvers
    getProcesses
}

proc Structural::xml::getSolutionStrategies { } {
    Model::InitVariables SolutionStrategyFileName strategydefinition.xml
    Model::getSolutionStrategies
}

proc Structural::xml::getElements { } {
    Model::InitVariables ElementsFileName Elements.xml
    Model::getElements
}

proc Structural::xml::getConditions { } {
    Model::InitVariables ConditionsFileName Conditions.xml
    Model::getConditions
}

proc Structural::xml::getConstitutiveLaws { } {
    Model::InitVariables ConstitutiveLawsFileName ConstitutiveLaws.xml
    Model::getConstitutiveLaws
}

proc Structural::xml::getSolvers { } {
    Model::InitVariables SolversFileName Solvers.xml
    Model::getSolvers
}

proc Structural::xml::getProcesses { } {
    Model::InitVariables ProcessesFileName Processes.xml
    Model::getProcesses
}


proc Structural::xml::getUniqueName {name} {
    return SM$name
}

Structural::xml::Init
