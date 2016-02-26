namespace eval Solid::xml {
     variable dir
}

proc Solid::xml::Init { } {
    variable dir
    Model::InitVariables dir $Solid::dir
    
    getSolutionStrategies
    getElements
    getConditions
    getConstitutiveLaws
    getSolvers
    getProcesses
}

proc Solid::xml::getSolutionStrategies { } {
    Model::InitVariables SolutionStrategyFileName strategydefinition.xml
    Model::getSolutionStrategies
}

proc Solid::xml::getElements { } {
    Model::InitVariables ElementsFileName Elements.xml
    Model::getElements
}

proc Solid::xml::getConditions { } {
    Model::InitVariables ConditionsFileName Conditions.xml
    Model::getConditions
}

proc Solid::xml::getConstitutiveLaws { } {
    Model::InitVariables ConstitutiveLawsFileName ConstitutiveLaws.xml
    Model::getConstitutiveLaws
}

proc Solid::xml::getSolvers { } {
    Model::InitVariables SolversFileName Solvers.xml
    Model::getSolvers
}

proc Solid::xml::getProcesses { } {
    Model::InitVariables ProcessesFileName Processes.xml
    Model::getProcesses
}


proc Solid::xml::getUniqueName {name} {
    return SM$name
}

Solid::xml::Init
