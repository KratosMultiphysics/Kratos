namespace eval SolidMechanics::xml {
     variable dir
}

proc SolidMechanics::xml::Init { } {
    variable dir
    Model::InitVariables dir $SolidMechanics::dir
    
    getSolutionStrategies
    getElements
    getConditions
    getConstitutiveLaws
    getSolvers
    getProcesses
}

proc SolidMechanics::xml::getSolutionStrategies { } {
    Model::InitVariables SolutionStrategyFileName strategydefinition.xml
    Model::getSolutionStrategies
}

proc SolidMechanics::xml::getElements { } {
    Model::InitVariables ElementsFileName Elements.xml
    Model::getElements
}

proc SolidMechanics::xml::getConditions { } {
    Model::InitVariables ConditionsFileName Conditions.xml
    Model::getConditions
}

proc SolidMechanics::xml::getConstitutiveLaws { } {
    Model::InitVariables ConstitutiveLawsFileName ConstitutiveLaws.xml
    Model::getConstitutiveLaws
}

proc SolidMechanics::xml::getSolvers { } {
    Model::InitVariables SolversFileName Solvers.xml
    Model::getSolvers
}

proc SolidMechanics::xml::getProcesses { } {
    Model::InitVariables ProcessesFileName Processes.xml
    Model::getProcesses
}


proc SolidMechanics::xml::getUniqueName {name} {
    return SM$name
}

SolidMechanics::xml::Init
