namespace eval BaseApp::xml {
     variable dir
}

proc BaseApp::xml::Init { } {
    variable dir
    Model::InitVariables dir $BaseApp::dir
    
    getSolutionStrategies
    getElements
    getConditions
    getConstitutiveLaws
    getSolvers
    getProcesses
}

proc BaseApp::xml::getSolutionStrategies { } {
    Model::InitVariables SolutionStrategyFileName strategydefinition.xml
    Model::getSolutionStrategies
}

proc BaseApp::xml::getElements { } {
    Model::InitVariables ElementsFileName Elements.xml
    Model::getElements
}

proc BaseApp::xml::getConditions { } {
    Model::InitVariables ConditionsFileName Conditions.xml
    Model::getConditions
}

proc BaseApp::xml::getConstitutiveLaws { } {
    Model::InitVariables ConstitutiveLawsFileName ConstitutiveLaws.xml
    Model::getConstitutiveLaws
}

proc BaseApp::xml::getSolvers { } {
    Model::InitVariables SolversFileName Solvers.xml
    Model::getSolvers
}

proc BaseApp::xml::getProcesses { } {
    Model::InitVariables ProcessesFileName Processes.xml
    Model::getProcesses
}


proc BaseApp::xml::getUniqueName {name} {
    return BA$name
}

BaseApp::xml::Init
