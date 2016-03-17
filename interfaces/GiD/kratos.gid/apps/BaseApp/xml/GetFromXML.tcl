namespace eval BaseApp::xml {
     variable dir
}

proc BaseApp::xml::Init { } {
    variable dir
    Model::InitVariables dir $BaseApp::dir
    
    Model::getSolutionStrategies strategydefinition.xml
    Model::getElements Elements.xml
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    Model::getProcesses Processes.xml
    Model::getConditions Conditions.xml
    Model::getSolvers Solvers.xml
}

proc BaseApp::xml::getUniqueName {name} {
    return BA$name
}

BaseApp::xml::Init
