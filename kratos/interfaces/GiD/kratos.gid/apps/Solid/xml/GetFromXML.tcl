namespace eval Solid::xml {
     variable dir
}

proc Solid::xml::Init { } {
     variable dir
     Model::InitVariables dir $Solid::dir

     Model::getSolutionStrategies Strategies.xml
     Model::getElements Elements.xml
     Model::getNodalConditions NodalConditions.xml
     Model::getConstitutiveLaws ConstitutiveLaws.xml
     Model::getProcesses Processes.xml
     Model::getConditions Conditions.xml
     Model::getSolvers "../../Common/xml/Solvers.xml"
     #Model::getSolvers Solvers.xml
}

proc Solid::xml::getUniqueName {name} {
    return SL$name
}

Solid::xml::Init
