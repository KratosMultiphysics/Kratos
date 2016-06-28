namespace eval FSI::xml {
    # Namespace variables declaration
    variable dir
}

proc FSI::xml::Init { } {
    # Namespace variables initialization
    variable dir
    Model::InitVariables dir $FSI::dir
    
    Model::ForgetSolutionStrategies
    Model::getSolutionStrategies "../../Fluid/xml/Strategies.xml"
    Model::getSolutionStrategies "../../Structural/xml/Strategies.xml"
    Model::getSolutionStrategies Strategies.xml
    Model::getProcesses Processes.xml
    Model::getConditions Conditions.xml
    
    Model::ForgetSolvers
    Model::getSolvers "../../Common/xml/Solvers.xml"
    Model::getSolvers Coupling_solvers.xml
}

proc FSI::xml::getUniqueName {name} {
    return ${::FSI::prefix}${name}
}


proc ::FSI::xml::MultiAppEvent {args} {
   if {$args eq "init"} {
        spdAux::parseRoutes
        spdAux::ConvertAllUniqueNames SL ST
   }
}


FSI::xml::Init
