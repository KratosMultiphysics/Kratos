namespace eval PotentialFluid::xml {
    # Namespace variables declaration
    variable dir
    variable lastImportMeshSize
    variable export_dir

}

proc PotentialFluid::xml::Init { } {
    # Namespace variables inicialization
    variable dir
    Model::DestroyEverything
    Model::InitVariables dir $PotentialFluid::dir

    Model::getSolutionStrategies Strategies.xml
    Model::getElements Elements.xml
    Model::getMaterials Materials.xml
    Model::getNodalConditions NodalConditions.xml
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    Model::getProcesses "../../Common/xml/Processes.xml"
    Model::getProcesses "../../Fluid/xml/Processes.xml"
    Model::getProcesses Processes.xml
    Model::getConditions Conditions.xml
    Model::getSolvers "../../Common/xml/Solvers.xml"
}


proc PotentialFluid::xml::MultiAppEvent {args} {
   if {$args eq "init"} {
     spdAux::parseRoutes
     spdAux::ConvertAllUniqueNames FL ${::PotentialFluid::prefix}
   }
}

proc PotentialFluid::xml::getUniqueName {name} {
    return ${::PotentialFluid::prefix}${name}
}

proc PotentialFluid::xml::CustomTree { args } {
    # Hide Results Cut planes
    #Fluid::xml::CustomTree {*}$args
    
    spdAux::SetValueOnTreeItem state hidden PTFLGravity
    
}

proc spdAux::injectConditions { basenode args} {
    set conditions [::Model::GetConditions [list ImplementedInApplication {FluidApplication CompressiblePotentialFlowApplication}]]
    spdAux::_injectCondsToTree $basenode $conditions
    $basenode delete
}

PotentialFluid::xml::Init
