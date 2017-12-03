namespace eval Fluid::xml {
    # Namespace variables declaration
    variable dir
}

proc Fluid::xml::Init { } {
    # Namespace variables inicialization
    variable dir
    Model::InitVariables dir $Fluid::dir
    
    Model::getSolutionStrategies Strategies.xml
    Model::getElements Elements.xml
    Model::getMaterials Materials.xml
    Model::getNodalConditions NodalConditions.xml
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    Model::getProcesses "../../Common/xml/Processes.xml"
    Model::getProcesses Processes.xml
    Model::getConditions Conditions.xml
    Model::getSolvers "../../Common/xml/Solvers.xml"
}

proc Fluid::xml::getUniqueName {name} {
    return ${::Fluid::prefix}${name}
}

proc Fluid::xml::CustomTree { args } {
    set root [customlib::GetBaseRoot]

    # Output control in output settings
    spdAux::SetValueOnTreeItem v time Results FileLabel
    spdAux::SetValueOnTreeItem v time Results OutputControlType

    # Drag in output settings
    if {[$root selectNodes "[spdAux::getRoute Results]/condition\[@n='Drag'\]"] eq ""} {
        gid_groups_conds::addF [spdAux::getRoute Results] include [list n Drag active 1 path {apps/Fluid/xml/Drag.spd}]
    }
    
    customlib::ProcessIncludes $::Kratos::kratos_private(Path)
    spdAux::parseRoutes

    # Nodal reactions in output settings
    if {[$root selectNodes "[spdAux::getRoute Results]/container\[@n='OnNodes'\]"] ne ""} {
        gid_groups_conds::addF "[spdAux::getRoute Results]/container\[@n='OnNodes'\]" value [list n REACTION pn "Reaction" v No values "Yes,No"]
    }
}

Fluid::xml::Init
