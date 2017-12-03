namespace eval Structural::xml {
     variable dir
}

proc Structural::xml::Init { } {
     variable dir
     Model::InitVariables dir $Structural::dir

    Model::getSolutionStrategies Strategies.xml
    Model::getElements Elements.xml
    Model::getMaterials Materials.xml
    Model::getNodalConditions NodalConditions.xml
    Model::getConstitutiveLaws ConstitutiveLaws.xml
    Model::getProcesses "../../Common/xml/Processes.xml"
    Model::getProcesses DeprecatedProcesses.xml
    Model::getProcesses Processes.xml
    Model::getConditions Conditions.xml
    Model::getSolvers "../../Common/xml/Solvers.xml"
}

proc Structural::xml::getUniqueName {name} {
    return ST$name
}

proc ::Structural::xml::MultiAppEvent {args} {

}

proc Structural::xml::CustomTree { args } {
    spdAux::SetValueOnTreeItem state hidden Results CutPlanes
    spdAux::SetValueOnTreeItem v SingleFile GiDOptions GiDMultiFileFlag
    spdAux::SetValueOnTreeItem v 1 GiDOptions EchoLevel
    
    set result_node [[customlib::GetBaseRoot] selectNodes "[spdAux::getRoute NodalResults]/value\[@n = 'CONTACT'\]"]
    if {$result_node ne "" } {$result_node delete}
}

proc Structural::xml::ProcCheckGeometryStructural {domNode args} {
    set ret "line,surface"
    if {$::Model::SpatialDimension eq "3D"} {
        set ret "line,surface,volume"
    }
    return $ret
}


proc Structural::xml::ProcGetSolutionStrategiesSolid { domNode args } {
    set names ""
    set pnames ""
    set solutionType [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute STSoluType]] v]
    set Sols [::Model::GetSolutionStrategies [list "SolutionType" $solutionType] ]
    set ids [list ]
    foreach ss $Sols {
        lappend ids [$ss getName]
        append names [$ss getName] ","
        append pnames [$ss getName] "," [$ss getPublicName] ","
    }
    set names [string range $names 0 end-1]
    set pnames [string range $pnames 0 end-1]

    $domNode setAttribute values $names
    set dv [lindex $ids 0]
    if {[$domNode getAttribute v] eq ""} {$domNode setAttribute v $dv}
    if {[$domNode getAttribute v] ni $ids} {$domNode setAttribute v $dv}
    #spdAux::RequestRefresh
    return $pnames
}

proc Structural::xml::ProcCheckNodalConditionStateSolid {domNode args} {
    # Overwritten the base function to add Solution Type restrictions
    set parts_un STParts
    if {[spdAux::getRoute $parts_un] ne ""} {
        set conditionId [$domNode @n]
        set elems [$domNode selectNodes "[spdAux::getRoute $parts_un]/group/value\[@n='Element'\]"]
        set elemnames [list ]
        foreach elem $elems { lappend elemnames [$elem @v]}
        set elemnames [lsort -unique $elemnames]
        
        set solutionType [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute STSoluType]] v]
        set params [list analysis_type $solutionType]
        if {[::Model::CheckElementsNodalCondition $conditionId $elemnames $params]} {return "normal"} else {return "hidden"}
    } {return "normal"}
}

proc Structural::xml::ProcCheckGeometrySolid {domNode args} {
    set ret "surface"
    if {$::Model::SpatialDimension eq "3D"} {
        set ret "surface,volume"
    }
    return $ret
}


proc Structural::xml::UpdateParts {domNode args} {
    set current [lindex [$domNode selectNodes "./group"] end]
    set element_name [get_domnode_attribute [$current selectNodes "./value\[@n = 'Element'\]"] v]
    set element [Model::getElement $element_name]
    set LocalAxesAutomaticFunction [$element getAttribute "LocalAxesAutomaticFunction"]
    if {$LocalAxesAutomaticFunction ne ""} {
        $LocalAxesAutomaticFunction $current
    }
}

proc Structural::xml::AddLocalAxesToBeamElement { current } {
    # set y_axis_deviation [get_domnode_attribute [$current selectNodes "./value\[@n = 'LOCAL_AXIS_ROTATION'\]"] v]
    # W $y_axis_deviation
    set group [get_domnode_attribute $current n]
    if {[GiD_EntitiesGroups get $group lines -count]} {
        foreach line [GiD_EntitiesGroups get $group lines] {
            GiD_Process MEscape Data Conditions AssignCond line_Local_axes change -Automatic- $line escape escape
            #set raw [lindex [lindex [GiD_Info conditions -localaxesmat line_Local_axes mesh $line] 0] 3]
        }
    }
}


############# procs #################
proc Structural::xml::ProcGetElementsStructural { domNode args } {
    set nodeApp [spdAux::GetAppIdFromNode $domNode]
    set sol_stratUN [apps::getAppUniqueName $nodeApp SolStrat]
    set schemeUN [apps::getAppUniqueName $nodeApp Scheme]
    if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] v] eq ""} {
        get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $sol_stratUN]] dict
    }
    if {[get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] v] eq ""} {
        get_domnode_attribute [$domNode selectNodes [spdAux::getRoute $schemeUN]] dict
    }
    
    #W "solStrat $sol_stratUN sch $schemeUN"
    set solStratName [::write::getValue $sol_stratUN]
    set schemeName [write::getValue $schemeUN]
    #W "$solStratName $schemeName"
    #W "************************************************************************"
    #W "$nodeApp $solStratName $schemeName"
    set elems [::Model::GetAvailableElements $solStratName $schemeName]
    #W "************************************************************************"

    set analysis_type [get_domnode_attribute [$domNode selectNodes [spdAux::getRoute STAnalysisType]] v]
    set params [list AnalysisType $analysis_type]

    set names [list ]
    set pnames [list ]
    foreach elem $elems {
        if {[$elem cumple {*}$args]} {
            set available {*}[split [$elem getAttribute AnalysisType] {,}]
            if {$analysis_type in $available} {
                lappend names [$elem getName]
                lappend pnames [$elem getName] 
                lappend pnames [$elem getPublicName]
            }
        }
    }
    set diction [join $pnames ","]
    set values [join $names ","]
    #W "[get_domnode_attribute $domNode v] $names"
    $domNode setAttribute values $values
    if {[get_domnode_attribute $domNode v] eq ""} {$domNode setAttribute v [lindex $names 0]}
    if {[get_domnode_attribute $domNode v] ni $names} {$domNode setAttribute v [lindex $names 0]; spdAux::RequestRefresh}
    #spdAux::RequestRefresh
    return $diction
}

Structural::xml::Init
