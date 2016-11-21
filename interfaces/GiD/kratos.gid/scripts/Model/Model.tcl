##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Model {
    variable SpatialDimension
    variable ValidSpatialDimensions
    variable SolutionStrategies
    variable Materials
    variable Elements
    variable Conditions
    variable NodalConditions
    variable ConstitutiveLaws
    variable Solvers
    variable Processes
    
    variable dir
}

proc Model::Init { } {
    variable SpatialDimension
    variable ValidSpatialDimensions
    variable dir
    variable SolutionStrategies
    variable Elements
    variable Materials
    variable Conditions
    variable NodalConditions
    variable ConstitutiveLaws
    variable Solvers
    variable Processes
    
    set dir $::Kratos::kratos_private(Path)
    
    set SolutionStrategies [list ]
    set Elements [list ]
    set Materials [list ]
    set Conditions [list ]
    set NodalConditions [list ]
    set ConstitutiveLaws [list ]
    set Solvers [list ]
    set Processes [list ]
    
    set SpatialDimension "3D"
    set ValidSpatialDimensions [list 2D 3D]
}

proc Model::InitVariables {varName varValue} {
    catch {
        set ::Model::$varName $varValue
    }
}

proc Model::getSolutionStrategies { SolutionStrategyFileName } {
    #variable SolutionStrategies
    variable dir
    
    #set SolutionStrategies [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $SolutionStrategyFileName]] doc
    
    ParseSolutionStrategies $doc
}

proc Model::getElements { ElementsFileName } {
    #variable Elements
    variable dir
    
    #set Elements [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $ElementsFileName]] doc
    
    ParseElements $doc
}
proc Model::getConditions { ConditionsFileName } {
    #variable Conditions
    variable dir
    
    #set Conditions [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $ConditionsFileName]] doc
    
    ParseConditions $doc
}
proc Model::getNodalConditions { NodalConditionsFileName } {
    #variable NodalConditions
    variable dir
    
    #set Conditions [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $NodalConditionsFileName]] doc
    #W [$doc asXML]
    ParseNodalConditions $doc
}

proc Model::getConstitutiveLaws { ConstitutiveLawsFileName } {
    #variable ConstitutiveLaws
    variable dir
    
    #set ConstitutiveLaws [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $ConstitutiveLawsFileName]] doc
    
    ParseConstitutiveLaws $doc
}

proc Model::getSolvers { SolversFileName } {
    #variable Solvers
    variable dir
    
    #set Solvers [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $SolversFileName]] doc
    
    ParseSolvers $doc
}

proc Model::getProcesses { ProcessesFileName } {
    #variable Processes
    variable dir
    
    #set Processes [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $ProcessesFileName]] doc
    
    ParseProcesses $doc
}
proc Model::getMaterials { MaterialsFileName } {
    variable dir
    dom parse [tDOM::xmlReadFile [file join $dir xml $MaterialsFileName]] doc
    ParseMaterials $doc
}

proc Model::DestroyEverything { } {
    Init
}

Model::Init