##################################################################################
#   This file is common for all Kratos Applications.
#   Do not change anything here unless it's strictly necessary.
##################################################################################

namespace eval Model {
    variable SpatialDimension
    variable SolutionStrategies
    variable Elements
    variable Conditions
    variable ConstitutiveLaws
    variable Solvers
    variable Processes
    
    variable dir
}

proc Model::Init { } {
    variable SpatialDimension
    variable dir
    variable SolutionStrategies
    variable Elements
    variable Conditions
    variable ConstitutiveLaws
    variable Solvers
    variable Processes
    
    set dir $::Kratos::kratos_private(Path)
    
    set SolutionStrategies [list ]
    set Elements [list ]
    set Conditions [list ]
    set ConstitutiveLaws [list ]
    set Solvers [list ]
    set Processes [list ]
    
    set SpatialDimension "3D"
}

proc Model::InitVariables {varName varValue} {
    catch {
        set ::Model::$varName $varValue
    }
}

proc Model::getSolutionStrategies { SolutionStrategyFileName } {
    variable SolutionStrategies
    variable dir
    
    set SolutionStrategies [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $SolutionStrategyFileName]] doc
    
    ParseSolutionStrategies $doc
}

proc Model::getElements { ElementsFileName } {
    variable Elements
    variable dir
    
    set Elements [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $ElementsFileName]] doc
    
    ParseElements $doc
}
proc Model::getConditions { ConditionsFileName } {
    variable Conditions
    variable dir
    
    set Conditions [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $ConditionsFileName]] doc
    
    ParseConditions $doc
}

proc Model::getConstitutiveLaws { ConstitutiveLawsFileName } {
    variable ConstitutiveLaws
    variable dir
    
    set ConstitutiveLaws [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $ConstitutiveLawsFileName]] doc
    
    ParseConstitutiveLaws $doc
}

proc Model::getSolvers { SolversFileName } {
    variable Solvers
    variable dir
    
    set Solvers [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $SolversFileName]] doc
    
    ParseSolvers $doc
}

proc Model::getProcesses { ProcessesFileName } {
    variable Processes
    variable dir
    
    set Processes [list ]
    dom parse [tDOM::xmlReadFile [file join $dir xml $ProcessesFileName]] doc
    
    ParseProcesses $doc
}

Model::Init