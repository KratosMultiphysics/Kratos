namespace eval Fluid::xml {
    # Namespace variables declaration
    variable dir
}

proc Fluid::xml::Init { } {
    # Namespace variables inicialization
    variable dir
    Model::InitVariables dir $Fluid::dir
    
    getSolutionStrategies
    getElements
    getConstitutiveLaws
    getConditions
    #getSolvers
}


proc Fluid::xml::getSolutionStrategies { } {
    Model::InitVariables SolutionStrategyFileName strategydefinition.xml
    Model::getSolutionStrategies
}

proc Fluid::xml::getElements { } {
    Model::InitVariables ElementsFileName Elements.xml
    Model::getElements
}
proc Fluid::xml::getConstitutiveLaws { } {
    Model::InitVariables ConstitutiveLawsFileName ConstitutiveLaws.xml
    Model::getConstitutiveLaws
}

proc Fluid::xml::getConditions { } {
    Model::InitVariables ConditionsFileName Conditions.xml
    Model::getConditions
}

proc Fluid::xml::getUniqueName {name} {
    return FL$name
}

Fluid::xml::Init
