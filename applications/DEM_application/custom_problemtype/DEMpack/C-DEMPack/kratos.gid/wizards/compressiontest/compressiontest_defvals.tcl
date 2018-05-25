namespace eval ::defvals {

}

namespace eval ::defvals::dem::wizdata {


}
proc ::defvals::dem::wizdata::init {} {
   
    variable TestType 
    set TestType 1
    
    variable ExperimentType 
    set ExperimentType "UCS"
    
    variable MeshType 
    set MeshType 1
    
    variable CalculationTime
    set CalculationTime "3.0e-3"
    
    variable DeltaTimeSelection 
    set DeltaTimeSelection "Automatic"

    variable DeltaTime 
    set DeltaTime "1e-7"
    
    variable LoadingVelocity 
    set LoadingVelocity "0.1"
    
    variable Pressure
    set Pressure "0.0"
    
    variable Displacement 
    set Displacement "Yes"
    
    variable Velocity 
    set Velocity "Yes"
    
    variable TotalForces 
    set TotalForces "No"
    
    variable Rhs 
    set Rhs "No"
    
    variable DampForces 
    set DampForces "No"
    
    variable AppliedForces 
    set AppliedForces "No"
    
    variable ContactSigma 
    set ContactSigma "No"
    
    variable ContactTau 
    set ContactTau "No"
    
    variable LocalContactForce 
    set LocalContactForce "No"
    
    variable FailureCriterionState 
    set FailureCriterionState "No"
    
    variable SinglePostprocessFile 
    set SinglePostprocessFile "Single"
    
    variable PrintOutputFile 
    set PrintOutputFile "2e-5"
    
    variable PrintGraphData 
    set PrintGraphData "1e-6"
    
    variable NumberOfThreads 
    set NumberOfThreads 1
    
    variable SelectOMPMPI 
    set SelectOMPMPI "OpenMP"
    
}

::defvals::dem::wizdata::init


namespace eval ::defvals::dem::wizmaterial {


}

proc ::defvals::dem::wizmaterial::init {} {
    
variable Density 
set Density 2300.0

variable YoungModulus 
set YoungModulus 2.8e10

variable PoissonRatio 
set PoissonRatio 0.20

variable ParticleFrictionAngle 
set ParticleFrictionAngle 0.25

variable LCS1 
set LCS1 20

variable LCS2 
set LCS2 30

variable LCS3 
set LCS3 60

variable YRC1 
set YRC1 4

variable YRC2 
set YRC2 6

variable YRC3 
set YRC3 22

variable PlasticYoungModulus 
set PlasticYoungModulus 2.8e10

variable PlasticYieldStress 
set PlasticYieldStress 20.0

variable DamageDeformationFactor 
set DamageDeformationFactor 0.2

variable TangentialStrength 
set TangentialStrength 18.0

variable NormalTensileStrength 
set NormalTensileStrength 4.0

variable InternalFrictionAngleCoeff 
set InternalFrictionAngleCoeff 0.9

}

::defvals::dem::wizmaterial::init