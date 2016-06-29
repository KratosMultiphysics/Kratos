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

proc ::FSI::xml::MokChannelFlexibleWall {args} {
    Kratos::ResetModel
    GiD_Process Mescape 'Layers ChangeName Layer0 Fluid escape
    
    GiD_Process Mescape Geometry Create Line 0 0 0 0 5 0 17.5 5 0 17.5 3 0 12.5 3 0 5 0 0 5 2.5 0 4.95 2.5 0 4.95 0 0 Close escape escape
    GiD_Process Mescape Geometry Create NurbsSurface 1 2 3 4 5 6 7 8 9 escape escape 

    GiD_Process 'Zoom Frame 
    
    GiD_Process 'Layers New Structure escape 
    GiD_Process 'Layers Off Fluid escape
    GiD_Process 'Layers ToUse Structure escape 

    GiD_Process Mescape Geometry Create Line 5 0 0 5 2.5 0 4.95 2.5 0 4.95 0 0  Close escape escape
    GiD_Process Mescape Geometry Create NurbsSurface 10 11 12 13 escape escape
    
    GiD_Process 'Layers Color Fluid 047186223 Transparent Fluid 255 escape 'Layers Color Structure 187119038 Transparent Structure 255 escape 

    GiD_Process 'Layers On Fluid escape
    
    GiD_Groups create Fluid
    GiD_Groups create Structure
    GiD_Groups create Inlet
    GiD_Groups create Outlet
    GiD_Groups create NoSlip
    GiD_Groups create FluidInterface
    GiD_Groups create FixedDisplacement
    GiD_Groups create StructureInterface

    GiD_EntitiesGroups assign Fluid surfaces 1
    GiD_EntitiesGroups assign Structure surfaces 2
    GiD_EntitiesGroups assign Inlet lines 1
    GiD_EntitiesGroups assign Outlet lines 3
    GiD_EntitiesGroups assign NoSlip lines {2 4 5 9}
    GiD_EntitiesGroups assign FluidInterface lines {6 7 8}
    GiD_EntitiesGroups assign FixedDisplacement lines 13
    GiD_EntitiesGroups assign StructureInterface lines {10 11 12}
    
    GidUtils::UpdateWindow GROUPS
    
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']} group {n Fluid}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n Element pn Element dict {[GetElements]} actualize_tree 1 values FractionalStep2D state hidden v FractionalStep2D}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n ConstitutiveLaw pn {Fluid type} actualize_tree 1 values Newtonian,HerschelBulkley dict {[GetConstitutiveLaws]} state normal v Newtonian}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n DENSITY pn Density state {[PartParamState]} unit_magnitude Density help {} v 956.0 units kg/m^3}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n VISCOSITY pn {Kinematic viscosity} state {[PartParamState]} unit_magnitude L^2/T help {Fluidized viscosity.} v 0.145 units m^2/s}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n YIELD_STRESS pn {Yield Stress} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n POWER_LAW_K pn {Consistency index (k)} state {[PartParamState]} unit_magnitude {} help {} v 1}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n POWER_LAW_N pn {Flow index (n)} state {[PartParamState]} unit_magnitude {} help {} v 1}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n YOUNG_MODULUS pn {Young Modulus} state {[PartParamState]} unit_magnitude P help {} v 206.9e9 units Pa}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n POISSON_RATIO pn {Poisson Ratio} state {[PartParamState]} unit_magnitude {} help {} v 0.29}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n THICKNESS pn Thickness state {[PartParamState]} unit_magnitude L help {} v 1.0 units m}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n KINEMATIC_HARDENING_MODULUS pn {Kinematic Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n REFERENCE_HARDENING_MODULUS pn {Reference Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n INFINITY_HARDENING_MODULUS pn {Infinity Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n HARDENING_EXPONENT pn {Hardening Exponent} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n DAMAGE_THRESHOLD pn {Damage Threshold} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n STRENGTH_RATIO pn {Strength Ratio} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']/group[@n='Fluid']} value {n FRACTURE_ENERGY pn {Fracture Energy} state {[PartParamState]} unit_magnitude {} help {} v 0}
    
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='Inlet2D']} group {n Inlet}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='Inlet2D']/group[@n='Inlet']} value {n factor pn Modulus unit_magnitude Velocity help {} state {} v 0.6067 units m/s}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='Inlet2D']/group[@n='Inlet']} value {n directionX wn {Inlet2D _X} pn {Direction X} help {} state {} v 1.0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='Inlet2D']/group[@n='Inlet']} value {n directionY wn {Inlet2D _Y} pn {Direction Y} help {} state {} v 0.0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='Inlet2D']/group[@n='Inlet']} value {n directionZ wn {Inlet2D _Z} pn {Direction Z} help {} state {[CheckDimension 3D]} v 0.0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='Outlet2D']} group {n Outlet}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='Outlet2D']/group[@n='Outlet']} value {n value pn Value unit_magnitude P help {} state {} v 0.0 units Pa}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='NoSlip2D']} group {n NoSlip}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']/condition[@n='FluidNoSlipInterface2D']} group {n FluidInterface}

    
    gid_groups_conds::setAttributesF {container[@n='FSI']/container[@n='Structural']/container[@n='StageInfo']/value[@n='SolutionType']} {v Dynamic}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']} group {n Structure}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n Element pn Element actualize_tree 1 values SmallDisplacementElement2D,TotalLagrangianElement2D,UpdatedLagrangianElement2D,UpdatedLagrangianElementUP2D dict {[GetElements]} state normal v SmallDisplacementElement2D}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n ConstitutiveLaw pn {Constitutive law} actualize_tree 1 values LinearElasticPlaneStrain2DLaw,LinearElasticPlaneStress2DLaw,IsotropicDamageSimoJuPlaneStrain2DLaw,IsotropicDamageSimoJuPlaneStress2DLaw dict {LinearElasticPlaneStrain2DLaw,Linear Elastic Plane Strain,LinearElasticPlaneStress2DLaw,Linear Elastic Plane Stress,IsotropicDamageSimoJuPlaneStrain2DLaw,Isotropic Damage Simo-Ju Plane Strain,IsotropicDamageSimoJuPlaneStress2DLaw,Isotropic Damage Simo-Ju Plane Stress} state {[UpdateDictAndReturnState]} v LinearElasticPlaneStrain2DLaw}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n SECTION_TYPE pn {Section type} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n THICKNESS pn Thickness state {[PartParamState]} unit_magnitude L help {} v 1.0 units m}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n DENSITY pn Density state {[PartParamState]} unit_magnitude Density help {} v 1500.0 units kg/m^3}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n VISCOSITY pn {Kinematic viscosity} state {[PartParamState]} unit_magnitude L^2/T help {Fluidized viscosity.} v 1e-6 units m^2/s}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n YIELD_STRESS pn {Yield Stress} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n POWER_LAW_K pn {Consistency index (k)} state {[PartParamState]} unit_magnitude {} help {} v 1}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n POWER_LAW_N pn {Flow index (n)} state {[PartParamState]} unit_magnitude {} help {} v 1}
    gid_groups_conds::addF -resolve_parametric 1 {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n YOUNG_MODULUS pn {Young Modulus} state {[PartParamState]} unit_magnitude P help {} v 2.3e6 units Pa}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n POISSON_RATIO pn {Poisson Ratio} state {[PartParamState]} unit_magnitude {} help {} v 0.45}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n KINEMATIC_HARDENING_MODULUS pn {Kinematic Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n REFERENCE_HARDENING_MODULUS pn {Reference Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n INFINITY_HARDENING_MODULUS pn {Infinity Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n HARDENING_EXPONENT pn {Hardening Exponent} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n DAMAGE_THRESHOLD pn {Damage Threshold} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n STRENGTH_RATIO pn {Strength Ratio} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']/group[@n='Structure']} value {n FRACTURE_ENERGY pn {Fracture Energy} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/container[@n='Boundary Conditions']/condition[@n='DISPLACEMENT']} group {n FixedDisplacement ov line}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/container[@n='Boundary Conditions']/condition[@n='DISPLACEMENT']/group[@n='FixedDisplacement']} value {n FixX pn {X Imposed} values 1,0 help {} state {} v 1}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/container[@n='Boundary Conditions']/condition[@n='DISPLACEMENT']/group[@n='FixedDisplacement']} value {n FixY pn {Y Imposed} values 1,0 help {} state {} v 1}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/container[@n='Boundary Conditions']/condition[@n='DISPLACEMENT']/group[@n='FixedDisplacement']} value {n FixZ pn {Z Imposed} values 1,0 help {} state {[CheckDimension 3D]} v 1}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/container[@n='Boundary Conditions']/condition[@n='DISPLACEMENT']/group[@n='FixedDisplacement']} value {n valueX wn {DISPLACEMENT _X} pn {Value X} help {} state {} v 0.0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/container[@n='Boundary Conditions']/condition[@n='DISPLACEMENT']/group[@n='FixedDisplacement']} value {n valueY wn {DISPLACEMENT _Y} pn {Value Y} help {} state {} v 0.0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/container[@n='Boundary Conditions']/condition[@n='DISPLACEMENT']/group[@n='FixedDisplacement']} value {n valueZ wn {DISPLACEMENT _Z} pn {Value Z} help {} state {[CheckDimension 3D]} v 0.0}
    gid_groups_conds::addF {container[@n='FSI']/container[@n='Structural']/container[@n='Loads']/condition[@n='Interface2D']} group {n StructureInterface}
    
    spdAux::RequestRefresh
}


FSI::xml::Init
