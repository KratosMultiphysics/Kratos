
proc ::FSI::examples::MokChannelFlexibleWall {args} {
    DrawMokChannelFlexibleWallGeometry
    TreeAssignationMokChannelFlexibleWall
}

proc FSI::examples::DrawMokChannelFlexibleWallGeometry {args} {
    Kratos::ResetModel
    GiD_Process Mescape 'Layers ChangeName Layer0 Fluid escape
    
    # Geometry creation
    set coordinates [list 0.5 0 0 0 0 0 0 0.5 0 1.75 0.5 0 1.75 0.3 0 1.35 0.3 0 1.3103 0.29881 0 1.2709 0.29453 0 1.232 0.28693 0 1.1937 0.27644 0 1.1565 0.26288 0 0.505 0 0]
    set fluidPoints [list ]
    foreach {x y z} $coordinates {
        lappend fluidPoints [GiD_Geometry create point append Fluid $x $y $z]
    }
    set coordinates [list 0.505 0.25 0 0.5 0.25 0]
    set fluidnterfacePoints [list ]
    foreach {x y z} $coordinates {
        lappend fluidnterfacePoints [GiD_Geometry create point append Fluid $x $y $z]
    }
    
    set fluidLines [list ]
    set initial [lindex $fluidPoints 0]
    foreach point [lrange $fluidPoints 1 end] {
        lappend fluidLines [GiD_Geometry create line append stline Fluid $initial $point]
        set initial $point
    }
    set fluidinteractionLines [list ]
    foreach point $fluidnterfacePoints {
        lappend fluidinteractionLines [GiD_Geometry create line append stline Fluid $initial $point]
        set initial $point
    }
    lappend fluidinteractionLines [GiD_Geometry create line append stline Fluid $initial [lindex $fluidPoints 0]]
    #set fluidSurface [GiD_Geometry create surface append plsurface Fluid [llength $fluidLines] {*}$fluidLines]
    set fluidalllines $fluidLines
    lappend fluidalllines {*}$fluidinteractionLines
    GiD_Process Mescape Geometry Create NurbsSurface {*}$fluidalllines escape escape 

    
    GiD_Process 'Layers New Structure escape 
    GiD_Process 'Layers Off Fluid escape
    GiD_Process 'Layers ToUse Structure escape
    
    
    set coordinates [list 0.505 0 0 0.505 0.25 0 0.5 0.25 0 0.5 0 0 ]
    set strucPoints [list ]
    foreach {x y z} $coordinates {
        lappend strucPoints [GiD_Geometry create point append Structure $x $y $z]
    }
    set strucLines [list ]
    set initial [lindex $strucPoints 0]
    foreach point [lrange $strucPoints 1 end] {
        lappend strucLines [GiD_Geometry create line append stline Structure $initial $point]
        set initial $point
    }
    lappend strucLines [GiD_Geometry create line append stline Structure $initial [lindex $strucPoints 0]]
    GiD_Process Mescape Geometry Create NurbsSurface {*}$strucLines escape escape 

    GiD_Process 'Layers Color Fluid 047186223 Transparent Fluid 255 escape 'Layers Color Structure 187119038 Transparent Structure 255 escape 
    GiD_Process 'Layers On Fluid escape
    
    if {$::Model::SpatialDimension eq "3D"} {
        GiD_Process 'Layers Off Structure escape Mescape
        GiD_Process Utilities Copy Surfaces Duplicate DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,0.5 1 escape Mescape
        GiD_Process 'Layers On Structure escape 'Layers Off Fluid escape Mescape
        GiD_Process Utilities Copy Surfaces Duplicate DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,0.5 2 escape Mescape
        GiD_Process 'Layers On Fluid escape 
        GiD_Process 'Layers Transparent Fluid 127 escape 
    }
    GiD_Process 'Zoom Frame
    GiD_Process 'Render Flat escape 
    
    # Group creation
    GiD_Groups create Fluid
    GiD_Groups create Structure
    GiD_Groups create Inlet
    GiD_Groups create Outlet
    GiD_Groups create NoSlip
    GiD_Groups create Slip
    GiD_Groups create FluidInterface
    GiD_Groups create FixedDisplacement
    GiD_Groups create StructureInterface

    # Group entities
    if {$::Model::SpatialDimension eq "3D"} {
        GiD_Groups create FluidFixedDisplacement_full
        GiD_Groups create FluidFixedDisplacement_lat
        GiD_EntitiesGroups assign Fluid volumes 1
        GiD_EntitiesGroups assign Structure volumes 2
        GiD_EntitiesGroups assign Inlet surfaces 4
        GiD_EntitiesGroups assign Outlet surfaces 6
        GiD_EntitiesGroups assign NoSlip surfaces {3 7 8 9 10 11 12 13}
        GiD_EntitiesGroups assign Slip surfaces {1 5 17} 
        GiD_EntitiesGroups assign FluidFixedDisplacement_full surfaces {3 4 5 6 7 8 9 10 11 12 13}
        GiD_EntitiesGroups assign FluidFixedDisplacement_lat surfaces {1 17}
        GiD_EntitiesGroups assign FluidInterface surfaces {14 15 16}
        GiD_EntitiesGroups assign FixedDisplacement surfaces {21}
        GiD_EntitiesGroups assign StructureInterface surfaces {18 19 20}
        
    } {
        GiD_Groups create FluidALEMeshBC
        GiD_EntitiesGroups assign Fluid surfaces 1
        GiD_EntitiesGroups assign Structure surfaces 2
        GiD_EntitiesGroups assign Inlet lines 2
        GiD_EntitiesGroups assign Outlet lines 4
        GiD_EntitiesGroups assign NoSlip lines $fluidLines
        GiD_EntitiesGroups unassign NoSlip lines {2 3 4}
        GiD_EntitiesGroups assign Slip lines 3
        GiD_EntitiesGroups assign FluidInterface lines $fluidinteractionLines
        GiD_EntitiesGroups assign FluidALEMeshBC lines $fluidLines
        GiD_EntitiesGroups assign FixedDisplacement lines [lindex $strucLines end]
        GiD_EntitiesGroups assign StructureInterface lines [lrange $strucLines 0 end-1]
    }
    GidUtils::UpdateWindow GROUPS
}

proc FSI::examples::TreeAssignationMokChannelFlexibleWall {args} {
    set nd $::Model::SpatialDimension
    
    set condtype line
    if {$::Model::SpatialDimension eq "3D"} { set condtype surface }
    # Fluid Parts
    set fluidParts {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']}
    gid_groups_conds::addF $fluidParts group {n Fluid}
    set fluidGroup "$fluidParts/group\[@n='Fluid'\]"
    gid_groups_conds::addF $fluidGroup value "n Element pn Element dict {\[GetElements\]} actualize_tree 1 values FractionalStep$nd state hidden v FractionalStep$nd"
    gid_groups_conds::addF $fluidGroup value {n ConstitutiveLaw pn {Fluid type} actualize_tree 1 values Newtonian,HerschelBulkley dict {[GetConstitutiveLaws]} state normal v Newtonian}
    gid_groups_conds::addF -resolve_parametric 1 $fluidGroup value {n DENSITY pn Density state {[PartParamState]} unit_magnitude Density help {} v 956.0 units kg/m^3}
    gid_groups_conds::addF -resolve_parametric 1 $fluidGroup value {n VISCOSITY pn {Kinematic viscosity} state {[PartParamState]} unit_magnitude L^2/T help {Fluidized viscosity.} v 0.145 units m^2/s}
    gid_groups_conds::addF $fluidGroup value {n YIELD_STRESS pn {Yield Stress} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF $fluidGroup value {n POWER_LAW_K pn {Consistency index (k)} state {[PartParamState]} unit_magnitude {} help {} v 1}
    gid_groups_conds::addF $fluidGroup value {n POWER_LAW_N pn {Flow index (n)} state {[PartParamState]} unit_magnitude {} help {} v 1}

    set fluidConditions {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']}
    # Fluid Interface
    set fluidInlet "$fluidConditions/condition\[@n='Inlet$nd'\]"
    
    # Fluid Inlet
    gid_groups_conds::addF $fluidInlet group "n Inlet ov $condtype"
    set fluidInletGroup "$fluidInlet/group\[@n='Inlet'\]"
    gid_groups_conds::addF -resolve_parametric 1 $fluidInletGroup value {n factor pn Modulus unit_magnitude Velocity help {} state {} v 0.6067 units m/s}
    gid_groups_conds::addF $fluidInletGroup value "n directionX wn {Inlet$nd _X} pn {Direction X} help {} state {} v 1.0"
    gid_groups_conds::addF $fluidInletGroup value "n directionY wn {Inlet$nd _Y} pn {Direction Y} help {} state {} v 0.0"
    gid_groups_conds::addF $fluidInletGroup value "n directionZ wn {Inlet$nd _Z} pn {Direction Z} help {} state {\[CheckDimension 3D\]} v 0.0"
        
    # Fluid Outlet
    set fluidOutlet "$fluidConditions/condition\[@n='Outlet$nd'\]"
    gid_groups_conds::addF $fluidOutlet group "n Outlet ov $condtype"
    gid_groups_conds::addF -resolve_parametric 1 "$fluidOutlet/group\[@n='Outlet'\]" value {n value pn Value unit_magnitude P help {} state {} v 0.0 units Pa}
    
    # Fluid Conditions
    gid_groups_conds::addF "$fluidConditions/condition\[@n='NoSlip$nd'\]" group "n NoSlip ov $condtype"
    gid_groups_conds::addF "$fluidConditions/condition\[@n='Slip$nd'\]" group "n Slip ov $condtype"
    gid_groups_conds::addF "$fluidConditions/condition\[@n='FluidNoSlipInterface$nd'\]" group "n FluidInterface ov $condtype"
    
    # Displacement 3D
    if {$nd eq "3D"} {
        set fluidDisplacement "$fluidConditions/condition\[@n='ALEMeshDisplacementBC3D'\]"
        gid_groups_conds::addF $fluidDisplacement group {n FluidFixedDisplacement_full ov surface}
        set fluidDisplacementGroup_full "$fluidDisplacement/group\[@n='FluidFixedDisplacement_full'\]"
        gid_groups_conds::addF $fluidDisplacementGroup_full value {n FixX pn {X Imposed} values 1,0 help {} state {} v 1}
        gid_groups_conds::addF $fluidDisplacementGroup_full value {n FixY pn {Y Imposed} values 1,0 help {} state {} v 1}
        gid_groups_conds::addF $fluidDisplacementGroup_full value {n FixZ pn {Z Imposed} values 1,0 help {} state {[CheckDimension 3D]} v 1}
        gid_groups_conds::addF $fluidDisplacementGroup_full value {n valueX wn {DISPLACEMENT _X} pn {Value X} help {} state {} v 0.0}
        gid_groups_conds::addF $fluidDisplacementGroup_full value {n valueY wn {DISPLACEMENT _Y} pn {Value Y} help {} state {} v 0.0}
        gid_groups_conds::addF $fluidDisplacementGroup_full value {n valueZ wn {DISPLACEMENT _Z} pn {Value Z} help {} state {[CheckDimension 3D]} v 0.0}
        gid_groups_conds::addF $fluidDisplacement group {n FluidFixedDisplacement_lat ov surface}
        set fluidDisplacementGroup_lat "$fluidDisplacement/group\[@n='FluidFixedDisplacement_lat'\]"
        gid_groups_conds::addF $fluidDisplacementGroup_lat value {n FixX pn {X Imposed} values 1,0 help {} state {} v 0}
        gid_groups_conds::addF $fluidDisplacementGroup_lat value {n FixY pn {Y Imposed} values 1,0 help {} state {} v 0}
        gid_groups_conds::addF $fluidDisplacementGroup_lat value {n FixZ pn {Z Imposed} values 1,0 help {} state {[CheckDimension 3D]} v 1}
        gid_groups_conds::addF $fluidDisplacementGroup_lat value {n valueX wn {DISPLACEMENT _X} pn {Value X} help {} state {} v 0.0}
        gid_groups_conds::addF $fluidDisplacementGroup_lat value {n valueY wn {DISPLACEMENT _Y} pn {Value Y} help {} state {} v 0.0}
        gid_groups_conds::addF $fluidDisplacementGroup_lat value {n valueZ wn {DISPLACEMENT _Z} pn {Value Z} help {} state {[CheckDimension 3D]} v 0.0}
    } {
        set fluidDisplacement "$fluidConditions/condition\[@n='ALEMeshDisplacementBC2D'\]"
        gid_groups_conds::addF $fluidDisplacement group {n FluidALEMeshBC ov line}
        set fluidDisplacementGroup "$fluidDisplacement/group\[@n='FluidALEMeshBC'\]"
        gid_groups_conds::addF $fluidDisplacementGroup value {n FixX pn {X Imposed} values 1,0 help {} state {} v 1}
        gid_groups_conds::addF $fluidDisplacementGroup value {n FixY pn {Y Imposed} values 1,0 help {} state {} v 1}
        gid_groups_conds::addF $fluidDisplacementGroup value {n FixZ pn {Z Imposed} values 1,0 help {} state {[CheckDimension 3D]} v 1}
        gid_groups_conds::addF $fluidDisplacementGroup value {n valueX wn {DISPLACEMENT _X} pn {Value X} help {} state {} v 0.0}
        gid_groups_conds::addF $fluidDisplacementGroup value {n valueY wn {DISPLACEMENT _Y} pn {Value Y} help {} state {} v 0.0}
        gid_groups_conds::addF $fluidDisplacementGroup value {n valueZ wn {DISPLACEMENT _Z} pn {Value Z} help {} state {[CheckDimension 3D]} v 0.0}
    }

    # Structural
    gid_groups_conds::setAttributesF {container[@n='FSI']/container[@n='Structural']/container[@n='StageInfo']/value[@n='SolutionType']} {v Dynamic}
    
    # Structural Parts
    set structParts {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']}
    if {$nd eq "3D"} {
        gid_groups_conds::addF $structParts group {n Structure ov volume}
    } {
        gid_groups_conds::addF $structParts group {n Structure ov surface}
    }
    
    set structPartsGroup "$structParts/group\[@n='Structure'\]"
    gid_groups_conds::addF $structPartsGroup value "n Element pn Element actualize_tree 1 dict {\[GetElements\]} state normal v SmallDisplacementElement$nd"
    gid_groups_conds::addF $structPartsGroup value "n ConstitutiveLaw pn {Constitutive law} actualize_tree 1 dict \[GetConstitutiveLaws\] state normal v LinearElasticPlaneStrain${nd}Law"
    gid_groups_conds::addF $structPartsGroup value {n SECTION_TYPE pn {Section type} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF -resolve_parametric 1 $structPartsGroup value {n THICKNESS pn Thickness state {[PartParamState]} unit_magnitude L help {} v 1.0 units m}
    gid_groups_conds::addF -resolve_parametric 1 $structPartsGroup value {n DENSITY pn Density state {[PartParamState]} unit_magnitude Density help {} v 1500.0 units kg/m^3}
    gid_groups_conds::addF -resolve_parametric 1 $structPartsGroup value {n VISCOSITY pn {Kinematic viscosity} state {[PartParamState]} unit_magnitude L^2/T help {Fluidized viscosity.} v 1e-6 units m^2/s}
    gid_groups_conds::addF $structPartsGroup value {n YIELD_STRESS pn {Yield Stress} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF -resolve_parametric 1 $structPartsGroup value {n YOUNG_MODULUS pn {Young Modulus} state {[PartParamState]} unit_magnitude P help {} v 2.3e6 units Pa}
    gid_groups_conds::addF $structPartsGroup value {n POISSON_RATIO pn {Poisson Ratio} state {[PartParamState]} unit_magnitude {} help {} v 0.45}
    gid_groups_conds::addF $structPartsGroup value {n KINEMATIC_HARDENING_MODULUS pn {Kinematic Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF $structPartsGroup value {n REFERENCE_HARDENING_MODULUS pn {Reference Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF $structPartsGroup value {n INFINITY_HARDENING_MODULUS pn {Infinity Hardening Modulus} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF $structPartsGroup value {n HARDENING_EXPONENT pn {Hardening Exponent} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF $structPartsGroup value {n DAMAGE_THRESHOLD pn {Damage Threshold} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF $structPartsGroup value {n STRENGTH_RATIO pn {Strength Ratio} state {[PartParamState]} unit_magnitude {} help {} v 0}
    gid_groups_conds::addF $structPartsGroup value {n FRACTURE_ENERGY pn {Fracture Energy} state {[PartParamState]} unit_magnitude {} help {} v 0}
    
    # Structural Displacement
    set structDisplacement {container[@n='FSI']/container[@n='Structural']/container[@n='Boundary Conditions']/condition[@n='DISPLACEMENT']}
    gid_groups_conds::addF $structDisplacement group {n FixedDisplacement ov line}
    set structDisplacementGroup "$structDisplacement/group\[@n='FixedDisplacement'\]"
    gid_groups_conds::addF $structDisplacementGroup value {n FixX pn {X Imposed} values 1,0 help {} state {} v 1}
    gid_groups_conds::addF $structDisplacementGroup value {n FixY pn {Y Imposed} values 1,0 help {} state {} v 1}
    gid_groups_conds::addF $structDisplacementGroup value {n FixZ pn {Z Imposed} values 1,0 help {} state {[CheckDimension 3D]} v 1}
    gid_groups_conds::addF $structDisplacementGroup value {n valueX wn {DISPLACEMENT _X} pn {Value X} help {} state {} v 0.0}
    gid_groups_conds::addF $structDisplacementGroup value {n valueY wn {DISPLACEMENT _Y} pn {Value Y} help {} state {} v 0.0}
    gid_groups_conds::addF $structDisplacementGroup value {n valueZ wn {DISPLACEMENT _Z} pn {Value Z} help {} state {[CheckDimension 3D]} v 0.0}
    
    # Structural Interface
    gid_groups_conds::addF "container\[@n='FSI'\]/container\[@n='Structural'\]/container\[@n='Loads'\]/condition\[@n='StructureInterface$nd'\]" group {n StructureInterface}

    spdAux::RequestRefresh
}
