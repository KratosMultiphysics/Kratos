
proc ::FSI::examples::HorizontalFlexibleBar {args} {
    DrawHorizontalFlexibleBarGeometry
    TreeAssignationHorizontalFlexibleBar
}

proc FSI::examples::DrawHorizontalFlexibleBarGeometry {args} {
    Kratos::ResetModel
    GiD_Process Mescape 'Layers ChangeName Layer0 Fluid escape
    
    # Geometry creation
    set coordinates [list 0 0 0 0 2 0 5 2 0 5 0 0]
    set fluidPoints [list ]
    foreach {x y z} $coordinates {
        lappend fluidPoints [GiD_Geometry create point append Fluid $x $y $z]
    }
    set fluidLines [list ]
    set initial [lindex $fluidPoints 0]
    foreach point [lrange $fluidPoints 1 end] {
        lappend fluidLines [GiD_Geometry create line append stline Fluid $initial $point]
        set initial $point
    }
    lappend fluidLines [GiD_Geometry create line append stline Fluid $initial [lindex $fluidPoints 0]]
    
    # Hole in the middle
    set coordinates [list 1 0.99 0 1 1.01 0 1.5 1.01 0 1.5 0.99 0]
    set fluidnterfacePoints [list ]
    foreach {x y z} $coordinates {
        lappend fluidnterfacePoints [GiD_Geometry create point append Fluid $x $y $z]
    }
    set fluidinteractionLines [list ]
    
    set initial [lindex $fluidnterfacePoints 0]
    foreach point [lrange $fluidnterfacePoints 1 end] {
        lappend fluidinteractionLines [GiD_Geometry create line append stline Fluid $initial $point]
        set initial $point
    }
    lappend fluidinteractionLines [GiD_Geometry create line append stline Fluid $initial [lindex $fluidnterfacePoints 0]]
    #set fluidSurface [GiD_Geometry create surface append plsurface Fluid [llength $fluidLines] {*}$fluidLines]

    GiD_Process Mescape Geometry Create NurbsSurface {*}$fluidLines escape escape
    GiD_Process MEscape Geometry Edit HoleNurb 1 {*}$fluidinteractionLines escape escape


    GiD_Process 'Layers New Structure escape 
    GiD_Process 'Layers Off Fluid escape
    GiD_Process 'Layers ToUse Structure escape
    
    
    set coordinates [list 1 0.99 0 1 1.01 0 1.5 1.01 0 1.5 0.99 0]
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
        GiD_Groups create FluidFixedDisplacement
        GiD_EntitiesGroups assign Fluid volumes 1
        GiD_EntitiesGroups assign Structure volumes 2
        GiD_EntitiesGroups assign Inlet surfaces 3
        GiD_EntitiesGroups assign Outlet surfaces 5
        GiD_EntitiesGroups assign NoSlip surfaces {1 11}
        GiD_EntitiesGroups assign Slip surfaces {4 6}
        GiD_EntitiesGroups assign FluidFixedDisplacement surfaces {1 11}
        GiD_EntitiesGroups assign FluidInterface surfaces {7 8 9 10}
        GiD_EntitiesGroups assign FixedDisplacement surfaces 12
        GiD_EntitiesGroups assign StructureInterface surfaces {2 12 13 14 15 16}
        
    } {
        GiD_EntitiesGroups assign Fluid surfaces 1
        GiD_EntitiesGroups assign Structure surfaces 2
        GiD_EntitiesGroups assign Inlet lines 1
        GiD_EntitiesGroups assign Outlet lines 3
        GiD_EntitiesGroups assign Slip lines {2 4}
        GiD_EntitiesGroups assign FluidInterface lines $fluidinteractionLines
        GiD_EntitiesGroups assign FixedDisplacement lines 9
        GiD_EntitiesGroups assign StructureInterface lines $strucLines
    }
    GidUtils::UpdateWindow GROUPS
}

proc FSI::examples::TreeAssignationHorizontalFlexibleBar {args} {
    set nd $::Model::SpatialDimension
    # Fluid Parts
    
    set condtype line
    if {$nd eq "3D"} { set condtype surface }
    set fluidParts {container[@n='FSI']/container[@n='Fluid']/condition[@n='Parts']}
    set fluidNode [spdAux::AddConditionGroupOnXPath $fluidParts Fluid]
    set props [list Element FractionalStep$nd ConstitutiveLaw Newtonian DENSITY 956.0 VISCOSITY 0.145 YIELD_STRESS 0 POWER_LAW_K 1 POWER_LAW_N 1]
    foreach {prop val} $props {
        set propnode [$fluidNode selectNodes "./value\[@n = '$prop'\]"]
        if {$propnode ne "" } {
            $propnode setAttribute v $val
        } else {
            W "Warning - Couldn't find property Fluid $prop"
        }
    }

    set fluidConditions {container[@n='FSI']/container[@n='Fluid']/container[@n='BoundaryConditions']}
    # Fluid Interface
    set fluidInlet "$fluidConditions/condition\[@n='AutomaticInlet$nd'\]"
    
    # Fluid Inlet
    set inletNode [spdAux::AddConditionGroupOnXPath $fluidInlet Inlet]
    $inletNode setAttribute ov $condtype
    set props [list ByFunction No modulus 0.6067 direction automatic_inwards_normal Interval Total]
    foreach {prop val} $props {
         set propnode [$inletNode selectNodes "./value\[@n = '$prop'\]"]
         if {$propnode ne "" } {
              $propnode setAttribute v $val
         } else {
            W "Warning - Couldn't find property Inlet $prop"
        }
    }  

    # Fluid Outlet
    set fluidOutlet "$fluidConditions/condition\[@n='Outlet$nd'\]"
    set outletNode [spdAux::AddConditionGroupOnXPath $fluidOutlet Outlet]
    $outletNode setAttribute ov $condtype
    set props [list value 0.0]
    foreach {prop val} $props {
         set propnode [$outletNode selectNodes "./value\[@n = '$prop'\]"]
         if {$propnode ne "" } {
              $propnode setAttribute v $val
         } else {
            W "Warning - Couldn't find property Outlet $prop"
        }
    }
    # Fluid Conditions
    if {$nd eq "3D"} {
        [spdAux::AddConditionGroupOnXPath "$fluidConditions/condition\[@n='NoSlip$nd'\]" NoSlip] setAttribute ov $condtype
    }
    [spdAux::AddConditionGroupOnXPath "$fluidConditions/condition\[@n='Slip$nd'\]" Slip] setAttribute ov $condtype
    [spdAux::AddConditionGroupOnXPath "$fluidConditions/condition\[@n='FluidNoSlipInterface$nd'\]" FluidInterface] setAttribute ov $condtype  
    
    # Displacement 3D
    if {$nd eq "3D"} {
        set fluidDisplacement "$fluidConditions/condition\[@n='DISPLACEMENT'\]"
        gid_groups_conds::addF $fluidDisplacement group {n FluidFixedDisplacement ov surface}
        set fluidDisplacementGroup "$fluidDisplacement/group\[@n='FluidFixedDisplacement'\]"
        gid_groups_conds::addF $fluidDisplacementGroup value {n FixX pn {X Imposed} values 1,0 help {} state {} v 0}
        gid_groups_conds::addF $fluidDisplacementGroup value {n FixY pn {Y Imposed} values 1,0 help {} state {} v 0}
        gid_groups_conds::addF $fluidDisplacementGroup value {n FixZ pn {Z Imposed} values 1,0 help {} state {[CheckDimension 3D]} v 1}
        gid_groups_conds::addF $fluidDisplacementGroup value {n valueX wn {DISPLACEMENT _X} pn {Value X} help {} state {} v 0.0}
        gid_groups_conds::addF $fluidDisplacementGroup value {n valueY wn {DISPLACEMENT _Y} pn {Value Y} help {} state {} v 0.0}
        gid_groups_conds::addF $fluidDisplacementGroup value {n valueZ wn {DISPLACEMENT _Z} pn {Value Z} help {} state {[CheckDimension 3D]} v 0.0}
    }
    
    

    # Structural
    gid_groups_conds::setAttributesF {container[@n='FSI']/container[@n='Structural']/container[@n='StageInfo']/value[@n='SolutionType']} {v Dynamic}
    
    # Structural Parts
    set structParts {container[@n='FSI']/container[@n='Structural']/condition[@n='Parts']}
    set structPartsNode [spdAux::AddConditionGroupOnXPath $structParts Structure]
    $structPartsNode setAttribute ov [expr {$nd == "3D" ? "volume" : "surface"}]
    set constLawNameStruc [expr {$nd == "3D" ? "LinearElastic3DLaw" : "LinearElasticPlaneStrain2DLaw"}]
    set props [list Element SmallDisplacementElement$nd ConstitutiveLaw $constLawNameStruc SECTION_TYPE 0 THICKNESS 1.0 DENSITY 1500.0 VISCOSITY 1e-6]
    lappend props YIELD_STRESS 0 YOUNG_MODULUS 2.3e6 POISSON_RATIO 0.45 
    foreach {prop val} $props {
         set propnode [$structPartsNode selectNodes "./value\[@n = '$prop'\]"]
         if {$propnode ne "" } {
              $propnode setAttribute v $val
         } else {
            W "Warning - Couldn't find property Structure $prop"
         }
    }
    
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
    gid_groups_conds::addF "container\[@n='FSI'\]/container\[@n='Structural'\]/container\[@n='Loads'\]/condition\[@n='Interface$nd'\]" group {n StructureInterface}

    spdAux::RequestRefresh
}