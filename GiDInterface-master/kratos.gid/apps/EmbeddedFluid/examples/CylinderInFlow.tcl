
proc ::EmbeddedFluid::examples::CylinderInFlow {args} {
    InitVariables
    DrawCylinderInFlowGeometry3D
    AssignGroupsCylinderInFlow3D
    AssignCylinderInFlowMeshSizes3D
    TreeAssignationCylinderInFlow3D

    AddMeshOptimizationPoints

    GiD_Process 'Zoom Frame
    GiD_Process 'Redraw
    GidUtils::UpdateWindow GROUPS
    GidUtils::UpdateWindow LAYER
}

proc EmbeddedFluid::examples::InitVariables { } {
    variable CylinderInFlow_Data
    set CylinderInFlow_Data(circle_center_x) 0.75
    set CylinderInFlow_Data(circle_center_y) 0.5
    set CylinderInFlow_Data(circle_center_z) 0.0
    set CylinderInFlow_Data(circle_radius)   0.1
}


# Draw Geometry
proc EmbeddedFluid::examples::DrawCylinderInFlowGeometry3D {args} {
    DrawCylinderInFlowGeometry2D
    GiD_Process Mescape Utilities Copy Surfaces Duplicate DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,1.0 1 escape escape escape
    GiD_Process Mescape Utilities Copy Surfaces Duplicate DoExtrude Surfaces MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,1.0 2 escape escape escape
    GiD_Layers edit opaque Fluid 0

}
proc EmbeddedFluid::examples::DrawCylinderInFlowGeometry2D {args} {
    Kratos::ResetModel
    GiD_Layers create Fluid
    GiD_Layers edit to_use Fluid

    # Geometry creation
    ## Points ##
    set coordinates [list 0 1 0 3.5 1 0 3.5 0 0 0 0 0]
    set fluidPoints [list ]
    foreach {x y z} $coordinates {
        lappend fluidPoints [GiD_Geometry create point append Fluid $x $y $z]
    }

    ## Lines ##
    set fluidLines [list ]
    set initial [lindex $fluidPoints 0]
    foreach point [lrange $fluidPoints 1 end] {
        lappend fluidLines [GiD_Geometry create line append stline Fluid $initial $point]
        set initial $point
    }
    lappend fluidLines [GiD_Geometry create line append stline Fluid $initial [lindex $fluidPoints 0]]

    ## Surface ##
    GiD_Process Mescape Geometry Create NurbsSurface {*}$fluidLines escape escape

    # Body #
    GiD_Layers create Body
    GiD_Layers edit to_use Body
    variable CylinderInFlow_Data
    set circle_center_x $CylinderInFlow_Data(circle_center_x)
    set circle_center_y $CylinderInFlow_Data(circle_center_y)
    set circle_center_z $CylinderInFlow_Data(circle_center_z)
    set circle_radius   $CylinderInFlow_Data(circle_radius)
    GiD_Process Mescape Geometry Create Object CirclePNR $circle_center_x $circle_center_y $circle_center_z 0.0 0.0 1.0 $circle_radius escape
    GiD_Process escape MEscape Geometry Edit DivideLine Multiple NumDivisions 2 5 escape escape 

    # GiD_Geometry delete surface 2

    # Create the hole
    GiD_Layers edit to_use Fluid
    # GiD_Process MEscape Geometry Edit HoleNurb 1 5 escape escape

}

# Group assign
proc EmbeddedFluid::examples::AssignGroupsCylinderInFlow3D {args} {
    # Create the groups
    GiD_Groups create Fluid
    GiD_Groups edit color Fluid "#26d1a8ff"
    GiD_EntitiesGroups assign Fluid volumes 1

    GiD_Groups create Inlet
    GiD_Groups edit color Inlet "#e0210fff"
    GiD_EntitiesGroups assign Inlet surfaces 6

    GiD_Groups create Outlet
    GiD_Groups edit color Outlet "#42eb71ff"
    GiD_EntitiesGroups assign Outlet surfaces 4

    GiD_Groups create No_Slip_Walls
    GiD_Groups edit color No_Slip_Walls "#3b3b3bff"
    GiD_EntitiesGroups assign No_Slip_Walls surfaces {1 3 5 7}

    GiD_Groups create No_Slip_Cylinder
    GiD_Groups edit color No_Slip_Cylinder "#3b3b3bff"
    GiD_EntitiesGroups assign No_Slip_Cylinder surfaces {8 9}
}

# Mesh sizes
proc EmbeddedFluid::examples::AssignCylinderInFlowMeshSizes3D {args} {
    set cylinder_mesh_size 0.005
    set walls_mesh_size 0.05
    set fluid_mesh_size 0.05
    GiD_Process Mescape Utilities Variables SizeTransitionsFactor 0.4 escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $cylinder_mesh_size {*}[GiD_EntitiesGroups get No_Slip_Cylinder surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $walls_mesh_size {*}[GiD_EntitiesGroups get Inlet surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $walls_mesh_size {*}[GiD_EntitiesGroups get Outlet surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $walls_mesh_size {*}[GiD_EntitiesGroups get No_Slip_Walls surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Volumes $fluid_mesh_size [GiD_EntitiesGroups get Fluid volumes] escape escape
    # Kratos::BeforeMeshGeneration $fluid_mesh_size
}
proc EmbeddedFluid::examples::AssignCylinderInFlowMeshSizes2D {args} {
    set cylinder_mesh_size 0.005
    set fluid_mesh_size 0.05
    GiD_Process Mescape Utilities Variables SizeTransitionsFactor 0.4 escape escape
    GiD_Process Mescape Meshing AssignSizes Lines $cylinder_mesh_size {*}[GiD_EntitiesGroups get No_Slip_Cylinder lines] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $fluid_mesh_size [GiD_EntitiesGroups get Fluid surfaces] escape escape
    # Kratos::BeforeMeshGeneration $fluid_mesh_size
}

# Tree assign
proc EmbeddedFluid::examples::TreeAssignationCylinderInFlow3D {args} {
    TreeAssignationCylinderInFlow2D
}
proc EmbeddedFluid::examples::TreeAssignationCylinderInFlow2D {args} {
    set nd $::Model::SpatialDimension
    set root [customlib::GetBaseRoot]

    set condtype line
    if {$nd eq "3D"} { set condtype surface }

    # Monolithic solution strategy set
    spdAux::SetValueOnTreeItem v "Monolithic" FLSolStrat

    # Fluid Parts
    set fluidParts [spdAux::getRoute "FLParts"]
    set fluidNode [spdAux::AddConditionGroupOnXPath $fluidParts Fluid]
    set props [list Element Monolithic$nd ConstitutiveLaw Newtonian DENSITY 1.0 DYNAMIC_VISCOSITY 0.003 YIELD_STRESS 0 POWER_LAW_K 1 POWER_LAW_N 1]
    foreach {prop val} $props {
        set propnode [$fluidNode selectNodes "./value\[@n = '$prop'\]"]
        if {$propnode ne "" } {
            $propnode setAttribute v $val
        } else {
            W "Warning - Couldn't find property Fluid $prop"
        }
    }

    set fluidConditions [spdAux::getRoute "FLBC"]

    # Fluid Inlet
    set fluidInlet "$fluidConditions/condition\[@n='AutomaticInlet$nd'\]"
    set inlets [list inlet1 0 1 "6*y*(1-y)*sin(pi*t*0.5)" inlet2 1 End "6*y*(1-y)"]
    ErasePreviousIntervals
    foreach {inlet_name ini end function} $inlets {
        spdAux::CreateInterval $inlet_name $ini $end
        GiD_Groups create "Inlet//$inlet_name"
        spdAux::AddIntervalGroup Inlet "Inlet//$inlet_name"
        set inletNode [spdAux::AddConditionGroupOnXPath $fluidInlet "Inlet//$inlet_name"]
        $inletNode setAttribute ov $condtype
        set props [list ByFunction Yes function_modulus $function direction automatic_inwards_normal Interval $inlet_name]
        foreach {prop val} $props {
             set propnode [$inletNode selectNodes "./value\[@n = '$prop'\]"]
             if {$propnode ne "" } {
                  $propnode setAttribute v $val
             } else {
                W "Warning - Couldn't find property Inlet $prop"
            }
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
    [spdAux::AddConditionGroupOnXPath "$fluidConditions/condition\[@n='NoSlip$nd'\]" No_Slip_Walls] setAttribute ov $condtype
    # [spdAux::AddConditionGroupOnXPath "$fluidConditions/condition\[@n='NoSlip$nd'\]" No_Slip_Cylinder] setAttribute ov $condtype

    # Time parameters
    set time_parameters [list EndTime 45 DeltaTime 0.1]
    set time_params_path [spdAux::getRoute "FLTimeParameters"]
    foreach {n v} $time_parameters {
        [$root selectNodes "$time_params_path/value\[@n = '$n'\]"] setAttribute v $v
    }
    # Output
    set time_parameters [list OutputControlType step OutputDeltaStep 1]
    set time_params_path [spdAux::getRoute "Results"]
    foreach {n v} $time_parameters {
        [$root selectNodes "$time_params_path/value\[@n = '$n'\]"] setAttribute v $v
    }
    # Parallelism
    set time_parameters [list ParallelSolutionType OpenMP OpenMPNumberOfThreads 4]
    set time_params_path [spdAux::getRoute "Parallelization"]
    foreach {n v} $time_parameters {
        [$root selectNodes "$time_params_path/value\[@n = '$n'\]"] setAttribute v $v
    }

    spdAux::RequestRefresh
}

proc EmbeddedFluid::examples::ErasePreviousIntervals { } {
    set root [customlib::GetBaseRoot]
    set interval_base [spdAux::getRoute "Intervals"]
    foreach int [$root selectNodes "$interval_base/blockdata\[@n='Interval'\]"] {
        if {[$int @name] ni [list Initial Total Custom1]} {$int delete}
    }
}

proc EmbeddedFluid::examples::AddMeshOptimizationPoints { } {
    set optimized_group "Optimized mesh"
    
    GiD_Layers create Mesh_Optimization
    GiD_Layers edit to_use Mesh_Optimization
    
    variable CylinderInFlow_Data
    set radius   $CylinderInFlow_Data(circle_radius)
    set center_x [expr $CylinderInFlow_Data(circle_center_x) + $radius +0.025]
    set center_y $CylinderInFlow_Data(circle_center_y)
    set center_z [expr $CylinderInFlow_Data(circle_center_z) + 0.025]
    set origin_point [GiD_Geometry create point append Mesh_Optimization $center_x $center_y $center_z]

    GiD_Groups create $optimized_group
    GiD_EntitiesGroups assign $optimized_group points $origin_point
    GiD_Process Mescape Utilities Copy Points Duplicate MaintainLayers MCopy 19 Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,0.05 $origin_point escape Mescape escape 
    
    set original_points [GiD_EntitiesGroups get $optimized_group points]
    GiD_Process Mescape Utilities Copy Points Duplicate MaintainLayers MCopy 15 Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.166,0.0,0.0 {*}$original_points escape Mescape escape 
    
    set original_points [GiD_EntitiesGroups get $optimized_group points]
    GiD_Process Mescape Meshing MeshCriteria ForcePointsTo VolumeMesh 1 escape {*}$original_points escape 
    GiD_Process Mescape Meshing AssignSizes Points 0.005 {*}$original_points escape escape

    GiD_Process Mescape Utilities Copy Points Duplicate MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.1,0.0  {*}$original_points escape 
    GiD_Process Mescape Utilities Copy Points Duplicate MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,-0.1,0.0 {*}$original_points escape 

}