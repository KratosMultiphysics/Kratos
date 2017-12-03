
proc ::Fluid::examples::CylinderInFlow {args} {
    DrawCylinderInFlowGeometry$::Model::SpatialDimension
    AssignGroupsCylinderInFlow$::Model::SpatialDimension
    AssignCylinderInFlowMeshSizes$::Model::SpatialDimension
    TreeAssignationCylinderInFlow$::Model::SpatialDimension

    GiD_Process 'Redraw
    GidUtils::UpdateWindow GROUPS
    GidUtils::UpdateWindow LAYER
}


# Draw Geometry
proc Fluid::examples::DrawCylinderInFlowGeometry3D {args} {
    DrawCylinderInFlowGeometry2D
    GiD_Process Mescape Utilities Copy Surfaces Duplicate DoExtrude Volumes MaintainLayers Translation FNoJoin 0.0,0.0,0.0 FNoJoin 0.0,0.0,1.0 1 escape escape escape
    GiD_Layers edit opaque Fluid 0

}
proc Fluid::examples::DrawCylinderInFlowGeometry2D {args} {
    Kratos::ResetModel
    GiD_Layers create Fluid
    GiD_Layers edit to_use Fluid

    # Geometry creation
    ## Points ##
    set coordinates [list 0 1 0 5 1 0 5 0 0 0 0 0]
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
    set circle_center_x 1.25
    set circle_center_y 0.5
    set circle_center_z 0.0
    set center_radius 0.1
    GiD_Process Mescape Geometry Create Object CirclePNR $circle_center_x $circle_center_y $circle_center_z 0.0 0.0 1.0 $center_radius escape
    GiD_Geometry delete surface 2

    # Create the hole
    GiD_Layers edit to_use Fluid
    GiD_Process MEscape Geometry Edit HoleNurb 1 5 escape escape

}


# Group assign
proc Fluid::examples::AssignGroupsCylinderInFlow2D {args} {
    # Create the groups
    GiD_Groups create Fluid
    GiD_Groups edit color Fluid "#26d1a8ff"
    GiD_EntitiesGroups assign Fluid surfaces 1

    GiD_Groups create Inlet
    GiD_Groups edit color Inlet "#e0210fff"
    GiD_EntitiesGroups assign Inlet lines 4

    GiD_Groups create Outlet
    GiD_Groups edit color Outlet "#42eb71ff"
    GiD_EntitiesGroups assign Outlet lines 2

    GiD_Groups create No_Slip_Walls
    GiD_Groups edit color No_Slip_Walls "#3b3b3bff"
    GiD_EntitiesGroups assign No_Slip_Walls lines {1 3}

    GiD_Groups create No_Slip_Cylinder
    GiD_Groups edit color No_Slip_Cylinder "#3b3b3bff"
    GiD_EntitiesGroups assign No_Slip_Cylinder lines 5
}
proc Fluid::examples::AssignGroupsCylinderInFlow3D {args} {
    # Create the groups
    GiD_Groups create Fluid
    GiD_Groups edit color Fluid "#26d1a8ff"
    GiD_EntitiesGroups assign Fluid volumes 1

    GiD_Groups create Inlet
    GiD_Groups edit color Inlet "#e0210fff"
    GiD_EntitiesGroups assign Inlet surfaces 5

    GiD_Groups create Outlet
    GiD_Groups edit color Outlet "#42eb71ff"
    GiD_EntitiesGroups assign Outlet surfaces 3

    GiD_Groups create No_Slip_Walls
    GiD_Groups edit color No_Slip_Walls "#3b3b3bff"
    GiD_EntitiesGroups assign No_Slip_Walls surfaces {1 2 4 7}

    GiD_Groups create No_Slip_Cylinder
    GiD_Groups edit color No_Slip_Cylinder "#3b3b3bff"
    GiD_EntitiesGroups assign No_Slip_Cylinder surfaces 6
}


# Mesh sizes
proc Fluid::examples::AssignCylinderInFlowMeshSizes3D {args} {
    set cylinder_mesh_size 0.005
    set walls_mesh_size 0.05
    set fluid_mesh_size 0.05
    GiD_Process Mescape Utilities Variables SizeTransitionsFactor 0.4 escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $cylinder_mesh_size {*}[GiD_EntitiesGroups get No_Slip_Cylinder surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $walls_mesh_size {*}[GiD_EntitiesGroups get Inlet surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $walls_mesh_size {*}[GiD_EntitiesGroups get Outlet surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $walls_mesh_size {*}[GiD_EntitiesGroups get No_Slip_Walls surfaces] escape escape
    GiD_Process Mescape Meshing AssignSizes Volumes $fluid_mesh_size [GiD_EntitiesGroups get Fluid volumes] escape escape
    Kratos::BeforeMeshGeneration $fluid_mesh_size
}
proc Fluid::examples::AssignCylinderInFlowMeshSizes2D {args} {
    set cylinder_mesh_size 0.005
    set fluid_mesh_size 0.05
    GiD_Process Mescape Utilities Variables SizeTransitionsFactor 0.4 escape escape
    GiD_Process Mescape Meshing AssignSizes Lines $cylinder_mesh_size {*}[GiD_EntitiesGroups get No_Slip_Cylinder lines] escape escape
    GiD_Process Mescape Meshing AssignSizes Surfaces $fluid_mesh_size [GiD_EntitiesGroups get Fluid surfaces] escape escape
    Kratos::BeforeMeshGeneration $fluid_mesh_size
}


# Tree assign
proc Fluid::examples::TreeAssignationCylinderInFlow3D {args} {
    TreeAssignationCylinderInFlow2D
}
proc Fluid::examples::TreeAssignationCylinderInFlow2D {args} {
    set nd $::Model::SpatialDimension
    set root [customlib::GetBaseRoot]

    set condtype line
    if {$nd eq "3D"} { set condtype surface }

    # Monolithic solution strategy set
    spdAux::SetValueOnTreeItem v "Monolithic" FLSolStrat

    # Fluid Parts
    set fluidParts [spdAux::getRoute "FLParts"]
    set fluidNode [spdAux::AddConditionGroupOnXPath $fluidParts Fluid]
    set props [list Element Monolithic$nd ConstitutiveLaw Newtonian DENSITY 1.0 VISCOSITY 0.003 YIELD_STRESS 0 POWER_LAW_K 1 POWER_LAW_N 1]
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
    [spdAux::AddConditionGroupOnXPath "$fluidConditions/condition\[@n='NoSlip$nd'\]" No_Slip_Cylinder] setAttribute ov $condtype

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

proc Fluid::examples::ErasePreviousIntervals { } {
    set root [customlib::GetBaseRoot]
    set interval_base [spdAux::getRoute "Intervals"]
    foreach int [$root selectNodes "$interval_base/blockdata\[@n='Interval'\]"] {
        if {[$int @name] ni [list Initial Total Custom1]} {$int delete}
    }
}
