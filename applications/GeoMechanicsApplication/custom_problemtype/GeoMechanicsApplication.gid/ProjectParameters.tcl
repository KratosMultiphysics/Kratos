proc WriteProjectParameters { basename dir problemtypedir TableDict} {

    ## Source auxiliar procedures
    source [file join $problemtypedir ProjectParametersAuxProcs.tcl]

    ## Start ProjectParameters.json file
    set filename [file join $dir ProjectParameters.json]
    set FileVar [open $filename w]

    set IsK0 [GiD_AccessValue get gendata Solution_Type]

    puts $FileVar "\{"

    ## problem_data
    puts $FileVar "    \"problem_data\": \{"
    puts $FileVar "        \"problem_name\":         \"$basename\","
    puts $FileVar "        \"start_time\":           [GiD_AccessValue get gendata Start_Time],"
    puts $FileVar "        \"end_time\":             [GiD_AccessValue get gendata End_Time],"
    puts $FileVar "        \"echo_level\":           [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "        \"parallel_type\":        \"[GiD_AccessValue get gendata Parallel_Configuration]\","
    puts $FileVar "        \"number_of_threads\":    [GiD_AccessValue get gendata Number_of_threads]"
    puts $FileVar "    \},"

    ## solver_settings
    puts $FileVar "    \"solver_settings\": \{"
    puts $FileVar "        \"solver_type\":                        \"U_Pw\","
    puts $FileVar "        \"model_part_name\":                    \"PorousDomain\","
    puts $FileVar "        \"domain_size\":                        [GiD_AccessValue get gendata Domain_Size],"
    puts $FileVar "        \"start_time\":                         [GiD_AccessValue get gendata Start_Time],"
    puts $FileVar "        \"model_import_settings\":              \{"
    puts $FileVar "            \"input_type\":       \"mdpa\","
    puts $FileVar "            \"input_filename\":   \"$basename\""
    puts $FileVar "        \},"
    puts $FileVar "        \"material_import_settings\":              \{"
    puts $FileVar "            \"materials_filename\":       \"MaterialParameters.json\""
    puts $FileVar "        \},"
    puts $FileVar "        \"time_stepping\":              \{"
    puts $FileVar "            \"time_step\":                [GiD_AccessValue get gendata Delta_Time],"
	puts $FileVar "            \"max_delta_time_factor\":    [GiD_AccessValue get gendata Max_Delta_Time_Factor]"
    puts $FileVar "        \},"
    puts $FileVar "        \"buffer_size\":                        2,"
    puts $FileVar "        \"echo_level\":                         [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "        \"clear_storage\":                      false,"
    puts $FileVar "        \"compute_reactions\":                  [GiD_AccessValue get gendata Write_Reactions],"
    puts $FileVar "        \"move_mesh_flag\":                     [GiD_AccessValue get gendata Move_Mesh],"

    set IsGapClosure [GiD_AccessValue get gendata Gap_Closure_Interface_Conditions]
    if {$IsGapClosure eq true} {
        puts $FileVar "        \"reform_dofs_at_each_step\":           true,"
    } else {
        puts $FileVar "        \"reform_dofs_at_each_step\":           [GiD_AccessValue get gendata Reform_Dofs_At_Each_Step],"
    }

    puts $FileVar "        \"nodal_smoothing\":                    [GiD_AccessValue get gendata Nodal_Smoothing],"
    puts $FileVar "        \"block_builder\":                      [GiD_AccessValue get gendata Block_Builder],"
    puts $FileVar "        \"solution_type\":                      \"[GiD_AccessValue get gendata Solution_Type]\","
    puts $FileVar "        \"scheme_type\":                        \"[GiD_AccessValue get gendata Scheme_Type]\","
    puts $FileVar "        \"reset_displacements\":                [GiD_AccessValue get gendata Reset_Displacements],"
    puts $FileVar "        \"newmark_beta\":                       [GiD_AccessValue get gendata Newmark_Beta],"
    puts $FileVar "        \"newmark_gamma\":                      [GiD_AccessValue get gendata Newmark_Gamma],"
    puts $FileVar "        \"newmark_theta\":                      [GiD_AccessValue get gendata Newmark_Theta],"
    puts $FileVar "        \"rayleigh_m\":                         [GiD_AccessValue get gendata Rayleigh_Mass],"
    puts $FileVar "        \"rayleigh_k\":                         [GiD_AccessValue get gendata Rayleigh_Stiffness],"
    puts $FileVar "        \"strategy_type\":                      \"[GiD_AccessValue get gendata Strategy_Type]\","
    puts $FileVar "        \"convergence_criterion\":              \"[GiD_AccessValue get gendata Convergence_Criterion]\","
    puts $FileVar "        \"displacement_relative_tolerance\":    [GiD_AccessValue get gendata Displacement_Relative_Tolerance],"
    puts $FileVar "        \"displacement_absolute_tolerance\":    [GiD_AccessValue get gendata Displacement_Absolute_Tolerance],"
    puts $FileVar "        \"residual_relative_tolerance\":        [GiD_AccessValue get gendata Residual_Relative_Tolerance],"
    puts $FileVar "        \"residual_absolute_tolerance\":        [GiD_AccessValue get gendata Residual_Absolute_Tolerance],"
    puts $FileVar "        \"min_iterations\":                      [GiD_AccessValue get gendata Min_Iterations],"   
    puts $FileVar "        \"max_iterations\":                      [GiD_AccessValue get gendata Max_Iterations],"
    puts $FileVar "        \"number_cycles\":                      [GiD_AccessValue get gendata Number_Of_Cycles],"
    puts $FileVar "        \"reduction_factor\":                   [GiD_AccessValue get gendata Reduction_Factor],"   
    puts $FileVar "        \"increase_factor\":                    [GiD_AccessValue get gendata Increase_Factor],"   
    puts $FileVar "        \"realised_factor\":                    [GiD_AccessValue get gendata Realised_Factor],"   
    puts $FileVar "        \"desired_iterations\":                 [GiD_AccessValue get gendata Desired_Iterations],"
    puts $FileVar "        \"max_radius_factor\":                  [GiD_AccessValue get gendata Max_Radius_Factor],"
    puts $FileVar "        \"min_radius_factor\":                  [GiD_AccessValue get gendata Min_Radius_Factor],"
    puts $FileVar "        \"calculate_reactions\":                [GiD_AccessValue get gendata Calculate_Reactions],"
    puts $FileVar "        \"max_line_search_iterations\":         [GiD_AccessValue get gendata Max_Line_Search_Iterations],"
    puts $FileVar "        \"first_alpha_value\":                  [GiD_AccessValue get gendata First_Alpha_Value],"
    puts $FileVar "        \"second_alpha_value\":                 [GiD_AccessValue get gendata Second_Alpha_Value],"
    puts $FileVar "        \"min_alpha\":                          [GiD_AccessValue get gendata Min_Alpha],"
    puts $FileVar "        \"max_alpha\":                          [GiD_AccessValue get gendata Max_Alpha],"
    puts $FileVar "        \"line_search_tolerance\":              [GiD_AccessValue get gendata Line_Search_Tolerance],"
    puts $FileVar "        \"rotation_dofs\":                      true,"
    ## linear_solver_settings
    puts $FileVar "        \"linear_solver_settings\":             \{"
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        if {[GiD_AccessValue get gendata Solver_Type] eq "AmgclMPISolver"} {
            puts $FileVar "            \"solver_type\":   \"AmgclMPISolver\","
            puts $FileVar "            \"krylov_type\":   \"fgmres\","
            puts $FileVar "            \"max_iteration\": 100,"
            puts $FileVar "            \"verbosity\":     [GiD_AccessValue get gendata Verbosity],"
            puts $FileVar "            \"tolerance\":     1.0e-6,"
            puts $FileVar "            \"scaling\":       [GiD_AccessValue get gendata Scaling]"
        } elseif {[GiD_AccessValue get gendata Solver_Type] eq "AztecSolver"} {
            puts $FileVar "            \"solver_type\":         \"AztecSolver\","
            puts $FileVar "            \"tolerance\":           1.0e-6,"
            puts $FileVar "            \"max_iteration\":       200,"
            puts $FileVar "            \"scaling\":             [GiD_AccessValue get gendata Scaling],"
            puts $FileVar "            \"preconditioner_type\": \"None\""
        } elseif {([GiD_AccessValue get gendata Solver_Type] eq "Klu") || ([GiD_AccessValue get gendata Solver_Type] eq "MultiLevelSolver")} {
            puts $FileVar "            \"solver_type\": \"[GiD_AccessValue get gendata Solver_Type]\","
            puts $FileVar "            \"scaling\":     [GiD_AccessValue get gendata Scaling]"
        } else {
            puts $FileVar "            \"solver_type\": \"Klu\","
            puts $FileVar "            \"scaling\":     false"
        }
    } else {
        if {([GiD_AccessValue get gendata Solver_Type] eq "amgcl") || ([GiD_AccessValue get gendata Solver_Type] eq "amgcl_ns")} {
            puts $FileVar "            \"solver_type\":     \"[GiD_AccessValue get gendata Solver_Type]\","
            puts $FileVar "            \"smoother_type\":   \"ilu0\","
            puts $FileVar "            \"krylov_type\":     \"gmres\","
            puts $FileVar "            \"coarsening_type\": \"aggregation\","
            puts $FileVar "            \"max_iteration\":   1000,"
            puts $FileVar "            \"verbosity\":       [GiD_AccessValue get gendata Verbosity],"
            puts $FileVar "            \"tolerance\":       1.0e-6,"
            puts $FileVar "            \"scaling\":         [GiD_AccessValue get gendata Scaling]"
        } elseif {[GiD_AccessValue get gendata Solver_Type] eq "bicgstab"} {
            puts $FileVar "            \"solver_type\":         \"[GiD_AccessValue get gendata Solver_Type]\","
            puts $FileVar "            \"tolerance\":           1.0e-6,"
            puts $FileVar "            \"max_iteration\":       1000,"
            puts $FileVar "            \"scaling\":             [GiD_AccessValue get gendata Scaling],"
            puts $FileVar "            \"preconditioner_type\": \"ilu0\""
        } else {
            puts $FileVar "            \"solver_type\":   \"[GiD_AccessValue get gendata Solver_Type]\","
            puts $FileVar "            \"scaling\":       [GiD_AccessValue get gendata Scaling]"
        }
    }
    puts $FileVar "        \},"
    ## problem_domain_sub_model_part_list
    set PutStrings \[
    # Soil_two_phase part
    AppendGroupNames PutStrings Soil_two_phase
    # Soil_drained part
    AppendGroupNames PutStrings Soil_drained
    # Soil_undrained part
    AppendGroupNames PutStrings Soil_undrained
    # Non_porous part
    AppendGroupNames PutStrings Non_porous
    # Beam part
    AppendGroupNames PutStrings Beam
    # Shell_thin_corotational part
    AppendGroupNames PutStrings Shell_thin_corotational
    # Shell_thick_corotational part
    AppendGroupNames PutStrings Shell_thick_corotational
    # Truss part
    AppendGroupNames PutStrings Truss
    # Anchor part
    AppendGroupNames PutStrings Anchor
    # Interface_drained Part
    AppendGroupNames PutStrings Interface_drained
    # Interface_undrained Part
    AppendGroupNames PutStrings Interface_undrained
    if {[GiD_Groups exists PropagationUnion_3d_6] eq 1} {
        append PutStrings \" PropagationUnion_3d_6 \" \]
    } else {
        set PutStrings [string trimright $PutStrings ,]
        append PutStrings \]
    }
    puts $FileVar "        \"problem_domain_sub_model_part_list\": $PutStrings,"
    ## processes_sub_model_part_list
    set PutStrings \[
    # Solid_Displacement
    AppendGroupNames PutStrings Solid_Displacement
    # Structural_Rotation
    AppendGroupNames PutStrings Structural_Rotation
    # Fluid_Pressure
    AppendGroupNames PutStrings Fluid_Pressure
    # Excavation
    AppendGroupNames PutStrings Excavation
    # Point_Load
    AppendGroupNames PutStrings Point_Load
    # Line_Load
    AppendGroupNames PutStrings Line_Load
    # Surface_Load
    AppendGroupNames PutStrings Surface_Load
    # Normal_Load
    AppendGroupNames PutStrings Normal_Load
    # Normal_Fluid_Flux
    AppendGroupNames PutStrings Normal_Fluid_Flux
    # Interface_Face_Load
    AppendGroupNames PutStrings Interface_Face_Load
    # Interface_Normal_Fluid_Flux
    AppendGroupNames PutStrings Interface_Normal_Fluid_Flux
    # Body_Acceleration
    AppendGroupNames PutStrings Body_Acceleration
    # Record_DISPLACEMENT
    AppendGroupNames PutStrings Record_DISPLACEMENT
    # Record_VELOCITY
    AppendGroupNames PutStrings Record_VELOCITY
    # Record_ACCELERATION
    AppendGroupNames PutStrings Record_ACCELERATION
    # Record_VOLUME_ACCELERATION
    AppendGroupNames PutStrings Record_VOLUME_ACCELERATION
    # Record_POINT_LOAD
    AppendGroupNames PutStrings Record_POINT_LOAD
    # Record_LINE_LOAD
    AppendGroupNames PutStrings Record_LINE_LOAD
    # Record_SURFACE_LOAD
    AppendGroupNames PutStrings Record_SURFACE_LOAD
    # Gap_Closure_Bars
    if {$IsGapClosure eq true} {
        set interface_Groups [list [GiD_Info conditions Interface_two_phase groups] [GiD_Info conditions Interface_drained groups] [GiD_Info conditions Interface_undrained groups]]
        foreach Groups $interface_Groups {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 135] eq true} {
                    append PutStrings \" Gap_Closure_Bars_[lindex [lindex $Groups $i] 1] \" ,
                }
            }
        }
    }

    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar "        \"processes_sub_model_part_list\":      $PutStrings,"
    ## body_domain_sub_model_part_list
    set PutStrings \[
    AppendGroupNames PutStrings Soil_two_phase
    AppendGroupNames PutStrings Soil_drained
    AppendGroupNames PutStrings Soil_undrained
    AppendGroupNames PutStrings Non_porous
    AppendGroupNames PutStrings Beam
    AppendGroupNames PutStrings Shell_thin_corotational
    AppendGroupNames PutStrings Shell_thick_corotational
    AppendGroupNames PutStrings Truss
    AppendGroupNames PutStrings Anchor
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    if {[GiD_AccessValue get gendata Strategy_Type] eq "Arc-Length"} {
        puts $FileVar "        \"body_domain_sub_model_part_list\":    $PutStrings,"
        ## loads_sub_model_part_list
        set PutStrings \[
        set iGroup 0
        # Point_Load
        AppendGroupNamesWithNum PutStrings iGroup Point_Load
        # Line_Load
        AppendGroupNamesWithNum PutStrings iGroup Line_Load
        # Surface_Load
        AppendGroupNamesWithNum PutStrings iGroup Surface_Load
        # Normal_Load
        AppendGroupNamesWithNum PutStrings iGroup Normal_Load
        # Normal_Fluid_Flux
        AppendGroupNamesWithNum PutStrings iGroup Normal_Fluid_Flux
        # Interface_Face_Load
        AppendGroupNamesWithNum PutStrings iGroup Interface_Face_Load
        # Interface_Normal_Fluid_Flux
        AppendGroupNamesWithNum PutStrings iGroup Interface_Normal_Fluid_Flux
        # Body_Acceleration
        AppendGroupNamesWithNum PutStrings iGroup Body_Acceleration
        if {$iGroup > 0} {
            set PutStrings [string trimright $PutStrings ,]
        }
        append PutStrings \]
        puts $FileVar "        \"loads_sub_model_part_list\":          $PutStrings,"
        ## loads_variable_list
        set PutStrings \[
        # Point_Load
        AppendGroupVariables PutStrings Point_Load POINT_LOAD
        # Line_Load
        AppendGroupVariables PutStrings Line_Load LINE_LOAD
        # Surface_Load
        AppendGroupVariables PutStrings Surface_Load SURFACE_LOAD
        # Normal_Load
        AppendGroupVariables PutStrings Normal_Load NORMAL_CONTACT_STRESS
        # Normal_Fluid_Flux
        AppendGroupVariables PutStrings Normal_Fluid_Flux NORMAL_FLUID_FLUX
        # Interface_Face_Load
        AppendGroupVariables PutStrings Interface_Face_Load FACE_LOAD
        # Interface_Normal_Fluid_Flux
        AppendGroupVariables PutStrings Interface_Normal_Fluid_Flux NORMAL_FLUID_FLUX
        # Body_Acceleration
        AppendGroupVariables PutStrings Body_Acceleration VOLUME_ACCELERATION
        if {$iGroup > 0} {
            set PutStrings [string trimright $PutStrings ,]
        }
        append PutStrings \]
        puts $FileVar "        \"loads_variable_list\":                $PutStrings"
        puts $FileVar "    \},"
    } else {
        puts $FileVar "        \"body_domain_sub_model_part_list\":    $PutStrings"
        puts $FileVar "    \},"
    }

    ## Output processes
    puts $FileVar "    \"output_processes\": \{"
    puts $FileVar "        \"gid_output\": \[\{"
    puts $FileVar "            \"python_module\": \"gid_output_process\","
    puts $FileVar "            \"kratos_module\": \"KratosMultiphysics\","
    puts $FileVar "            \"process_name\": \"GiDOutputProcess\","
    puts $FileVar "            \"Parameters\":    \{"
    puts $FileVar "                \"model_part_name\": \"PorousDomain.porous_computational_model_part\","
    puts $FileVar "                \"output_name\": \"$basename\","
    puts $FileVar "                \"postprocess_parameters\": \{"
    puts $FileVar "                    \"result_file_configuration\": \{"
    puts $FileVar "                        \"gidpost_flags\":       \{"
    puts $FileVar "                            \"WriteDeformedMeshFlag\": \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $FileVar "                            \"WriteConditionsFlag\":   \"[GiD_AccessValue get gendata Write_conditions]\","
    puts $FileVar "                            \"GiDPostMode\":           \"[GiD_AccessValue get gendata GiD_post_mode]\","
    puts $FileVar "                            \"MultiFileFlag\":         \"[GiD_AccessValue get gendata Multi_file_flag]\""
    puts $FileVar "                        \},"
    puts $FileVar "                        \"file_label\":          \"[GiD_AccessValue get gendata File_label]\","
    puts $FileVar "                        \"output_control_type\": \"[GiD_AccessValue get gendata Output_control_type]\","
    puts $FileVar "                        \"output_interval\":     [GiD_AccessValue get gendata Output_Interval],"
    puts $FileVar "                        \"body_output\":         [GiD_AccessValue get gendata Body_output],"
    puts $FileVar "                        \"node_output\":         [GiD_AccessValue get gendata Node_output],"
    puts $FileVar "                        \"skin_output\":         [GiD_AccessValue get gendata Skin_output],"
    puts $FileVar "                        \"plane_output\":        \[\],"
    # nodal_results
    set PutStrings \[
    set iGroup 0
    AppendOutputVariables PutStrings iGroup Write_Solid_Displacement DISPLACEMENT
    AppendOutputVariables PutStrings iGroup Write_Solid_Displacement TOTAL_DISPLACEMENT
    AppendOutputVariables PutStrings iGroup Write_Structural_Rotation ROTATION
    AppendOutputVariables PutStrings iGroup Write_Fluid_Pressure WATER_PRESSURE
    if {[GiD_AccessValue get gendata Write_Reactions] eq true} {
        incr iGroup
        append PutStrings \" REACTION \" , \" REACTION_WATER_PRESSURE \" ,
    }
    AppendOutputVariables PutStrings iGroup Write_Point_Load POINT_LOAD
    AppendOutputVariables PutStrings iGroup Write_Line_Load LINE_LOAD
    AppendOutputVariables PutStrings iGroup Write_Surface_Load SURFACE_LOAD
    AppendOutputVariables PutStrings iGroup Write_Normal_Load NORMAL_CONTACT_STRESS
    AppendOutputVariables PutStrings iGroup Write_Tangential_Load TANGENTIAL_CONTACT_STRESS
    AppendOutputVariables PutStrings iGroup Write_Normal_Fluid_Flux NORMAL_FLUID_FLUX
    AppendOutputVariables PutStrings iGroup Write_Body_Acceleration VOLUME_ACCELERATION
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        incr iGroup
        append PutStrings \" PARTITION_INDEX \" ,
    }
    # Nodal smoothed variables
    if {[GiD_AccessValue get gendata Nodal_Smoothing] eq true} {
        AppendOutputVariables PutStrings iGroup Write_Effective_Stress NODAL_CAUCHY_STRESS_TENSOR
        AppendOutputVariables PutStrings iGroup Write_Damage NODAL_DAMAGE_VARIABLE
        AppendOutputVariables PutStrings iGroup Write_Joint_Width NODAL_JOINT_WIDTH
        AppendOutputVariables PutStrings iGroup Write_Damage NODAL_JOINT_DAMAGE
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $FileVar "                        \"nodal_results\":       $PutStrings,"
    # gauss_point_results
    set PutStrings \[
    set iGroup 0
    AppendOutputVariables PutStrings iGroup Write_Structural_Moment MOMENT
    AppendOutputVariables PutStrings iGroup Write_Structural_Force FORCE
    AppendOutputVariables PutStrings iGroup Write_Strain GREEN_LAGRANGE_STRAIN_TENSOR
    AppendOutputVariables PutStrings iGroup Write_Effective_Stress CAUCHY_STRESS_TENSOR
    AppendOutputVariables PutStrings iGroup Write_Total_Stress TOTAL_STRESS_TENSOR
    AppendOutputVariables PutStrings iGroup Write_Von_Mises_Stress VON_MISES_STRESS
    AppendOutputVariables PutStrings iGroup Write_Fluid_Flux FLUID_FLUX_VECTOR
    AppendOutputVariables PutStrings iGroup Write_Permeability PERMEABILITY_MATRIX
    AppendOutputVariables PutStrings iGroup Write_Damage DAMAGE_VARIABLE
    AppendOutputVariables PutStrings iGroup Write_Joint_Width JOINT_WIDTH
    AppendOutputVariables PutStrings iGroup Write_Local_Stress_Vector LOCAL_STRESS_VECTOR
    AppendOutputVariables PutStrings iGroup Write_Local_Relative_Displacement LOCAL_RELATIVE_DISPLACEMENT_VECTOR
    AppendOutputVariables PutStrings iGroup Write_Local_Fluid_Flux LOCAL_FLUID_FLUX_VECTOR
    AppendOutputVariables PutStrings iGroup Write_Local_Permeability LOCAL_PERMEABILITY_MATRIX
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $FileVar "                        \"gauss_point_results\": $PutStrings"
    puts $FileVar "                    \},"
    puts $FileVar "                    \"point_data_configuration\":  \[\]"
    puts $FileVar "                \}"
    puts $FileVar "            \}"
    puts $FileVar "        \}\]"
    puts $FileVar "    \},"

    ## Processes
    puts $FileVar "    \"processes\": \{"
    ## constraints_process_list
    set Groups [GiD_Info conditions Solid_Displacement groups]
    set NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Structural_Rotation groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Fluid_Pressure groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Excavation groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Record_DISPLACEMENT groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Record_VELOCITY groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Record_ACCELERATION groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Record_VOLUME_ACCELERATION groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Record_POINT_LOAD groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Record_LINE_LOAD groups]
    incr NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Record_SURFACE_LOAD groups]
    incr NumGroups [llength $Groups]

    set iGroup 0
    puts $FileVar "        \"constraints_process_list\": \[\{"
    # Solid_Displacement
    set Groups [GiD_Info conditions Solid_Displacement groups]
    WriteConstraintVectorProcess FileVar iGroup $Groups volumes DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups surfaces DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups lines DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups points DISPLACEMENT $TableDict $NumGroups

    # Structural_Rotation
    set Groups [GiD_Info conditions Structural_Rotation groups]
    WriteConstraintVectorProcess FileVar iGroup $Groups volumes ROTATION $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups surfaces ROTATION $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups lines ROTATION $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups points ROTATION $TableDict $NumGroups

    # Note: it is important to write processes in the following order to account for intersections between conditions
    # Fluid_Pressure
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    WritePressureConstraintProcess FileVar iGroup $Groups volumes WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups surfaces WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups lines WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups points WATER_PRESSURE $TableDict $NumGroups

    # Excavation
    set Groups [GiD_Info conditions Excavation groups]
    WriteExcavationConstraintProcess FileVar iGroup $Groups volumes EXCAVATION $NumGroups
    WriteExcavationConstraintProcess FileVar iGroup $Groups surfaces EXCAVATION $NumGroups
    WriteExcavationConstraintProcess FileVar iGroup $Groups lines EXCAVATION $NumGroups
    WriteExcavationConstraintProcess FileVar iGroup $Groups points EXCAVATION $NumGroups

    # Record_DISPLACEMENT
    set Groups [GiD_Info conditions Record_DISPLACEMENT groups]
    WriteResultVectorProcess FileVar iGroup $Groups points DISPLACEMENT $NumGroups

    # Record_VELOCITY
    set Groups [GiD_Info conditions Record_VELOCITY groups]
    WriteResultVectorProcess FileVar iGroup $Groups points DISPLACEMENT $NumGroups

    # Record_ACCELERATION
    set Groups [GiD_Info conditions Record_ACCELERATION groups]
    WriteResultVectorProcess FileVar iGroup $Groups points DISPLACEMENT $NumGroups

    # Record_VOLUME_ACCELERATION
    set Groups [GiD_Info conditions Record_VOLUME_ACCELERATION groups]
    WriteResultVectorProcess FileVar iGroup $Groups points DISPLACEMENT $NumGroups

    # Record_POINT_LOAD
    set Groups [GiD_Info conditions Record_POINT_LOAD groups]
    WriteResultVectorProcess FileVar iGroup $Groups points DISPLACEMENT $NumGroups

    # Record_LINE_LOAD
    set Groups [GiD_Info conditions Record_LINE_LOAD groups]
    WriteResultVectorProcess FileVar iGroup $Groups points DISPLACEMENT $NumGroups

    # Record_SURFACE_LOAD
    set Groups [GiD_Info conditions Record_SURFACE_LOAD groups]
    WriteResultVectorProcess FileVar iGroup $Groups points DISPLACEMENT $NumGroups


    puts $FileVar "    \}\],"
                
    ## loads_process_list
    set Groups [GiD_Info conditions Point_Load groups]
    set NumGroups [llength $Groups]

    set Groups [GiD_Info conditions Line_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Surface_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Normal_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Body_Acceleration groups]
    incr NumGroups [llength $Groups]

    if {$NumGroups > 0} {
        set iGroup 0
        puts $FileVar "    \"loads_process_list\": \[\{"
        # Point_Load
        set Groups [GiD_Info conditions Point_Load groups]
        WriteLoadVectorProcess FileVar iGroup $Groups POINT_LOAD $TableDict $NumGroups
        # Line_Load
        set Groups [GiD_Info conditions Line_Load groups]
        WriteLoadVectorProcess FileVar iGroup $Groups LINE_LOAD $TableDict $NumGroups
        # Surface_Load
        set Groups [GiD_Info conditions Surface_Load groups]
        WriteLoadVectorProcess FileVar iGroup $Groups SURFACE_LOAD $TableDict $NumGroups
        # Normal_Load
        set Groups [GiD_Info conditions Normal_Load groups]
        WriteNormalLoadProcess FileVar iGroup $Groups NORMAL_CONTACT_STRESS $TableDict $NumGroups
        # Normal_Fluid_Flux
        set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
        WriteLoadScalarProcess FileVar iGroup $Groups NORMAL_FLUID_FLUX $TableDict $NumGroups
        # Interface_Face_Load
        set Groups [GiD_Info conditions Interface_Face_Load groups]
        WriteLoadVectorProcess FileVar iGroup $Groups FACE_LOAD $TableDict $NumGroups
        # Interface_Normal_Fluid_Flux
        set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
        WriteLoadScalarProcess FileVar iGroup $Groups NORMAL_FLUID_FLUX $TableDict $NumGroups
        # Body_Acceleration
        set Groups [GiD_Info conditions Body_Acceleration groups]
        WriteLoadVectorProcess FileVar iGroup $Groups VOLUME_ACCELERATION $TableDict $NumGroups
    } else {
        puts $FileVar "    \"loads_process_list\":       \[\],"
    }

    ## auxiliar_process_list
    set NumGroups 0
    if {$IsGapClosure eq true} {
        set interface_Groups [list [GiD_Info conditions Interface_two_phase groups] [GiD_Info conditions Interface_drained groups] [GiD_Info conditions Interface_undrained groups]]
        foreach Groups $interface_Groups {
            for {set i 0} {$i < [llength $Groups]} {incr i} {
                if {[lindex [lindex $Groups $i] 135] eq true} {
                    incr NumGroups
                }
            }
        }
    }
    
    if {$NumGroups > 0} {
        set iGroup 0
        puts $FileVar "        \"auxiliar_process_list\": \[\{"
        # Gap_Closure_Bars
        if {$IsGapClosure eq true} {
            set interface_Groups [list [GiD_Info conditions Interface_two_phase groups] [GiD_Info conditions Interface_drained groups] [GiD_Info conditions Interface_undrained groups]]
            foreach Groups $interface_Groups {
                WriteGapClosureInterfaceProcess FileVar iGroup $Groups $NumGroups
            }
        }
    } else {
        puts $FileVar "        \"auxiliar_process_list\": \[\]"
    }

    puts $FileVar "    \}"
    puts $FileVar "\}"

    close $FileVar
}
