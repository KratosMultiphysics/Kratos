proc WriteProjectParameters { basename dir problemtypedir TableDict} {

    ## Source auxiliar procedures
    source [file join $problemtypedir ProjectParametersAuxProcs.tcl]

    ## Start ProjectParameters.json file
    set filename [file join $dir ProjectParameters.json]
    set FileVar [open $filename w]

    puts $FileVar "\{"

    ## problem_data
    puts $FileVar "    \"problem_data\": \{"
    puts $FileVar "        \"problem_name\":         \"$basename\","
    puts $FileVar "        \"start_time\":           [GiD_AccessValue get gendata Start_Time],"
    puts $FileVar "        \"end_time\":             [GiD_AccessValue get gendata End_Time],"
    puts $FileVar "        \"echo_level\":           [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "        \"parallel_type\":        \"[GiD_AccessValue get gendata Parallel_Configuration]\","
    puts $FileVar "        \"number_of_threads\":    [GiD_AccessValue get gendata Number_of_threads],"
    if {[GiD_AccessValue get gendata Initial_Stresses] eq false} {
        puts $FileVar "        \"fracture_utility\":     [GiD_AccessValue get gendata Fracture_Propagation]"
    } else {
        puts $FileVar "        \"fracture_utility\":     [GiD_AccessValue get gendata Fracture_Propagation],"
        puts $FileVar "        \"initial_stress_utility_settings\":   \{"
        puts $FileVar "            \"mode\":       \"[GiD_AccessValue get gendata Mode]\","
        puts $FileVar "            \"initial_input_filename\":   \"initial_$basename\""
        puts $FileVar "        \}"
    }
    puts $FileVar "    \},"

    ## solver_settings
    puts $FileVar "    \"solver_settings\": \{"
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        puts $FileVar "        \"solver_type\":                        \"poromechanics_MPI_U_Pw_solver\","
    } else {
        puts $FileVar "        \"solver_type\":                        \"poromechanics_U_Pw_solver\","
    }
    puts $FileVar "        \"model_part_name\":                    \"PorousModelPart\","
    puts $FileVar "        \"domain_size\":                        [GiD_AccessValue get gendata Domain_Size],"
    puts $FileVar "        \"start_time\":                         [GiD_AccessValue get gendata Start_Time],"
    puts $FileVar "        \"time_step\":                          [GiD_AccessValue get gendata Delta_Time],"
    puts $FileVar "        \"model_import_settings\":              \{"
    puts $FileVar "            \"input_type\":    \"mdpa\","
    puts $FileVar "            \"input_filename\":    \"$basename\""
    puts $FileVar "        \},"
    puts $FileVar "        \"material_import_settings\": \{"
    puts $FileVar "            \"materials_filename\":    \"PoroMaterials.json\""
    puts $FileVar "        \},"
    puts $FileVar "        \"buffer_size\":                        2,"
    puts $FileVar "        \"echo_level\":                         [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "        \"clear_storage\":                      false,"
    puts $FileVar "        \"compute_reactions\":                  [GiD_AccessValue get gendata Write_Reactions],"
    puts $FileVar "        \"move_mesh_flag\":                     [GiD_AccessValue get gendata Move_Mesh],"
    set IsPeriodic [GiD_AccessValue get gendata Periodic_Interface_Conditions]
    if {$IsPeriodic eq true} {
        puts $FileVar "        \"periodic_interface_conditions\":      true,"
        puts $FileVar "        \"reform_dofs_at_each_step\":           true,"
        puts $FileVar "        \"nodal_smoothing\":                    true,"
    } else {
        puts $FileVar "        \"periodic_interface_conditions\":      false,"
        puts $FileVar "        \"reform_dofs_at_each_step\":           [GiD_AccessValue get gendata Reform_Dofs_At_Each_Step],"
        puts $FileVar "        \"nodal_smoothing\":                    [GiD_AccessValue get gendata Nodal_Smoothing],"
    }
    puts $FileVar "        \"block_builder\":                      [GiD_AccessValue get gendata Block_Builder],"
    puts $FileVar "        \"solution_type\":                      \"[GiD_AccessValue get gendata Solution_Type]\","
    puts $FileVar "        \"scheme_type\":                        \"[GiD_AccessValue get gendata Scheme_Type]\","
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
    puts $FileVar "        \"max_iteration\":                      [GiD_AccessValue get gendata Max_Iterations],"
    puts $FileVar "        \"desired_iterations\":                 [GiD_AccessValue get gendata Desired_Iterations],"
    puts $FileVar "        \"max_radius_factor\":                  [GiD_AccessValue get gendata Max_Radius_Factor],"
    puts $FileVar "        \"min_radius_factor\":                  [GiD_AccessValue get gendata Min_Radius_Factor],"
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        puts $FileVar "        \"nonlocal_damage\":                    false,"
    } else {
        puts $FileVar "        \"nonlocal_damage\":                    [GiD_AccessValue get gendata Non-local_Damage],"
    }
    puts $FileVar "        \"characteristic_length\":              [GiD_AccessValue get gendata Characteristic_Length],"
    ## linear_solver_settings
    puts $FileVar "        \"linear_solver_settings\":             \{"
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        if {[GiD_AccessValue get gendata Solver_Type] eq "amgcl"} {
            puts $FileVar "            \"solver_type\":   \"amgcl\","
            puts $FileVar "            \"krylov_type\":   \"fgmres\","
            puts $FileVar "            \"max_iteration\": 100,"
            puts $FileVar "            \"verbosity\":     [GiD_AccessValue get gendata Verbosity],"
            puts $FileVar "            \"tolerance\":     1.0e-6,"
            puts $FileVar "            \"scaling\":       [GiD_AccessValue get gendata Scaling]"
        } elseif {[GiD_AccessValue get gendata Solver_Type] eq "aztec"} {
            puts $FileVar "            \"solver_type\":         \"aztec\","
            puts $FileVar "            \"tolerance\":           1.0e-6,"
            puts $FileVar "            \"max_iteration\":       200,"
            puts $FileVar "            \"scaling\":             [GiD_AccessValue get gendata Scaling],"
            puts $FileVar "            \"preconditioner_type\": \"None\""
        } elseif {([GiD_AccessValue get gendata Solver_Type] eq "klu") || ([GiD_AccessValue get gendata Solver_Type] eq "multi_level")} {
            puts $FileVar "            \"solver_type\": \"[GiD_AccessValue get gendata Solver_Type]\","
            puts $FileVar "            \"scaling\":     [GiD_AccessValue get gendata Scaling]"
        } else {
            puts $FileVar "            \"solver_type\": \"klu\","
            puts $FileVar "            \"scaling\":     false"
        }
    } else {
        if {[GiD_AccessValue get gendata Solver_Type] eq "amgcl"} {
            puts $FileVar "            \"solver_type\":     \"amgcl\","
            puts $FileVar "            \"smoother_type\":   \"ilu0\","
            puts $FileVar "            \"krylov_type\":     \"gmres\","
            puts $FileVar "            \"coarsening_type\": \"aggregation\","
            puts $FileVar "            \"max_iteration\":   100,"
            puts $FileVar "            \"verbosity\":       [GiD_AccessValue get gendata Verbosity],"
            puts $FileVar "            \"tolerance\":       1.0e-6,"
            puts $FileVar "            \"scaling\":         [GiD_AccessValue get gendata Scaling]"
        } elseif {[GiD_AccessValue get gendata Solver_Type] eq "bicgstab"} {
            puts $FileVar "            \"solver_type\":         \"bicgstab\","
            puts $FileVar "            \"tolerance\":           1.0e-6,"
            puts $FileVar "            \"max_iteration\":       100,"
            puts $FileVar "            \"scaling\":             [GiD_AccessValue get gendata Scaling],"
            puts $FileVar "            \"preconditioner_type\": \"ilu0\""
        } elseif {([GiD_AccessValue get gendata Solver_Type] eq "skyline_lu_factorization") || ([GiD_AccessValue get gendata Solver_Type] eq "ExternalSolversApplication.super_lu")} {
            puts $FileVar "            \"solver_type\":   \"[GiD_AccessValue get gendata Solver_Type]\""
        } else {
            puts $FileVar "            \"solver_type\":   \"ExternalSolversApplication.super_lu\""
        }
    }
    puts $FileVar "        \},"
    ## problem_domain_sub_model_part_list
    set PutStrings \[
    # Body_Part
    AppendGroupNames PutStrings Body_Part
    # Interface_Part
    AppendGroupNames PutStrings Interface_Part
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
    # Fluid_Pressure
    AppendGroupNames PutStrings Fluid_Pressure
    # Force
    AppendGroupNames PutStrings Force
    # Face_Load
    AppendGroupNames PutStrings Face_Load
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
    # Periodic_Bars
    if {$IsPeriodic eq true} {
        set Groups [GiD_Info conditions Interface_Part groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            if {[lindex [lindex $Groups $i] 20] eq true} {
                append PutStrings \" Periodic_Bars_[lindex [lindex $Groups $i] 1] \" ,
            }
        }
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar "        \"processes_sub_model_part_list\":      $PutStrings,"
    ## body_domain_sub_model_part_list
    set PutStrings \[
    AppendGroupNames PutStrings Body_Part
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    if {[GiD_AccessValue get gendata Strategy_Type] eq "arc_length"} {
        puts $FileVar "        \"body_domain_sub_model_part_list\":    $PutStrings,"
        ## loads_sub_model_part_list
        set PutStrings \[
        set iGroup 0
        # Force
        AppendGroupNamesWithNum PutStrings iGroup Force
        # Face_Load
        AppendGroupNamesWithNum PutStrings iGroup Face_Load
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
        # Force
        AppendGroupVariables PutStrings Force FORCE
        # Face_Load
        AppendGroupVariables PutStrings Face_Load FACE_LOAD
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
    puts $FileVar "                \"model_part_name\": \"PorousModelPart.porous_computational_model_part\","
    puts $FileVar "                \"output_name\": \"$basename\","
    puts $FileVar "                \"postprocess_parameters\": \{"
    puts $FileVar "                    \"result_file_configuration\": \{"
    puts $FileVar "                        \"gidpost_flags\":       \{"
    puts $FileVar "                            \"WriteDeformedMeshFlag\": \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $FileVar "                            \"WriteConditionsFlag\":   \"[GiD_AccessValue get gendata Write_conditions]\","
    if { ([GiD_AccessValue get gendata Fracture_Propagation] eq true) || ($IsPeriodic eq true) } {
        puts $FileVar "                            \"GiDPostMode\":           \"GiD_PostAscii\","
        puts $FileVar "                            \"MultiFileFlag\":         \"MultipleFiles\""
        puts $FileVar "                        \},"
        puts $FileVar "                        \"file_label\":          \"time\","
    } else {
        puts $FileVar "                            \"GiDPostMode\":           \"[GiD_AccessValue get gendata GiD_post_mode]\","
        puts $FileVar "                            \"MultiFileFlag\":         \"[GiD_AccessValue get gendata Multi_file_flag]\""
        puts $FileVar "                        \},"
        puts $FileVar "                        \"file_label\":          \"[GiD_AccessValue get gendata File_label]\","
    }
    puts $FileVar "                        \"output_control_type\": \"[GiD_AccessValue get gendata Output_control_type]\","
    puts $FileVar "                        \"output_frequency\":    [GiD_AccessValue get gendata Output_frequency],"
    puts $FileVar "                        \"body_output\":         [GiD_AccessValue get gendata Body_output],"
    puts $FileVar "                        \"node_output\":         [GiD_AccessValue get gendata Node_output],"
    puts $FileVar "                        \"skin_output\":         [GiD_AccessValue get gendata Skin_output],"
    puts $FileVar "                        \"plane_output\":        \[\],"
    # nodal_results
    set PutStrings \[
    set iGroup 0
    AppendOutputVariables PutStrings iGroup Write_Solid_Displacement DISPLACEMENT
    AppendOutputVariables PutStrings iGroup Write_Fluid_Pressure WATER_PRESSURE
    if {[GiD_AccessValue get gendata Write_Reactions] eq true} {
        incr iGroup
        append PutStrings \" REACTION \" , \" REACTION_WATER_PRESSURE \" ,
    }
    AppendOutputVariables PutStrings iGroup Write_Force FORCE
    AppendOutputVariables PutStrings iGroup Write_Face_Load FACE_LOAD
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
    AppendOutputVariables PutStrings iGroup Write_Initial_Stress INITIAL_STRESS_TENSOR
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $FileVar "                        \"nodal_results\":       $PutStrings,"
    # gauss_point_results
    set PutStrings \[
    set iGroup 0
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
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    incr NumGroups [llength $Groups]
    set iGroup 0
    puts $FileVar "        \"constraints_process_list\": \[\{"
    # Solid_Displacement
    set Groups [GiD_Info conditions Solid_Displacement groups]
    WriteConstraintVectorProcess FileVar iGroup $Groups volumes DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups surfaces DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups lines DISPLACEMENT $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups points DISPLACEMENT $TableDict $NumGroups
    # Note: it is important to write processes in the following order to account for intersections between conditions
    # Fluid_Pressure
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    WritePressureConstraintProcess FileVar iGroup $Groups volumes WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups surfaces WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups lines WATER_PRESSURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups points WATER_PRESSURE $TableDict $NumGroups
    ## loads_process_list
    set Groups [GiD_Info conditions Force groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Face_Load groups]
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
        puts $FileVar "        \"loads_process_list\": \[\{"
        # Force
        set Groups [GiD_Info conditions Force groups]
        WriteLoadVectorProcess FileVar iGroup $Groups FORCE $TableDict $NumGroups
        # Face_Load
        set Groups [GiD_Info conditions Face_Load groups]
        WriteLoadVectorProcess FileVar iGroup $Groups FACE_LOAD $TableDict $NumGroups
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
        puts $FileVar "        \"loads_process_list\":       \[\],"
    }
    ## auxiliar_process_list
    set NumGroups 0
    if {$IsPeriodic eq true} {
        set Groups [GiD_Info conditions Interface_Part groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            if {[lindex [lindex $Groups $i] 20] eq true} {
                incr NumGroups
            }
        }
    }
    if {$NumGroups > 0} {
        set iGroup 0
        puts $FileVar "        \"auxiliar_process_list\": \[\{"
        # Periodic_Bars
        if {$IsPeriodic eq true} {
            set Groups [GiD_Info conditions Interface_Part groups]
            WritePeriodicInterfaceProcess FileVar iGroup $Groups $NumGroups
        }
    } else {
        puts $FileVar "        \"auxiliar_process_list\": \[\]"
    }

    puts $FileVar "    \}"
    puts $FileVar "\}"

    close $FileVar
}
