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
    puts $FileVar "        \"model_part_name\":      \"FluidTransportDomain\","
    puts $FileVar "        \"domain_size\":          [GiD_AccessValue get gendata Domain_Size],"
    puts $FileVar "        \"start_time\":           [GiD_AccessValue get gendata Start_Time],"
    puts $FileVar "        \"end_time\":             [GiD_AccessValue get gendata End_Time],"
    puts $FileVar "        \"time_step\":            [GiD_AccessValue get gendata Delta_Time],"
    puts $FileVar "        \"parallel_type\":        \"[GiD_AccessValue get gendata Parallel_Configuration]\","
    puts $FileVar "        \"number_of_threads\":    [GiD_AccessValue get gendata Number_of_threads]"
    puts $FileVar "    \},"

    ## solver_settings
    puts $FileVar "    \"solver_settings\": \{"
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        puts $FileVar "        \"solver_type\":                        \"fluid_transport_solver\","
    } else {
        puts $FileVar "        \"solver_type\":                        \"fluid_transport_solver\","
    }
    puts $FileVar "        \"model_import_settings\":              \{"
    puts $FileVar "            \"input_type\":       \"mdpa\","
    puts $FileVar "            \"input_filename\":   \"$basename\","
    puts $FileVar "            \"input_file_label\": 0"
    puts $FileVar "        \},"
    puts $FileVar "        \"buffer_size\":                        3,"
    puts $FileVar "        \"echo_level\":                         [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "        \"clear_storage\":                      false,"
    puts $FileVar "        \"compute_reactions\":                  [GiD_AccessValue get gendata Write_Reactions],"
    puts $FileVar "        \"move_mesh_flag\":                     [GiD_AccessValue get gendata Move_Mesh],"
    puts $FileVar "        \"reform_dofs_at_each_step\":           [GiD_AccessValue get gendata Reform_Dofs_At_Each_Step],"
    puts $FileVar "        \"block_builder\":                      [GiD_AccessValue get gendata Block_Builder],"
    puts $FileVar "        \"solution_type\":                      \"[GiD_AccessValue get gendata Solution_Type]\","
    puts $FileVar "        \"scheme_type\":                        \"[GiD_AccessValue get gendata Scheme]\","
    puts $FileVar "        \"newmark_theta\":                      [GiD_AccessValue get gendata Newmark_Theta],"
    puts $FileVar "        \"strategy_type\":                      \"[GiD_AccessValue get gendata Strategy_Type]\","
    puts $FileVar "        \"convergence_criterion\":              \"[GiD_AccessValue get gendata Convergence_Criterion]\","
    puts $FileVar "        \"displacement_relative_tolerance\":    [GiD_AccessValue get gendata Displacement_Relative_Tolerance],"
    puts $FileVar "        \"displacement_absolute_tolerance\":    [GiD_AccessValue get gendata Displacement_Absolute_Tolerance],"
    puts $FileVar "        \"residual_relative_tolerance\":        [GiD_AccessValue get gendata Residual_Relative_Tolerance],"
    puts $FileVar "        \"residual_absolute_tolerance\":        [GiD_AccessValue get gendata Residual_Absolute_Tolerance],"
    puts $FileVar "        \"max_iteration\":                      [GiD_AccessValue get gendata Max_Iterations],"
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
        if {[GiD_AccessValue get gendata Solver_Type] eq "AMGCL"} {
            puts $FileVar "            \"solver_type\":     \"AMGCL\","
            puts $FileVar "            \"smoother_type\":   \"ilu0\","
            puts $FileVar "            \"krylov_type\":     \"gmres\","
            puts $FileVar "            \"coarsening_type\": \"aggregation\","
            puts $FileVar "            \"max_iteration\":   100,"
            puts $FileVar "            \"verbosity\":       [GiD_AccessValue get gendata Verbosity],"
            puts $FileVar "            \"tolerance\":       1.0e-6,"
            puts $FileVar "            \"scaling\":         [GiD_AccessValue get gendata Scaling]"
        } elseif {[GiD_AccessValue get gendata Solver_Type] eq "BICGSTABSolver"} {
            puts $FileVar "            \"solver_type\":         \"BICGSTABSolver\","
            puts $FileVar "            \"tolerance\":           1.0e-6,"
            puts $FileVar "            \"max_iteration\":       100,"
            puts $FileVar "            \"scaling\":             [GiD_AccessValue get gendata Scaling],"
            puts $FileVar "            \"preconditioner_type\": \"ILU0Preconditioner\""
        } elseif {([GiD_AccessValue get gendata Solver_Type] eq "SkylineLUFactorizationSolver") || ([GiD_AccessValue get gendata Solver_Type] eq "SuperLUSolver")} {
            puts $FileVar "            \"solver_type\":   \"[GiD_AccessValue get gendata Solver_Type]\""
        } else {
            puts $FileVar "            \"solver_type\":   \"SuperLUSolver\""
        }
    }
    puts $FileVar "        \},"
    ## problem_domain_sub_model_part_list
    set PutStrings \[
    # Body_Part
    AppendGroupNames PutStrings Body_Part
    set PutStrings [string trimright $PutStrings ,]
    puts $FileVar "        \"problem_domain_sub_model_part_list\": $PutStrings\],"
    ## processes_sub_model_part_list
    set PutStrings \[
    # Velocity
    AppendGroupNames PutStrings Velocity
    # Phi_Value
    AppendGroupNames PutStrings Phi_Value
    # Face_Heat_Flux
    AppendGroupNames PutStrings Face_Heat_Flux
    # Q_Source
    AppendGroupNames PutStrings Q_Source
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar "        \"processes_sub_model_part_list\":      $PutStrings"

    puts $FileVar "    \},"

    ## output_configuration
    puts $FileVar "    \"output_configuration\": \{"
    puts $FileVar "        \"result_file_configuration\": \{"
    puts $FileVar "            \"gidpost_flags\":       \{"
    puts $FileVar "                \"WriteDeformedMeshFlag\": \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $FileVar "                \"WriteConditionsFlag\":   \"[GiD_AccessValue get gendata Write_conditions]\","
    puts $FileVar "                \"GiDPostMode\":           \"[GiD_AccessValue get gendata GiD_post_mode]\","
    puts $FileVar "                \"MultiFileFlag\":         \"[GiD_AccessValue get gendata Multi_file_flag]\""
    puts $FileVar "            \},"
    puts $FileVar "            \"file_label\":          \"[GiD_AccessValue get gendata File_label]\","
    puts $FileVar "            \"output_control_type\": \"[GiD_AccessValue get gendata Output_control_type]\","
    puts $FileVar "            \"output_frequency\":    [GiD_AccessValue get gendata Output_frequency],"
    puts $FileVar "            \"body_output\":         [GiD_AccessValue get gendata Body_output],"
    puts $FileVar "            \"node_output\":         [GiD_AccessValue get gendata Node_output],"
    puts $FileVar "            \"skin_output\":         [GiD_AccessValue get gendata Skin_output],"
    puts $FileVar "            \"plane_output\":        \[\],"
    # nodal_results
    set PutStrings \[
    set iGroup 0
    AppendOutputVariables PutStrings iGroup Write_Velocity VELOCITY

    set SolutionType [GiD_AccessValue get gendata Solution_Type]
    set SchemeType [GiD_AccessValue get gendata Scheme]

    if {$SolutionType eq "Steady"} {

    AppendOutputVariables PutStrings iGroup Write_Phi_Value TEMPERATURE

    } else {
        if {$SchemeType eq "Implicit"} {

            AppendOutputVariables PutStrings iGroup Write_Phi_Value TEMPERATURE
            AppendOutputVariables PutStrings iGroup Write_Phi_Value PHI_THETA

        } else {

            AppendOutputVariables PutStrings iGroup Write_Phi_Value TEMPERATURE

        }

    }

    AppendOutputVariables PutStrings iGroup Write_Normals_Value NORMAL
    if {[GiD_AccessValue get gendata Write_Reactions] eq true} {
        incr iGroup
        append PutStrings \" REACTION_FLUX \" ,
    }
    AppendOutputVariables PutStrings iGroup Write_Face_Heat_Flux FACE_HEAT_FLUX
    AppendOutputVariables PutStrings iGroup Write_Q_Source HEAT_FLUX
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        incr iGroup
        append PutStrings \" PARTITION_INDEX \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $FileVar "            \"nodal_results\":       $PutStrings,"
    # gauss_point_results
    set PutStrings \[
    set iGroup 0
    AppendOutputVariables PutStrings iGroup Write_Peclet PECLET
    AppendOutputVariables PutStrings iGroup Write_FIC_Beta FIC_BETA
    AppendOutputVariables PutStrings iGroup Write_Phi_Gradient PHI_GRADIENT

    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $FileVar "            \"gauss_point_results\": $PutStrings"
    puts $FileVar "        \},"
    puts $FileVar "        \"point_data_configuration\":  \[\]"

    puts $FileVar "    \},"

    ## constraints_process_list
    set Groups [GiD_Info conditions Velocity groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Phi_Value groups]
    incr NumGroups [llength $Groups]
    set iGroup 0
    puts $FileVar "    \"constraints_process_list\": \[\{"
    # Velocity
    set Groups [GiD_Info conditions Velocity groups]
    WriteConstraintVectorProcess FileVar iGroup $Groups volumes VELOCITY $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups surfaces VELOCITY $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups lines VELOCITY $TableDict $NumGroups
    WriteConstraintVectorProcess FileVar iGroup $Groups points VELOCITY $TableDict $NumGroups
    # Note: it is important to write processes in the following order to account for intersections between conditions
    # Phi_Value
    set Groups [GiD_Info conditions Phi_Value groups]

    if {$SolutionType eq "Steady"} {

    WritePressureConstraintProcess FileVar iGroup $Groups volumes TEMPERATURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups surfaces TEMPERATURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups lines TEMPERATURE $TableDict $NumGroups
    WritePressureConstraintProcess FileVar iGroup $Groups points TEMPERATURE $TableDict $NumGroups

    } else {

        if {$SchemeType eq "Implicit"} {

            incr NumGroups [llength $Groups]

            WritePressureConstraintProcess FileVar iGroup $Groups volumes PHI_THETA $TableDict $NumGroups
            WritePressureConstraintProcess FileVar iGroup $Groups surfaces PHI_THETA $TableDict $NumGroups
            WritePressureConstraintProcess FileVar iGroup $Groups lines PHI_THETA $TableDict $NumGroups
            WritePressureConstraintProcess FileVar iGroup $Groups points PHI_THETA $TableDict $NumGroups

            WriteTempConstraintProcess FileVar iGroup $Groups TEMPERATURE $TableDict $NumGroups

        } else {

            WritePressureConstraintProcess FileVar iGroup $Groups volumes TEMPERATURE $TableDict $NumGroups
            WritePressureConstraintProcess FileVar iGroup $Groups surfaces TEMPERATURE $TableDict $NumGroups
            WritePressureConstraintProcess FileVar iGroup $Groups lines TEMPERATURE $TableDict $NumGroups
            WritePressureConstraintProcess FileVar iGroup $Groups points TEMPERATURE $TableDict $NumGroups

        }

    }

    ## loads_process_list
    set Groups [GiD_Info conditions Face_Heat_Flux groups]
    incr NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Q_Source groups]
    set NumGroups [llength $Groups]

    if {$NumGroups > 0} {
        set iGroup 0
        puts $FileVar "    \"loads_process_list\": \[\{"
        # Face_Heat_Flux
        set Groups [GiD_Info conditions Face_Heat_Flux groups]
        WriteLoadScalarProcess FileVar iGroup $Groups FACE_HEAT_FLUX $TableDict $NumGroups
        # Q_Source
        set Groups [GiD_Info conditions Q_Source groups]
        WriteLoadScalarProcess FileVar iGroup $Groups HEAT_FLUX $TableDict $NumGroups
    } else {
        puts $FileVar "    \"loads_process_list\":       \[\]"
    }

    puts $FileVar "\}"

    close $FileVar
}
