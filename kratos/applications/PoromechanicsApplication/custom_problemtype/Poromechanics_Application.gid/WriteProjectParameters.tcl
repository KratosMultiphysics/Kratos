proc WriteProjectParameters { basename dir TableList} {
    
    ## Start ProjectParameters.json file
    set filename [file join $dir ProjectParameters.json]
    set varfile [open $filename w]
    
    puts $varfile "\{"
        
    ## problem_data
    puts $varfile "    \"problem_data\": \{"
    puts $varfile "        \"problem_name\":         \"$basename\","
    puts $varfile "        \"model_part_name\":      \"PorousDomain\","
    puts $varfile "        \"domain_size\":          [GiD_AccessValue get gendata Domain_Size],"
    puts $varfile "        \"start_time\":           [GiD_AccessValue get gendata Start_Time],"
    puts $varfile "        \"end_time\":             [GiD_AccessValue get gendata End_Time],"
    puts $varfile "        \"time_step\":            [GiD_AccessValue get gendata Delta_Time],"
    puts $varfile "        \"OMP_threads\":          [GiD_AccessValue get gendata OMP_Threads]"
    puts $varfile "    \},"
    
    ## solver_settings
    puts $varfile "    \"solver_settings\": \{"
    puts $varfile "        \"solver_type\":                        \"poromechanics_U_Pw_solver\","
    puts $varfile "        \"model_import_settings\":              \{"
    puts $varfile "            \"input_type\":     \"mdpa\","
    puts $varfile "            \"input_filename\": \"$basename\""
    puts $varfile "        \},"
    puts $varfile "        \"buffer_size\":                        2,"
    puts $varfile "        \"echo_level\":                         [GiD_AccessValue get gendata Echo_Level],"
    puts $varfile "        \"reform_dofs_at_each_step\":           false,"
    puts $varfile "        \"compute_reactions\":                  [GiD_AccessValue get gendata Write_Reactions],"
    puts $varfile "        \"move_mesh_flag\":                     true,"
    puts $varfile "        \"solution_type\":                      \"[GiD_AccessValue get gendata Solution_Type]\","
    puts $varfile "        \"scheme_type\":                        \"[GiD_AccessValue get gendata Scheme_Type]\","
    puts $varfile "        \"newmark_beta\":                       [GiD_AccessValue get gendata Newmark_Beta],"
    puts $varfile "        \"newmark_gamma\":                      [GiD_AccessValue get gendata Newmark_Gamma],"
    puts $varfile "        \"newmark_theta\":                      [GiD_AccessValue get gendata Newmark_Theta],"
    puts $varfile "        \"rayleigh_m\":                         [GiD_AccessValue get gendata Rayleigh_Mass],"
    puts $varfile "        \"rayleigh_k\":                         [GiD_AccessValue get gendata Rayleigh_Stiffness],"
    puts $varfile "        \"strategy_type\":                      \"[GiD_AccessValue get gendata Strategy_Type]\","
    puts $varfile "        \"fracture_propagation\":               [GiD_AccessValue get gendata Fracture_Propagation],"
    puts $varfile "        \"convergence_criterion\":              \"[GiD_AccessValue get gendata Convergence_Criterion]\","
    puts $varfile "        \"displacement_relative_tolerance\":    [GiD_AccessValue get gendata Displacement_Relative_Tolerance],"
    puts $varfile "        \"displacement_absolute_tolerance\":    [GiD_AccessValue get gendata Displacement_Absolute_Tolerance],"
    puts $varfile "        \"residual_relative_tolerance\":        [GiD_AccessValue get gendata Residual_Relative_Tolerance],"
    puts $varfile "        \"residual_absolute_tolerance\":        [GiD_AccessValue get gendata Residual_Absolute_Tolerance],"
    puts $varfile "        \"max_iteration\":                      [GiD_AccessValue get gendata Max_Iterations],"
    puts $varfile "        \"desired_iterations\":                 [GiD_AccessValue get gendata Desired_Iterations],"
    puts $varfile "        \"max_radius_factor\":                  [GiD_AccessValue get gendata Max_Radius_Factor],"
    puts $varfile "        \"min_radius_factor\":                  [GiD_AccessValue get gendata Min_Radius_Factor],"
    puts $varfile "        \"builder\":                            \"[GiD_AccessValue get gendata Builder]\","
    puts $varfile "        \"nonlocal_damage\":                    [GiD_AccessValue get gendata Non-local_Damage],"
    puts $varfile "        \"characteristic_length\":              [GiD_AccessValue get gendata Characteristic_Length],"
    ## linear_solver_settings
    puts $varfile "        \"linear_solver_settings\":             \{"
    if {[GiD_AccessValue get gendata Solver_Type]=="AMGCL"} {
        puts $varfile "            \"solver_type\":                    \"AMGCL\","
        puts $varfile "            \"smoother_type\":                  \"ilu0\","
        puts $varfile "            \"krylov_type\":                    \"gmres\","
        puts $varfile "            \"coarsening_type\":                \"aggregation\","
        puts $varfile "            \"max_iteration\":                  100,"
        puts $varfile "            \"provide_coordinates\":            false,"
        puts $varfile "            \"gmres_krylov_space_dimension\":   100,"
        puts $varfile "            \"verbosity\":                      [GiD_AccessValue get gendata Verbosity],"
        puts $varfile "            \"tolerance\":                      1.0e-6,"
        puts $varfile "            \"scaling\":                        [GiD_AccessValue get gendata Scaling],"
        puts $varfile "            \"block_size\":                     1,"
        puts $varfile "            \"use_block_matrices_if_possible\": true,"
        puts $varfile "            \"coarse_enough\":                  5000"
    } elseif {[GiD_AccessValue get gendata Solver_Type]=="BICGSTABSolver"} {
        puts $varfile "            \"solver_type\":                    \"[GiD_AccessValue get gendata Solver_Type]\","
        puts $varfile "            \"tolerance\":                      1.0e-5,"
        puts $varfile "            \"max_iteration\":                  100,"
        puts $varfile "            \"scaling\":                        [GiD_AccessValue get gendata Scaling],"
        puts $varfile "            \"preconditioner_type\":            \"ILU0Preconditioner\""
    } else {
        puts $varfile "            \"solver_type\":                    \"[GiD_AccessValue get gendata Solver_Type]\","
        puts $varfile "            \"tolerance\":                      1.0e-5,"
        puts $varfile "            \"max_iteration\":                  100"
    }
    puts $varfile "        \},"
    ## problem_domain_sub_model_part_list
    set PutStrings \[
    # Body_Part
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Part
    set Groups [GiD_Info conditions Interface_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"problem_domain_sub_model_part_list\": $PutStrings,"
    ## body_domain_sub_model_part_list
    set PutStrings \[
    set Groups [GiD_Info conditions Body_Part groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"body_domain_sub_model_part_list\":      $PutStrings,"
    ## processes_sub_model_part_list
    set PutStrings \[
    # Solid_Displacement
    set Groups [GiD_Info conditions Solid_Displacement groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Fluid_Pressure
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Force
    set Groups [GiD_Info conditions Force groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Body_Acceleration
    set Groups [GiD_Info conditions Body_Acceleration groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $varfile "        \"processes_sub_model_part_list\":      $PutStrings,"
    ## loads_sub_model_part_list
    set PutStrings \[
    set iGroup 0
    # Force
    set Groups [GiD_Info conditions Force groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    # Body_Acceleration
    set Groups [GiD_Info conditions Body_Acceleration groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        incr iGroup
        append PutStrings \" [lindex [lindex $Groups $i] 1] \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $varfile "        \"loads_sub_model_part_list\":          $PutStrings,"
    ## loads_variable_list
    set PutStrings \[
    # Force
    set Groups [GiD_Info conditions Force groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" FORCE \" ,
    }
    # Face_Load
    set Groups [GiD_Info conditions Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" FACE_LOAD \" ,
    }
    # Normal_Load
    set Groups [GiD_Info conditions Normal_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" NORMAL_CONTACT_STRESS \" ,
    }
    # Normal_Fluid_Flux
    set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" NORMAL_FLUID_FLUX \" ,
    }
    # Interface_Face_Load
    set Groups [GiD_Info conditions Interface_Face_Load groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" FACE_LOAD \" ,
    }
    # Interface_Normal_Fluid_Flux
    set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" NORMAL_FLUID_FLUX \" ,
    }
    # Body_Acceleration
    set Groups [GiD_Info conditions Body_Acceleration groups]
    for {set i 0} {$i < [llength $Groups]} {incr i} {
        append PutStrings \" VOLUME_ACCELERATION \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $varfile "        \"loads_variable_list\":                $PutStrings"
    puts $varfile "    \},"
    
    ## output_configuration
    puts $varfile "    \"output_configuration\": \{"
    puts $varfile "        \"result_file_configuration\": \{"
    puts $varfile "            \"gidpost_flags\":       \{"
    puts $varfile "                \"GiDPostMode\":           \"[GiD_AccessValue get gendata GiD_post_mode]\","
    puts $varfile "                \"WriteDeformedMeshFlag\": \"[GiD_AccessValue get gendata Write_deformed_mesh]\","
    puts $varfile "                \"WriteConditionsFlag\":   \"[GiD_AccessValue get gendata Write_conditions]\","
    if {[GiD_AccessValue get gendata Fracture_Propagation]==true} {
        puts $varfile "                \"MultiFileFlag\":         \"MultipleFiles\""
        puts $varfile "            \},"
        puts $varfile "            \"file_label\":          \"time\","
        puts $varfile "            \"output_control_type\": \"time\","
    } else {
        puts $varfile "                \"MultiFileFlag\":         \"[GiD_AccessValue get gendata Multi_file_flag]\""
        puts $varfile "            \},"
        puts $varfile "            \"file_label\":          \"[GiD_AccessValue get gendata File_label]\","
        puts $varfile "            \"output_control_type\": \"[GiD_AccessValue get gendata Output_control_type]\","
    }
    puts $varfile "            \"output_frequency\":    [GiD_AccessValue get gendata Output_frequency],"
    puts $varfile "            \"body_output\":         [GiD_AccessValue get gendata Body_output],"
    puts $varfile "            \"node_output\":         [GiD_AccessValue get gendata Node_output],"
    puts $varfile "            \"skin_output\":         [GiD_AccessValue get gendata Skin_output],"
    puts $varfile "            \"plane_output\":        \[\],"
    # nodal_results
    set PutStrings \[
    set iGroup 0
    if {[GiD_AccessValue get gendata Write_Solid_Displacement]==true} {
        incr iGroup
        append PutStrings \" DISPLACEMENT \" ,
    }
    if {[GiD_AccessValue get gendata Write_Fluid_Pressure]==true} {
        incr iGroup
        append PutStrings \" WATER_PRESSURE \" ,
    }
    if {[GiD_AccessValue get gendata Write_Reactions]==true} {
        incr iGroup
        append PutStrings \" REACTION \" , \" REACTION_WATER_PRESSURE \" ,
    }
    if {[GiD_AccessValue get gendata Write_Force]==true} {
        incr iGroup
        append PutStrings \" FORCE \" ,
    }
    if {[GiD_AccessValue get gendata Write_Face_Load]==true} {
        incr iGroup
        append PutStrings \" FACE_LOAD \" ,
    }
    if {[GiD_AccessValue get gendata Write_Normal_Load]==true} {
        incr iGroup
        append PutStrings \" NORMAL_CONTACT_STRESS \" ,
    }
    if {[GiD_AccessValue get gendata Write_Tangential_Load]==true} {
        incr iGroup
        append PutStrings \" TANGENTIAL_CONTACT_STRESS \" ,
    }
    if {[GiD_AccessValue get gendata Write_Normal_Fluid_Flux]==true} {
        incr iGroup
        append PutStrings \" NORMAL_FLUID_FLUX \" ,
    }
    if {[GiD_AccessValue get gendata Write_Body_Acceleration]==true} {
        incr iGroup
        append PutStrings \" VOLUME_ACCELERATION \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $varfile "            \"nodal_results\":       $PutStrings,"
    # gauss_point_results
    set PutStrings \[
    set iGroup 0
    if {[GiD_AccessValue get gendata Write_Strain]==true} {
        incr iGroup
        append PutStrings \" GREEN_LAGRANGE_STRAIN_TENSOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Effective_Stress]==true} {
        incr iGroup
        append PutStrings \" CAUCHY_STRESS_TENSOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Total_Stress]==true} {
        incr iGroup
        append PutStrings \" TOTAL_STRESS_TENSOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Von_Mises_Stress]==true} {
        incr iGroup
        append PutStrings \" VON_MISES_STRESS \" ,
    }
    if {[GiD_AccessValue get gendata Write_Fluid_Flux]==true} {
        incr iGroup
        append PutStrings \" FLUID_FLUX_VECTOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Permeability]==true} {
        incr iGroup
        append PutStrings \" PERMEABILITY_MATRIX \" ,
    }
    if {[GiD_AccessValue get gendata Write_Damage]==true} {
        incr iGroup
        append PutStrings \" DAMAGE_VARIABLE \" ,
    }
    if {[GiD_AccessValue get gendata Write_Local_Stress_Vector]==true} {
        incr iGroup
        append PutStrings \" LOCAL_STRESS_VECTOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Local_Relative_Displacement]==true} {
        incr iGroup
        append PutStrings \" LOCAL_RELATIVE_DISPLACEMENT_VECTOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Local_Fluid_Flux]==true} {
        incr iGroup
        append PutStrings \" LOCAL_FLUID_FLUX_VECTOR \" ,
    }
    if {[GiD_AccessValue get gendata Write_Local_Permeability]==true} {
        incr iGroup
        append PutStrings \" LOCAL_PERMEABILITY_MATRIX \" ,
    }
    if {$iGroup > 0} {
        set PutStrings [string trimright $PutStrings ,]
    }
    append PutStrings \]
    puts $varfile "            \"gauss_point_results\": $PutStrings"
    puts $varfile "        \},"
    puts $varfile "        \"point_data_configuration\":  \[\]"
    puts $varfile "    \},"
    
    ## restart_options
    puts $varfile "    \"restart_options\": \{"
    puts $varfile "        \"SaveRestart\":      false,"
    puts $varfile "        \"RestartFrequency\": 0,"
    puts $varfile "        \"LoadRestart\":      false,"
    puts $varfile "        \"Restart_Step\":     0"
    puts $varfile "    \},"

    ## constraints_process_list
    set Groups [GiD_Info conditions Solid_Displacement groups]
    set NumGroups [llength $Groups]
    set Groups [GiD_Info conditions Fluid_Pressure groups]
    incr NumGroups [llength $Groups]
    if {$NumGroups > 0} {
        set iGroup 0
        puts $varfile "    \"constraints_process_list\": \[\{"
        # Solid_Displacement
        set Groups [GiD_Info conditions Solid_Displacement groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] volumes]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"python_module\": \"apply_constraint_vector_table_process\","
                puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":  \"ApplyConstraintVectorTableProcess\","
                puts $varfile "        \"Parameters\":    \{"
                puts $varfile "            \"mesh_id\":         0,"
                puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":   \"DISPLACEMENT\","
                puts $varfile "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 13]\],"
                puts $varfile "            \"is_fixed\":        \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 10],[lindex [lindex $Groups $i] 15]\],"
                puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":           \[0,0,0\]"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
                }
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] surfaces]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"python_module\": \"apply_constraint_vector_table_process\","
                puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":  \"ApplyConstraintVectorTableProcess\","
                puts $varfile "        \"Parameters\":    \{"
                puts $varfile "            \"mesh_id\":         0,"
                puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":   \"DISPLACEMENT\","
                puts $varfile "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 13]\],"
                puts $varfile "            \"is_fixed\":        \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 10],[lindex [lindex $Groups $i] 15]\],"
                puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":           \[0,0,0\]"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
                }
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] lines]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"python_module\": \"apply_constraint_vector_table_process\","
                puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":  \"ApplyConstraintVectorTableProcess\","
                puts $varfile "        \"Parameters\":    \{"
                puts $varfile "            \"mesh_id\":         0,"
                puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":   \"DISPLACEMENT\","
                puts $varfile "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 13]\],"
                puts $varfile "            \"is_fixed\":        \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 10],[lindex [lindex $Groups $i] 15]\],"
                puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":           \[0,0,0\]"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
                }
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] points]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"python_module\": \"apply_constraint_vector_table_process\","
                puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":  \"ApplyConstraintVectorTableProcess\","
                puts $varfile "        \"Parameters\":    \{"
                puts $varfile "            \"mesh_id\":         0,"
                puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":   \"DISPLACEMENT\","
                puts $varfile "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 13]\],"
                puts $varfile "            \"is_fixed\":        \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 10],[lindex [lindex $Groups $i] 15]\],"
                puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 9],[lindex [lindex $Groups $i] 14]\],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":           \[0,0,0\]"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
                }
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        # Fluid_Pressure
        set Groups [GiD_Info conditions Fluid_Pressure groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] volumes]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"python_module\": \"apply_pore_pressure_table_process\","
                puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":  \"ApplyPorePressureTableProcess\","
                puts $varfile "        \"Parameters\":    \{"
                puts $varfile "            \"mesh_id\":              0,"
                puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":        \"WATER_PRESSURE\","
                puts $varfile "            \"is_fixed\":             [lindex [lindex $Groups $i] 8],"
                puts $varfile "            \"value\":                [lindex [lindex $Groups $i] 4],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":                0,"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":                $AuxList,"
                }
                if {[lindex [lindex $Groups $i] 3] == "Hydrostatic"} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"hydrostatic\":          $PutStrings,"
                if {[lindex [lindex $Groups $i] 5] == "Y"} {
                    set PutStrings 2
                } elseif {[lindex [lindex $Groups $i] 5] == "Z"} {
                    set PutStrings 3
                } else {
                    set PutStrings 1
                }
                puts $varfile "            \"gravity_direction\":    $PutStrings,"
                puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] surfaces]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"python_module\": \"apply_pore_pressure_table_process\","
                puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":  \"ApplyPorePressureTableProcess\","
                puts $varfile "        \"Parameters\":    \{"
                puts $varfile "            \"mesh_id\":              0,"
                puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":        \"WATER_PRESSURE\","
                puts $varfile "            \"is_fixed\":             [lindex [lindex $Groups $i] 8],"
                puts $varfile "            \"value\":                [lindex [lindex $Groups $i] 4],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":                0,"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":                $AuxList,"
                }
                if {[lindex [lindex $Groups $i] 3] == "Hydrostatic"} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"hydrostatic\":          $PutStrings,"
                if {[lindex [lindex $Groups $i] 5] == "Y"} {
                    set PutStrings 2
                } elseif {[lindex [lindex $Groups $i] 5] == "Z"} {
                    set PutStrings 3
                } else {
                    set PutStrings 1
                }
                puts $varfile "            \"gravity_direction\":    $PutStrings,"
                puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] lines]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"python_module\": \"apply_pore_pressure_table_process\","
                puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":  \"ApplyPorePressureTableProcess\","
                puts $varfile "        \"Parameters\":    \{"
                puts $varfile "            \"mesh_id\":              0,"
                puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":        \"WATER_PRESSURE\","
                puts $varfile "            \"is_fixed\":             [lindex [lindex $Groups $i] 8],"
                puts $varfile "            \"value\":                [lindex [lindex $Groups $i] 4],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":                0,"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":                $AuxList,"
                }
                if {[lindex [lindex $Groups $i] 3] == "Hydrostatic"} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"hydrostatic\":          $PutStrings,"
                if {[lindex [lindex $Groups $i] 5] == "Y"} {
                    set PutStrings 2
                } elseif {[lindex [lindex $Groups $i] 5] == "Z"} {
                    set PutStrings 3
                } else {
                    set PutStrings 1
                }
                puts $varfile "            \"gravity_direction\":    $PutStrings,"
                puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] points]
            if {[llength $Entities] > 0} {
                incr iGroup
                puts $varfile "        \"python_module\": \"apply_pore_pressure_table_process\","
                puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
                puts $varfile "        \"process_name\":  \"ApplyPorePressureTableProcess\","
                puts $varfile "        \"Parameters\":    \{"
                puts $varfile "            \"mesh_id\":              0,"
                puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
                puts $varfile "            \"variable_name\":        \"WATER_PRESSURE\","
                puts $varfile "            \"is_fixed\":             [lindex [lindex $Groups $i] 8],"
                puts $varfile "            \"value\":                [lindex [lindex $Groups $i] 4],"
                if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                    puts $varfile "            \"table\":                0,"
                } else {
                    set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                    set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                    puts $varfile "            \"table\":                $AuxList,"
                }
                if {[lindex [lindex $Groups $i] 3] == "Hydrostatic"} {
                    set PutStrings true
                } else {
                    set PutStrings false
                }
                puts $varfile "            \"hydrostatic\":          $PutStrings,"
                if {[lindex [lindex $Groups $i] 5] == "Y"} {
                    set PutStrings 2
                } elseif {[lindex [lindex $Groups $i] 5] == "Z"} {
                    set PutStrings 3
                } else {
                    set PutStrings 1
                }
                puts $varfile "            \"gravity_direction\":    $PutStrings,"
                puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 6],"
                puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 7]"
                puts $varfile "        \}"
                if {$iGroup < $NumGroups} {
                    puts $varfile "    \},\{"
                } else {
                    puts $varfile "    \}\],"
                }
            }
        }

    } else {
        puts $varfile "    \"constraints_process_list\": \[\],"
    }
    
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
        puts $varfile "    \"loads_process_list\": \[\{"
        # Force
        set Groups [GiD_Info conditions Force groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"python_module\": \"apply_load_vector_table_process\","
            puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":  \"ApplyLoadVectorTableProcess\","
            puts $varfile "        \"Parameters\":    \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"FORCE\","
            puts $varfile "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 7],[lindex [lindex $Groups $i] 11]\],"
            puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           \[0,0,0\]"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\]"
            }
        }
        # Face_Load
        set Groups [GiD_Info conditions Face_Load groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"python_module\": \"apply_load_vector_table_process\","
            puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":  \"ApplyLoadVectorTableProcess\","
            puts $varfile "        \"Parameters\":    \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"FACE_LOAD\","
            puts $varfile "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 7],[lindex [lindex $Groups $i] 11]\],"
            puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           \[0,0,0\]"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\]"
            }
        }
        # Normal_Load
        set Groups [GiD_Info conditions Normal_Load groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"python_module\": \"apply_normal_load_table_process\","
            puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":  \"ApplyNormalLoadTableProcess\","
            puts $varfile "        \"Parameters\":    \{"
            puts $varfile "            \"mesh_id\":              0,"
            puts $varfile "            \"model_part_name\":      \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":        \"NORMAL_CONTACT_STRESS\","
            puts $varfile "            \"active\":               \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 11]\],"
            puts $varfile "            \"value\":                \[[lindex [lindex $Groups $i] 5],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":                \[0,0\],"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":                \[[lindex $AuxList 0],[lindex $AuxList 1]\],"
            }
            if {[lindex [lindex $Groups $i] 4] == "Hydrostatic"} {
                set PutStrings true
            } else {
                set PutStrings false
            }
            puts $varfile "            \"hydrostatic\":          $PutStrings,"
            if {[lindex [lindex $Groups $i] 6] == "Y"} {
                set PutStrings 2
            } elseif {[lindex [lindex $Groups $i] 6] == "Z"} {
                set PutStrings 3
            } else {
                set PutStrings 1
            }
            puts $varfile "            \"gravity_direction\":    $PutStrings,"
            puts $varfile "            \"reference_coordinate\": [lindex [lindex $Groups $i] 7],"
            puts $varfile "            \"specific_weight\":      [lindex [lindex $Groups $i] 8]"
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\]"
            }
        }
        # Normal_Fluid_Flux
        set Groups [GiD_Info conditions Normal_Fluid_Flux groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"python_module\": \"apply_load_scalar_table_process\","
            puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":  \"ApplyLoadScalarTableProcess\","
            puts $varfile "        \"Parameters\":    \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"NORMAL_FLUID_FLUX\","
            puts $varfile "            \"value\":           [lindex [lindex $Groups $i] 3],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           0"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           $AuxList"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\]"
            }
        }
        # Interface_Face_Load
        set Groups [GiD_Info conditions Interface_Face_Load groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"python_module\": \"apply_load_vector_table_process\","
            puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":  \"ApplyLoadVectorTableProcess\","
            puts $varfile "        \"Parameters\":    \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"FACE_LOAD\","
            puts $varfile "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 7],[lindex [lindex $Groups $i] 11]\],"
            puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           \[0,0,0\]"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\]"
            }
        }
        # Interface_Normal_Fluid_Flux
        set Groups [GiD_Info conditions Interface_Normal_Fluid_Flux groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"python_module\": \"apply_load_scalar_table_process\","
            puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":  \"ApplyLoadScalarTableProcess\","
            puts $varfile "        \"Parameters\":    \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"NORMAL_FLUID_FLUX\","
            puts $varfile "            \"value\":           [lindex [lindex $Groups $i] 3],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           0"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           $AuxList"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\]"
            }
        }
        # Body_Acceleration
        set Groups [GiD_Info conditions Body_Acceleration groups]
        for {set i 0} {$i < [llength $Groups]} {incr i} {
            incr iGroup
            puts $varfile "        \"python_module\": \"apply_load_vector_table_process\","
            puts $varfile "        \"kratos_module\": \"KratosMultiphysics.PoromechanicsApplication\","
            puts $varfile "        \"process_name\":  \"ApplyLoadVectorTableProcess\","
            puts $varfile "        \"Parameters\":    \{"
            puts $varfile "            \"mesh_id\":         0,"
            puts $varfile "            \"model_part_name\": \"[lindex [lindex $Groups $i] 1]\","
            puts $varfile "            \"variable_name\":   \"VOLUME_ACCELERATION\","
            puts $varfile "            \"active\":          \[[lindex [lindex $Groups $i] 3],[lindex [lindex $Groups $i] 7],[lindex [lindex $Groups $i] 11]\],"
            puts $varfile "            \"value\":           \[[lindex [lindex $Groups $i] 4],[lindex [lindex $Groups $i] 8],[lindex [lindex $Groups $i] 12]\],"
            if {[GiD_AccessValue get gendata Strategy_Type]=="Arc-Length"} {
                puts $varfile "            \"table\":           \[0,0,0\]"
            } else {
                set SearchInList [lsearch $TableList [lindex [lindex $Groups $i] 1]*]
                set AuxList [lindex $TableList [expr { $SearchInList+1 }]]
                puts $varfile "            \"table\":           \[[lindex $AuxList 0],[lindex $AuxList 1],[lindex $AuxList 2]\]"
            }
            puts $varfile "        \}"
            if {$iGroup < $NumGroups} {
                puts $varfile "    \},\{"
            } else {
                puts $varfile "    \}\]"
            }
        }
    } else {
        puts $varfile "    \"loads_process_list\":       \[\]"
    }

    puts $varfile "\}"
    
    close $varfile
}
