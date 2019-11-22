MODULE solver_globals
    IMPLICIT NONE
    INTEGER :: num_nodes
    INTEGER :: echo_level
    REAL :: delta_time = 0.1
    REAL :: end_time = 0.5

END MODULE solver_globals

program dummy_solver_fortran
    use co_sim_io
    use solver_globals, only: num_nodes, echo_level, delta_time, end_time
    implicit none

    INTEGER :: size
    REAL(8), DIMENSION(:), pointer :: f_ptr_data_int
    INTEGER, DIMENSION(:), pointer :: f_ptr_data_real
    type (c_ptr) :: c_ptr_data


    CHARACTER(100) :: input_char
    INTEGER :: level_of_coupling

    ! ! INTEGER :: size
    ! ! REAL(8), DIMENSION(:), pointer :: f_ptr_data_int
    ! ! INTEGER, DIMENSION(:), pointer :: f_ptr_data_real
    ! ! type (c_ptr) :: c_ptr_data
    ! ! call CoSimIO_ExportMesh("FortranSolver"//c_null_char, "temperature"//c_null_char, size, c_ptr_data)

    ! call c_f_pointer(c_ptr_data, f_ptr_data_real, [size])

    ! call FreeCMemory(c_ptr_data)


    ! Make sure the right number of command line arguments have been provided
    IF(COMMAND_ARGUMENT_COUNT().NE.3)THEN
      WRITE(*,*)'ERROR, THREE COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
      STOP
    ENDIF

    ! Read the command line arguments
    CALL GET_COMMAND_ARGUMENT(1,input_char)
    READ(input_char,*)level_of_coupling
    CALL GET_COMMAND_ARGUMENT(2,input_char)
    READ(input_char,*)num_nodes
    CALL GET_COMMAND_ARGUMENT(3,input_char)
    READ(input_char,*)echo_level

    IF (level_of_coupling == 0) THEN
        call run_solution_loop()
    ELSE IF (level_of_coupling == 1) THEN
        call run_solution_loop_weak_coupling()
    ELSE IF (level_of_coupling == 2) THEN
        call run_solution_loop_strong_coupling()
    ELSE IF (level_of_coupling == 3) THEN
        call run_solution_co_simulation_orchestrated()
    ELSE
        WRITE(*,*)'ERROR, WRONG LEVEL OF COUPLING; CAN ONLY BE 0, 1, 2, 3, STOPPING'
        STOP
    END IF

    CONTAINS

    subroutine solver_print(string_to_print)
        character(*), INTENT(IN) :: string_to_print
        IF (echo_level > 0) THEN
            write(*,*) "Solver [F]: ", string_to_print
        End IF
    end subroutine solver_print

    !!! solving-functions !!!
    SUBROUTINE advance_in_time (time)
        REAL, intent(inout):: time
        character(len=1024) :: print_string
        time = time + delta_time
        write(print_string,*) "AdvanceInTime; new time: ", time
        call solver_print(adjustl(trim(print_string)))
    END SUBROUTINE advance_in_time

    SUBROUTINE initialize_solution_step ()
        call solver_print("InitializeSolutionStep")
    END SUBROUTINE initialize_solution_step

    SUBROUTINE solve_solution_step ()
        call solver_print("SolveSolutionStep")
    END SUBROUTINE solve_solution_step

    SUBROUTINE finalize_solution_step ()
        call solver_print("FinalizeSolutionStep")
    END SUBROUTINE finalize_solution_step

    !!! data-exchange functions !!!
    SUBROUTINE import_data_from_co_sim (connection_name, identifier)
        character(*), INTENT(IN) :: connection_name
        character(*), INTENT(IN) :: identifier
        call solver_print("    Before Importing Data from CoSim")
        ! call CoSimIO_ImportData(connection_name//c_null_char, identifier//c_null_char, size, c_ptr_data)
        call solver_print("    After Importing Data from CoSim")
    END SUBROUTINE import_data_from_co_sim

    SUBROUTINE export_data_to_co_sim (connection_name, identifier)
        character(*), INTENT(IN) :: connection_name
        character(*), INTENT(IN) :: identifier
        call solver_print("    Before Exporting Data to CoSim")
        ! call CoSimIO_ExportData(connection_name//c_null_char, identifier//c_null_char, size, c_ptr_data)
        call solver_print("    After Exporting Data to CoSim")
    END SUBROUTINE export_data_to_co_sim

    SUBROUTINE import_mesh_from_co_sim (connection_name, identifier)
        character(*), INTENT(IN) :: connection_name
        character(*), INTENT(IN) :: identifier
        call solver_print("    Before Importing Mesh from CoSim")
        ! call CoSimIO_ImportMesh(connection_name//c_null_char, identifier//c_null_char, ...)
        call solver_print("    After Importing Mesh from CoSim")
    END SUBROUTINE import_mesh_from_co_sim

    SUBROUTINE export_mesh_to_co_sim (connection_name, identifier)
        character(*), INTENT(IN) :: connection_name
        character(*), INTENT(IN) :: identifier
        call solver_print("    Before Exporting Mesh to CoSim")
        ! call CoSimIO_ExportMesh(connection_name//c_null_char, identifier//c_null_char, ...)
        call solver_print("    After Exporting Mesh to CoSim")
    END SUBROUTINE export_mesh_to_co_sim


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Functions for executing the different versions of CoSimulation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE run_solution_loop ()
        real :: current_time = 0.0

        call solver_print("Doing STANDALONE simulation")

        do while (current_time < end_time)
            write(*,*)
            call advance_in_time(current_time)
            call initialize_solution_step()
            call solve_solution_step()
            call finalize_solution_step()
        end do

        call solver_print("Exiting")
    END SUBROUTINE run_solution_loop


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE run_solution_loop_weak_coupling ()
        real :: current_time = 0.0

        call solver_print("Doing COUPLED simulation (weakly coupled)")

        call CoSimIO_Connect("FortranSolver"//c_null_char, "FortranSolverInputFile"//c_null_char)

        call import_mesh_from_co_sim("FortranSolver", "interface_1")
        call export_mesh_to_co_sim("FortranSolver", "interface_2")

        do while (current_time < end_time)
            write(*,*)
            call advance_in_time(current_time)
            call initialize_solution_step()
            call import_data_from_co_sim("FortranSolver", "pressure")
            call solve_solution_step()
            call export_data_to_co_sim("FortranSolver", "temperature")
            call finalize_solution_step()
        end do

        call CoSimIO_Disconnect("FortranSolver"//c_null_char)

        call solver_print("Exiting")
    END SUBROUTINE run_solution_loop_weak_coupling


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE run_solution_loop_strong_coupling ()
        real :: current_time = 0.0
        integer :: convergence_signal = 0

        call solver_print("Doing COUPLED simulation (strongly coupled)")

        call CoSimIO_Connect("FortranSolver"//c_null_char, "FortranSolverInputFile"//c_null_char)

        call import_mesh_from_co_sim("FortranSolver", "interface_1")
        call export_mesh_to_co_sim("FortranSolver", "interface_2")

        do while (current_time < end_time)
            write(*,*)
            call advance_in_time(current_time)
            call initialize_solution_step()
            do ! convergence loop, execute until convergence is achieved
                call import_data_from_co_sim("FortranSolver", "pressure")
                call solve_solution_step()
                call export_data_to_co_sim("FortranSolver", "temperature")
                call CoSimIO_IsConverged("FortranSolver"//c_null_char, convergence_signal)
                if(convergence_signal.NE.0) exit
            end do

            call finalize_solution_step()

        end do

        call CoSimIO_Disconnect("FortranSolver"//c_null_char)

        call solver_print("Exiting")
    END SUBROUTINE run_solution_loop_strong_coupling


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE run_solution_co_simulation_orchestrated ()
        call solver_print("Doing COUPLED simulation (orchestrated by CoSimulation)")

        call CoSimIO_Connect("FortranSolver"//c_null_char, "FortranSolverInputFile"//c_null_char)

        call CoSimIO_RegisterAdvanceInTime("FortranSolver"//c_null_char, C_FUNLOC(advance_in_time))
        call CoSimIO_RegisterSolvingFunction("FortranSolver"//c_null_char,&
            "InitializeSolutionStep"//c_null_char, C_FUNLOC(initialize_solution_step))
        call CoSimIO_RegisterSolvingFunction("FortranSolver"//c_null_char,&
            "SolveSolutionStep"//c_null_char, C_FUNLOC(solve_solution_step))

        call CoSimIO_RegisterDataExchangeFunction("FortranSolver"//c_null_char,&
            "ImportData"//c_null_char, C_FUNLOC(import_data_from_co_sim))
        call CoSimIO_RegisterDataExchangeFunction("FortranSolver"//c_null_char,&
            "ExportData"//c_null_char, C_FUNLOC(export_data_to_co_sim))

        call CoSimIO_RegisterDataExchangeFunction("FortranSolver"//c_null_char,&
            "ImportMesh"//c_null_char, C_FUNLOC(import_mesh_from_co_sim))
        call CoSimIO_RegisterDataExchangeFunction("FortranSolver"//c_null_char,&
            "ExportMesh"//c_null_char, C_FUNLOC(export_mesh_to_co_sim))

        call CoSimIO_Run("FortranSolver"//c_null_char)

        call CoSimIO_Disconnect("FortranSolver"//c_null_char)

        call solver_print("Exiting")
    END SUBROUTINE run_solution_co_simulation_orchestrated

end program