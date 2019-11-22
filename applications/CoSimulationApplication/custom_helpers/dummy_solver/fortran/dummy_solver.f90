program dummy_solver_fortran
    use co_sim_io

    INTEGER :: size
    REAL(8), DIMENSION(:), pointer :: f_ptr_data_int
    INTEGER, DIMENSION(:), pointer :: f_ptr_data_real
    type (c_ptr) :: c_ptr_data

    call CoSimIO_Connect("FortranSolver"//c_null_char, "FortranSolverInputFile"//c_null_char)

    call CoSimIO_ImportData("FortranSolver"//c_null_char, "pressure"//c_null_char, size, c_ptr_data)
    call CoSimIO_ExportData("FortranSolver"//c_null_char, "temperature"//c_null_char, size, c_ptr_data)


    ! INTEGER :: size
    ! REAL(8), DIMENSION(:), pointer :: f_ptr_data_int
    ! INTEGER, DIMENSION(:), pointer :: f_ptr_data_real
    ! type (c_ptr) :: c_ptr_data
    ! call CoSimIO_ExportMesh("FortranSolver"//c_null_char, "temperature"//c_null_char, size, c_ptr_data)

    call c_f_pointer(c_ptr_data, f_ptr_data_real, [size])

    call FreeCMemory(c_ptr_data)

    call CoSimIO_RegisterAdvanceInTime("FortranSolver"//c_null_char, C_FUNLOC(advance_in_time))
    call CoSimIO_RegisterSolvingFunction("FortranSolver"//c_null_char,&
        "InitializeSolutionStep"//c_null_char, C_FUNLOC(initialize_solution_step))
    call CoSimIO_RegisterSolvingFunction("FortranSolver"//c_null_char,&
        "SolveSolutionStep"//c_null_char, C_FUNLOC(solve_solution_step))
    call CoSimIO_RegisterSolvingFunction("FortranSolver"//c_null_char,&
        "FinalizeSolutionStep"//c_null_char, C_FUNLOC(finalize_solution_step))

    call CoSimIO_Run("FortranSolver"//c_null_char)
    call CoSimIO_Disconnect("FortranSolver"//c_null_char)

    CONTAINS

    SUBROUTINE advance_in_time (time)
        REAL, intent(inout):: time

        time = time + 0.1
        write(*,*) "AdvanceInTime"
    END SUBROUTINE advance_in_time

    SUBROUTINE initialize_solution_step ()
        write(*,*) "InitializeSolutionStep"
    END SUBROUTINE initialize_solution_step

    SUBROUTINE solve_solution_step ()
        write(*,*) "SolveSolutionStep"
    END SUBROUTINE solve_solution_step

    SUBROUTINE finalize_solution_step ()
        write(*,*) "FinalizeSolutionStep"
    END SUBROUTINE finalize_solution_step
end program