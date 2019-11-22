! PROGRAM dummy_solver_fortran
!     USE co_sim_io

!     IMPLICIT NONE

!     INTEGER :: num_nodes, num_elems
!     INTEGER, DIMENSION(4) :: node_ids, elems
!     INTEGER, DIMENSION(1) :: num_nodes_per_elem
!     REAL(8), DIMENSION(12) :: nodes, data_send, data_recv
!     ! CALL EMPIRE_API_Connect('fortranInput.xml')

!     ! num_nodes = 4
!     ! num_elems = 1
!     ! node_ids  = (/ 1, 2, 3, 4 /)
!     ! elems     = (/ 1, 2, 3, 4 /)
!     ! num_nodes_per_elem = (/ 4 /)
!     ! nodes     = (/ 0., 0., 0., &
!     !                0., 1., 0., &
!     !                1., 1., 0., &
!     !                1., 0., 0.  /)
!     ! data_recv = 0D0
!     ! data_send = 1D0

!     ! CALL EMPIRE_API_sendMesh(" ",num_nodes, num_elems, nodes, node_ids, &
!     !                          num_nodes_per_elem, elems)

!     ! CALL EMPIRE_API_sendDataField(" ",12, data_send)
!     ! CALL EMPIRE_API_recvDataField(" ",12, data_recv)

!     ! write(*,*) 'data_recv: ', data_recv

!     ! CALL EMPIRE_API_Disconnect()

! END PROGRAM dummy_solver_fortran


program dummy_solver_fortran
    use co_sim_io

    ! call cosim_register (C_FUNLOC(advance_in_time),          "AdvanceInTime"//c_null_char)
    ! call cosim_register (C_FUNLOC(initialize_solution_step), "InitializeSolutionStep"//c_null_char)
    ! call cosim_register (C_FUNLOC(solve_solution_step),      "SolveSolutionStep"//c_null_char)
    ! call cosim_register (C_FUNLOC(finalize_solution_step),   "FinalizeSolutionStep"//c_null_char)

    ! call cosim_run()







    INTEGER :: size
    REAL(8), DIMENSION(:), pointer :: f_ptr_data_int
    INTEGER, DIMENSION(:), pointer :: f_ptr_data_real
    type (c_ptr) :: c_ptr_data


    ! write(*,*) "F BEFORE calling c_sub"
    ! call c_sub(dim, W_ptr)
    ! convert c pointer to fortran pointer

    call CoSimIO_ImportData("FortranSolver"//c_null_char, "pressure"//c_null_char, size, c_ptr_data)
    call c_f_pointer(c_ptr_data, f_ptr_data_real, [size])

    call FreeCMemory(c_ptr_data)

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