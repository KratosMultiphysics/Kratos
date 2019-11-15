PROGRAM dummy_solver_fortran
    USE co_sim_io

    IMPLICIT NONE

    INTEGER :: num_nodes, num_elems
    INTEGER, DIMENSION(4) :: node_ids, elems
    INTEGER, DIMENSION(1) :: num_nodes_per_elem
    REAL(8), DIMENSION(12) :: nodes, data_send, data_recv
    ! CALL EMPIRE_API_Connect('fortranInput.xml')

    ! num_nodes = 4
    ! num_elems = 1
    ! node_ids  = (/ 1, 2, 3, 4 /)
    ! elems     = (/ 1, 2, 3, 4 /)
    ! num_nodes_per_elem = (/ 4 /)
    ! nodes     = (/ 0., 0., 0., &
    !                0., 1., 0., &
    !                1., 1., 0., &
    !                1., 0., 0.  /)
    ! data_recv = 0D0
    ! data_send = 1D0

    ! CALL EMPIRE_API_sendMesh(" ",num_nodes, num_elems, nodes, node_ids, &
    !                          num_nodes_per_elem, elems)

    ! CALL EMPIRE_API_sendDataField(" ",12, data_send)
    ! CALL EMPIRE_API_recvDataField(" ",12, data_recv)

    ! write(*,*) 'data_recv: ', data_recv

    ! CALL EMPIRE_API_Disconnect()

  END PROGRAM dummy_solver_fortran
