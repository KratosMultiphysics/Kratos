MODULE co_sim_io
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

    ! Define interface of C function.
    INTERFACE
        SUBROUTINE co_sim_run () BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
        END SUBROUTINE co_sim_run

        SUBROUTINE co_sim_register (func, name) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            TYPE(C_FUNPTR), INTENT(IN), VALUE :: func
            CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: name
        END SUBROUTINE co_sim_register

    END INTERFACE

END MODULE co_sim_io