MODULE co_sim_io
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

    ! Define interface of C function.
    INTERFACE

        SUBROUTINE CoSimIO_Connect (ConnectionName, SettingsFileName) BIND(C, NAME="CoSimIO_Connect")
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: SettingsFileName
        END SUBROUTINE CoSimIO_Connect

        SUBROUTINE CoSimIO_Disconnect (ConnectionName) BIND(C, NAME="CoSimIO_Disconnect")
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
        END SUBROUTINE CoSimIO_Disconnect


        SUBROUTINE CoSimIO_ImportData(ConnectionName, Identifier, Size, Data) BIND (C, NAME="CoSimIO_ImportData")
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: Identifier
            INTEGER(C_INT), INTENT(INOUT)  :: Size
            type(C_PTR) :: Data
        END SUBROUTINE CoSimIO_ImportData

        SUBROUTINE CoSimIO_ExportData(ConnectionName, Identifier, Size, Data) BIND (C, NAME="CoSimIO_ExportData")
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: Identifier
            INTEGER(C_INT), INTENT(IN)  :: Size
            type(C_PTR) :: Data
        END SUBROUTINE CoSimIO_ExportData


        SUBROUTINE CoSimIO_ImportMesh(ConnectionName, Identifier, NumberOfNodes, &
                                      NumberOfElements, NodalCoordinates, &
                                      ElementConnectivities, ElementTypes&
                                      ) BIND (C, NAME="CoSimIO_ImportMesh")
        USE, INTRINSIC :: ISO_C_BINDING
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
        CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: Identifier
        INTEGER(C_INT), INTENT(INOUT) :: NumberOfNodes
        INTEGER(C_INT), INTENT(INOUT) :: NumberOfElements
        type(C_PTR) :: NodalCoordinates
        type(C_PTR) :: ElementConnectivities
        type(C_PTR) :: ElementTypes
        END SUBROUTINE CoSimIO_ImportMesh

        SUBROUTINE CoSimIO_ExportMesh(ConnectionName, Identifier, NumberOfNodes, &
                                      NumberOfElements, NodalCoordinates, &
                                      ElementConnectivities, ElementTypes&
                                      ) BIND (C, NAME="CoSimIO_ExportMesh")
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: Identifier
            INTEGER(C_INT), INTENT(IN) :: NumberOfNodes
            INTEGER(C_INT), INTENT(IN) :: NumberOfElements
            type(C_PTR) :: NodalCoordinates
            type(C_PTR) :: ElementConnectivities
            type(C_PTR) :: ElementTypes
        END SUBROUTINE CoSimIO_ExportMesh


        SUBROUTINE CoSimIO_Run (ConnectionName) BIND(C, NAME="CoSimIO_Run")
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
        END SUBROUTINE CoSimIO_Run

        SUBROUTINE CoSimIO_Register (ConnectionName, FunctionName, FunctionPointer) BIND(C, NAME="CoSimIO_Register")
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
            CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: FunctionName
            TYPE(C_FUNPTR), INTENT(IN), VALUE :: FunctionPointer
        END SUBROUTINE CoSimIO_Register

        SUBROUTINE CoSimIO_IsConverged (ConnectionName) BIND(C, NAME="CoSimIO_IsConverged")
            USE, INTRINSIC :: ISO_C_BINDING
            CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: ConnectionName
        END SUBROUTINE CoSimIO_IsConverged

        SUBROUTINE FreeCMemory(Data) BIND (C,NAME='_FreeMemory')
            USE, INTRINSIC :: ISO_C_BINDING
            type(C_PTR) :: Data
        END SUBROUTINE FreeCMemory

    END INTERFACE

END MODULE co_sim_io