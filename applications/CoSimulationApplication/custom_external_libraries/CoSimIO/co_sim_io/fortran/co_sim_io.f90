!     ______     _____ _           ________
!    / ____/___ / ___/(_)___ ___  /  _/ __ |
!   / /   / __ \\__ \/ / __ `__ \ / // / / /
!  / /___/ /_/ /__/ / / / / / / // // /_/ /
!  \____/\____/____/_/_/ /_/ /_/___/\____/
!  Kratos CoSimulationApplication
!
!  License:         BSD License, see license.txt
!
!  Main authors:    Philipp Bucher (https://github.com/philbucher), Ankit Nanda(https://github.com/Ankit-Nanda-altair)
!

      MODULE co_sim_io
          USE, INTRINSIC :: ISO_C_BINDING
          IMPLICIT NONE
          
          !START: co_sim_io_c_info.h structures
          TYPE, BIND (C) ::  CoSimIO_Info
              TYPE(c_ptr)      :: PtrCppInfo              
          END TYPE CoSimIO_Info
          
          !END: co_sim_io_c_info.h structures
          
          !START: co_sim_io_c_model_part.h structures
          TYPE, BIND (C) ::  CoSimIO_Node
              TYPE(c_ptr)      :: PtrCppNode              
          END TYPE CoSimIO_Node
          
          TYPE, BIND (C) ::  CoSimIO_Element
              TYPE(c_ptr)      :: PtrCppElement              
          END TYPE CoSimIO_Element
          
          TYPE, BIND (C) ::  CoSimIO_ModelPart
              TYPE(c_ptr)      :: PtrCppModelPart              
          END TYPE CoSimIO_ModelPart
          
          ENUM, BIND(C) 
                ENUMERATOR:: CoSimIO_ConnectionStatus = -1
                ENUMERATOR:: CoSimIO_NotConnected
                ENUMERATOR:: CoSimIO_Connected
                ENUMERATOR:: CoSimIO_Disconnected
                ENUMERATOR:: CoSimIO_ConnectionError
                ENUMERATOR:: CoSimIO_DisconnectionError
          END ENUM

          ENUM, BIND(C)  
              ENUMERATOR:: CoSimIO_ElementType = -1
              ENUMERATOR:: CoSimIO_Hexahedra3D20
              ENUMERATOR:: CoSimIO_Hexahedra3D27
              ENUMERATOR:: CoSimIO_Hexahedra3D8
              ENUMERATOR:: CoSimIO_Prism3D15
              ENUMERATOR:: CoSimIO_Prism3D6
              ENUMERATOR:: CoSimIO_Pyramid3D13
              ENUMERATOR:: CoSimIO_Pyramid3D5
              ENUMERATOR:: CoSimIO_Quadrilateral2D4
              ENUMERATOR:: CoSimIO_Quadrilateral2D8
              ENUMERATOR:: CoSimIO_Quadrilateral2D9
              ENUMERATOR:: CoSimIO_Quadrilateral3D4
              ENUMERATOR:: CoSimIO_Quadrilateral3D8
              ENUMERATOR:: CoSimIO_Quadrilateral3D9
              ENUMERATOR:: CoSimIO_Tetrahedra3D10
              ENUMERATOR:: CoSimIO_Tetrahedra3D4
              ENUMERATOR:: CoSimIO_Triangle2D3
              ENUMERATOR:: CoSimIO_Triangle2D6
              ENUMERATOR:: CoSimIO_Triangle3D3
              ENUMERATOR:: CoSimIO_Triangle3D6
              ENUMERATOR:: CoSimIO_Line2D2
              ENUMERATOR:: CoSimIO_Line2D3
              ENUMERATOR:: CoSimIO_Line3D2
              ENUMERATOR:: CoSimIO_Line3D3
              ENUMERATOR:: CoSimIO_Point2D
              ENUMERATOR:: CoSimIO_Point3D
          END ENUM
          !END: co_sim_io_c_model_part.h structures
      
          INTERFACE
              
          ! START: co_sim_io_c.h interface functions and subroutines
              TYPE(CoSimIO_Info) FUNCTION CoSimIO_Hello() BIND(C, NAME="CoSimIO_Hello")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
              END FUNCTION CoSimIO_Hello
      
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_Connect (I_Settings) BIND(C, NAME="CoSimIO_Connect")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE      :: I_Settings
              END FUNCTION CoSimIO_Connect

              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_Disconnect (I_Info) BIND(C, NAME="CoSimIO_Disconnect")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE      :: I_Info
              END FUNCTION CoSimIO_Disconnect

              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_ImportData (I_Info, O_Size, O_Data) BIND(C, NAME="CoSimIO_ImportData")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE                                     :: I_Info
                  INTEGER(KIND=c_int),pointer                                   :: O_Size
                  REAL(KIND=c_double),pointer,DIMENSION(:)                      :: O_Data
              END FUNCTION CoSimIO_ImportData      
      
      
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_ImportData_fortran (I_Info, O_Size, O_Data) 
     &                                                        BIND(C, NAME="CoSimIO_ImportData_fortran")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE                                     :: I_Info
                  INTEGER(KIND=c_int),VALUE                                   :: O_Size
                  REAL(KIND=c_double),pointer,DIMENSION(:)                      :: O_Data
              END FUNCTION CoSimIO_ImportData_fortran    
          
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_ExportData (I_Info, I_Size, I_Data) BIND(C, NAME="CoSimIO_ExportData")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE         :: I_Info
                  INTEGER(KIND = c_int), VALUE      :: I_Size
                  REAL(KIND=c_double),DIMENSION(*)                :: I_Data
              END FUNCTION CoSimIO_ExportData              
          
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_ImportMesh (I_Info, O_ModelPart) BIND(C, NAME="CoSimIO_ImportMesh")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart, CoSimIO_Info 
                  TYPE(CoSimIO_Info), VALUE             :: I_Info
                  TYPE(CoSimIO_ModelPart), VALUE        :: O_ModelPart
              END FUNCTION CoSimIO_ImportMesh             
          
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_ExportMesh (I_Info, I_ModelPart) BIND(C, NAME="CoSimIO_ExportMesh")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart, CoSimIO_Info 
                  TYPE(CoSimIO_Info), VALUE             :: I_Info
                  TYPE(CoSimIO_ModelPart), VALUE        :: I_ModelPart
              END FUNCTION CoSimIO_ExportMesh
      
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_ImportInfo (I_Info) BIND(C, NAME="CoSimIO_ImportInfo")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE      :: I_Info
              END FUNCTION CoSimIO_ImportInfo
      
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_ExportInfo (I_Info) BIND(C, NAME="CoSimIO_ExportInfo")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE      :: I_Info
              END FUNCTION CoSimIO_ExportInfo
        
               INTEGER(KIND=C_INT) FUNCTION CoSimIO_Register (I_Info, I_FunctionPointer) BIND(C, NAME="CoSimIO_Register")
                   USE, INTRINSIC :: ISO_C_BINDING
                   IMPORT CoSimIO_Info
                   TYPE(CoSimIO_Info), VALUE    :: I_Info
                   TYPE(C_FUNPTR), VALUE        :: I_FunctionPointer
               END FUNCTION CoSimIO_Register
      
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_Run (I_Info) BIND(C, NAME="CoSimIO_Run")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE      :: I_Info
              END FUNCTION CoSimIO_Run
      
      
              TYPE(c_ptr) FUNCTION CoSimIO_Malloc (size)
                  USE, INTRINSIC :: ISO_C_BINDING
                  REAL(KIND=c_double), VALUE :: size
              END FUNCTION CoSimIO_Malloc
            
            
              SUBROUTINE CoSimIO_Free(ptr_h) BIND(C, NAME='CoSimIO_Free')
                  USE, INTRINSIC :: ISO_C_BINDING
                  !IMPORT CoSimIO_Info
                  TYPE(c_ptr),INTENT(IN), VALUE :: ptr_h              
              END SUBROUTINE CoSimIO_Free
      ! END: co_sim_io_c.h interface
      
      !START: co_sim_io_c_info.h interface
              TYPE(CoSimIO_Info) FUNCTION CoSimIO_CreateInfo() BIND(C, NAME="CoSimIO_CreateInfo")  
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
              END FUNCTION CoSimIO_CreateInfo
      
              TYPE(CoSimIO_Info)  FUNCTION CoSimIO_CopyInfo (I_Info) BIND(C, NAME="CoSimIO_CopyInfo")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE      :: I_Info
              END FUNCTION CoSimIO_CopyInfo
      
              INTEGER(KIND=c_int)  FUNCTION CoSimIO_FreeInfo (I_Info) BIND(C, NAME="CoSimIO_FreeInfo")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE      :: I_Info
              END FUNCTION CoSimIO_FreeInfo          
            
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_Info_Has (I_Info, I_Key) BIND(C, NAME="CoSimIO_Info_Has")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE           :: I_Info
                  CHARACTER(len=1, KIND=c_char), DIMENSION(*)       :: I_Key
              END FUNCTION CoSimIO_Info_Has
            
            
              SUBROUTINE CoSimIO_Info_Erase(I_Info, I_Key) BIND(C, NAME='CoSimIO_Info_Erase')
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), INTENT(IN), VALUE           :: I_Info
                  CHARACTER(len=1, KIND=c_char), DIMENSION(*) , INTENT(IN)      :: I_Key            
              END SUBROUTINE CoSimIO_Info_Erase
            
            
              SUBROUTINE CoSimIO_Info_Clear(I_Info) BIND(C, NAME='CoSimIO_Info_Clear')
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), INTENT(IN) , VALUE          :: I_Info          
              END SUBROUTINE CoSimIO_Info_Clear
      
              INTEGER(KIND=c_int)  FUNCTION CoSimIO_Info_Size (I_Info) BIND(C, NAME="CoSimIO_Info_Size")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE      :: I_Info
              END FUNCTION CoSimIO_Info_Size            
            
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_Info_GetInt (I_Info, I_Key) BIND(C, NAME="CoSimIO_Info_GetInt")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE                                :: I_Info
                  CHARACTER(len=1, KIND=c_char), DIMENSION(*)              :: I_Key
              END FUNCTION CoSimIO_Info_GetInt            
            
      
              REAL(KIND=c_double) FUNCTION CoSimIO_Info_GetDouble (I_Info, I_Key) BIND(C, NAME="CoSimIO_Info_GetDouble")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE                               :: I_Info
                  CHARACTER(len=1, KIND=c_char) , DIMENSION(*)            :: I_Key
              END FUNCTION CoSimIO_Info_GetDouble           
            
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_Info_GetBool (I_Info, I_Key) BIND(C, NAME="CoSimIO_Info_GetBool")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE                                :: I_Info
                  CHARACTER(len=1, KIND=c_char), DIMENSION(*)              :: I_Key
              END FUNCTION CoSimIO_Info_GetBool            
            
      
              type(c_ptr) FUNCTION CoSimIO_Info_GetString (I_Info, I_Key) BIND(C, NAME="CoSimIO_Info_GetString")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), VALUE                                   :: I_Info
                  CHARACTER(LEN = 1, KIND=c_char), DIMENSION(*)               :: I_Key
                  
              END FUNCTION CoSimIO_Info_GetString           
            
      
              TYPE(CoSimIO_Info) FUNCTION CoSimIO_Info_GetInfo (I_Info, I_Key) BIND(C, NAME="CoSimIO_Info_GetInfo")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info),VALUE                                  :: I_Info
                  CHARACTER(len=1, KIND=c_char) , DIMENSION(*)              :: I_Key
              END FUNCTION CoSimIO_Info_GetInfo
            
            
              SUBROUTINE CoSimIO_Info_SetInt(I_Info, I_Key, I_Value) BIND(C, NAME='CoSimIO_Info_SetInt')
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), INTENT(IN), VALUE                             :: I_Info
                  CHARACTER(len=1, KIND=c_char), DIMENSION(*) , INTENT(IN)          :: I_Key 
                  INTEGER(KIND = c_int), INTENT(IN), VALUE                          :: I_Value
              END SUBROUTINE CoSimIO_Info_SetInt
            
            
              SUBROUTINE CoSimIO_Info_SetDouble(I_Info, I_Key, I_Value) BIND(C, NAME='CoSimIO_Info_SetDouble')
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), INTENT(IN), VALUE                     :: I_Info
                  CHARACTER(len=1, KIND=c_char), DIMENSION(*) , INTENT(IN)                :: I_Key 
                  REAL(KIND = c_double), INTENT(IN), VALUE                  :: I_Value
              END SUBROUTINE CoSimIO_Info_SetDouble
            
            
              SUBROUTINE CoSimIO_Info_SetBool(I_Info, I_Key, I_Value) BIND(C, NAME='CoSimIO_Info_SetBool')
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), INTENT(IN), VALUE           :: I_Info
                  CHARACTER(len=1, KIND=c_char) , DIMENSION(*), INTENT(IN)      :: I_Key 
                  INTEGER(KIND = c_int), INTENT(IN), VALUE        :: I_Value
              END SUBROUTINE CoSimIO_Info_SetBool
            
            
              SUBROUTINE CoSimIO_Info_SetString(I_Info, I_Key, I_Value) BIND(C, NAME='CoSimIO_Info_SetString')
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), INTENT(IN),VALUE              :: I_Info
                  CHARACTER(KIND=c_char), INTENT(IN)       :: I_Key(*)
                  CHARACTER( KIND=c_char), INTENT(IN)       :: I_Value(20)
              END SUBROUTINE CoSimIO_Info_SetString
            
            
              SUBROUTINE CoSimIO_Info_SetInfo(I_Info, I_Key, I_Value) BIND(C, NAME='CoSimIO_Info_SetInfo')
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Info
                  TYPE(CoSimIO_Info), INTENT(IN), VALUE           :: I_Info
                  CHARACTER(len=1, KIND=c_char), DIMENSION(*) , INTENT(IN)      :: I_Key 
                  TYPE(CoSimIO_Info), INTENT(IN), VALUE           :: I_Value
              END SUBROUTINE CoSimIO_Info_SetInfo
              
      ! END: co_sim_io_c_info.h interface  
      
      ! START: co_sim_io_c_model_part.h interface          
            ! NODE FUNCTIONS
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_Node_Id (I_Node) BIND(C, NAME="CoSimIO_Node_Id")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node
                  TYPE(CoSimIO_Node), VALUE                                :: I_Node
              END FUNCTION CoSimIO_Node_Id

              REAL(KIND=c_double) FUNCTION CoSimIO_Node_X (I_Node) BIND(C, NAME="CoSimIO_Node_X")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node
                  TYPE(CoSimIO_Node), VALUE                                :: I_Node
              END FUNCTION CoSimIO_Node_X

              REAL(KIND=c_double) FUNCTION CoSimIO_Node_Y (I_Node) BIND(C, NAME="CoSimIO_Node_Y")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node
                  TYPE(CoSimIO_Node), VALUE                                :: I_Node
              END FUNCTION CoSimIO_Node_Y

              REAL(KIND=c_double) FUNCTION CoSimIO_Node_Z (I_Node) BIND(C, NAME="CoSimIO_Node_Z")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node
                  TYPE(CoSimIO_Node), VALUE                                :: I_Node
              END FUNCTION CoSimIO_Node_Z

              REAL(KIND=c_double) FUNCTION CoSimIO_Node_Coordinate (I_Node, I_Index) BIND(C, NAME="CoSimIO_Node_Coordinate")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node
                  TYPE(CoSimIO_Node), VALUE                               :: I_Node
                  INTEGER , VALUE                                         :: I_Index
              END FUNCTION CoSimIO_Node_Coordinate
      
              !ELEMENT FUNCTIONS
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_Element_Id (I_Element) BIND(C, NAME="CoSimIO_Element_Id")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Element
                  TYPE(CoSimIO_Element), VALUE                                :: I_Element
              END FUNCTION CoSimIO_Element_Id
      
              INTEGER(KIND(CoSimIO_ElementType)) FUNCTION CoSimIO_Element_Type (I_Element) BIND(C, NAME="CoSimIO_Element_Type")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Element , CoSimIO_ElementType
                  TYPE(CoSimIO_Element), VALUE                                :: I_Element
              END FUNCTION CoSimIO_Element_Type
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_Element_NumberOfNodes (I_Element) BIND(C, NAME="CoSimIO_Element_NumberOfNodes")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Element
                  TYPE(CoSimIO_Element), VALUE                                :: I_Element
              END FUNCTION CoSimIO_Element_NumberOfNodes
      
              TYPE(CoSimIO_Node) FUNCTION CoSimIO_Element_GetNodeByIndex (I_Element,I_Index) 
     &                                                BIND(C, NAME="CoSimIO_Element_GetNodeByIndex")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node, CoSimIO_Element
                  TYPE(CoSimIO_Element), VALUE                                :: I_Element
                  INTEGER , VALUE                                             :: I_Index
              END FUNCTION CoSimIO_Element_GetNodeByIndex
      
              !MODELPART FUNCTIONS
                
            
      
              TYPE(CoSimIO_ModelPart) FUNCTION CoSimIO_CreateModelPart (I_Name) BIND(C, NAME="CoSimIO_CreateModelPart")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart
                  CHARACTER(KIND=c_char), INTENT(IN)       :: I_Name(*)
              END FUNCTION CoSimIO_CreateModelPart
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_FreeModelPart (I_ModelPart) BIND(C, NAME="CoSimIO_FreeModelPart")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
              END FUNCTION CoSimIO_FreeModelPart
      
              type(c_ptr) FUNCTION CoSimIO_ModelPart_Name (I_ModelPart) BIND(C, NAME="CoSimIO_ModelPart_Name")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
              END FUNCTION CoSimIO_ModelPart_Name
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_ModelPart_NumberOfNodes (I_ModelPart) 
     &                                                    BIND(C, NAME="CoSimIO_ModelPart_NumberOfNodes")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
              END FUNCTION CoSimIO_ModelPart_NumberOfNodes
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_ModelPart_NumberOfLocalNodes (I_ModelPart) 
     &                                                    BIND(C, NAME="CoSimIO_ModelPart_NumberOfLocalNodes")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
              END FUNCTION CoSimIO_ModelPart_NumberOfLocalNodes
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_ModelPart_NumberOfGhostNodes (I_ModelPart) 
     &                                                    BIND(C, NAME="CoSimIO_ModelPart_NumberOfGhostNodes")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
              END FUNCTION CoSimIO_ModelPart_NumberOfGhostNodes
      
              INTEGER(KIND=c_int) FUNCTION CoSimIO_ModelPart_NumberOfElements (I_ModelPart) 
     &                                                    BIND(C, NAME="CoSimIO_ModelPart_NumberOfElements")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
              END FUNCTION CoSimIO_ModelPart_NumberOfElements
      
              TYPE(CoSimIO_Node) FUNCTION CoSimIO_ModelPart_GetNodeByIndex (I_ModelPart,I_Index) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_GetNodeByIndex")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node, CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
                  INTEGER , VALUE                                               :: I_Index
              END FUNCTION CoSimIO_ModelPart_GetNodeByIndex
      
              TYPE(CoSimIO_Node) FUNCTION CoSimIO_ModelPart_GetLocalNodeByIndex (I_ModelPart,I_Index) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_GetLocalNodeByIndex")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node, CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
                  INTEGER , VALUE                                               :: I_Index
              END FUNCTION CoSimIO_ModelPart_GetLocalNodeByIndex
      
              TYPE(CoSimIO_Node) FUNCTION CoSimIO_ModelPart_GetGhostNodeByIndex (I_ModelPart,I_Index) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_GetGhostNodeByIndex")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node, CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
                  INTEGER , VALUE                                               :: I_Index
              END FUNCTION CoSimIO_ModelPart_GetGhostNodeByIndex
      
              TYPE(CoSimIO_Node) FUNCTION CoSimIO_ModelPart_GetNodeById (I_ModelPart,I_Index) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_GetNodeById")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node, CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                :: I_ModelPart
                  INTEGER , VALUE                                               :: I_Index
              END FUNCTION CoSimIO_ModelPart_GetNodeById
      
              TYPE(CoSimIO_Element) FUNCTION CoSimIO_ModelPart_GetElementByIndex (I_ModelPart,I_Index) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_GetElementByIndex")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Element, CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                              :: I_ModelPart
                  INTEGER , VALUE                                             :: I_Index
              END FUNCTION CoSimIO_ModelPart_GetElementByIndex
      
              TYPE(CoSimIO_Element) FUNCTION CoSimIO_ModelPart_GetElementById (I_ModelPart,I_Index) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_GetElementById")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Element, CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                              :: I_ModelPart
                  INTEGER , VALUE                                             :: I_Index
              END FUNCTION CoSimIO_ModelPart_GetElementById
            
            
              SUBROUTINE CoSimIO_ModelPart_Clear(I_ModelPart) BIND(C, NAME='CoSimIO_ModelPart_Clear')
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), INTENT(IN) , VALUE          :: I_ModelPart          
              END SUBROUTINE CoSimIO_ModelPart_Clear
      
              TYPE(CoSimIO_Node) FUNCTION CoSimIO_ModelPart_CreateNewNode (I_ModelPart,I_Id,I_X,I_Y,I_Z) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_CreateNewNode")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node, CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                                         :: I_ModelPart
                  INTEGER , VALUE                                                        :: I_Id
                  REAL(kind=c_double), VALUE                                             :: I_X
                  REAL(kind=c_double), VALUE                                             :: I_Y
                  REAL(kind=c_double), VALUE                                             :: I_Z
              END FUNCTION CoSimIO_ModelPart_CreateNewNode
      
              TYPE(CoSimIO_Node) FUNCTION CoSimIO_ModelPart_CreateNewGhostNode (I_ModelPart,I_Id,I_X,I_Y,I_Z,PartitionIndex) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_CreateNewGhostNode")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Node, CoSimIO_ModelPart
                  TYPE(CoSimIO_ModelPart), VALUE                             :: I_ModelPart
                  INTEGER , VALUE                                            :: I_Id
                  REAL(kind=c_double), VALUE                                 :: I_X
                  REAL(kind=c_double), VALUE                                 :: I_Y
                  REAL(kind=c_double), VALUE                                 :: I_Z
                  INTEGER , VALUE                                            :: PartitionIndex
              END FUNCTION CoSimIO_ModelPart_CreateNewGhostNode
      
              TYPE(CoSimIO_Element) FUNCTION CoSimIO_ModelPart_CreateNewElement (I_ModelPart,I_Id, I_Type, I_Connectivities) 
     &                                                BIND(C, NAME="CoSimIO_ModelPart_CreateNewElement")
                  USE, INTRINSIC :: ISO_C_BINDING
                  IMPORT CoSimIO_Element, CoSimIO_ModelPart, CoSimIO_ElementType
                  TYPE(CoSimIO_ModelPart), VALUE                             :: I_ModelPart
                  INTEGER , VALUE                                            :: I_Id
                  INTEGER(KIND(CoSimIO_ElementType)), VALUE                  :: I_Type
                  INTEGER(KIND=c_int),DIMENSION(*)                           :: I_Connectivities
              END FUNCTION CoSimIO_ModelPart_CreateNewElement
      
      
      
  
      ! END : co_sim_io_c_model_part.h interface 
      

      
      END INTERFACE
      
      CONTAINS
      
        FUNCTION CoSimIO_Info_GetString_Wrapper(I_Info, I_Key) RESULT(f_string)
        
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr,c_f_pointer,c_char,c_null_char
            
            TYPE(CoSimIO_Info), VALUE                                   :: I_Info
            CHARACTER(LEN = 1, KIND=c_char), DIMENSION(*)               :: I_Key
            TYPE(c_ptr)                                                 :: c_string_pointer
            CHARACTER(len=:), ALLOCATABLE                               :: f_string
            CHARACTER(kind=c_char), DIMENSION(:), POINTER               :: char_pointer => null()
            CHARACTER(len=255)                                          :: temp_string
            INTEGER                                                     :: i,length
            
            IF (ALLOCATED(f_string)) DEALLOCATE (f_string)
            
            c_string_pointer = CoSimIO_Info_GetString(I_Info, I_Key)
            
            CALL c_f_pointer(c_string_pointer,char_pointer,[255])

            temp_string = " "
            DO i=1,255
              IF (char_pointer(i)==c_null_char) THEN
                length=i-1;EXIT
              END IF
              temp_string(i:i)=char_pointer(i)
            END DO
            
            ALLOCATE(CHARACTER(len=length)::f_string)
            f_string = TRIM(temp_string(1:length))//ACHAR(0)
        END FUNCTION CoSimIO_Info_GetString_Wrapper
        
      
        FUNCTION CoSimIO_ModelPart_Name_Wrapper(I_ModelPart) RESULT(f_string)
        
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr,c_f_pointer,c_char,c_null_char
            
            TYPE(CoSimIO_ModelPart), VALUE                              :: I_ModelPart
            TYPE(c_ptr)                                                 :: c_string_pointer
            CHARACTER(len=:), ALLOCATABLE                               :: f_string
            CHARACTER(kind=c_char), DIMENSION(:), POINTER               :: char_pointer => null()
            CHARACTER(len=255)                                          :: temp_string
            INTEGER                                                     :: i,length
            
            IF (ALLOCATED(f_string)) DEALLOCATE (f_string)
            
            c_string_pointer = CoSimIO_ModelPart_Name (I_ModelPart)
            
            CALL c_f_pointer(c_string_pointer,char_pointer,[255])

            temp_string = " "
            DO i=1,255
              IF (char_pointer(i)==c_null_char) THEN
                length=i-1;EXIT
              END IF
              temp_string(i:i)=char_pointer(i)
            END DO
            
            ALLOCATE(CHARACTER(len=length)::f_string)
            f_string = TRIM(temp_string(1:length))//ACHAR(0)
            
        END FUNCTION CoSimIO_ModelPart_Name_Wrapper
      
      END MODULE co_sim_io