//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2010-05-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ZOLTAN_PARTITION_UTILITY)
#define  KRATOS_ZOLTAN_PARTITION_UTILITY

#ifdef _OPENMP
#include <omp.h>
#endif

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

// Mpi
#include "mpi.h"
#include "zoltan_cpp.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/dof.h"
#include "includes/variables.h"

#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "processes/node_erase_process.h"
#include "custom_utilities/parallel_fill_communicator.h"


namespace Kratos
{
///This Function is designed to refine a mesh of triangles or tetrahedra in MPI
///This is achieved by using the trilinos epetra facilities. Please note that Trilinos has
///to be patched to work correctly, meaning that versions > 10.6 are needed to compile this

class ZoltanPartitionUtility
{
public:

    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef vector<Matrix> Matrix_Order_Tensor;
    typedef vector<Vector> Vector_Order_Tensor;
    typedef vector<Vector_Order_Tensor> Node_Vector_Order_Tensor;
    typedef Node < 3 > PointType;
    typedef Node < 3 > ::Pointer PointPointerType;
    typedef std::vector<PointType::Pointer> PointVector;
    typedef PointVector::iterator PointIterator;
    
    typedef struct{
      int           numGlobalPoints;
      int           numLocalPoints;
      ZOLTAN_ID_PTR myGlobalIDs;
      float         *x;
      float         *y;
      float         *z;
    } MESH_DATA;

    /**constructor:
     *@param ModelPart& the model part to be refined
     *@param Epetra_MpiComm the Epetra Communicator to be used
     */

    ZoltanPartitionUtility()
    {
    }

    ~ZoltanPartitionUtility()
    {
    }
    
    static int get_number_of_objects(void *data, 
                                     int *ierr)
    {
        MESH_DATA *mesh= (MESH_DATA *)data;
        *ierr = ZOLTAN_OK;
        
        return mesh->numLocalPoints;
    }

    static void get_object_list(void *data, 
                                int sizeGID, 
                                int sizeLID,
                                ZOLTAN_ID_PTR globalID, 
                                ZOLTAN_ID_PTR localID,
                                int wgt_dim, 
                                float *obj_wgts, 
                                int *ierr)
    {
        int i;
        MESH_DATA *mesh= (MESH_DATA *)data;
        *ierr = ZOLTAN_OK;

        /* In this example, return the IDs of our objects, but no weights.
        * Zoltan will assume equally weighted objects.
        */

        for (i=0; i<mesh->numLocalPoints; i++)
        {
            globalID[i] = mesh->myGlobalIDs[i];
            localID[i] = i;
        }
    }

    static int get_num_geometry(void *data, 
                                int *ierr)
    {
        *ierr = ZOLTAN_OK;
        
        return 2;
    }

    static void get_geometry_list(void *data, 
                                  int sizeGID, 
                                  int sizeLID,
                                  int num_obj,
                                  ZOLTAN_ID_PTR globalID, 
                                  ZOLTAN_ID_PTR localID,
                                  int num_dim, 
                                  double *geom_vec, 
                                  int *ierr)
    {
        MESH_DATA *mesh= (MESH_DATA *)data;

        if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 2))
        {
            *ierr = ZOLTAN_FATAL;
            
            return;
        }

        *ierr = ZOLTAN_OK;

        for (int i=0;  i < num_obj ; i++)
        {
            geom_vec[3*i + 0] = (double)mesh->x[i];
            geom_vec[3*i + 1] = (double)mesh->y[i];
            geom_vec[3*i + 2] = (double)mesh->z[i];
        }

        return;
    }
    
    void read_model_part(ModelPart& rModelPart, MESH_DATA *myMesh)
    {
        myMesh->numGlobalPoints = rModelPart.GetCommunicator().LocalMesh().Nodes().size();
        myMesh->numLocalPoints  = rModelPart.GetCommunicator().LocalMesh().Nodes().size();
        
        for(ModelPart::NodesContainerType::iterator i_node = rModelPart.GetCommunicator().LocalMesh().Nodes().begin();
            i_node != rModelPart.GetCommunicator().LocalMesh().Nodes().end();
            i_node++
           )
        {
            unsigned int index = i_node - rModelPart.GetCommunicator().LocalMesh().Nodes().begin();
          
            myMesh->myGlobalIDs[index] = i_node->Id();
            myMesh->x[index] = i_node->X();
            myMesh->y[index] = i_node->Y();
            myMesh->z[index] = i_node->Z();
        }
    }

    void CalculatePartition(ModelPart& rModelPart)
    {
        KRATOS_TRY
        
        int changes, numGidEntries, numLidEntries, numImport, numExport;
        ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
        int *importProcs, *importToPart, *exportProcs, *exportToPart;
        
        float version;
        struct Zoltan_Struct *zz;
        MESH_DATA myMesh;
                
        read_model_part(rModelPart, &myMesh);
        
        Zoltan_Initialize(0, NULL, &version);    // argc & argv can be null as they are used for the MPI_init call.
        
        // Dynamically create Zoltan object.
        zz = Zoltan_Create(MPI_COMM_WORLD);
        
        // Set general parameters
        Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
        Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
        Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
        Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
        Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
        
        // Set RCB parameters
        Zoltan_Set_Param(zz, "KEEP_CUTS", "1");
        Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
        Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1");
        
        // Query functions, to provide geometry to Zoltan
        Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &myMesh);
        Zoltan_Set_Obj_List_Fn(zz, get_object_list, &myMesh);
        Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &myMesh);
        Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, &myMesh);
        
        // Perform partitioning 
        Zoltan_LB_Partition(zz,// input (all remaining fields are output)
          &changes,             // 1 if partitioning was changed, 0 otherwise
          &numGidEntries,       // Number of integers used for a global ID
          &numLidEntries,       // Number of integers used for a local ID
          &numImport,           // Number of vertices to be sent to me
          &importGlobalGids,    // Global IDs of vertices to be sent to me
          &importLocalGids,     // Local IDs of vertices to be sent to me
          &importProcs,         // Process rank for source of each incoming vertex
          &importToPart,        // New partition for each incoming vertex
          &numExport,           // Number of vertices I must send to other processes
          &exportGlobalGids,    // Global IDs of the vertices I must send
          &exportLocalGids,     // Local IDs of the vertices I must send
          &exportProcs,         // Process to which I send each of the vertices
          &exportToPart);       // Partition to which each vertex will belong

        // Explicitly delete the Zoltan object
        Zoltan_Destroy(&zz);

        KRATOS_CATCH("")
    }
    ///************************************************************************************************
    ///************************************************************************************************
    ///function to print DETAILED mesh information. WARNING: to be used for debugging only as many informations
    ///are plotted

    void PrintDebugInfo()
    {
        KRATOS_TRY


        KRATOS_CATCH("");
    }

};

} // namespace Kratos.

#endif // KRATOS_ZOLTAN_PARTITION_UTILITY  defined 


