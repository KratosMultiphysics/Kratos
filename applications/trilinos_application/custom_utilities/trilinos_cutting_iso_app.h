
//   Project Name:        Kratos
//   Last Modified by:    $Author: pbecker $
//   Date:                $Date: 2010-05-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_TRILINOS_LOCAL_CUTTING_ISO_APP)
#define KRATOS_TRILINOS_LOCAL_CUTTING_ISO_APP 


#ifdef _OPENMP
#include <omp.h>
#endif




// System includes
#include <string>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <algorithm>


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/dof.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"


// #include "containers/array_1d.h"
// #include "processes/find_nodal_neighbours_process.h"
// #include "processes/find_elements_neighbours_process.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
//#include "utilities/split_triangle.h"
#include "utilities/split_triangle.c"
#include "utilities/split_tetrahedra.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/geometry_utilities.h"
//#include "processes/node_erase_process.h"
#include "custom_utilities/parallel_fill_communicator.h"
// #include "spatial_containers/spatial_containers.h"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"

#include "Epetra_FEVector.h"

#include "Epetra_Import.h"
#include "Epetra_MpiComm.h"
#include "parallel_fill_communicator.h"



namespace Kratos
{
/// This Function is designed to generate isosurfaces from a given variable.
/// using the origin model part (tetraedras) it creates the isosurface in a different model part.
/// this is useful to avoid printing the whole thetraedra mesh when models are too large, just as with the cutting application
/// It creates nodes and triangles that define the isosurfaces and it interpolates the variables from the original model part

class TrilinosCuttingIsosurfaceApplication
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

    /**constructor:
     *@param ModelPart& the model part to be refined
     *@param Epetra_MpiComm the Epetra Communicator to be used
     */

    TrilinosCuttingIsosurfaceApplication(Epetra_MpiComm& Comm) : mrComm(Comm)
    {
        //smallest_edge=1.0;
    } //

    ~TrilinosCuttingIsosurfaceApplication()
    {
    }


    ///************************************************************************************************
    ///************************************************************************************************




    /// ADDSKINCONDITIONS: THIS FUNCTION ADDS TO THE NEW MODEL PART THE DATA OF THE CONDITIONS BELONGING TO THE OLD MODEL PART, THE NODES COORDINATES ALREADY EXIST. WE ONLY NEED TO COPY THEM INTO THE NEW MODEL PART.
    /** this function adds the skin condtion.
        WARNING: They have to be triangles and it CAN'T be empty, otherwise a segmentation fault will appear
     * @param mr_model_part . original model part
     * @param mr__new_model_part . destinantion model part
     * @param plane_number . layer to add the conditions (integer)
     **/

    void AddSkinConditions(ModelPart& mr_model_part, ModelPart& mr_new_model_part, int plane_number )
    {
        ModelPart& this_model_part = mr_model_part;
        ModelPart& new_model_part = mr_new_model_part;
        if (mrComm.MyPID() == 0)
        {
            KRATOS_WATCH("Adding Skin Conditions to the new model part, added in layer:")
            KRATOS_WATCH(plane_number)
        }
        KRATOS_TRY


        NodesArrayType& rNodes_old = this_model_part.Nodes();        //needed to find the position of the node in the node array
        NodesArrayType::iterator it_begin_node_old = rNodes_old.ptr_begin();

        int nlocal_nodes = this_model_part.GetCommunicator().LocalMesh().Nodes().size(); //only local nodes
        int nnodes = this_model_part.Nodes().size(); //both local nodes and the ones belonging to other CPU

        //first we create a non overlapping map (only local nodes)
        int *local_ids = new int[nlocal_nodes];
        int k = 0;
        for (ModelPart::NodesContainerType::iterator it = this_model_part.GetCommunicator().LocalMesh().NodesBegin(); it != this_model_part.GetCommunicator().LocalMesh().NodesEnd(); it++)
        {
            local_ids[k++] = it->Id() - 1;
        }
        Kratos::shared_ptr<Epetra_Map> pmy_map = Kratos::make_shared<Epetra_Map>(-1, nlocal_nodes, local_ids, 0, mrComm);
        delete [] local_ids;

        //now create a map that has overlapping elements ... that is both local and ghosts
        //generate a map with the ids of the nodes
        int *ids = new int[nnodes];
        k = 0;
        for (ModelPart::NodesContainerType::iterator it = this_model_part.NodesBegin(); it != this_model_part.NodesEnd(); it++)
        {
            ids[k++] = it->Id() - 1;
        }
        Kratos::shared_ptr<Epetra_Map> pmy_ov_map = Kratos::make_shared<Epetra_Map>(-1, nnodes, ids, 0, mrComm);
        delete [] ids;


        //now we create a non overlapping vector:aux_non_overlapping_graph
        //this one will store which processor is the owner of each of the node's we'll be creating
        Kratos::shared_ptr<Epetra_FEVector > aux_non_overlapping_graph = Kratos::make_shared<Epetra_FEVector>(*pmy_map,1,false);

        aux_non_overlapping_graph->PutScalar(-1.0); //zero means this node will not have to be cloned into the new model part

        //the id of our processor:
        double this_partition_index = double(mrComm.MyPID());

        //here we'll store the nodes that will be cloned. NOT ids but position in the nodes array (full, not only owned)
        bool* used_nodes = new bool  [nnodes];
        for (int index = 0; index!=nnodes; ++index) used_nodes[index]=false; //starting as zero (no needed nodes)


        //now we have to search for the conditions. they can only be triangles, otherwise -> not stored
        int number_of_local_conditions=0;
        int aux_ids;
        for (ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); it++)
        {
            Geometry<Node < 3 > >& geom = it->GetGeometry();
            if (geom.size()==3)
            {
                ++number_of_local_conditions; //conditions are always owned
                for (unsigned int i = 0; i < geom.size(); i++)
                {
                    aux_ids = geom[i].Id() - 1;
                    int node_position = this_model_part.Nodes().find(geom[i].Id()) - it_begin_node_old; //probably there-s a better way to do this, i only need the position in the array, (not the ID)
                    used_nodes[node_position]=true; //we will have to clone this node into the new model part, no matter if owned or not.
                    int ierr = aux_non_overlapping_graph->ReplaceGlobalValues( 1 , &aux_ids, &this_partition_index,0); // saving that this processor owns this node
                    if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure ->ln 183", "");
                }
            }
        }


        int ierr = -1;
        ierr = aux_non_overlapping_graph->GlobalAssemble(Insert,true); //Epetra_CombineMode mode=Add);
        if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure --> ln 249", "");
        //now in our local graph we have also the nodes that are required by other processors


        double* local_non_ov = new double  [nlocal_nodes]; //a human readeable copy of the FEvector
        ierr = aux_non_overlapping_graph->ExtractCopy(local_non_ov,nlocal_nodes);
        if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure", "");


        int n_owned_nonzeros = 0;

        for (int index=0; index!=nlocal_nodes; ++index)
        {
            if (local_non_ov[index]>(-0.5))
            {
                //counting the number of new nodes this processor will own
                ++n_owned_nonzeros;
            }
        }

        // COUNTING THE NUMBER OF NODES AND CONDITIONS THAT WILL BE CREATED BY PREVIOUS PROCESSES OF THE CONDITIONS
        int nodes_before;
        mrComm.ScanSum(&n_owned_nonzeros, &nodes_before, 1);
        nodes_before -= n_owned_nonzeros;      //number of nodes created by other previous processors used for conditions

        int conditions_before = -1; // same but for conditions
        mrComm.ScanSum(&number_of_local_conditions, &conditions_before, 1);
        conditions_before -=number_of_local_conditions;

        //NOW WE COUNT THE NUMBER OF NODES AND TRIANGLES CREATED BY PREVIOUS CUTTING PLANES *from the new model part
        //find our the total number of nodes
        int number_of_local_nodes = new_model_part.Nodes().size(); // CHANGED from THIS_MODEL_PART to NEW_MODEL_PART. these are the nodes that i have from previous cuts
        int number_of_old_nodes = -1;                                //we're also counting ghost nodes! there will be 'holes' in the node ids between each cutting plane and the following one. ( = ghost nodes)
        mrComm.SumAll(&number_of_local_nodes, &number_of_old_nodes, 1);
        if (number_of_old_nodes < 0) number_of_old_nodes = 0;

        int number_of_local_old_conditions = new_model_part.Conditions().size(); // CHANGED from THIS_MODEL_PART to NEW_MODEL_PART. these are the nodes that i have from previous cuts
        int number_of_old_conditions = -1;                                //we're also counting ghost nodes! there will be 'holes' in the node ids between each cutting plane and the following one. ( = ghost nodes)
        mrComm.SumAll(&number_of_local_old_conditions, &number_of_old_conditions, 1);
        //KRATOS_WATCH(total_number_of_nodes)
        if (number_of_old_conditions < 0) number_of_old_conditions = 0;


        //adding original variables plus the 2 new ones that we need
        new_model_part.GetNodalSolutionStepVariablesList() = this_model_part.GetNodalSolutionStepVariablesList();
        new_model_part.AddNodalSolutionStepVariable(FATHER_NODES);
        new_model_part.AddNodalSolutionStepVariable(WEIGHT_FATHER_NODES);


        //info from the original model part
        ConditionsArrayType& rConditions = this_model_part.Conditions();
        //ConditionsArrayType::iterator cond_it_begin = rConditions.ptr_begin();
	// ConditionsArrayType::iterator cond_it_end = rConditions.ptr_end();

        //ConditionsArrayType& rConditions_new = new_model_part.Conditions();
        //ConditionsArrayType::iterator cond_it_end_new = rConditions_new.ptr_end();
	//        ConditionsArrayType::iterator cond_it_begin_new = rConditions_new.ptr_begin();

        //NodesArrayType& rNodes_new = new_model_part.Nodes();        //i need the model part just to check the id of the new nodes.
        //NodesArrayType::iterator it_end_node_new = rNodes_new.ptr_end();
        //NodesArrayType::iterator it_begin_node_new = rNodes_new.ptr_begin();


        Kratos::shared_ptr<Epetra_FEVector> IDs_non_overlapping_graph = Kratos::make_shared<Epetra_FEVector>(*pmy_map,1,false); //name self explaining
        //KRATOS_WATCH(number_of_old_nodes) ; KRATOS_WATCH(nodes_before);
        int node_id=number_of_old_nodes+nodes_before; //nodes we have previously.
        for (int index=0; index!=nlocal_nodes; ++index)
        {
            if (local_non_ov[index]>(-0.5))    //actually it can only belong to this processor, otherwise it has a -1, meaning it doesnt have to be cloned
            {
                ++node_id;
                double node_id_double=double(node_id);
                ierr = IDs_non_overlapping_graph->ReplaceMyValue ( index  ,  0 ,  node_id_double ); //saving the ID
                KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;
            }
        }

        ierr = -1;
        ierr = IDs_non_overlapping_graph->GlobalAssemble(Insert,true); //Epetra_CombineMode mode=Add);
        if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure --> ln 249", "");

        //KRATOS_WATCH('line333')

        //NOW WE MUST CREATE THE VECTOR CONTAINING BOTH OWNED AND NOT OWNED NODES.
        //ACTUALLY THEY'RE CREATED THE SAME WAY, BUT THE OWNER PROCESSOR IS THE ONE THAT SETS THE IDS
        //ONCE DONE, THESE IDs ARE SHARED AND NODES THAT ARE REAPETED ARE CREATED TWICE BY DIFFERENT PROCESSOR BUT WITH THE SAME ID number.
        Epetra_Import importer(*pmy_ov_map, *pmy_map);
        Kratos::shared_ptr<Epetra_FEVector> IDs_overlap = Kratos::make_shared<Epetra_FEVector>(*pmy_ov_map,1,false);

        ierr = IDs_overlap->Import(*IDs_non_overlapping_graph, importer, Insert);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

        double* local_IDs_ov = new double  [nnodes]; //copy in a human readeable format. the FEvector refuses to use the operator []
        ierr = IDs_overlap->ExtractCopy(local_IDs_ov,nnodes);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

        //we must now create an overlapping fevector with te partition index:
        Kratos::shared_ptr<Epetra_FEVector> Partition_overlap = Kratos::make_shared<Epetra_FEVector>(*pmy_ov_map,1,false);

        ierr = Partition_overlap->Import(*aux_non_overlapping_graph, importer, Insert);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

        double* local_Partition_ov = new double  [nnodes]; //copy in a human readeable format. the FEvector refuses to use the operator []
        ierr = Partition_overlap->ExtractCopy(local_Partition_ov,nnodes);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

        ///KRATOS_WATCH('line349')

        for (int index=0; index!=nnodes; ++index)  //now we have to save the nodes!.
        {
            // KRATOS_WATCH(auxiliar);
            if (used_nodes[index]==true)
            {
                //KRATOS_WATCH(index)
                int node_number= int(local_IDs_ov[index]);

                ModelPart::NodesContainerType::iterator it_node = this_model_part.Nodes().begin()+index; //CHANGE THIS! this only work if the nodes in the original model part are consecutive. or??

                //Node < 3 > ::Pointer pnode = new_model_part.CreateNewNode(node_number, it_node->X(), it_node->Y(), it_node->Z());  //recordar que es el nueevo model part!!
                Node < 3 >::Pointer pnode = Node < 3 > ::Pointer (new Node < 3 >(node_number, it_node->X(), it_node->Y(), it_node->Z()));
                pnode->SetSolutionStepVariablesList(&(new_model_part.GetNodalSolutionStepVariablesList()));

                pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());
                pnode->GetValue(FATHER_NODES).resize(0);
                pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node.base() ) );       // we keep the same size despite we only need one. to have everyhing with the same size
                pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node.base() ) );
                pnode-> GetValue(WEIGHT_FATHER_NODES) = 1.0;  //
                pnode->X0() = it_node->X0();
                pnode->Y0() = it_node->Y0();
                pnode->Z0() = it_node->Z0();
                pnode->FastGetSolutionStepValue(PARTITION_INDEX)=local_Partition_ov[index];
                new_model_part.Nodes().push_back(pnode);
                //std::cout <<  mrComm.MyPID() << " " << pnode->Id() <<" " << pnode->X0() << " " <<pnode->Y0() <<pnode->Z0() <<std::endl;
            }
        }
        new_model_part.Nodes().Sort();

        new_model_part.Nodes().Unique();

        //KRATOS_WATCH('line404')

        //NOW WE MUST COPY THE CONDITIONS
        int triangle_ID = number_of_old_conditions + conditions_before;
        vector<int>  triangle_nodes(3); //here we'll save the nodes' ids with the new node names
        Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");         //condition type
        Properties::Pointer properties = this_model_part.GetMesh().pGetProperties(plane_number); 		//this will allow us later to turn this layer on/off in GID

        for(ModelPart::ConditionsContainerType::iterator i_condition = rConditions.begin() ; i_condition != rConditions.end() ; i_condition++) //looping all the conditions
        {
            Geometry<Node<3> >&geom = i_condition->GetGeometry(); //current condition(nodes, etc)
            if (geom.size()==3)
            {
                for(unsigned int i = 0; i < i_condition->GetGeometry().size() ; i++)         //looping the nodes
                {
                    int position = (geom[i].Id()) - 1; //the id of the node minus one is the position in the Condition_Node array (position0 = node1)
                    int local_position = (*pmy_ov_map).LID(position); //changing from global id to local id
                    triangle_nodes[i]=local_IDs_ov[local_position] ; // saving the i nodeId
                } //nodes id saved. now we have to create the element.
                Triangle3D3<Node<3> > geometry(
                    new_model_part.Nodes()(triangle_nodes[0]),  //condition to be added
                    new_model_part.Nodes()(triangle_nodes[1]),
                    new_model_part.Nodes()(triangle_nodes[2])
                );
                ++triangle_ID;
                Condition::Pointer p_condition = rReferenceCondition.Create(triangle_ID, geometry, properties);
                new_model_part.Conditions().push_back(p_condition); //and done! added a new triangloe to the new model part
            }
        }

        Clear();
        ParallelFillCommunicator(new_model_part).Execute(); //changed from PrintDebugInfo to Execute
        if (mrComm.MyPID() == 0) std::cout << "copyng conditions and recalculation plan have been completed" << std::endl;
        KRATOS_CATCH("")
    }





    template <class TDataType>
    //this is the function to create isosurfaces from a scalar variable. the other (a component from a vectorial variable).
    void GenerateVariableCut(ModelPart& mr_model_part, ModelPart& mr_new_model_part, Variable<TDataType>& variable, double isovalue, int plane_number, float tolerance)
    {
        KRATOS_TRY
        if (mrComm.MyPID() == 0)
        {
	  std::cout <<"Generating Isosurface with the following data:"<<std::endl;
            KRATOS_WATCH(variable);
            KRATOS_WATCH(isovalue);
        }

        ModelPart& this_model_part = mr_model_part;
        //ModelPart& new_model_part = mr_new_model_part;

        vector<int> Elems_In_Plane(this_model_part.Elements().size()); //our (int) vector, where we write 1 when the element is cut by the cutting plane. when it is 2, it means we have 2 triangles (4 cutting points)
        int number_of_triangles = 0;

        Kratos::shared_ptr<Epetra_FECrsMatrix> p_edge_ids; //helper matrix to assign ids to the edges to be  refined
        Kratos::shared_ptr<Epetra_FECrsMatrix> p_partition_ids; //helper matrix to assign a partition to the edges
        vector<int> List_New_Nodes; ///* the news nodes
        vector<int> partition_new_nodes; ///* the news nodes
        vector<array_1d<int, 2 > > father_node_ids; ///* edges where are the news nodes
        vector< array_1d<double, 3 > > Coordinate_New_Node; ///* the coordinate of the new nodes
        Kratos::shared_ptr<Epetra_FECrsMatrix> used_nodes_matrix;

        PointerVector< Element > New_Elements;
        PointerVector< Condition > New_Conditions;

        CSR_Row_Matrix(mr_model_part, p_edge_ids, used_nodes_matrix);
        MPI_Barrier(MPI_COMM_WORLD);
        if (mrComm.MyPID() == 0) std::cout << "index matrix constructed" << std::endl;

        FirstLoop(mr_model_part, p_edge_ids, p_partition_ids, variable, isovalue , number_of_triangles, Elems_In_Plane, tolerance, used_nodes_matrix);
        MPI_Barrier(MPI_COMM_WORLD);
        if (mrComm.MyPID() == 0) std::cout << "Search_Edge_To_Be_Refined completed" << std::endl;

        Create_List_Of_New_Nodes(mr_model_part, mr_new_model_part, p_edge_ids, p_partition_ids, List_New_Nodes, partition_new_nodes, father_node_ids, used_nodes_matrix);
        MPI_Barrier(MPI_COMM_WORLD);
        if (mrComm.MyPID() == 0) std::cout << "Create_List_Of_New_Nodes completed" << std::endl;

        Calculate_Coordinate_And_Insert_New_Nodes(mr_model_part, mr_new_model_part, father_node_ids, List_New_Nodes, partition_new_nodes, variable, isovalue , tolerance);
        MPI_Barrier(MPI_COMM_WORLD);
        if (mrComm.MyPID() == 0) std::cout << "Calculate_Coordinate_And_Insert_New_Nodes completed" << std::endl;

        GenerateElements(mr_model_part, mr_new_model_part, Elems_In_Plane, p_edge_ids, plane_number, number_of_triangles, variable);
        MPI_Barrier(MPI_COMM_WORLD);
        if (mrComm.MyPID() == 0) std::cout << "finished generating elements" << std::endl;

        //fill the communicator
        ParallelFillCommunicator(mr_new_model_part).Execute();
        if (mrComm.MyPID() == 0) std::cout << "recalculation of communication plan completed" << std::endl;

        //clean up the data
        Clear();

        KRATOS_CATCH("")
    }



    void Clear()
    {
        KRATOS_TRY
        Kratos::shared_ptr<Epetra_Map> empty_map;
        empty_map.swap(mp_non_overlapping_map);

        mtotal_number_of_existing_nodes = 0;

        Kratos::shared_ptr<Epetra_FECrsGraph> empty1;
        mp_non_overlapping_graph.swap(empty1);

        Kratos::shared_ptr<Epetra_CrsGraph> empty2;
        mp_overlapping_graph.swap(empty2);


        KRATOS_CATCH("");
    }



    void CSR_Row_Matrix(
        ModelPart& this_model_part,
        Kratos::shared_ptr<Epetra_FECrsMatrix>& p_edge_ids,
        Kratos::shared_ptr<Epetra_FECrsMatrix>& used_nodes_matrix)
    {
        KRATOS_TRY

        //generate a map with the ids of the nodes
        int nlocal_nodes = this_model_part.GetCommunicator().LocalMesh().Nodes().size();
        int *local_ids = new int[nlocal_nodes];
        int k = 0;
        for (ModelPart::NodesContainerType::iterator it = this_model_part.GetCommunicator().LocalMesh().NodesBegin(); it != this_model_part.GetCommunicator().LocalMesh().NodesEnd(); it++)
        {
            local_ids[k++] = it->Id() - 1;
        }

        Kratos::shared_ptr<Epetra_Map> pmy_map = Kratos::make_shared<Epetra_Map>(-1, nlocal_nodes, local_ids, 0, mrComm);
        mp_non_overlapping_map.swap(pmy_map);
        delete [] local_ids;

        //now create a matrix that has overlapping elements ... that is both local and ghosts
        //generate a map with the ids of the nodes
        int nnodes = this_model_part.Nodes().size();
        int *ids = new int[nnodes];
        k = 0;
        for (ModelPart::NodesContainerType::iterator it = this_model_part.NodesBegin(); it != this_model_part.NodesEnd(); it++)
        {
            ids[k++] = it->Id() - 1;
        }

        Kratos::shared_ptr<Epetra_Map> pmy_ov_map = Kratos::make_shared<Epetra_Map>(-1, nnodes, ids, 0, mrComm);
        mp_overlapping_map.swap(pmy_ov_map);
        delete [] ids;

        //generate the graph
        int guess_row_size = 20;
        Kratos::shared_ptr<Epetra_FECrsGraph> aux_non_overlapping_graph = Kratos::make_shared<Epetra_FECrsGraph>(Copy, *mp_non_overlapping_map, guess_row_size);
        aux_non_overlapping_graph.swap(mp_non_overlapping_graph);

        int aux_ids[4];
        for (ModelPart::ElementsContainerType::iterator it = this_model_part.ElementsBegin(); it != this_model_part.ElementsEnd(); it++)
        {
            Geometry<Node < 3 > >& geom = it->GetGeometry();
            for (unsigned int i = 0; i < geom.size(); i++)
                aux_ids[i] = geom[i].Id() - 1;

            int ierr = mp_non_overlapping_graph->InsertGlobalIndices(geom.size(), aux_ids, geom.size(), aux_ids);
            KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;
        }
        for (ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); it++)
        {
            Geometry<Node < 3 > >& geom = it->GetGeometry();
            for (unsigned int i = 0; i < geom.size(); i++)
                aux_ids[i] = geom[i].Id() - 1;

            int ierr = mp_non_overlapping_graph->InsertGlobalIndices(geom.size(), aux_ids, geom.size(), aux_ids);
            KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;
        }
        mp_non_overlapping_graph->GlobalAssemble();


        //fill the edge_matrix
        Kratos::shared_ptr<Epetra_FECrsMatrix> pA = Kratos::make_shared<Epetra_FECrsMatrix>(Copy, *mp_non_overlapping_graph);
        pA->PutScalar(-1.0);

        pA.swap(p_edge_ids);


        Kratos::shared_ptr<Epetra_FECrsMatrix> pB = Kratos::make_shared<Epetra_FECrsMatrix>(Copy, *mp_overlapping_map, guess_row_size);
        pB->PutScalar(0.0);

        pB.swap(used_nodes_matrix);

        KRATOS_CATCH("")

    }
    ///************************************************************************************************
    ///************************************************************************************************

    void FirstLoop(ModelPart& this_model_part, Kratos::shared_ptr<Epetra_FECrsMatrix>& p_edge_ids,
                   Kratos::shared_ptr<Epetra_FECrsMatrix>& p_partition_ids, Variable<double>& variable,
                   double isovalue, int& number_of_triangles, vector<int>& Elems_In_Plane, double tolerance,
                   Kratos::shared_ptr<Epetra_FECrsMatrix>& used_nodes_matrix)//

    {
        KRATOS_TRY
        ElementsArrayType& rElements = this_model_part.Elements();

        Kratos::shared_ptr<Epetra_FECrsMatrix> p_nonoverlapping_partitions
        = Kratos::make_shared<Epetra_FECrsMatrix>(Copy, *mp_non_overlapping_graph);
        p_nonoverlapping_partitions->PutScalar(-1.0);

        double this_partition_index = double(mrComm.MyPID());

        //first of all create a matrix with no overlap
        ElementsArrayType::iterator it_begin = rElements.ptr_begin(); //+element_partition[k];
        ElementsArrayType::iterator it_end = rElements.ptr_end(); //+element_partition[k+1];

        double node_value; // node to closest point in the plane
        double neigh_value; //other node of the edge (neighbour) to closest point in the plane
        double diff_node_value;                // difference between the imposed value of the variable and its value in the node
        double diff_neigh_value;; //distance between the two nodes of the edge
        //double diff_node_neigh;
        //array_1d<double, 3 > temp_dist; //aux segment
        //array_1d<double, 3 > node_coord; //
        //array_1d<double, 3 > neigh_coord; //
        array_1d<unsigned int, 4 > list_matching_nodes; // used to save the new nodes that match exactly old nodes  (very unlikely, but might be 4 for very plane elements)
        unsigned int exact_nodes = 0;
        unsigned int outside_nodes = 0;
        number_of_triangles = 0;
        int current_element = 0; //current element. it's a position. NOT ID!
        int number_of_cuts = 0; //this is the counter explained in the following lines

        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
        {
            ++current_element;
            number_of_cuts = 0;
            exact_nodes = 0;
            outside_nodes = 0;
            Geometry<Node < 3 > >&geom = it->GetGeometry(); //geometry of the element
            for (unsigned int i = 0; i < it->GetGeometry().size(); i++) //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
                //when we have a triangle inside a thetraedra, its edges (or nodes) must be cut 3 times by the plane. if we loop all 2 times we can have a counter. when it's = 6 then we have a triangle. when tetraedras are cutted 8 times then we have 2 triangles (or a cuatrilateral, the same)
            {
                node_value= geom[i].FastGetSolutionStepValue(variable);
                diff_node_value = isovalue - node_value; // dist = (xnode-xp)*versor closest point-plane distance
                for (unsigned int j = 0; j < it->GetGeometry().size(); j++) //  looping on the neighbours
                {
                    //unsigned int this_node_cutted_edges=0;
                    if (i != j) //(cant link node with itself)
                    {
                        neigh_value= geom[j].FastGetSolutionStepValue(variable);
                        diff_neigh_value = isovalue - neigh_value;
                        //diff_node_neigh = node_value - neigh_value;
                        //now that we have the two points of the edge defined we can check whether it is cut by the plane or not
                        //bool isovernode = false; // if true, then it can't be between the nodes

                        if (fabs(diff_node_value) < (tolerance)  ) //then our node is part of the plane (this should have been done before the loop on neighbours, but this way it is easier to read .
                        {
                            int index_i = geom[i].Id() - 1; // i node id
                            double value = -1.0;
                            double true_value = -1.0;
                            int ierr = p_edge_ids->SumIntoGlobalValues(1, &index_i, 1, &index_i, &value);
                            if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure --> ln 235", "");
                            ierr = p_nonoverlapping_partitions->ReplaceGlobalValues(1, &index_i, 1, &index_i, &this_partition_index);
                            if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure --> ln 237", "");
                            ierr = used_nodes_matrix->ReplaceGlobalValues(1, &index_i, 1, &index_i, &true_value);
                            if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure --> ln 237", "");
                            //isovernode = true;
                            number_of_cuts += 2; //since its neighbour wont take this case as a cut, we must save 2 cuts instead of one. (to reach number_of_cuts=6),
                            ++exact_nodes;
                            list_matching_nodes[i] = geom[i].Id();
                            if ((diff_node_value * diff_neigh_value) > 0.0) //if true, it means that this node, despite being close, is actually outside the element. So if we have 2 or more ouside, we will not add this triangle in this element because we consider it property of another element
                                ++outside_nodes;
                            break; //we exit the j loop.
                        }
                        if ((diff_node_value * diff_neigh_value) < 0.0 && (fabs(diff_neigh_value)>(tolerance))) // this means one is on top of the plane and the other on the bottom, no need to do more checks, it's in between!
                        {
                            int index_i = geom[i].Id() - 1; //i node id
                            int index_j = geom[j].Id() - 1; //j node id
                            double value = -1.0;
                            double true_value = -1.0;
                            ++number_of_cuts;
                            if (index_j > index_i)
                            {
                                int ierr = p_edge_ids->SumIntoGlobalValues(1, &index_i, 1, &index_j, &value);
                                if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure --> ln 235", "");
                                ierr = p_nonoverlapping_partitions->ReplaceGlobalValues(1, &index_i, 1, &index_j, &this_partition_index);
                                if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure --> ln 237", "");
                                ierr = used_nodes_matrix->ReplaceGlobalValues(1, &index_i, 1, &index_j, &true_value);
                                if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure --> ln 237", "");
                            }
                        }
                    } //closing the i!=j if
                } //closing the neighbour (j) loop
            } //closing the nodes (i) loop


            //now we have to save the data. we should get a list with the elements that will genereate triangles and the total number of triangles
            Elems_In_Plane[current_element - 1] = 0; //we initialize as 0
            if (exact_nodes < 3 || outside_nodes<2 )   //this means at least one new node has to be generated
            {
                if (number_of_cuts == 6) //it can be 8, in that case we have 2 triangles (the cut generates a square)
                {
                    number_of_triangles += 1;
                    Elems_In_Plane[current_element - 1] = 1; //i still don't know the number of the node so i'll have to do another loop later to assign to define node id's of each triangular element
                }
                else if (number_of_cuts == 8)   // 2 triangles in the element!
                {
                    number_of_triangles += 2;
                    Elems_In_Plane[ current_element - 1] = 2;
                }
            }
        } //closing the elem loop



        int ierr = -1;
        ierr = p_edge_ids->GlobalAssemble(true, Add);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

        ierr = p_nonoverlapping_partitions->GlobalAssemble(true, Insert);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

        //import the non overlapping matrix into the overlapping one
        int MaxNumEntries = p_edge_ids->MaxNumEntries() + 5;
        Epetra_Import importer(*mp_overlapping_map, *mp_non_overlapping_map);

        Kratos::shared_ptr<Epetra_FECrsMatrix> pAoverlap = Kratos::make_shared<Epetra_FECrsMatrix>(Copy, *mp_overlapping_map, MaxNumEntries);
        Kratos::shared_ptr<Epetra_FECrsMatrix> paux = Kratos::make_shared<Epetra_FECrsMatrix>(Copy, *mp_overlapping_map, MaxNumEntries);

        ierr = pAoverlap->Import(*p_edge_ids, importer, Insert);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;
        ierr = pAoverlap->FillComplete();
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

        ierr = paux->Import(*p_nonoverlapping_partitions, importer, Insert);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;
        ierr = paux->FillComplete();
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

        //perform a pointer swap - after this we have overlapping matrices with the values syncronized
        pAoverlap.swap(p_edge_ids);
        paux.swap(p_partition_ids);

        //sace the overlapping graph
        Kratos::shared_ptr<Epetra_CrsGraph> pg = Kratos::make_shared<Epetra_CrsGraph>(p_edge_ids->Graph());
        mp_overlapping_graph.swap(pg);
        //KRATOS_WATCH(number_of_triangles)
        KRATOS_CATCH("")
    }

    ///************************************************************************************************
    ///************************************************************************************************



    void Create_List_Of_New_Nodes(ModelPart& this_model_part, ModelPart& new_model_part,
                                  Kratos::shared_ptr<Epetra_FECrsMatrix>& p_edge_ids,
                                  Kratos::shared_ptr<Epetra_FECrsMatrix>& p_partition_ids,
                                  vector<int> &List_New_Nodes,
                                  vector<int> &partition_new_nodes,
                                  vector<array_1d<int, 2 > >& father_node_ids,
                                  Kratos::shared_ptr<Epetra_FECrsMatrix>& used_nodes_matrix)
    {
        KRATOS_TRY
        //here we count the new nodes on the local mesh
        int NumMyRows = p_edge_ids->NumMyRows();
        int Row; // iterator on rows
        int Col; // iterator on cols
        int MaxNumEntries = p_edge_ids->MaxNumEntries();
        double * id_values = new double[MaxNumEntries];
        double * partition_values = new double[MaxNumEntries];
        int * Indices = new int[MaxNumEntries];
        int NumEntries;
        int GlobalRow;

        int NumEntries_aux;
        int * Indices_aux = new int[MaxNumEntries];
        double * id_values_aux = new double[MaxNumEntries];

        int n_owned_nonzeros = 0;
        int n_new_nonzeros = 0; //this are the nonzeros owned or not, but that must be known
        double this_partition_index = mrComm.MyPID();

        for (Row = 0; Row < NumMyRows; ++Row)
        {
            GlobalRow = p_edge_ids->GRID(Row);
            int ierr = p_edge_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, id_values, Indices);
            if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure", "");

            ierr = p_partition_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, partition_values, Indices);
            if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure", "");

            for (Col = 0; Col != NumEntries; Col++)
            {
                if (Indices[Col] >= Row)
                    if (id_values[Col] < -1.5 && mp_overlapping_map->MyGID(Indices[Col]) == true)
                    {
                        n_new_nonzeros = n_new_nonzeros + 1;
                        if (partition_values[Col] == this_partition_index)
                            n_owned_nonzeros = n_owned_nonzeros + 1;
                    }

            }
        }
        //use the scan sum
        int nodes_before = -1;
        mrComm.ScanSum(&n_owned_nonzeros, &nodes_before, 1);
        nodes_before = nodes_before - n_owned_nonzeros;
        if (nodes_before < 0)
            KRATOS_THROW_ERROR(std::logic_error, "problem with scan sum ... giving a negative number of nodes before", "")

            //find our the total number of nodes
            int number_of_local_nodes = new_model_part.Nodes().size(); // CHANGED from THIS_MODEL_PART to NEW_MODEL_PART. these are the nodes that i have from previous cuts
        int total_number_of_nodes = -1;
        mrComm.SumAll(&number_of_local_nodes, &total_number_of_nodes, 1);
        //KRATOS_WATCH(total_number_of_nodes)
        if (total_number_of_nodes < 0) total_number_of_nodes = 0; //this means we're in the first cutting plane.
        mtotal_number_of_existing_nodes = total_number_of_nodes;
        //the ids of the new nodes we will create AND OWN will be between
        //start_id and end_id;
        //non local nodes may have ids out of this limits
        int start_id = total_number_of_nodes + nodes_before + 1;
        int end_id = start_id + n_owned_nonzeros;

        //now distribute the ids of the new nodes so that they will be the same over all of the processors.
        Kratos::shared_ptr<Epetra_FECrsMatrix> plocal_ids = Kratos::make_shared<Epetra_FECrsMatrix>(Copy, *mp_non_overlapping_graph);
        plocal_ids->PutScalar(-1);

        int id = start_id;

        for (Row = 0; Row < NumMyRows; ++Row)
        {
            GlobalRow = p_edge_ids->GRID(Row);

            int ierr = p_edge_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, id_values, Indices);
            KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

            ierr = p_partition_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, partition_values, Indices);
            KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;

            for (Col = 0; Col < NumEntries; ++Col)
            {
                if (Indices[Col] >= Row)
                    if (id_values[Col] < -1.5 &&
                            partition_values[Col] == this_partition_index &&
                            mp_overlapping_map->MyGID(Indices[Col]) == true)
                    {
                        double did = double(id);
                        plocal_ids->ReplaceGlobalValues(1, &GlobalRow, 1, &(Indices[Col]), &did);
                        id++;
                    }
            }
        }
        plocal_ids->GlobalAssemble(true, Insert);

        KRATOS_ERROR_IF(id != end_id) << "the own node count is not verified...some error occurred" << std::endl;
        //now import the non local elements
        Epetra_Import importer(*mp_overlapping_map, *mp_non_overlapping_map);


        int ierr = p_edge_ids->Import(*plocal_ids, importer, Insert);
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;
        p_edge_ids->FillComplete();
        KRATOS_ERROR_IF(ierr < 0) << "epetra failure in line" << __LINE__ << std::endl;
        //            paux.swap(p_edge_ids);

        ///* New Id  of the nodes
        List_New_Nodes.resize(n_new_nonzeros);
        partition_new_nodes.resize(n_new_nonzeros);
        father_node_ids.resize(n_new_nonzeros);

        int k = 0;
        for (Row = 0; Row < NumMyRows; ++Row)
        {
            GlobalRow = p_edge_ids->GRID(Row);
            int num_id_entries = -1;
            int ierr = p_edge_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, num_id_entries, id_values, Indices);
            if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure ->ln420", "");

            ierr = p_partition_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, partition_values, Indices);
            if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure ->ln423", "");

            //ierr = used_nodes_matrix->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, used_nodes_matrix_row, Indices);
            //if (ierr < 0) KRATOS_THROW_ERROR(std::logic_error, "epetra failure ->ln423", "");

            if (NumEntries != num_id_entries) KRATOS_THROW_ERROR(std::logic_error, "we should have the same number of new_ids and of partition_values", "");
            for (Col = 0; Col < NumEntries; ++Col)
            {
                if (Indices[Col] >= Row)
                {
                    //nodes are to be created only if they are identified by a positive ID
                    if (id_values[Col] >= 0 && mp_overlapping_map->MyGID(Indices[Col]) == true)
                    {
                        List_New_Nodes[k] = id_values[Col];
                        partition_new_nodes[k] = partition_values[Col];

                        father_node_ids[k][0] = GlobalRow + 1;
                        father_node_ids[k][1] = Indices[Col] + 1; //+1;
                        int IsNeeded = GetUpperTriangularMatrixValue(used_nodes_matrix, GlobalRow , (Indices[Col]), MaxNumEntries, NumEntries_aux, Indices_aux, id_values_aux);
                        if (IsNeeded ==0)
                        {
                            father_node_ids[k][0] = -1;
                            father_node_ids[k][1] = -1;
                        }
                        k++;
                    }
                    if (id_values[Col] < -1.5) //at this point every edge to be refined should have an Id!!
                        KRATOS_THROW_ERROR(std::logic_error, "edge to be refined without id assigned", "")
                    }
            }
        }
        if (k != int(n_new_nonzeros)) KRATOS_THROW_ERROR(std::logic_error, "number of new nodes check failed", "")


            delete [] Indices;
        delete [] id_values;
        delete [] partition_values;

        //KRATOS_WATCH(n_new_nonzeros)
        KRATOS_CATCH("")

    }

    ///************************************************************************************************
    ///************************************************************************************************
    // insert the new nodes in the model part and interopolate the variables

    void Calculate_Coordinate_And_Insert_New_Nodes(ModelPart& this_model_part, ModelPart& new_model_part,
            const vector<array_1d<int, 2 > >& father_node_ids,
            const vector<int> &List_New_Nodes,
            const vector<int> &partition_new_nodes,
            Variable<double>& variable, double isovalue, float tolerance)
    {
        KRATOS_TRY
        array_1d<double, 3 > Coord_Node_1;
        array_1d<double, 3 > Coord_Node_2;
        array_1d<double, 3 > vector_distance;
        //array_1d<double, 3 > Xp_2;
        //array_1d<double, 3 > intersection;
        //array_1d<double, 3 > temp_dist;
        double node_value;
        double neigh_value;
        double diff_node_value;
        double diff_neigh_value;
        double diff_node_neigh;
        //double dist_node_point;
        //double dist_node_neigh;
        //double dist_node_intersect;
        double weight;
        vector< array_1d<double, 3 > > Coordinate_New_Node;
        Coordinate_New_Node.resize(father_node_ids.size());

        new_model_part.GetNodalSolutionStepVariablesList() = this_model_part.GetNodalSolutionStepVariablesList();
        new_model_part.AddNodalSolutionStepVariable(FATHER_NODES);
        new_model_part.AddNodalSolutionStepVariable(WEIGHT_FATHER_NODES);

        PointerVector< Node < 3 > > new_nodes;

        MPI_Barrier(MPI_COMM_WORLD);
        if ((father_node_ids.size())!=0)
        {
            for (unsigned int i = 0; i < father_node_ids.size(); i++)    /// calculating the coordinate of the new nodes
            {

                //getting father node 1
                const int& node_i = father_node_ids[i][0];
                const int& node_j = father_node_ids[i][1];
                if (node_i>0)
                {
                    ModelPart::NodesContainerType::iterator it_node1 = this_model_part.Nodes().find(node_i);
                    if (it_node1 == this_model_part.NodesEnd())
                    {
		      std::cout << "- father node 1 - looking for Id " << node_i << " " << node_j << std::endl;
                        KRATOS_THROW_ERROR(std::logic_error, "error inexisting node", "")
                    }
                    noalias(Coord_Node_1) = it_node1->Coordinates();

                    //getting father node 2
                    ModelPart::NodesContainerType::iterator it_node2 = this_model_part.Nodes().find(node_j);
                    if (it_node2 == this_model_part.NodesEnd())
                    {
		      std::cout << "- father node 2 - looking for Id " << node_i << " " << node_j << std::endl;
                        KRATOS_THROW_ERROR(std::logic_error, "error inexisting node", "")
                    }
                    noalias(Coord_Node_2) = it_node2->Coordinates();

                    //noalias(temp_dist) = Coord_Node_1;
                    //noalias(temp_dist) -= Xp; //temp_dist =node_coord-Xpoint
                    //dist_node_point = inner_prod(temp_dist, versor); // dist = (xnode-xp)*versor closest point-plane distance
                    //dist_node_point = fabs(dist_node_point);
                    node_value = it_node1->FastGetSolutionStepValue(variable);
                    diff_node_value = isovalue - node_value; // dist = (xnode-xp)*versor closest point-plane distance

                    neigh_value = it_node2->FastGetSolutionStepValue(variable);
                    diff_neigh_value = isovalue - neigh_value; // dist = (xnode-xp)*versor closest point-plane distance

                    diff_node_neigh= node_value - neigh_value;
                    //Xp_1 = Xp - Coord_Node_1;
                    vector_distance = Coord_Node_2 - Coord_Node_1;
                    //dist_node_intersect = (inner_prod(versor, Xp_1)) / (inner_prod(versor, Xp_2)); //line-plane interesection, this is a RELATIVE distance. ====>   point= Node1 + (Node2-Node1)*dist_node_intersect
                    //dist_node_neigh = sqrt(pow((Coord_Node_1[0] - Coord_Node_2[0]), 2) + pow((Coord_Node_1[1] - Coord_Node_2[1]), 2) + pow((Coord_Node_1[2] - Coord_Node_2[2]), 2)); // distance between node and neighbour
                    if (fabs(diff_node_value) <  tolerance) weight = 0.0;
                    else if (fabs(diff_neigh_value) <  tolerance) weight = 1.0;
                    else weight = fabs(diff_node_value / diff_node_neigh);

                    if (weight > 1.05) weight=1.0; //KRATOS_WATCH("**** something's wrong! weight higher than 1! ****");

                    for (unsigned int index = 0; index != 3; ++index) //we loop the 3 coordinates)
                        if (father_node_ids[i][0] != father_node_ids[i][1])
                            Coordinate_New_Node[i][index] = Coord_Node_1[index] + vector_distance[index] * weight;
                        else
                            Coordinate_New_Node[i][index] = Coord_Node_1[index]; //when both nodes are the same it doesnt make any sense to interpolate

                    //Node < 3 > ::Pointer pnode = new_model_part.CreateNewNode(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]);  //recordar que es el nueevo model part!!
                    Node < 3 >::Pointer pnode = Node < 3 > ::Pointer (new Node < 3 >(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]));
                    pnode->SetSolutionStepVariablesList(&(new_model_part.GetNodalSolutionStepVariablesList()));
                    pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());
                    pnode->GetValue(FATHER_NODES).resize(0);
                    pnode->GetValue(FATHER_NODES).push_back(Node < 3 > ::WeakPointer(*it_node1.base()));
                    pnode->GetValue(FATHER_NODES).push_back(Node < 3 > ::WeakPointer(*it_node2.base()));
                    pnode-> GetValue(WEIGHT_FATHER_NODES) = weight;

                    pnode->X0() = weight * (it_node2->X0()) + (1.0 - weight) * it_node1->X0();
                    pnode->Y0() = weight * (it_node2->Y0()) + (1.0 - weight) * it_node1->Y0();
                    pnode->Z0() = weight * (it_node2->Z0()) + (1.0 - weight) * it_node1->Z0();
                    pnode->FastGetSolutionStepValue(PARTITION_INDEX) =  double(partition_new_nodes[i]) ;
                    new_model_part.Nodes().push_back(pnode);
                }

            }
            new_model_part.Nodes().Sort();
            unsigned int current_size = new_model_part.Nodes().size();
            new_model_part.Nodes().Unique();

            if(current_size != new_model_part.Nodes().size())
                KRATOS_THROW_ERROR(std::logic_error,"some node was duplicated! error","");


        }//closing the if we have nodes
        // else {KRATOS_WATCH("no nodes here")}

        KRATOS_CATCH("")
    }


    /// ************************************************************************************************
    /// ************************************************************************************************
    /// as the name says...

    void GenerateElements(
        ModelPart& this_model_part, ModelPart& new_model_part, vector<int> Elems_In_Plane,
        const Kratos::shared_ptr<Epetra_FECrsMatrix> p_edge_ids,
        int plane_number, int& number_of_triangles,
        Variable<double>& variable
    )
    {
        KRATOS_TRY

        array_1d<double, 3 > temp_vector1;
        array_1d<double, 3 > temp_vector2;
        array_1d<double, 3 > temp_vector3;
        array_1d<double, 3 > temp_vector4;
        array_1d<double, 3 > temp_vector5;
        array_1d<int, 6 > nodes_for_2triang; //to be used when there are 2 triangles
        double dist2; //to be used when there are 2 triangles in the tetraedra
        double dist3;
        double control;
        unsigned int temp_int;

        DenseMatrix<int> new_conectivity;

        int total_existing_elements = -1; //warning, they're conditions, not elements!
        int local_existing_elements = new_model_part.Conditions().size();
        mrComm.SumAll(&local_existing_elements, &total_existing_elements, 1);

        if (total_existing_elements < 0) total_existing_elements = 0;

        Element const rReferenceElement;
        PointerVector< Element > Old_Elements;

        int MaxNumEntries = p_edge_ids->MaxNumEntries();
        double * id_values = new double[MaxNumEntries];
        int * Indices = new int[MaxNumEntries];
        int NumEntries;



        ElementsArrayType& rElements_old = this_model_part.Elements();
        ElementsArrayType::iterator it_begin_old = rElements_old.ptr_begin();
        ElementsArrayType::iterator it_end_old = rElements_old.ptr_end();

        Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");
        Properties::Pointer properties = this_model_part.GetMesh().pGetProperties(plane_number);

        int triangle_id_int = 0; //counter for the new elements generated
        int current_element = 0; // counter to know in which thetraedra we are in.
        unsigned int triangle_nodes = 0; // number of nodes already saved (of the current element)
        bool new_node = false; //used to check whether the current node has been saved or not

        array_1d<int, 4 > TriangleNodesArray; //nodes of the element to be generated. 4 in case they're 2 triangles
        for (unsigned int k = 0; k != 4; ++k)
        {
            TriangleNodesArray[k] = 0;
        } //initializing in 0, meaning we have no nodes yet

        int elements_before = -1;
        mrComm.ScanSum(&number_of_triangles, &elements_before, 1); //the elements that have already been created in previous threads. (but not in previous cuts, that is stored in total_existing_elements
        if (elements_before < 0) elements_before = 0;
        else elements_before-=number_of_triangles;

        array_1d<double, 3 > gradient;
        for (int j = 0; j < 3; j++) //we reset the value to zero to avoid warnings
            gradient(j) = 0.0;
        ///we enter the element loop
        for (ElementsArrayType::iterator it = it_begin_old; it != it_end_old; ++it)
        {
            /////////////
            if (Elems_In_Plane[current_element - 1] != 0) //meaning wi will need this element
            {
                for (int j = 0; j < 3; j++) //we reset the value to zero
                    gradient(j) = 0.0;
                double Volume;
                BoundedMatrix<double, 4, 3> DN_DX;
                array_1d<double, 4 > N;
                Geometry< Node < 3 > >& geom = it->GetGeometry();
                GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                for (int I = 0; I < 4; I++)
                {
                    double node_value = geom[I].FastGetSolutionStepValue(variable);
                    for (int j = 0; j < 3; j++)
                        gradient(j) += DN_DX(I, j) * node_value;
                }

            }

            ++current_element;
            triangle_nodes = 0; //starting, no nodes yet
            ///we enter in the if for only one triangle in the tetraedra
            if (Elems_In_Plane[current_element - 1] == 1) //do not forget than can be both 1 or 2 triangles per tetraedra. this is the simplest case. no need to check anything, we just create an element with the 3 nodes
            {
                for (int counter = 0; counter != 4; ++counter) TriangleNodesArray[counter] = 0;
                //checking element conectivities
                for (unsigned int i = 0; i < it->GetGeometry().size(); i++)
                {
                    Geometry<Node < 3 > >&geom = it->GetGeometry(); //i node of the element
                    for (unsigned int j = 0; j < it->GetGeometry().size(); j++) //j node of the element
                    {
                        new_node = true; //by default it's a new node
                        int index_i = geom[i].Id() - 1; //i node id
                        int index_j = geom[j].Id() - 1;
                        int NodeId = GetUpperTriangularMatrixValue(p_edge_ids, index_i, index_j, MaxNumEntries, NumEntries, Indices, id_values);

                        for (unsigned int l = 0; l != 3; ++l)
                        {
                            if (TriangleNodesArray[l] == NodeId //if we have already saved this node or it has not been cutted, then we have no new node to add (coord(i,j)=-1)
                                    || NodeId < 1)
                                new_node = false;
                        }
                        //if it's a new node and the indexes are correct:
                        if (new_node && index_i <= index_j)
                        {
                            TriangleNodesArray[triangle_nodes] = NodeId;
                            triangle_nodes++;
                        }
                        if (triangle_nodes == 3) break; //if we have already found 3 nodes then we can exit
                    } //closing j node loop
                    if (triangle_nodes == 3) break; //egal
                } //closing j node loop
                //now we have to check that the normal of the element matches the one of the plane (they could have opposite directions)
                ModelPart::NodesContainerType::iterator it_node1 = new_model_part.Nodes().find(TriangleNodesArray[0]);
                noalias(temp_vector1) = it_node1->Coordinates(); //node 1

                it_node1 = new_model_part.Nodes().find(TriangleNodesArray[1]);
                noalias(temp_vector2) = it_node1->Coordinates(); //node 2

                it_node1 = new_model_part.Nodes().find(TriangleNodesArray[2]);
                noalias(temp_vector3) = it_node1->Coordinates(); //nodo 3

                temp_vector3 -= temp_vector1; //first edge
                temp_vector2 -= temp_vector1; //second edge


                MathUtils<double>::CrossProduct(temp_vector4, temp_vector2, temp_vector3); //multiplying the 2 edges gives us a normal vector to the element

                if (inner_prod(temp_vector4, gradient) < 0.0) //if the signs do not match then they have opposite directions
                {
                    temp_int = TriangleNodesArray[2];
                    TriangleNodesArray[2] = TriangleNodesArray[1];
                    TriangleNodesArray[1] = temp_int; //we switch 2 nodes and ready
                }

                //generate new Elements
                Triangle3D3<Node < 3 > > geom(
                    new_model_part.Nodes()(TriangleNodesArray[0]),
                    new_model_part.Nodes()(TriangleNodesArray[1]),
                    new_model_part.Nodes()(TriangleNodesArray[2])
                );

                Condition::Pointer p_condition = rReferenceCondition.Create(triangle_id_int + 1 + elements_before + total_existing_elements, geom, properties); //creating the element using the reference element. notice we are using the first element to avoid overriting nodes created by other cutting planes
                new_model_part.Conditions().push_back(p_condition);
                ++triangle_id_int;

            }//closing if  (1 triangle)

            ///entering now the if for 2 triangles inside the tetraedra.
            if (Elems_In_Plane[current_element - 1] == 2) //now we have 2 elements. we cant just create 2 elements with a random node order because they might overlap and not cover the whole area defined by the trapezoid
            {
                //to fix this we'll first create a plane. see below
                for (int counter = 0; counter != 4; ++counter) TriangleNodesArray[counter] = 0;
                //checking conectivities to find nodes
                for (unsigned int i = 0; i < it->GetGeometry().size(); i++) //nodo i
                {
                    Geometry<Node < 3 > >&geom = it->GetGeometry();
                    for (unsigned int j = 0; j < it->GetGeometry().size(); j++) //nodo j
                    {
                        new_node = true;
                        int index_i = geom[i].Id() - 1;
                        int index_j = geom[j].Id() - 1;
                        int NodeId = GetUpperTriangularMatrixValue(p_edge_ids, index_i, index_j, MaxNumEntries, NumEntries, Indices, id_values);

                        for (unsigned int l = 0; l != 3; ++l)
                        {
                            if (TriangleNodesArray[l] == NodeId //same as the part with only one triangle (look above)
                                    || NodeId < 1)
                            {
                                new_node = false;
                            }
                        }
                        if (new_node && index_i < index_j)
                        {
                            TriangleNodesArray[triangle_nodes] = NodeId;
                            triangle_nodes++;
                        }
                        if (triangle_nodes == 4) break; //once we've found the 4 nodes we can exit
                    } //closing i loop
                    if (triangle_nodes == 4) break;
                }
                if (triangle_nodes!=4)
                {
                    for (int counter = 0; counter != 4; ++counter) KRATOS_WATCH(TriangleNodesArray[counter]);
                    KRATOS_THROW_ERROR(  std::logic_error, "ln 1289", "")
                }
                //now we have to start checking the angles. the easiest way (i think) is creating a new plane using the original plane and a segment created by 2 nodes
                // using crossproduct we get a perpendicular plane. (either point can be used as origin).
                //since the cuadrilateral is created by the cut of teatraedra,
                //none of its internal angles can exceed 180 degrees and hence our new plane divides the cuadrilateral into 2 triangles if the distances to the other points have different signs (one on top and the other on the bottom of this new plane)
                //otherwise this edge is just an edge of the cuadrilateral and we have to look for another.
                //so let's begin! (we'll keep an origin node and we'll loop different nodes as the end of the segment till we find one that satisfies our criteria)
                ModelPart::NodesContainerType::iterator it_node1 = new_model_part.Nodes().find(TriangleNodesArray[0]);
                noalias(temp_vector1) = it_node1->Coordinates(); //nodo 1 (origin)
                int jjj;
                int kkk;
                bool error_bool=true;

                for (int iii = 1; iii != 4; ++iii)   //end node of the segment that will be used to create the plane (will be contained too)
                {
                    it_node1 = new_model_part.Nodes().find(TriangleNodesArray[iii]); //i node. we always keep node 0 as origin
                    noalias(temp_vector2) = (it_node1->Coordinates()); //node2 (end)
                    noalias(temp_vector3) = temp_vector2 - temp_vector1; //segment 1-2
                    //now i have to create the new plane
                    MathUtils<double>::CrossProduct(temp_vector4, gradient, temp_vector3); //done. now temp_vector4 is the (normal to the) new plane, perpendicular to the one containing the triangles
                    //the origin of the plane is temp_vector1 (temp_vector2 could also be used)
                    //now we need to check distances to the other nodes (i+2 (let's call them jjj and i+3=kkk since we can't go futher than i=3)
                    if (iii == 1)
                    {
                        jjj = 2;
                        kkk = 3;
                    }
                    else if (iii == 2)
                    {
                        jjj = 3;
                        kkk = 1;
                    }
                    else
                    {
                        jjj = 1;
                        kkk = 2;
                    }

                    it_node1 = new_model_part.Nodes().find(TriangleNodesArray[jjj]);
                    noalias(temp_vector2) = it_node1->Coordinates(); //one of the remaining nodes;

                    it_node1 = new_model_part.Nodes().find(TriangleNodesArray[kkk]);
                    noalias(temp_vector3) = it_node1->Coordinates(); //the other remaining node;

                    noalias(temp_vector2) -= temp_vector1; // minus origin point of the plane
                    noalias(temp_vector3) -= temp_vector1;


                    dist2 = inner_prod(temp_vector2, temp_vector4); // dot product
                    dist3 = inner_prod(temp_vector3, temp_vector4);
                    control = dist2*dist3;
                    //and that's it. we now have to check if the distance have different signs. to do so we multiply :
                    if (control < 0.0) //we have the right one! one node on each side of the plane generated by nodes 0 and iii
                    {
                        nodes_for_2triang[0] = TriangleNodesArray[0];
                        nodes_for_2triang[1] = TriangleNodesArray[jjj];
                        nodes_for_2triang[2] = TriangleNodesArray[iii]; //finish first triangle
                        nodes_for_2triang[3] = TriangleNodesArray[iii];
                        nodes_for_2triang[4] = TriangleNodesArray[kkk];
                        nodes_for_2triang[5] = TriangleNodesArray[0]; //finish 2nd triangle
                        //	KRATOS_WATCH(nodes_for_2triang);
                        error_bool=false;
                        break; //no need to keep looking, we can exit the loop

                    } //closing the if

                }//by the time this finishes i should already have TriangleNodesArray
                if (error_bool)    //for (int counter = 0; counter != 4; ++counter) KRATOS_WATCH(TriangleNodesArray[counter]);
                {
                    nodes_for_2triang[0] = TriangleNodesArray[0];
                    nodes_for_2triang[1] = TriangleNodesArray[1];
                    nodes_for_2triang[2] = TriangleNodesArray[2]; //finish first triangle
                    nodes_for_2triang[3] = TriangleNodesArray[3];
                    nodes_for_2triang[4] = TriangleNodesArray[2];
                    nodes_for_2triang[5] = TriangleNodesArray[0];
                } //finish 2nd triangle } //KRATOS_THROW_ERROR(std::logic_error, "ln 1349", "")  }

                //checking if the normal to our element is oriented correctly, just as we did when we had only 1 triangle (not commented here)
                for (int index = 0; index != 2; ++index) //for triangle 1 and 2
                {
                    //KRATOS_WATCH(nodes_for_2triang[index * 3 + 0])
                    //KRATOS_WATCH(nodes_for_2triang[index * 3 + 1])
                    //KRATOS_WATCH(nodes_for_2triang[index * 3 + 2])
                    //if (error_bool)KRATOS_WATCH(nodes_for_2triang[index * 3 + 0]);
                    it_node1 = new_model_part.Nodes().find(nodes_for_2triang[index * 3 + 0]);
                    noalias(temp_vector1) = it_node1->Coordinates(); //node 1

                    it_node1 = new_model_part.Nodes().find(nodes_for_2triang[index * 3 + 1]);
                    noalias(temp_vector2) = it_node1->Coordinates(); //node 2

                    it_node1 = new_model_part.Nodes().find(nodes_for_2triang[index * 3 + 2]);
                    noalias(temp_vector3) = it_node1->Coordinates(); //node 3

                    temp_vector3 -= temp_vector1;
                    temp_vector2 -= temp_vector1;
                    MathUtils<double>::CrossProduct(temp_vector4, temp_vector2, temp_vector3);

                    if (inner_prod(temp_vector4, gradient) < 0.0)
                    {
                        temp_int = nodes_for_2triang[index * 3 + 2];
                        nodes_for_2triang[index * 3 + 2] = nodes_for_2triang[index * 3 + 1];
                        nodes_for_2triang[index * 3 + 1] = temp_int;
                    }

                    Triangle3D3<Node < 3 > > geom(
                        new_model_part.Nodes()(nodes_for_2triang[index * 3 + 0]),
                        new_model_part.Nodes()(nodes_for_2triang[index * 3 + 1]),
                        new_model_part.Nodes()(nodes_for_2triang[index * 3 + 2])
                    );

                    Condition::Pointer p_condition = rReferenceCondition.Create(triangle_id_int + 1 + elements_before + total_existing_elements, geom, properties);

                    new_model_part.Conditions().push_back(p_condition);
                    ++triangle_id_int;

                    for (int counter = 0; counter != 4; ++counter) TriangleNodesArray[counter] = 0; //resetting, just in case

                }//cierro el index

            }//closing if elems_in_plane=2

        }//closing element loops

        KRATOS_CATCH("")



    }

    /// ************************************************************************************************
    /// ************************************************************************************************

    void UpdateCutData(ModelPart& new_model_part, ModelPart& old_model_part)
    {
      if (mrComm.MyPID() == 0) std::cout <<"Updating cut data:"<<std::endl;
        int step_data_size = old_model_part.GetNodalSolutionStepDataSize();

        //looping the nodes, no data is assigned to elements
        for (ModelPart::NodesContainerType::iterator it = new_model_part.NodesBegin(); it != new_model_part.NodesEnd(); it++)
        {
            double* node0_data = it->GetValue(FATHER_NODES)[0].SolutionStepData().Data(0); //current step only, (since we'll call this every timestep
            double* node1_data = it->GetValue(FATHER_NODES)[1].SolutionStepData().Data(0);
            double weight = it->GetValue(WEIGHT_FATHER_NODES);
            double* step_data = (it)->SolutionStepData().Data(0);
            double partition_index= it->FastGetSolutionStepValue(PARTITION_INDEX);

            //now we only have to copy the information from node_data to step_data
            for (int j = 0; j < step_data_size; j++) //looping all the variables and interpolating using weight
            {
                step_data[j] = (1.0-weight) * node0_data[j] + ( weight) * node1_data[j];
            }
            it->FastGetSolutionStepValue(PARTITION_INDEX)=partition_index;
        }//closing node loop
    }//closing subroutine


    void DeleteCutData(ModelPart& new_model_part)
    {
      if (mrComm.MyPID() == 0) std::cout <<"Deleting cut data:"<<std::endl;
        new_model_part.Nodes().clear();
        new_model_part.Conditions().clear();
        new_model_part.Elements().clear();
    }


protected:

    double smallest_edge;

    Epetra_MpiComm& mrComm;
    Kratos::shared_ptr<Epetra_Map> mp_overlapping_map;
    Kratos::shared_ptr<Epetra_Map> mp_non_overlapping_map;
    int mtotal_number_of_existing_nodes;

    Kratos::shared_ptr<Epetra_FECrsGraph> mp_non_overlapping_graph;
    Kratos::shared_ptr<Epetra_CrsGraph> mp_overlapping_graph;

    ///this function transfers the Constitutive Law internal variables from the father to the child.
    ///note that this is done through the vector Variable INTERNAL_VARIABLES which should
    ///also contain the geometric data needed for this.


    double GetValueFromRow(int row, int j, int row_size, int* indices, double* values)
    {
        for (int i = 0; i < row_size; i++)
        {
            if (indices[i] == j)
            {
                return values[i];
            }
        }
        return -1;
        //KRATOS_THROW_ERROR(std::logic_error, "expected index not found", "")
    }

    double GetUpperTriangularMatrixValue(const Kratos::shared_ptr<Epetra_FECrsMatrix>& p_edge_ids, int index_0, int index_1, int& MaxNumEntries, int& NumEntries, int* Indices, double* values)
    {
        double value;
        if (index_0 > index_1)
        {
            p_edge_ids->ExtractGlobalRowCopy(index_1, MaxNumEntries, NumEntries, values, Indices); //Coord(index_1, index_0);
            value = this->GetValueFromRow(index_1, index_0, NumEntries, Indices, values);
        }
        else
        {
            p_edge_ids->ExtractGlobalRowCopy(index_0, MaxNumEntries, NumEntries, values, Indices);
            value = this->GetValueFromRow(index_0, index_1, NumEntries, Indices, values);
        }
        return value;

    }





};



} // namespace Kratos.

#endif // KRATOS_TRILINOS_LOCAL_CUTTING_ISO_APP  defined 


