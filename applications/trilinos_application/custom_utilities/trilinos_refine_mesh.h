
///WARNING
/// it is needed compute he neirgbourg

//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2010-05-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_TRILINOS_LOCAL_REFINE_SIMPLEX_MESH)
#define  KRATOS_TRILINOS_LOCAL_REFINE_SIMPLEX_MESH


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
#include "utilities/split_tetrahedra.c"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "processes/node_erase_process.h"
#include "custom_utilities/parallel_fill_communicator.h"
// #include "spatial_containers/spatial_containers.h"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_MpiComm.h"
#include "parallel_fill_communicator.h"



namespace Kratos
{
    ///This Function is designed to refine a mesh of triangles or tetrahedra in MPI
    ///This is achieved by using the trilinos epetra facilities. Please note that Trilinos has
    ///to be patched to work correctly, meaning that versions > 10.6 are needed to compile this

    class TrilinosRefineMesh
    {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ModelPart::ConditionsContainerType ConditionsArrayType;
        typedef boost::numeric::ublas::vector<Matrix> Matrix_Order_Tensor;
        typedef boost::numeric::ublas::vector<Vector> Vector_Order_Tensor;
        typedef boost::numeric::ublas::vector<Vector_Order_Tensor> Node_Vector_Order_Tensor;
        typedef Node < 3 > PointType;
        typedef Node < 3 > ::Pointer PointPointerType;
        typedef std::vector<PointType::Pointer> PointVector;
        typedef PointVector::iterator PointIterator;

        /**constructor:
         *@param ModelPart& the model part to be refined
         *@param Epetra_MpiComm the Epetra Communicator to be used
         */

        TrilinosRefineMesh(ModelPart& model_part, Epetra_MpiComm& Comm) : mr_model_part(model_part), mrComm(Comm)
        {
        }

        ~TrilinosRefineMesh()
        {
        }


        ///************************************************************************************************
        ///************************************************************************************************

        /**Local_Refine_Mesh is the umbrella function that performs the refinement.
         * The elements to be refined are identified by the flag SPLIT_ELEMENT=true
         * current level of refinement is stored element by element in the flag
         * REFINEMENT_LEVEL of type integer
         * NOTE: the refinement is performed with the ONLY AIM of delivering a conformant mesh,
         * that is, no step is performed to ensure the quality of the resultign mesh.
         * The user is expected to select the elements to be refined so to guarantee that some
         * quality level is retained. No effort is made to keep the resulting mesh balanced, consequently
         * a load balancing step should be performed after the refinement step.
         * @param refine_on_reference --> it controls if the interpolation of data (both internal and nodal)
         *        should be performed on the initial domain (if set it to true) OR on the deformed one.
         *        this is very important for example for refining correctly total lagrangian elements.
         *        If the parameter is set to two, the model_part is expected to contain the DISPLACEMENT
         *        variable
         * @param interpolate_internal_variables --> this flag controls if the internal variables (Constitutive Law vars)
         *        should be interpolated or not
         */
        void Local_Refine_Mesh(bool refine_on_reference, bool interpolate_internal_variables, int domain_size)
        {
            KRATOS_TRY

            if (refine_on_reference == true)
                if (!(mr_model_part.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT)))
                    KRATOS_ERROR(std::logic_error, "DISPLACEMENT Variable is not in the model part -- needed if refine_on_reference = true", "");

            if (!(mr_model_part.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX)))
                KRATOS_ERROR(std::logic_error, "PARTITION_INDEX Variable is not in the model part -- mpi not possible", "");

            boost::shared_ptr<Epetra_FECrsMatrix> p_edge_ids; //helper matrix to assign ids to the edges to be  refined
            boost::shared_ptr<Epetra_FECrsMatrix> p_partition_ids; //helper matrix to assign a partition to the edges
            boost::numeric::ublas::vector<int> List_New_Nodes; ///* the news nodes
            boost::numeric::ublas::vector<int> partition_new_nodes; ///* the news nodes
            boost::numeric::ublas::vector<array_1d<int, 2 > > father_node_ids; ///* edges where are the news nodes
            boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node; ///* the coordinate of the new nodes

            PointerVector< Element > New_Elements;
            PointerVector< Condition > New_Conditions;
            if (refine_on_reference == true)
            {
                for (ModelPart::NodesContainerType::iterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); it++)
                {
                    it->X() = it->X0();
                    it->Y() = it->Y0();
                    it->Z() = it->Z0();
                }
            }

            if (mrComm.MyPID() == 0) cout << "beginning the refinement" << endl;

            CSR_Row_Matrix(mr_model_part, p_edge_ids);
            if (mrComm.MyPID() == 0) cout << "index matrix constructed" << endl;

            Search_Edge_To_Be_Refined(mr_model_part, p_edge_ids, p_partition_ids);
            if (mrComm.MyPID() == 0) cout << "Search_Edge_To_Be_Refined completed" << endl;

            Create_List_Of_New_Nodes(mr_model_part, p_edge_ids, p_partition_ids, List_New_Nodes, partition_new_nodes, father_node_ids);
            if (mrComm.MyPID() == 0) cout << "Create_List_Of_New_Nodes completed" << endl;

            Calculate_Coordinate_And_Insert_New_Nodes(mr_model_part, father_node_ids, List_New_Nodes, partition_new_nodes);
            if (mrComm.MyPID() == 0) cout << "Calculate_Coordinate_And_Insert_New_Nodes completed" << endl;

            if (domain_size == 2)
                Erase_Old_Element_And_Create_New_Triangle_Element(mr_model_part, p_edge_ids, New_Elements, interpolate_internal_variables);
            else if (domain_size == 3)
            {
                Erase_Old_Element_And_Create_New_Tetra_Element(mr_model_part, p_edge_ids, New_Elements, interpolate_internal_variables);
                Erase_Old_Condition_And_Create_New_Triangle_Conditions(mr_model_part, p_edge_ids, New_Conditions, interpolate_internal_variables);
            } else
                KRATOS_ERROR(std::logic_error, "domain size can be either 2 or 3!", "");

            
            if (mrComm.MyPID() == 0) cout << "Erase_Old_Element_And_Create_New completed" << endl;


            Renumbering_Elements(mr_model_part, New_Elements);
            Renumbering_Conditions(mr_model_part, New_Conditions);
            if (mrComm.MyPID() == 0) cout << "Renumbering completed" << endl;






            if (refine_on_reference == true)
            {
                for (ModelPart::NodesContainerType::iterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); it++)
                {
                    const array_1d<double, 3 > & disp = it->FastGetSolutionStepValue(DISPLACEMENT);
                    it->X() = it->X0() + disp[0];
                    it->Y() = it->Y0() + disp[1];
                    it->Z() = it->Z0() + disp[2];
                }
            }

            //fill the communicator
            ParallelFillCommunicator(mr_model_part).Execute();
            if (mrComm.MyPID() == 0) cout << "recalculation of communication plan completed" << endl;

            //clean up the data
            Clear();

            KRATOS_CATCH("")
        }
        ///************************************************************************************************
        ///************************************************************************************************
        ///function to print DETAILED mesh information. WARNING: to be used for debugging only as many informations
        ///are plotted

        void PrintDebugInfo()
        {
            KRATOS_TRY

            //print ghost mesh
            cout << "proc = " << mrComm.MyPID() << "ghost mesh nodes" << endl;
            for (ModelPart::NodesContainerType::iterator it = mr_model_part.GetCommunicator().GhostMesh().NodesBegin();
                    it != mr_model_part.GetCommunicator().GhostMesh().NodesEnd();
                    it++)
                cout << it->Id() << " " << it->Coordinates() << endl;

            KRATOS_CATCH("");

        }

        ///This function frees all of the memory used
        void Clear()
        {
            KRATOS_TRY
            boost::shared_ptr<Epetra_Map> empty_map;
            empty_map.swap(mp_non_overlapping_map);

            mtotal_number_of_existing_nodes = 0;

            boost::shared_ptr<Epetra_FECrsGraph> empty1;
            mp_non_overlapping_graph.swap(empty1);

            boost::shared_ptr<Epetra_CrsGraph> empty2;
            mp_overlapping_graph.swap(empty2);


            KRATOS_CATCH("");
        }



        ///************************************************************************************************
        ///************************************************************************************************
        ///helper function to resize the CSR_Row_Matrix used internally to store the ids of the newly
        ///created nodes

        void CSR_Row_Matrix(
                ModelPart& this_model_part,
                boost::shared_ptr<Epetra_FECrsMatrix>& p_edge_ids)
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

            boost::shared_ptr<Epetra_Map> pmy_map = boost::shared_ptr<Epetra_Map > (new Epetra_Map(-1, nlocal_nodes, local_ids, 0, mrComm));
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

            boost::shared_ptr<Epetra_Map> pmy_ov_map = boost::shared_ptr<Epetra_Map > (new Epetra_Map(-1, nnodes, ids, 0, mrComm));
            mp_overlapping_map.swap(pmy_ov_map);
            delete [] ids;

            //generate the graph
            int guess_row_size = 20;
            boost::shared_ptr<Epetra_FECrsGraph > aux_non_overlapping_graph = boost::shared_ptr<Epetra_FECrsGraph > (new Epetra_FECrsGraph(Copy, *mp_non_overlapping_map, guess_row_size));
            aux_non_overlapping_graph.swap(mp_non_overlapping_graph);
            //            Epetra_FECrsGraph Agraph(Copy, *mp_non_overlapping_map, guess_row_size);

            int aux_ids[4]; //this works both for 2d and 3d
            for (ModelPart::ElementsContainerType::iterator it = this_model_part.ElementsBegin(); it != this_model_part.ElementsEnd(); it++)
            {
                Geometry<Node < 3 > >& geom = it->GetGeometry();
                for (unsigned int i = 0; i < geom.size(); i++)
                    aux_ids[i] = geom[i].Id() - 1;

                int ierr = mp_non_overlapping_graph->InsertGlobalIndices(geom.size(), aux_ids, geom.size(), aux_ids);
                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure ->ln174", "");
            }
            for (ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); it++)
            {
                Geometry<Node < 3 > >& geom = it->GetGeometry();
                for (unsigned int i = 0; i < geom.size(); i++)
                    aux_ids[i] = geom[i].Id() - 1;

                int ierr = mp_non_overlapping_graph->InsertGlobalIndices(geom.size(), aux_ids, geom.size(), aux_ids);
                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure ->ln 183", "");
            }
            mp_non_overlapping_graph->GlobalAssemble();



            //fill the edge_matrix
            boost::shared_ptr<Epetra_FECrsMatrix> pA = boost::shared_ptr<Epetra_FECrsMatrix > (new Epetra_FECrsMatrix(Copy, *mp_non_overlapping_graph));
            pA->PutScalar(-1.0);

            pA.swap(p_edge_ids);


            KRATOS_CATCH("")

        }

        ///************************************************************************************************
        ///************************************************************************************************
        ///Auxiliary function to identify the edges to be refined. It performs a loop over the elements
        ///and if the element is marked for being splitted, it marks for splitting all of its edges

        void Search_Edge_To_Be_Refined(
                ModelPart& this_model_part,
                boost::shared_ptr<Epetra_FECrsMatrix>& p_edge_ids,
                boost::shared_ptr<Epetra_FECrsMatrix>& p_partition_ids)
        {
            KRATOS_TRY
            ElementsArrayType& rElements = this_model_part.Elements();

            boost::shared_ptr<Epetra_FECrsMatrix> p_nonoverlapping_partitions
                    = boost::shared_ptr<Epetra_FECrsMatrix > (new Epetra_FECrsMatrix(Copy, *mp_non_overlapping_graph));
            p_nonoverlapping_partitions->PutScalar(-1.0);

            double this_partition_index = double(mrComm.MyPID());
            //            KRATOS_ERROR(std::logic_error,"aaaaaaaaaaa","")

            //first of all create a matrix with no overlap
            ElementsArrayType::iterator it_begin = rElements.ptr_begin(); //+element_partition[k];
            ElementsArrayType::iterator it_end = rElements.ptr_end(); //+element_partition[k+1];
            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                //unsigned int level = it->GetValue(REFINEMENT_LEVEL);
                if (it->GetValue(SPLIT_ELEMENT) == true)
                {
                    Element::GeometryType& geom = it->GetGeometry(); // Nodos del elemento
                    for (unsigned int i = 0; i < geom.size(); i++)
                    {
                        int index_i = geom[i].Id() - 1;
                        for (unsigned int j = 0; j < geom.size(); j++)
                        {
                            int index_j = geom[j].Id() - 1;
                            double value = -1.0;
                            if (index_j > index_i)
                            {
                                int ierr = p_edge_ids->SumIntoGlobalValues(1, &index_i, 1, &index_j, &value);
                                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure --> ln 235", "");
                                ierr = p_nonoverlapping_partitions->ReplaceGlobalValues(1, &index_i, 1, &index_j, &this_partition_index);
                                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure --> ln 237", "");
                                //                        Coord(index_i, index_j) = -2;
                            }
                        }
                    }
                }
            }
            int ierr = -1;
            ierr = p_edge_ids->GlobalAssemble(true, Add);
            if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure --> ln 246", "");

            ierr = p_nonoverlapping_partitions->GlobalAssemble(true, Insert);
            if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure --> ln 249", "");

            //import the non overlapping matrix into the overlapping one
            int MaxNumEntries = p_edge_ids->MaxNumEntries() + 5;
            Epetra_Import importer(*mp_overlapping_map, *mp_non_overlapping_map);

            boost::shared_ptr<Epetra_FECrsMatrix> pAoverlap = boost::shared_ptr<Epetra_FECrsMatrix > (new Epetra_FECrsMatrix(Copy, *mp_overlapping_map, MaxNumEntries));
            boost::shared_ptr<Epetra_FECrsMatrix> paux = boost::shared_ptr<Epetra_FECrsMatrix > (new Epetra_FECrsMatrix(Copy, *mp_overlapping_map, MaxNumEntries));

            ierr = pAoverlap->Import(*p_edge_ids, importer, Insert);
            if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");
            ierr = pAoverlap->FillComplete();
            if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");

            ierr = paux->Import(*p_nonoverlapping_partitions, importer, Insert);
            if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");
            ierr = paux->FillComplete();
            if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");

            //perform a pointer swap - after this we have overlapping matrices with the values syncronized
            pAoverlap.swap(p_edge_ids);
            paux.swap(p_partition_ids);

            //sace the overlapping graph
            boost::shared_ptr<Epetra_CrsGraph > pg = boost::shared_ptr<Epetra_CrsGraph > (new Epetra_CrsGraph(p_edge_ids->Graph()));
            mp_overlapping_graph.swap(pg);

            KRATOS_CATCH("")
        }

        ///************************************************************************************************
        ///************************************************************************************************
        ///It defines the list of new nodes to be created basing on the data marked
        ///in the auxiliary CSR matrix. It also assigns the Id to the new nodes

        void Create_List_Of_New_Nodes(ModelPart& this_model_part,
                boost::shared_ptr<Epetra_FECrsMatrix>& p_edge_ids,
                boost::shared_ptr<Epetra_FECrsMatrix>& p_partition_ids,
                boost::numeric::ublas::vector<int> &List_New_Nodes,
                boost::numeric::ublas::vector<int> &partition_new_nodes,
                boost::numeric::ublas::vector<array_1d<int, 2 > >& father_node_ids)
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

            int n_owned_nonzeros = 0;
            int n_new_nonzeros = 0; //this are the nonzeros owned or not, but that must be known
            double this_partition_index = mrComm.MyPID();

            for (Row = 0; Row < NumMyRows; ++Row)
            {
                GlobalRow = p_edge_ids->GRID(Row);
                int ierr = p_edge_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, id_values, Indices);
                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");

                ierr = p_partition_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, partition_values, Indices);
                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");

                for (Col = 0; Col != NumEntries; Col++)
                {
                    if (Indices[Col] > Row)
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
                KRATOS_ERROR(std::logic_error, "problem with scan sum ... giving a negative number of nodes before", "")

                //find our the total number of nodes
                int number_of_local_nodes = this_model_part.GetCommunicator().LocalMesh().Nodes().size();
            int total_number_of_nodes = -1;
            mrComm.SumAll(&number_of_local_nodes, &total_number_of_nodes, 1);

//            int counter=0;
//            for(ModelPart::NodesContainerType::iterator iii=this_model_part.NodesBegin(); iii!=this_model_part.NodesEnd(); iii++)
//                if(iii->FastGetSolutionStepValue(PARTITION_INDEX)==mrComm.MyPID())
//                    counter++;
//
//            if(counter!=number_of_local_nodes )
//                KRATOS_ERROR(std::logic_error,"local mesh is not up to date!","");

            mtotal_number_of_existing_nodes = total_number_of_nodes;

//for(ModelPart::NodesContainerType::iterator iii=this_model_part.NodesBegin(); iii!=this_model_part.NodesEnd(); iii++)
//    if(iii->Id() > mtotal_number_of_existing_nodes)
//        KRATOS_ERROR(std::logic_error,"node ids are not contiguous","");


            //the ids of the new nodes we will create AND OWN will be between
            //start_id and end_id;
            //non local nodes may have ids out of this limits
            int start_id = total_number_of_nodes + nodes_before + 1;
            int end_id = start_id + n_owned_nonzeros;

            //now distribute the ids of the new nodes so that they will be the same over all of the processors.
            boost::shared_ptr<Epetra_FECrsMatrix> plocal_ids = boost::shared_ptr<Epetra_FECrsMatrix > (new Epetra_FECrsMatrix(Copy, *mp_non_overlapping_graph));
            plocal_ids->PutScalar(-1);

            int id = start_id;
            //            Epetra_FECrsGraph Agraph(Copy, *mp_non_overlapping_map, MaxNumEntries);
            for (Row = 0; Row < NumMyRows; ++Row)
            {
                GlobalRow = p_edge_ids->GRID(Row);

                int ierr = p_edge_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, id_values, Indices);
                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");

                ierr = p_partition_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, partition_values, Indices);
                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");

                for (Col = 0; Col < NumEntries; ++Col)
                {
                    if (Indices[Col] > Row)
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

            if (id != end_id) KRATOS_ERROR(std::logic_error, "the own node count is not verified...some error occurred", "");
            //now import the non local elements
            Epetra_Import importer(*mp_overlapping_map, *mp_non_overlapping_map);

            //            boost::shared_ptr<Epetra_FECrsMatrix> paux = boost::shared_ptr<Epetra_FECrsMatrix > (new Epetra_FECrsMatrix(Copy, *mp_overlapping_graph));


            int ierr = p_edge_ids->Import(*plocal_ids, importer, Insert);
            if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");
            p_edge_ids->FillComplete();
            if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure", "");
            //            paux.swap(p_edge_ids);


            ///* New Id de los Nodos
            List_New_Nodes.resize(n_new_nonzeros);
            partition_new_nodes.resize(n_new_nonzeros);
            father_node_ids.resize(n_new_nonzeros);



            int k = 0;
            for (Row = 0; Row < NumMyRows; ++Row)
            {
                GlobalRow = p_edge_ids->GRID(Row);
                int num_id_entries = -1;
                int ierr = p_edge_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, num_id_entries, id_values, Indices);
                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure ->ln420", "");

                ierr = p_partition_ids->ExtractGlobalRowCopy(GlobalRow, MaxNumEntries, NumEntries, partition_values, Indices);
                if (ierr < 0) KRATOS_ERROR(std::logic_error, "epetra failure ->ln423", "");

                if (NumEntries != num_id_entries) KRATOS_ERROR(std::logic_error, "we should have the same number of new_ids and of partition_values", "")
                    for (Col = 0; Col < NumEntries; ++Col)
                    {
                        if (Indices[Col] > Row)
                        {
                            //nodes are to be created only if they are identified by a positive ID
                            //AND if they both appear in the overlapping map for the current processor
                            if (id_values[Col] >= 0 && mp_overlapping_map->MyGID(Indices[Col]) == true)
                            {
//if(this_model_part.Nodes().find(id_values[Col])!=this_model_part.Nodes().end())
//    KRATOS_ERROR(std::logic_error,"node already existing","");

                                List_New_Nodes[k] = id_values[Col];
                                partition_new_nodes[k] = partition_values[Col];

                                father_node_ids[k][0] = GlobalRow + 1;
                                father_node_ids[k][1] = Indices[Col] + 1; //+1;

                                k++;
                            }
                            if (id_values[Col] < -1.5) //at this point every edge to be refined should have an Id!!
                                KRATOS_ERROR(std::logic_error, "edge to be refined without id assigned", "")
                            }
                    }
            }
            if (k != int(n_new_nonzeros)) KRATOS_ERROR(std::logic_error, "number of new nodes check failed", "")


                delete [] Indices;
            delete [] id_values;
            delete [] partition_values;

            KRATOS_CATCH("")

        }

        ///************************************************************************************************
        ///************************************************************************************************
        // insert the news nodes in the model part and interopolate the variables

        void Calculate_Coordinate_And_Insert_New_Nodes(ModelPart& this_model_part,
                const boost::numeric::ublas::vector<array_1d<int, 2 > >& father_node_ids,
                const boost::numeric::ublas::vector<int> &List_New_Nodes,
                const boost::numeric::ublas::vector<int> &partition_new_nodes)
        {
            KRATOS_TRY
            array_1d<double, 3 > Coord_Node_1;
            array_1d<double, 3 > Coord_Node_2;
            boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
            Coordinate_New_Node.resize(father_node_ids.size());
            unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
            Node < 3 > ::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();

            //check that all of the nodes have the same number of dofs
            for (ModelPart::NodesContainerType::iterator it = this_model_part.NodesBegin(); it != this_model_part.NodesEnd(); it++)
                if ((it->GetDofs()).size() != reference_dofs.size())
                {
                    cout << "reference_dof_list" << *(this_model_part.NodesBegin()) << endl;
                    cout << "inconsistent dof list found" << *it << endl;
                    KRATOS_ERROR(std::logic_error, "list of dofs is not the same on all of the nodes!", "")
                }

            PointerVector< Node < 3 > > new_nodes;

            for (unsigned int i = 0; i < father_node_ids.size(); i++)
            {
                /// calculating the coordinate of the news nodes

                //getting father node 1
                const int& node_i = father_node_ids[i][0];
                const int& node_j = father_node_ids[i][1];
                ModelPart::NodesContainerType::iterator it_node1 = this_model_part.Nodes().find(node_i);
                if (it_node1 == this_model_part.NodesEnd())
                {
                    cout << "- father node 1 - looking for Id " << node_i << " " << node_j << endl;
                    KRATOS_ERROR(std::logic_error, "error inexisting node", "")
                }
                std::size_t pos1 = it_node1 - this_model_part.NodesBegin();
                noalias(Coord_Node_1) = it_node1->Coordinates();

                //getting father node 2
                ModelPart::NodesContainerType::iterator it_node2 = this_model_part.Nodes().find(node_j);
                if (it_node2 == this_model_part.NodesEnd())
                {
                    cout << "- father node 2 - looking for Id " << node_i << " " << node_j << endl;
                    KRATOS_ERROR(std::logic_error, "error inexisting node", "")
                }
                std::size_t pos2 = it_node2 - this_model_part.NodesBegin();
                noalias(Coord_Node_2) = it_node2->Coordinates();


                noalias(Coordinate_New_Node[i]) = 0.50 * (Coord_Node_1 + Coord_Node_2);
                /// inserting the news node in the model part

//                if(this_model_part.Nodes().find(List_New_Nodes[i]) != this_model_part.Nodes().end())
//                {
//                    cout << this_model_part.Nodes().find(List_New_Nodes[i])->Id() << " " << this_model_part.Nodes().find(List_New_Nodes[i])->Coordinates() << "new coords " << Coordinate_New_Node[i] ;
//                    cout << "position_found = " << this_model_part.Nodes().find(List_New_Nodes[i]) - this_model_part.Nodes().begin() << " number of existing nodes " << this_model_part.Nodes().size();
//                    cout << "last id in model_part = " << (this_model_part.Nodes().end()-1)->Id() << "problematic Id() = " << List_New_Nodes[i] << endl;
//                    KRATOS_ERROR(std::logic_error,"attempting to create a duplicated node","")
//                }


                //                Node < 3 > ::Pointer pnode = this_model_part.CreateNewNode(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]);
                Node < 3 > ::Pointer pnode = AuxCreateNewNode(this_model_part, List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]);
                new_nodes.push_back(pnode);



                it_node1 = this_model_part.NodesBegin() + pos1;
                it_node2 = this_model_part.NodesBegin() + pos2;

                pnode->GetValue(FATHER_NODES).resize(0);
                pnode->GetValue(FATHER_NODES).push_back(Node < 3 > ::WeakPointer(*it_node1.base()));
                pnode->GetValue(FATHER_NODES).push_back(Node < 3 > ::WeakPointer(*it_node2.base()));

                pnode->X0() = 0.5 * (it_node1->X0() + it_node2->X0());
                pnode->Y0() = 0.5 * (it_node1->Y0() + it_node2->Y0());
                pnode->Z0() = 0.5 * (it_node1->Z0() + it_node2->Z0());


                for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
                {
                    Node < 3 > ::DofType& rDof = *iii;
                    Node < 3 > ::DofType::Pointer p_new_dof = pnode->pAddDof(rDof);
                    if (it_node1->IsFixed(iii->GetVariable()) == true && it_node2->IsFixed(iii->GetVariable()) == true)
                        (p_new_dof)->FixDof();
                    else
                    {
                        (p_new_dof)->FreeDof();
                    }

                }

                ///* intepolating the data
                unsigned int buffer_size = pnode->GetBufferSize();
                for (unsigned int step = 0; step < buffer_size; step++)
                {
                    double* new_step_data = pnode->SolutionStepData().Data(step);
                    double* step_data1 = it_node1->SolutionStepData().Data(step);
                    double* step_data2 = it_node2->SolutionStepData().Data(step);
                    ///*copying this data in the position of the vector we are interested in
                    for (unsigned int j = 0; j < step_data_size; j++)
                    {
                        new_step_data[j] = 0.5 * (step_data1[j] + step_data2[j]);
                    }
                }

                pnode->FastGetSolutionStepValue(PARTITION_INDEX) = partition_new_nodes[i];


                /// WARNING =  only for reactions;
                //                const double zero = 0.00;
                //                for (Node < 3 > ::DofsContainerType::iterator iii = pnode->GetDofs().begin(); iii != pnode->GetDofs().end(); iii++)
                //                {
                //                    if (pnode->IsFixed(iii->GetVariable()) == false)
                //                    {
                //
                //                        iii->GetSolutionStepReactionValue() = zero;
                //
                //                    }
                //                }

            }

            for (PointerVector<Node < 3 > >::iterator it = new_nodes.begin(); it != new_nodes.end(); it++)
            {
//                cout << (*it.base())->Id() << endl;

//               here soemthing wrong
               this_model_part.GetMesh(0).Nodes().push_back(*it.base());
            }

            unsigned int current_size = this_model_part.Nodes().size();

            this_model_part.Nodes().Unique();

            if(current_size != this_model_part.Nodes().size())
                KRATOS_ERROR(std::logic_error,"some node was duplicated! error","");


            KRATOS_CATCH("")
        }


        /// ************************************************************************************************
        /// ************************************************************************************************
        /// as the name says...

        void Erase_Old_Element_And_Create_New_Triangle_Element(
                ModelPart& this_model_part,
                const boost::shared_ptr<Epetra_FECrsMatrix> p_edge_ids,
                PointerVector< Element >& New_Elements,
                bool interpolate_internal_variables
                )
        {
            KRATOS_TRY
            boost::numeric::ublas::matrix<int> new_conectivity;

            int total_existing_elements = -1;
            int local_existing_elements = this_model_part.Elements().size();
            mrComm.SumAll(&local_existing_elements, &total_existing_elements, 1);

            Element const rReferenceElement;
            unsigned int to_be_deleted = 0;
            unsigned int large_id = total_existing_elements * 15;
            bool create_element = false;
            int edge_ids[3];
            int t[12];
            int nel = 0;
            int splitted_edges = 0;
            int nint = 0;
            int aux[6];

            ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();
            PointerVector< Element > Old_Elements;

            unsigned int current_id = (this_model_part.Elements().end() - 1)->Id() + 1;
            for (ModelPart::ElementsContainerType::iterator it = this_model_part.Elements().begin();
                    it != this_model_part.Elements().end(); ++it)
            {
                for (unsigned int i = 0; i < 12; i++)
                {
                    t[i] = -1;
                }
                Element::GeometryType& geom = it->GetGeometry();
                Calculate_Triangle_Edges(geom, p_edge_ids, edge_ids, aux);

                ///* crea las nuevas conectividades
                create_element = Split_Triangle(edge_ids, t, &nel, &splitted_edges, &nint);

                ///* crea los nuevos elementos
                if (create_element == true)
                {
                    to_be_deleted++;
                    for (int i = 0; i < nel; i++)
                    {
                        int i0, i1, i2;
                        TriangleGetNewConnectivityGID(i, t, aux, &i0, &i1, &i2);

                        Triangle2D3<Node < 3 > > geom(
                                this_model_part.Nodes()(i0),
                                this_model_part.Nodes()(i1),
                                this_model_part.Nodes()(i2)
                                );


                        Element::Pointer p_element;
                        p_element = it->Create(current_id, geom, it->pGetProperties());
                        p_element->Initialize();
                        p_element->InitializeSolutionStep(rCurrentProcessInfo);
                        p_element->FinalizeSolutionStep(rCurrentProcessInfo);

                        /// setting the internal variables in the child elem
                        if (interpolate_internal_variables == true)
                            InterpolateInternalVariables(nel, *it.base(), p_element, rCurrentProcessInfo);

                        // Transfer elemental variables
                        p_element->Data() = it->Data();
                        //const unsigned int& level = it->GetValue(REFINEMENT_LEVEL);
                        p_element->GetValue(SPLIT_ELEMENT) = false;
                        //p_element->SetValue(REFINEMENT_LEVEL, 1);
                        New_Elements.push_back(p_element);
                        current_id++;

                    }
                    it->SetId(large_id);
                    large_id++;
                }

            }

            ///* adding news elements to the model part
            for (PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new != New_Elements.end(); it_new++)
            {
                this_model_part.Elements().push_back(*(it_new.base()));
            }

            ///* all of the elements to be erased are at the end
            this_model_part.Elements().Sort();

            ///*now remove all of the "old" elements
            this_model_part.Elements().erase(this_model_part.Elements().end() - to_be_deleted, this_model_part.Elements().end());

            int number_of_elements_before = -1;
            int number_of_own_elements = this_model_part.Elements().size();
            mrComm.ScanSum(&number_of_own_elements, &number_of_elements_before, 1);
            number_of_elements_before = number_of_elements_before - number_of_own_elements;
            if (number_of_elements_before < 0)
                KRATOS_ERROR(std::logic_error, "problem with scan sum ... giving a negative number of nodes before", "");

            int start_elem_id = number_of_elements_before + 1;
            //            int end_elem_id = number_of_elements_before + number_of_own_elements;

            int id_counter = start_elem_id;
            for (ModelPart::ElementsContainerType::iterator it = this_model_part.Elements().begin();
                    it != this_model_part.Elements().end(); ++it)
            {
                it->SetId(id_counter++);
            }

            //                //find our the total number of nodes
            //                int number_of_local_nodes = this_model_part.GetCommunicator().LocalMesh().Nodes().size();
            //            int total_number_of_nodes = -1;
            //            mrComm.SumAll(&number_of_local_nodes, &total_number_of_nodes, 1);


            KRATOS_CATCH("")
        }

        /// ************************************************************************************************
        /// ************************************************************************************************

        void Calculate_Triangle_Edges(Element::GeometryType& geom,
                const boost::shared_ptr<Epetra_FECrsMatrix> p_edge_ids,
                int* edge_ids,
                int aux[6]
                )
        {
            KRATOS_TRY
                    int index_0 = geom[0].Id() - 1;
            int index_1 = geom[1].Id() - 1;
            int index_2 = geom[2].Id() - 1;

            aux[0] = geom[0].Id();
            aux[1] = geom[1].Id();
            aux[2] = geom[2].Id();

            //here we count the new nodes on the local mesh
            int MaxNumEntries = p_edge_ids->MaxNumEntries();
            double * id_values = new double[MaxNumEntries];
            int * Indices = new int[MaxNumEntries];
            int NumEntries;

            aux[3] = GetUpperTriangularMatrixValue(p_edge_ids, index_0, index_1, MaxNumEntries, NumEntries, Indices, id_values);
            aux[4] = GetUpperTriangularMatrixValue(p_edge_ids, index_1, index_2, MaxNumEntries, NumEntries, Indices, id_values);
            aux[5] = GetUpperTriangularMatrixValue(p_edge_ids, index_2, index_0, MaxNumEntries, NumEntries, Indices, id_values);

            TriangleSplitMode(aux, edge_ids);


            KRATOS_CATCH("")
        }


        /// ************************************************************************************************
        /// ************************************************************************************************
        /// function to ensure that the elements and conditions are numbered consecutively

        void Renumbering_Elements(ModelPart& this_model_part,
                PointerVector< Element >& New_Elements
                )
        {
            KRATOS_TRY
                    //compute element ids
                    int nel = this_model_part.Elements().size();
            int el_before = -1;
            mrComm.ScanSum(&nel, &el_before, 1);
            //            int el_end_id = el_before + 1;
            int el_start_id = el_before - nel + 1;
            if (el_start_id < 0) KRATOS_ERROR(std::logic_error, "wrong id in renumbering of the elements", "")

                unsigned int id_elem = el_start_id;
            for (ModelPart::ElementsContainerType::iterator it = this_model_part.ElementsBegin();
                    it != this_model_part.ElementsEnd(); ++it)
            {
                if (it->Id() != id_elem)
                {
                    it->SetId(id_elem);
                }
                id_elem++;
            }

            //verify that the nodes are numbered consecutively
            int local_max_id = (mr_model_part.GetCommunicator().LocalMesh().NodesEnd() - 1)->Id();
            int local_number_of_nodes = mr_model_part.GetCommunicator().LocalMesh().Nodes().size();
            int max_id = -1;
            int tot_nnodes = -1;
            mrComm.MaxAll(&local_max_id, &max_id, 1);
            mrComm.SumAll(&local_number_of_nodes, &tot_nnodes, 1);
            if (max_id != tot_nnodes)
                KRATOS_ERROR(std::logic_error, "node ids are not consecutive", "");

            int tot_elements;
            int local_max_elem_id = (mr_model_part.Elements().end() - 1)->Id();
            int local_number_of_elem = mr_model_part.Elements().size();
            mrComm.MaxAll(&local_max_elem_id, &max_id, 1);
            mrComm.SumAll(&local_number_of_elem, &tot_elements, 1);
            if (max_id != tot_elements)
                KRATOS_ERROR(std::logic_error, "element ids are not consecutive", "");

            KRATOS_CATCH("")

        }

        /// ************************************************************************************************
        /// ************************************************************************************************
        /// function to ensure that the Conditions and conditions are numbered consecutively

        void Renumbering_Conditions(ModelPart& this_model_part,
                PointerVector< Condition >& New_Conditions
                )
        {
            KRATOS_TRY


                    //compute Condition ids
                    int nel = this_model_part.Conditions().size();

            if (nel != 0)
            {
                int el_before = -1;
                mrComm.ScanSum(&nel, &el_before, 1);
                //            int el_end_id = el_before + 1;
                int el_start_id = el_before - nel + 1;
                if (el_start_id < 0) KRATOS_ERROR(std::logic_error, "wrong id in renumbering of the Conditions", "")

                    unsigned int id_elem = el_start_id;
                for (ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin();
                        it != this_model_part.ConditionsEnd(); ++it)
                {
                    if (it->Id() != id_elem)
                    {
                        it->SetId(id_elem);
                    }
                    id_elem++;
                }

                //verify that the nodes are numbered consecutively
                int local_max_id = (mr_model_part.GetCommunicator().LocalMesh().NodesEnd() - 1)->Id();
                int local_number_of_nodes = mr_model_part.GetCommunicator().LocalMesh().Nodes().size();
                int max_id = -1;
                int tot_nnodes = -1;
                mrComm.MaxAll(&local_max_id, &max_id, 1);
                mrComm.SumAll(&local_number_of_nodes, &tot_nnodes, 1);
                if (max_id != tot_nnodes)
                    KRATOS_ERROR(std::logic_error, "node ids are not consecutive", "");

                int tot_Conditions;
                int local_max_elem_id = (mr_model_part.Conditions().end() - 1)->Id();
                int local_number_of_elem = mr_model_part.Conditions().size();
                mrComm.MaxAll(&local_max_elem_id, &max_id, 1);
                mrComm.SumAll(&local_number_of_elem, &tot_Conditions, 1);
                if (max_id != tot_Conditions)
                    KRATOS_ERROR(std::logic_error, "Condition ids are not consecutive", "");

            }
            KRATOS_CATCH("")

        }


    protected:
        ModelPart& mr_model_part;
        Epetra_MpiComm& mrComm;
        boost::shared_ptr<Epetra_Map> mp_overlapping_map;
        boost::shared_ptr<Epetra_Map> mp_non_overlapping_map;
        int mtotal_number_of_existing_nodes;

        boost::shared_ptr<Epetra_FECrsGraph> mp_non_overlapping_graph;
        boost::shared_ptr<Epetra_CrsGraph> mp_overlapping_graph;

        ///this function transfers the Constitutive Law internal variables from the father to the child.
        ///note that this is done through the vector Variable INTERNAL_VARIABLES which should
        ///also contain the geometric data needed for this.

        void InterpolateInternalVariables(const int& nel,
                const Element::Pointer father_elem,
                Element::Pointer child_elem,
                ProcessInfo& rCurrentProcessInfo)
        {
            std::vector<Vector> values;
            father_elem->GetValueOnIntegrationPoints(INTERNAL_VARIABLES, values, rCurrentProcessInfo);
            /* /// WARNING =  Calculando la longitud ponderada de fractura del elemento. Solo valido para isotropic_damage
             Element::GeometryType& geom_father = father_elem->GetGeometry();
             Element::GeometryType& geom_child  = child_elem->GetGeometry();
             double area_father = geom_father.Area();
             double area_child  = geom_child.Area();
             values[0][4]       = (area_child/area_father) * values[0][4];
             */
            child_elem->SetValueOnIntegrationPoints(INTERNAL_VARIABLES, values, rCurrentProcessInfo);

        }

        double GetValueFromRow(int row, int j, int row_size, int* indices, double* values)
        {
            for (int i = 0; i < row_size; i++)
            {
                if (indices[i] == j)
                {
                    return values[i];
                }
            }

            KRATOS_ERROR(std::logic_error, "expected index not found", "")
        }

        double GetUpperTriangularMatrixValue(const boost::shared_ptr<Epetra_FECrsMatrix>& p_edge_ids, int index_0, int index_1, int& MaxNumEntries, int& NumEntries, int* Indices, double* values)
        {
            double value;
            if (index_0 > index_1)
            {
                p_edge_ids->ExtractGlobalRowCopy(index_1, MaxNumEntries, NumEntries, values, Indices); //Coord(index_1, index_0);
                value = this->GetValueFromRow(index_1, index_0, NumEntries, Indices, values);
            } else
            {
                p_edge_ids->ExtractGlobalRowCopy(index_0, MaxNumEntries, NumEntries, values, Indices);
                value = this->GetValueFromRow(index_0, index_1, NumEntries, Indices, values);
            }
            return value;

        }

        void Erase_Old_Element_And_Create_New_Tetra_Element(
                ModelPart& this_model_part,
                const boost::shared_ptr<Epetra_FECrsMatrix> p_edge_ids,
                PointerVector< Element >& New_Elements,
                bool interpolate_internal_variables
                )
        {
            if (this_model_part.Elements().size() > 0)
            {
                //		  boost::numeric::ublas::matrix<int> new_conectivity;

                //		  Element const rReferenceElement;
                unsigned int to_be_deleted = 0;

                int total_existing_elements = -1;
                int local_existing_elements = this_model_part.Elements().size();
                mrComm.SumAll(&local_existing_elements, &total_existing_elements, 1);

                unsigned int large_id = total_existing_elements * 15;
                unsigned int current_id = local_existing_elements + 1;
                bool create_element = false;

                ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();
                int edge_ids[6];
                int t[56];
                for (unsigned int i = 0; i < 56; i++)
                {
                    t[i] = -1;
                }
                int nel = 0;
                int splitted_edges = 0;
                int internal_node = 0;

                ModelPart::NodesContainerType center_nodes;

                for (ModelPart::ElementsContainerType::iterator it = this_model_part.ElementsBegin(); it != this_model_part.ElementsEnd(); ++it)
                {
                    //KRATOS_WATCH(it->Id())
                    Element::GeometryType& geom = it->GetGeometry();

                    int index_0 = geom[0].Id() - 1;
                    int index_1 = geom[1].Id() - 1;
                    int index_2 = geom[2].Id() - 1;
                    int index_3 = geom[3].Id() - 1;

                    //put the global ids in aux
                    int aux[11];
                    aux[0] = geom[0].Id();
                    aux[1] = geom[1].Id();
                    aux[2] = geom[2].Id();
                    aux[3] = geom[3].Id();

                    //here we count the new nodes on the local mesh
                    int MaxNumEntries = p_edge_ids->MaxNumEntries();
                    double * id_values = new double[MaxNumEntries];
                    int * Indices = new int[MaxNumEntries];
                    int NumEntries;

                    aux[4] = GetUpperTriangularMatrixValue(p_edge_ids, index_0, index_1, MaxNumEntries, NumEntries, Indices, id_values);
                    aux[5] = GetUpperTriangularMatrixValue(p_edge_ids, index_0, index_2, MaxNumEntries, NumEntries, Indices, id_values);
                    aux[6] = GetUpperTriangularMatrixValue(p_edge_ids, index_0, index_3, MaxNumEntries, NumEntries, Indices, id_values);
                    aux[7] = GetUpperTriangularMatrixValue(p_edge_ids, index_1, index_2, MaxNumEntries, NumEntries, Indices, id_values);
                    aux[8] = GetUpperTriangularMatrixValue(p_edge_ids, index_1, index_3, MaxNumEntries, NumEntries, Indices, id_values);
                    aux[9] = GetUpperTriangularMatrixValue(p_edge_ids, index_2, index_3, MaxNumEntries, NumEntries, Indices, id_values);

                    TetrahedraSplitMode(aux, edge_ids);

                    create_element = Split_Tetrahedra(edge_ids, t, &nel, &splitted_edges, &internal_node);



                    if (create_element == true)
                    {
                        if (internal_node == 1)
                        {
                            aux[10] = CreateCenterNode(geom, this_model_part, center_nodes);
                        }

                        to_be_deleted++;
                        //create the new connectivity
                        for (int i = 0; i < nel; i++)
                        {

                            unsigned int base = i * 4;
                            unsigned int i0 = aux[t[base]];
                            unsigned int i1 = aux[t[base + 1]];
                            unsigned int i2 = aux[t[base + 2]];
                            unsigned int i3 = aux[t[base + 3]];
                            // KRATOS_WATCH(i0)
                            // KRATOS_WATCH(i1)
                            // KRATOS_WATCH(i2)
                            // KRATOS_WATCH(i3)
                            // KRATOS_WATCH("-------------------------------------------" )


                            Tetrahedra3D4<Node < 3 > > geom(
                                    this_model_part.Nodes()(i0),
                                    this_model_part.Nodes()(i1),
                                    this_model_part.Nodes()(i2),
                                    this_model_part.Nodes()(i3)
                                    );

                            //generate new element by cloning the base one
                            Element::Pointer p_element = it->Create(current_id, geom, it->pGetProperties());
                            p_element->Initialize();
                            p_element->InitializeSolutionStep(rCurrentProcessInfo);
                            p_element->FinalizeSolutionStep(rCurrentProcessInfo);
                            New_Elements.push_back(p_element);

                            /// setting the internal variables in the child elem
                            if (interpolate_internal_variables == true)
                                InterpolateInternalVariables(nel, *it.base(), p_element, rCurrentProcessInfo);

                            // Transfer elemental variables
                            p_element->Data() = it->Data();
                            p_element->GetValue(SPLIT_ELEMENT) = false;

                            current_id++;

                        }
                        it->SetId(large_id);
                        large_id++;
                    }

                    //		      for (unsigned int i = 0; i < 32; i++)
                    //		      {
                    //			  t[i] = -1;
                    //		      }
                }

                ///* all of the elements to be erased are at the end
                this_model_part.Elements().Sort();

                ///*now remove all of the "old" elements
                this_model_part.Elements().erase(this_model_part.Elements().end() - to_be_deleted, this_model_part.Elements().end());

                unsigned int total_size = this_model_part.Elements().size() + New_Elements.size();
                this_model_part.Elements().reserve(total_size);

                ///* adding news elements to the model part
                for (PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new != New_Elements.end(); it_new++)
                {

                    this_model_part.Elements().push_back(*(it_new.base()));
                }

                int number_of_elements_before = -1;
                int number_of_own_elements = this_model_part.Elements().size();
                mrComm.ScanSum(&number_of_own_elements, &number_of_elements_before, 1);
                number_of_elements_before = number_of_elements_before - number_of_own_elements;
                if (number_of_elements_before < 0)
                    KRATOS_ERROR(std::logic_error, "problem with scan sum ... giving a negative number of nodes before", "");

                int start_elem_id = number_of_elements_before + 1;

                int id_counter = start_elem_id;
                for (ModelPart::ElementsContainerType::iterator it = this_model_part.Elements().begin();
                        it != this_model_part.Elements().end(); ++it)
                {
                    it->SetId(id_counter++);
                }

                AssignIdToCenterNode(center_nodes);





            }



        }

        unsigned int CreateCenterNode(Geometry<Node < 3 > >& geom, ModelPart& model_part, ModelPart::NodesContainerType& center_nodes)
        {
            //determine a new unique id
            unsigned int new_id = (model_part.NodesEnd() - 1)->Id() + 1;

            if (model_part.Nodes().find(new_id) != model_part.NodesEnd())
                KRATOS_ERROR(std::logic_error, "id is already being used", "");


            //determine the coordinates of the new node
            double X = (geom[0].X() + geom[1].X() + geom[2].X() + geom[3].X()) / 4.0;
            double Y = (geom[0].Y() + geom[1].Y() + geom[2].Y() + geom[3].Y()) / 4.0;
            double Z = (geom[0].Z() + geom[1].Z() + geom[2].Z() + geom[3].Z()) / 4.0;

            double X0 = (geom[0].X0() + geom[1].X0() + geom[2].X0() + geom[3].X0()) / 4.0;
            double Y0 = (geom[0].Y0() + geom[1].Y0() + geom[2].Y0() + geom[3].Y0()) / 4.0;
            double Z0 = (geom[0].Z0() + geom[1].Z0() + geom[2].Z0() + geom[3].Z0()) / 4.0;

            //generate the new node
            Node < 3 > ::Pointer pnode = model_part.CreateNewNode(new_id, X, Y, Z);

            unsigned int buffer_size = model_part.NodesBegin()->GetBufferSize();
            pnode->SetBufferSize(buffer_size);

            pnode->X0() = X0;
            pnode->Y0() = Y0;
            pnode->Z0() = Z0;

            //add the dofs
            Node < 3 > ::DofsContainerType& reference_dofs = (model_part.NodesBegin())->GetDofs();
            unsigned int step_data_size = model_part.GetNodalSolutionStepDataSize();

            for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
            {
                Node < 3 > ::DofType& rDof = *iii;
                Node < 3 > ::DofType::Pointer p_new_dof = pnode->pAddDof(rDof);

                //the variables are left as free for the internal node
                (p_new_dof)->FreeDof();

                /*       if(it_node1->IsFixed(iii->GetVariable()) == true && it_node2->IsFixed(iii->GetVariable()) == true)
                         (p_new_dof)->FixDof();
                       else
                          { (p_new_dof)->FreeDof();}       */

            }

            ///* intepolating the data

            for (unsigned int step = 0; step < buffer_size; step++)
            {
                double* new_step_data = pnode->SolutionStepData().Data(step);
                double* step_data1 = geom[0].SolutionStepData().Data(step);
                double* step_data2 = geom[1].SolutionStepData().Data(step);
                double* step_data3 = geom[2].SolutionStepData().Data(step);
                double* step_data4 = geom[3].SolutionStepData().Data(step);
                ///*copying this data in the position of the vector we are interested in
                for (unsigned int j = 0; j < step_data_size; j++)
                {
                    new_step_data[j] = 0.25 * (step_data1[j] + step_data2[j] + step_data3[j] + step_data4[j]);
                }
            }

            pnode->FastGetSolutionStepValue(PARTITION_INDEX) = mrComm.MyPID();

            center_nodes.push_back(pnode);

            return new_id;
        }

        void AssignIdToCenterNode(ModelPart::NodesContainerType& center_nodes)
        {
            int total_number_of_nodes_nocenter = 1;
            int local_number_of_center_nodes = center_nodes.size();

            int rank = mrComm.MyPID();

            int number_of_nodes_before = -1;
            mrComm.ScanSum(&local_number_of_center_nodes, &number_of_nodes_before, 1);
            number_of_nodes_before = number_of_nodes_before - local_number_of_center_nodes;
            if (number_of_nodes_before < 0)
                KRATOS_ERROR(std::logic_error, "problem with scan sum ... giving a negative number of nodes before", "");

            //find our the total number of nodes (excluding center nodes)
            int number_of_own_nodes_nocenter = 0;
            for (ModelPart::NodesContainerType::iterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); it++)
                if (it->FastGetSolutionStepValue(PARTITION_INDEX) == rank)
                    number_of_own_nodes_nocenter++;
            number_of_own_nodes_nocenter -= local_number_of_center_nodes;
            mrComm.SumAll(&number_of_own_nodes_nocenter, &total_number_of_nodes_nocenter, 1);

            int counter = total_number_of_nodes_nocenter + number_of_nodes_before + 1;
            for (ModelPart::NodesContainerType::iterator it = center_nodes.begin(); it != center_nodes.end(); it++)
                it->SetId(counter++);

        }


        /// ************************************************************************************************
        /// ************************************************************************************************
        /// as the name says...

        void Erase_Old_Condition_And_Create_New_Triangle_Conditions(
                ModelPart& this_model_part,
                const boost::shared_ptr<Epetra_FECrsMatrix> p_edge_ids,
                PointerVector< Condition >& New_Conditions,
                bool interpolate_internal_variables
                )
        {
            KRATOS_TRY

            if (this_model_part.Conditions().size() != 0)
            {
                boost::numeric::ublas::matrix<int> new_conectivity;

                int total_existing_Conditions = -1;
                int local_existing_Conditions = this_model_part.Conditions().size();
                mrComm.SumAll(&local_existing_Conditions, &total_existing_Conditions, 1);

                unsigned int to_be_deleted = 0;
                unsigned int large_id = total_existing_Conditions * 4;
                bool create_Condition = false;
                int edge_ids[3];
                int t[12];
                int nel = 0;
                int splitted_edges = 0;
                int nint = 0;
                int aux[6];

                ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();
                PointerVector< Condition > Old_Conditions;

                unsigned int current_id = (this_model_part.Conditions().end() - 1)->Id() + 1;
                for (ModelPart::ConditionsContainerType::iterator it = this_model_part.Conditions().begin();
                        it != this_model_part.Conditions().end(); ++it)
                {
                    for (unsigned int i = 0; i < 12; i++)
                    {
                        t[i] = -1;
                    }
                    Condition::GeometryType& geom = it->GetGeometry();
                    Calculate_Triangle_Edges(geom, p_edge_ids, edge_ids, aux);

                    ///* crea las nuevas conectividades
                    create_Condition = Split_Triangle(edge_ids, t, &nel, &splitted_edges, &nint);

                    ///* crea los nuevos Conditionos
                    if (create_Condition == true)
                    {
                        to_be_deleted++;
                        for (int i = 0; i < nel; i++)
                        {
                            int i0, i1, i2;
                            TriangleGetNewConnectivityGID(i, t, aux, &i0, &i1, &i2);

//                            cout << i0 << " " << i1 << " " << i2 << endl;

                            Triangle3D3<Node < 3 > > geom(
                                    this_model_part.Nodes()(i0),
                                    this_model_part.Nodes()(i1),
                                    this_model_part.Nodes()(i2)
                                    );
                            KRATOS_WATCH("ln1304");

                            Condition::Pointer p_Condition;
                            p_Condition = it->Create(current_id, geom, it->pGetProperties());
                            p_Condition->Initialize();
                            p_Condition->InitializeSolutionStep(rCurrentProcessInfo);
                            p_Condition->FinalizeSolutionStep(rCurrentProcessInfo);

                            /// setting the internal variables in the child elem
                            //                        if (interpolate_internal_variables == true)
                            //                            InterpolateInternalVariables(nel, *it.base(), p_Condition, rCurrentProcessInfo);

                            // Transfer Conditional variables
                            p_Condition->Data() = it->Data();
                            //p_Condition->SetValue(REFINEMENT_LEVEL, 1);
                            New_Conditions.push_back(p_Condition);
                            current_id++;

                        }
                        it->SetId(large_id);
                        large_id++;
                    }

                }



                ///* all of the Conditions to be erased are at the end
                this_model_part.Conditions().Sort();

                ///*now remove all of the "old" Conditions
                this_model_part.Conditions().erase(this_model_part.Conditions().end() - to_be_deleted, this_model_part.Conditions().end());

                this_model_part.Conditions().reserve(this_model_part.Conditions().size() + New_Conditions.size());
                ///* adding news Conditions to the model part
                for (PointerVector< Condition >::iterator it_new = New_Conditions.begin(); it_new != New_Conditions.end(); it_new++)
                {
                    this_model_part.Conditions().push_back(*(it_new.base()));
                }

                int number_of_Conditions_before = -1;
                int number_of_own_Conditions = this_model_part.Conditions().size();
                mrComm.ScanSum(&number_of_own_Conditions, &number_of_Conditions_before, 1);
                number_of_Conditions_before = number_of_Conditions_before - number_of_own_Conditions;
                if (number_of_Conditions_before < 0)
                    KRATOS_ERROR(std::logic_error, "problem with scan sum ... giving a negative number of nodes before", "");

                int start_elem_id = number_of_Conditions_before + 1;
                //            int end_elem_id = number_of_Conditions_before + number_of_own_Conditions;

                int id_counter = start_elem_id;
                for (ModelPart::ConditionsContainerType::iterator it = this_model_part.Conditions().begin();
                        it != this_model_part.Conditions().end(); ++it)
                {
                    it->SetId(id_counter++);
                }

                //                //find our the total number of nodes
                //                int number_of_local_nodes = this_model_part.GetCommunicator().LocalMesh().Nodes().size();
                //            int total_number_of_nodes = -1;
                //            mrComm.SumAll(&number_of_local_nodes, &total_number_of_nodes, 1);

            }

            KRATOS_CATCH("")
        }

        Node<3>::Pointer AuxCreateNewNode(ModelPart& r_model_part, int Id, double x, double y, double z)
        {
            //create a new node
            Node<3>::Pointer p_new_node = Node<3>::Pointer(new Node<3>(Id, x, y, z));

            // Giving model part's variables list to the node
            p_new_node->SetSolutionStepVariablesList(&(r_model_part.GetNodalSolutionStepVariablesList()));

            p_new_node->SetBufferSize(r_model_part.NodesBegin()->GetBufferSize());
 
            return p_new_node;
        }



    };








} // namespace Kratos.

#endif // KRATOS_TRILINOS_LOCAL_REFINE_SIMPLEX_MESH  defined 


