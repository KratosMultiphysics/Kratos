//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2010-05-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_LOCAL_REFINE_TETRAHEDRA_MESH)
#define  KRATOS_LOCAL_REFINE_TETRAHEDRA_MESH


#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"
#include <boost/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>


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
#include "containers/array_1d.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_elements_neighbours_process.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "utilities/split_tetrahedra.c"
#include "utilities/split_triangle.c"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "processes/node_erase_process.h" 
#include "spatial_containers/spatial_containers.h"


namespace Kratos
{

    class Local_Refine_Tetrahedra_Mesh
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

        Local_Refine_Tetrahedra_Mesh(ModelPart& model_part) : mr_model_part(model_part)
        {
        }

        ~Local_Refine_Tetrahedra_Mesh()
        {
        }

	///This function performs a local refinement of a given tetrahedral mesh.
	///the "skin" of the mesh is also refined accordingly
	///The resulting mesh is guaranteed to be conformant, and (in principle) termination of the algorithm is guaranteed
	///to identify an element to be splitted do:
	///elementpointer->SetValue(SPLIT_ELEMENT,true)
	///all of the internal variables are interpolated for the newly created nodes. 
	///if a degree of freedom is fixed at both ends of a given edge, than the new node created on that edge is also fixed
	///@param refine_on_reference the interpolation of the variables is performed on the undeformed domain (requires DISPLACEMENT)
	///@param interpolate_internal_variables flag to specify if constitutive law variables should be interpolated or not
	///WARNING: nodal neighbours are assumed in this function and need to be updated on exit.
        void Local_Refine_Mesh(bool refine_on_reference, bool interpolate_internal_variables)
        {

            KRATOS_TRY
            
	    if(refine_on_reference==true)
	      if(!(mr_model_part.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT)) )
		  KRATOS_ERROR(std::logic_error,"DISPLACEMENT Variable is not in the model part -- needed if refine_on_reference = true","")

            boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;
            boost::numeric::ublas::vector<int> List_New_Nodes;
            compressed_matrix<int> Coord;
            WeakPointerVector< Node < 3 > > new_nodes;
            ModelPart& this_model_part = mr_model_part;

            ElementsArrayType& rElements = this_model_part.Elements();
            ElementsArrayType::iterator it_begin = rElements.ptr_begin();
            ElementsArrayType::iterator it_end = rElements.ptr_end();
            PointerVector< Element > New_Elements;
            New_Elements.reserve(20);



            ///WARNING = La numeracion local viene de 0 a 9

            if (refine_on_reference == true)
            {
                for (ModelPart::NodesContainerType::iterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); it++)
                {
                    it->X() = it->X0();
                    it->Y() = it->Y0();
                    it->Z() = it->Z0();
                }
            }
            // KRATOS_WATCH("line 106");
            CSR_Row_Matrix(this_model_part, Coord);
            // KRATOS_WATCH("line 107");
            Search_Edge_To_Be_Refined(this_model_part, Coord);
            // KRATOS_WATCH("line 108");
            Create_List_Of_New_Nodes(this_model_part, Coord, List_New_Nodes, Position_Node);
            // KRATOS_WATCH("line 110");
            Calculate_Coordinate_And_Insert_New_Nodes(this_model_part, new_nodes, Position_Node, List_New_Nodes);
            // KRATOS_WATCH("line 111");
            // KRATOS_WATCH(List_New_Nodes);
            Erase_Old_Element_And_Create_New_Tetra_Element(this_model_part, Coord, New_Elements);
            Erase_Old_Conditions_And_Create_New(this_model_part, Coord);
            Renumbering_Elements_And_Nodes(this_model_part);


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



            KRATOS_CATCH("")

        }


        ///************************************************************************************************
        ///************************************************************************************************

        void CSR_Row_Matrix(ModelPart& this_model_part, compressed_matrix<int>& Coord)
        {
            NodesArrayType& pNodes = this_model_part.Nodes();
            Coord.resize(pNodes.size(), pNodes.size());
            NodesArrayType::iterator i_begin = pNodes.ptr_begin();
            NodesArrayType::iterator i_end = pNodes.ptr_end();

	    std::vector<unsigned int> aux(10000);
	    
            for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i)
            {
                int index_i = i->Id() - 1;
                WeakPointerVector< Node < 3 > >& neighb_nodes = i->GetValue(NEIGHBOUR_NODES);
                
                unsigned int active = 0;
                for (WeakPointerVector< Node < 3 > >::iterator inode = neighb_nodes.begin();
                        inode != neighb_nodes.end(); inode++)
                {
                    int index_j = inode->Id() - 1;
                    if (index_j > index_i)
                    {
                        aux[active] = index_j;
                        active++;
                    }
                }
                std::sort(aux.begin(), aux.begin() + active);

                for (unsigned int k = 0; k < active; k++)
                {
                    Coord.push_back(index_i, aux[k], -1);

                }

            }

        }

        ///************************************************************************************************
        ///************************************************************************************************

        void Search_Edge_To_Be_Refined(ModelPart& this_model_part, compressed_matrix<int>& Coord)
        {

            ElementsArrayType& rElements = this_model_part.Elements();

            ElementsArrayType::iterator it_begin = rElements.ptr_begin(); //+element_partition[k];
            ElementsArrayType::iterator it_end = rElements.ptr_end(); //+element_partition[k+1];
            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                /*        unsigned int level = it->GetValue(REFINEMENT_LEVEL);
                         if( level > 0 )*/
                if (it->GetValue(SPLIT_ELEMENT) == true)
                {
                    Element::GeometryType& geom = it->GetGeometry(); // Nodos del elemento
                    for (unsigned int i = 0; i < geom.size(); i++)
                    {
                        int index_i = geom[i].Id() - 1;
                        for (unsigned int j = 0; j < geom.size(); j++)
                        {
                            int index_j = geom[j].Id() - 1;
                            if (index_j > index_i)
                            {
                                Coord(index_i, index_j) = -2;
                            }
                        }
                    }
                }

            }
        }

        ///************************************************************************************************
        ///************************************************************************************************

        void Create_List_Of_New_Nodes(ModelPart& this_model_part, compressed_matrix<int>& Coord, boost::numeric::ublas::vector<int> &List_New_Nodes,
                boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node)
        {

            unsigned int number_of_new_nodes = 0;

            NodesArrayType& pNodes = this_model_part.Nodes();

            typedef compressed_matrix<int>::iterator1 i1_t;
            typedef compressed_matrix<int>::iterator2 i2_t;

            ///*WARNING
            for (i1_t i1 = Coord.begin1(); i1 != Coord.end1(); ++i1)
            {
                for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
                {
                    if (Coord(i2.index1(), i2.index2()) == -2)
                    {
                        number_of_new_nodes++;
                    }
                }
            }


            ///* New Id de los Nodos
            List_New_Nodes.resize(number_of_new_nodes);
            int total_node = pNodes.size();
            for (unsigned int i = 0; i < number_of_new_nodes; i++)
            {
                List_New_Nodes[i] = total_node + i + 1;
            }


            ///* setting edges -2 to the new id of the new node
            ///* WARNING
            Position_Node.resize(number_of_new_nodes);
            unsigned int index = 0;
            for (i1_t i1 = Coord.begin1(); i1 != Coord.end1(); ++i1)
            {
                for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
                {
                    if (Coord(i2.index1(), i2.index2()) == -2)
                    {
                        Coord(i2.index1(), i2.index2()) = List_New_Nodes[index];
                        Position_Node[index][0] = i2.index1() + 1;
                        Position_Node[index][1] = i2.index2() + 1;
                        index++;
                    }
                }
            }

        }


        ///************************************************************************************************
        ///************************************************************************************************

        void Calculate_Coordinate_And_Insert_New_Nodes(ModelPart& this_model_part,
                WeakPointerVector< Node < 3 > >& new_nodes,
                const boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node,
                const boost::numeric::ublas::vector<int> &List_New_Nodes)
        {

            array_1d<double, 3 > Coord_Node_1;
            array_1d<double, 3 > Coord_Node_2;
            boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
            Coordinate_New_Node.resize(Position_Node.size());
            unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
            Node < 3 > ::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();

            new_nodes.reserve(10);


            for (unsigned int i = 0; i < Position_Node.size(); i++)
            {

                /// calculating the coordinate of the news nodes
                const int& node_i = Position_Node[i][0];
                const int& node_j = Position_Node[i][1];
                ModelPart::NodesContainerType::iterator it_node1 = this_model_part.Nodes().find(node_i);
                std::size_t pos1 = it_node1 - this_model_part.NodesBegin();
                noalias(Coord_Node_1) = it_node1->Coordinates();
                ModelPart::NodesContainerType::iterator it_node2 = this_model_part.Nodes().find(node_j);
                std::size_t pos2 = it_node2 - this_model_part.NodesBegin();
                noalias(Coord_Node_2) = it_node2->Coordinates();
                noalias(Coordinate_New_Node[i]) = 0.50 * (Coord_Node_1 + Coord_Node_2);

                //      KRATOS_WATCH("line 294");

                /// inserting the news node in the model part
                Node < 3 > ::Pointer pnode = this_model_part.CreateNewNode(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]);
                pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());

                it_node1 = this_model_part.NodesBegin() + pos1;
                it_node2 = this_model_part.NodesBegin() + pos2;

                pnode->GetValue(FATHER_NODES).resize(0);
                pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node1.base() ) );
                pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node2.base() ) );

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
                //       KRATOS_WATCH("line 318");

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

                //      KRATOS_WATCH("line 334");

                /// WARNING =  only for reactions;
                /*      const double zero = 0.00;
                      for(Node<3>::DofsContainerType::iterator iii = pnode->GetDofs().begin();    iii != pnode->GetDofs().end(); iii++)
                       {
                          if(pnode->IsFixed(iii->GetVariable())==false)
                                 {
                                    iii->GetSolutionStepReactionValue() = zero;
                                   }
                       }  */

                //     KRATOS_WATCH("line 346");
                new_nodes.push_back(pnode);
            }
        }

        void Renumbering_Elements_And_Nodes(ModelPart& this_model_part)
        {

            unsigned int id_node = 1;
            unsigned int id_elem = 1;
            unsigned int id_cond = 1;
            NodesArrayType& pNodes = this_model_part.Nodes();
            NodesArrayType::iterator i_begin = pNodes.ptr_begin();
            NodesArrayType::iterator i_end = pNodes.ptr_end();
            //ProcessInfo& rCurrentProcessInfo  = this_model_part.GetProcessInfo();


            for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i)
            {
                if (i->Id() != id_node)
                { //std::cout<< "Setting Id of Node  " << i->Id() << " by " <<  id_node << std::endl;
                    i->SetId(id_node);
                }
                id_node++;
            }

            ElementsArrayType& rElements = this_model_part.Elements();
            ElementsArrayType::iterator it_begin = rElements.ptr_begin();
            ElementsArrayType::iterator it_end = rElements.ptr_end();

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                if (it->Id() != id_elem)
                {
                    it->SetId(id_elem);
                }
                id_elem++;
            }

            for (ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); ++it)
            {
                if (it->Id() != id_cond)
                {
                    it->SetId(id_cond);
                }
                id_cond++;
            }

        }

        inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
        {
            partitions.resize(number_of_threads + 1);
            int partition_size = number_of_rows / number_of_threads;
            partitions[0] = 0;
            partitions[number_of_threads] = number_of_rows;
            for (unsigned int i = 1; i < number_of_threads; i++)
                partitions[i] = partitions[i - 1] + partition_size;
        }

        void Erase_Old_Element_And_Create_New_Tetra_Element(
                ModelPart& this_model_part,
                const compressed_matrix<int>& Coord,
                PointerVector< Element >& New_Elements
                )
        {

            boost::numeric::ublas::matrix<int> new_conectivity;
            ElementsArrayType& rElements = this_model_part.Elements();
            ElementsArrayType::iterator it_begin = rElements.ptr_begin();
            ElementsArrayType::iterator it_end = rElements.ptr_end();
            Element const rReferenceElement;
            unsigned int to_be_deleted = 0;
            unsigned int large_id = (rElements.end() - 1)->Id() * 15;
            unsigned int current_id = (rElements.end() - 1)->Id() + 1;
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

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                //KRATOS_WATCH(it->Id())
                Element::GeometryType& geom = it->GetGeometry();

                int index_0 = geom[0].Id() - 1;
                int index_1 = geom[1].Id() - 1;
                int index_2 = geom[2].Id() - 1;
                int index_3 = geom[3].Id() - 1;

                //put the global ids in aux
                array_1d<int, 11 > aux;
                aux[0] = geom[0].Id();
                aux[1] = geom[1].Id();
                aux[2] = geom[2].Id();
                aux[3] = geom[3].Id();

                if (index_0 > index_1)
                    aux[4] = Coord(index_1, index_0);
                else
                    aux[4] = Coord(index_0, index_1);


                if (index_0 > index_2)
                    aux[5] = Coord(index_2, index_0);
                else
                    aux[5] = Coord(index_0, index_2);


                if (index_0 > index_3)
                    aux[6] = Coord(index_3, index_0);
                else
                    aux[6] = Coord(index_0, index_3);


                if (index_1 > index_2)
                    aux[7] = Coord(index_2, index_1);
                else
                    aux[7] = Coord(index_1, index_2);


                if (index_1 > index_3)
                    aux[8] = Coord(index_3, index_1);
                else
                    aux[8] = Coord(index_1, index_3);


                if (index_2 > index_3)
                    aux[9] = Coord(index_3, index_2);
                else
                    aux[9] = Coord(index_2, index_3);

                // KRATOS_WATCH(   it->Id() );
                // KRATOS_WATCH(index_0);
                // KRATOS_WATCH(index_1);
                // KRATOS_WATCH(index_2);
                // KRATOS_WATCH(index_3);
                // KRATOS_WATCH(   aux );

                ///**************************************************

                //edge 01
                if (aux[4] < 0)
                    if (index_0 > index_1) edge_ids[0] = 0;
                    else edge_ids[0] = 1;
                else
                    edge_ids[0] = 4;

                //edge 02
                if (aux[5] < 0)
                    if (index_0 > index_2) edge_ids[1] = 0;
                    else edge_ids[1] = 2;
                else
                    edge_ids[1] = 5;

                //edge 03
                if (aux[6] < 0)
                    if (index_0 > index_3) edge_ids[2] = 0;
                    else edge_ids[2] = 3;
                else
                    edge_ids[2] = 6;

                //edge 12
                if (aux[7] < 0)
                    if (index_1 > index_2) edge_ids[3] = 1;
                    else edge_ids[3] = 2;
                else
                    edge_ids[3] = 7;

                //edge 13
                if (aux[8] < 0)
                    if (index_1 > index_3) edge_ids[4] = 1;
                    else edge_ids[4] = 3;
                else
                    edge_ids[4] = 8;

                //edge 23
                if (aux[9] < 0)
                    if (index_2 > index_3) edge_ids[5] = 2;
                    else edge_ids[5] = 3;
                else
                    edge_ids[5] = 9;

                //do the split
                //           for (unsigned int i = 0; i<6; i++) { std::cout<< "edges[" << i << "] = " << edge_ids[i] << std::endl;}
                //           KRATOS_WATCH("-------------------------------------------" )


                create_element = Split(edge_ids, t, &nel, &splitted_edges, &internal_node);

                if (internal_node == 1)
                {
                    // 	    std::cout << "creating internal node" << std::endl;
                    //generate new internal node
                    aux[10] = CreateCenterNode(geom, this_model_part);
		    
		    bool verified = false;
		    for(int iii=0; iii<nel*4; iii++)
		      if(t[iii] == 10)
			verified = true;
		      
		    if(verified == false)
		    {
		      KRATOS_WATCH(nel);
			for(int iii=0; iii<nel*4; iii++)
			  std::cout << t[iii] << std::endl;
			
		      KRATOS_ERROR(std::logic_error,"internal node is created but not used","");
		    }
		    
                }

                /*		KRATOS_ERROR(std::logic_error,"case not handled","");   */
                /*KRATOS_WATCH(splitted_edges);
                KRATOS_WATCH(internal_node);*/
                //           KRATOS_WATCH(aux)
                //KRATOS_WATCH(splitted_edges)

                //           for (unsigned int i = 0; i<56; i++) { std::cout<< "t[" << i << "] = " << t[i] << std::endl;}
                //           KRATOS_WATCH("-------------------------------------------" )


                if (create_element == true)
                {

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
                        Element::Pointer p_element;
                        p_element = it->Create(current_id, geom, it->pGetProperties());
                        New_Elements.push_back(p_element);

                        //pass the REFINEMENT_LEVEL and ensure that no further refinement is performed
                        p_element->GetValue(REFINEMENT_LEVEL) = it->GetValue(REFINEMENT_LEVEL);
                        p_element->GetValue(SPLIT_ELEMENT) = false;

                        current_id++;

                    }
                    it->SetId(large_id);
                    large_id++;
                }

                for (unsigned int i = 0; i < 32; i++)
                {
                    t[i] = -1;
                }
            }
            
            ///* all of the elements to be erased are at the end
            rElements.Sort();

            ///*now remove all of the "old" elements
            rElements.erase(this_model_part.Elements().end() - to_be_deleted, this_model_part.Elements().end());
	    
	    unsigned int total_size = this_model_part.Elements().size()+ New_Elements.size();
	    this_model_part.Elements().reserve(total_size);

            ///* adding news elements to the model part
            for (PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new != New_Elements.end(); it_new++)
            {
                it_new->Initialize();
                it_new->InitializeSolutionStep(rCurrentProcessInfo);
                it_new->FinalizeSolutionStep(rCurrentProcessInfo);
                rElements.push_back(*(it_new.base()));
            }
            
//             //renumber
// 	    unsigned int my_index = 1;
// 	    for(ModelPart::ElementsContainerType::iterator it=this_model_part.ElementsBegin(); it!=this_model_part.ElementsEnd(); it++)
// 	      it->SetId(my_index++);

            

        }

        unsigned int CreateCenterNode(Geometry<Node < 3 > >& geom, ModelPart& model_part)
        {
            //determine a new unique id
            unsigned int new_id = (model_part.NodesEnd() - 1)->Id() + 1;

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

            /// WARNING =  only for reactions; (commented as the new node is free on all dofs)
            /*      const double zero = 0.00;
                  for(Node<3>::DofsContainerType::iterator iii = pnode->GetDofs().begin();    iii != pnode->GetDofs().end(); iii++)
                   {
                      if(pnode->IsFixed(iii->GetVariable())==false)
                             {
                                iii->GetSolutionStepReactionValue() = zero;
                               }
                   } */

            /*KRATOS_WATCH(*(model_part.NodesEnd()-1));
                 model_part.Nodes().push_back(pnode);*/

            //      KRATOS_WATCH(*(model_part.NodesEnd()-1));

            return new_id;
        }

        void Erase_Old_Conditions_And_Create_New(
                ModelPart& this_model_part,
                const compressed_matrix<int>& Coord
                )

        {
            PointerVector< Condition > New_Conditions;

	    ConditionsArrayType& rConditions = this_model_part.Conditions();
            ConditionsArrayType::iterator it_begin = rConditions.ptr_begin();
            ConditionsArrayType::iterator it_end = rConditions.ptr_end();
            unsigned int to_be_deleted = 0;
            unsigned int large_id = (rConditions.end() - 1)->Id() * 7;
            int  edge_ids[3];       
	    int  t[12];       
	    int  nel             = 0;
	    int  splitted_edges  = 0; 
	    int  nint            = 0;
	    array_1d<int,6> aux;

            ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();


            unsigned int current_id = (rConditions.end() - 1)->Id() + 1;
            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it)
            {

                Condition::GeometryType& geom = it->GetGeometry();

                if (geom.size() == 3)
                {
		    Calculate_Edges(geom, Coord, edge_ids, aux);
		    
                    ///* crea las nuevas conectividades
		    bool create_condition =  Split_Triangle(edge_ids, t, &nel, &splitted_edges, &nint); 

		    ///* crea los nuevos elementos           
		    if(create_condition==true)
		    {
		      to_be_deleted++;  
		      for(int i=0; i<nel; i++)
			{
			  
			  unsigned int base = i * 3;
			  unsigned int i0   = aux[t[base]];
			  unsigned int i1   = aux[t[base+1]];
			  unsigned int i2   = aux[t[base+2]];

			  Triangle3D3<Node<3> > newgeom(
							this_model_part.Nodes()(i0),
							this_model_part.Nodes()(i1),
							this_model_part.Nodes()(i2) 
							);

			      
			      Condition::Pointer pcond = it->Create(current_id, newgeom, it->pGetProperties());
			      
			      const int is_slip = it->GetValue(IS_STRUCTURE);
			      pcond->GetValue(IS_STRUCTURE) = is_slip;
		
			      New_Conditions.push_back(pcond);
			      current_id++; 
						
		      }
			it->SetId(large_id);   
			large_id++;
		    }
		}
	    }
            	    
            ///* all of the elements to be erased are at the end
            this_model_part.Conditions().Sort();

            ///*now remove all of the "old" elements
            this_model_part.Conditions().erase(this_model_part.Conditions().end() - to_be_deleted, this_model_part.Conditions().end());

	    unsigned int total_size = this_model_part.Conditions().size()+ New_Conditions.size();
	    this_model_part.Conditions().reserve(total_size);
	    

            ///* adding news elements to the model part
            //Vector temp;
            for (PointerVector< Condition >::iterator it_new = New_Conditions.begin(); it_new != New_Conditions.end(); it_new++)
            {

                it_new->Initialize();
                it_new->InitializeSolutionStep(rCurrentProcessInfo);
                 it_new->FinalizeSolutionStep(rCurrentProcessInfo);
                this_model_part.Conditions().push_back(*(it_new.base()));
            }
            
            //renumber
	    unsigned int my_index = 1;
	    for(ModelPart::ConditionsContainerType::iterator it=this_model_part.ConditionsBegin(); it!=this_model_part.ConditionsEnd(); it++)
	      it->SetId(my_index++);

        }

	/// ************************************************************************************************          
	/// ************************************************************************************************ 
	void  Calculate_Edges(Element::GeometryType& geom,
			      const compressed_matrix<int>& Coord,
			      int*  edge_ids,
			      array_1d<int,6>& aux
			      )
	{           
		    int index_0 = geom[0].Id()-1;
		    int index_1 = geom[1].Id()-1;
		    int index_2 = geom[2].Id()-1; 
		    
		    aux[0] = geom[0].Id();
		    aux[1] = geom[1].Id();
		    aux[2] = geom[2].Id();
		    //------------------------------------------------------------------------- 
		    if(index_0 > index_1)
		      aux[3] = Coord(index_1, index_0);
		    else
		      aux[3] = Coord(index_0, index_1); 

	  
		    if(index_1 > index_2)
		      aux[4] = Coord(index_2, index_1);
		    else 
		      aux[4] = Coord(index_1, index_2 );

		    
		    if(index_2 > index_0)
		      aux[5] = Coord(index_0, index_2);
		    else 
		      aux[5] = Coord(index_2, index_0 );           
		    //-------------------------------------------------------------------------
		    
		    //edge 01
		    if(aux[3] < 0)
			if(index_0 > index_1) edge_ids[0] = 0;                           
			else edge_ids[0] = 1;
		    else
			edge_ids[0] = 3;

		    //edge 12
		    if(aux[4] < 0)
			if(index_1 > index_2) edge_ids[1] = 1;
			else edge_ids[1] = 2;
		    else
			edge_ids[1] = 4;

		    //edge 20
		    if(aux[5] < 0)
			if(index_2 > index_0) edge_ids[2] = 2;
			else edge_ids[2] = 0;
		    else
			edge_ids[2] = 5;
	  
	}      

    protected:
        ModelPart& mr_model_part;

    };
}



#endif