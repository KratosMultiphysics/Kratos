//   Project Name:        KratosMeshingApplication
//   Last Modified by:    $Author:  VMataix $
//   Date:                $Date: 20-04-2016 $
//   Revision:            $Revision:    1.0 $
//

#if !defined(KRATOS_LOCAL_REFINE_PRISM_MESH)
#define  KRATOS_LOCAL_REFINE_PRISM_MESH

#ifdef _OPENMP
#include <omp.h>
#endif

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

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

#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"

#include "utilities/split_prism.hpp"
#include "geometries/prism_3d_6.h"
#include "processes/node_erase_process.h"

namespace Kratos
{

class Local_Refine_Prism_Mesh
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

    /************************************* CONSTRUCTOR *********************************/
    /***********************************************************************************/
    
    Local_Refine_Prism_Mesh(ModelPart& model_part) : mr_model_part(model_part)
    {

    }

    /************************************* DESTRUCTOR **********************************/
    /***********************************************************************************/
    
    ~Local_Refine_Prism_Mesh()
    {
      
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * Refine the mesh locally, call all the commands necessaries to compute the remeshing
    * @param refine_on_reference: Boolean that defines if refine or not the mesh according to the reference
    * @param interpolate_internal_variables: Boolean that defines if to interpolate or not the internal variables
    */
    
    void Local_Refine_Mesh(
            bool refine_on_reference,
            bool interpolate_internal_variables
    )
    {
        KRATOS_TRY;

        if (refine_on_reference == true)
        {
            if (!(mr_model_part.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT)))
            {
                KRATOS_THROW_ERROR(std::logic_error, "DISPLACEMENT Variable is not in the model part -- needed if refine_on_reference = true", "");
            }
        }

        compressed_matrix<int> Coord;                                              // The matrix that stores all the index of the geometry
        boost::numeric::ublas::vector<int> List_New_Nodes;                         // The news nodes
        boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;           // Edges where are the news nodes
        boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node; // The coordinate of the new nodes

        PointerVector< Element > New_Elements;
	
        if (refine_on_reference == true)
        {
            for (ModelPart::NodesContainerType::iterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); it++)
            {
                it->X() = it->X0();
                it->Y() = it->Y0();
                it->Z() = it->Z0();
            }
        }

        /* Calling all the functions necessaries to refine the mesh */
        CSR_Row_Matrix(mr_model_part, Coord);
        Search_Edge_To_Be_Refined(mr_model_part, Coord);
        Create_List_Of_New_Nodes(mr_model_part, Coord, List_New_Nodes, Position_Node);
        Calculate_Coordinate_And_Insert_New_Nodes(mr_model_part, Position_Node, List_New_Nodes);
        Erase_Old_Element_And_Create_New_Prism_Element(mr_model_part, Coord, New_Elements, interpolate_internal_variables);
        Erase_Old_Conditions_And_Create_New(mr_model_part, Coord);
        Renumering_Elements_And_Nodes(mr_model_part, New_Elements);

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

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * This function initialises the matrix Cord
    * @return Coord: The matrix that stores all the index of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */
    
    void CSR_Row_Matrix(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord
    )
    {
        NodesArrayType& pNodes = this_model_part.Nodes();
        NodesArrayType::iterator it_begin = pNodes.ptr_begin();
        NodesArrayType::iterator it_end   = pNodes.ptr_end();

        Coord.resize(pNodes.size(), pNodes.size(), false);

        for(NodesArrayType::iterator i = it_begin; i!=it_end; i++)
        {
            int index_i = i->Id() - 1; // WARNING: MESH MUST BE IN ORDER
            WeakPointerVector< Node < 3 > >& neighb_nodes = i->GetValue(NEIGHBOUR_NODES);

            std::vector<unsigned int> aux(neighb_nodes.size());
            unsigned int active = 0;
            for (WeakPointerVector< Node < 3 > >::iterator inode = neighb_nodes.begin(); inode != neighb_nodes.end(); inode++)
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

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * This functions looks for potential edges that could be refined
    * @return Coord: The matrix that stores all the index of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */
    
    void Search_Edge_To_Be_Refined(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord
    )
    {
        ElementsArrayType& rElements = this_model_part.Elements();
        ElementsArrayType::iterator it_begin = rElements.ptr_begin();
        ElementsArrayType::iterator it_end   = rElements.ptr_end();

        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
        {
            if (it->GetValue(SPLIT_ELEMENT) == true)
            {
                Element::GeometryType& geom = it->GetGeometry(); // Nodes of the element
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

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * It creates the list of new nodes
    * @return this_model_part: The model part of the model (it is the input too)
    * @return Coord: The matrix that stores all the index of the geometry
    * @return List_New_Nodes: List that contents the index of the new nodes to be created
    * @return Position_Node: The vector that contents the position in the edge of the new nodes
    */
    
    void Create_List_Of_New_Nodes(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord,
            boost::numeric::ublas::vector<int> &List_New_Nodes,
            boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node
    )
    {
        unsigned int number_of_new_nodes = 0;

        NodesArrayType& pNodes = this_model_part.Nodes();

        typedef compressed_matrix<int>::iterator1 i1_t;
        typedef compressed_matrix<int>::iterator2 i2_t;

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

        // New ID of the nodes
        List_New_Nodes.resize(number_of_new_nodes, false);
        int total_node = pNodes.size();
        for (unsigned int i = 0; i < number_of_new_nodes; i++)
        {
            List_New_Nodes[i] = total_node + i + 1;
        }

        // Setting edges -2 to the new id of the new node
        Position_Node.resize(number_of_new_nodes, false);
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

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * Computes the coordinates of the new nodes in the center of the edges's triangle.
    * Insert the news nodes in the model part and interopolate the variables
    * @param List_New_Nodes: List that contents the index of the new nodes to be created
    * @param Position_Node: The vector that contents the position in the edge of the new nodes
    * @return this_model_part: The model part of the model (it is the input too)
    */
    
    void Calculate_Coordinate_And_Insert_New_Nodes(
            ModelPart& this_model_part,
            const boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node,
            const boost::numeric::ublas::vector<int> &List_New_Nodes
    )
    {
        array_1d<double, 3 > Coord_Node_1;
        array_1d<double, 3 > Coord_Node_2;
        boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
        Coordinate_New_Node.resize(Position_Node.size(), false);
        unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
        Node < 3 > ::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();

        for (unsigned int i = 0; i < Position_Node.size(); i++)
        {
            /* Calculating the coordinate of the news nodes */
            const int& node_i = Position_Node[i][0];
            const int& node_j = Position_Node[i][1];
	    
            ModelPart::NodesContainerType::iterator it_node1 = this_model_part.Nodes().find(node_i);
            std::size_t pos1 = it_node1 - this_model_part.NodesBegin();
            noalias(Coord_Node_1) = it_node1->Coordinates();
	    
            ModelPart::NodesContainerType::iterator it_node2 = this_model_part.Nodes().find(node_j);
            std::size_t pos2 = it_node2 - this_model_part.NodesBegin();
            noalias(Coord_Node_2) = it_node2->Coordinates();
	    
            noalias(Coordinate_New_Node[i]) = 0.50 * (Coord_Node_1 + Coord_Node_2);
            
            /* Inserting the news node in the model part */
            Node < 3 >::Pointer pnode = Node < 3 >::Pointer(new Node < 3 >(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]));
            pnode->SetSolutionStepVariablesList( this_model_part.NodesBegin()->pGetVariablesList() );
            pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());

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
                {
                    (p_new_dof)->FixDof();
                }
                else
                {
                    (p_new_dof)->FreeDof();
                }
            }

            // Interpolating the data
            unsigned int buffer_size = pnode->GetBufferSize();
            for (unsigned int step = 0; step < buffer_size; step++)
            {
                double* new_step_data = pnode->SolutionStepData().Data(step);
                double* step_data1 = it_node1->SolutionStepData().Data(step);
                double* step_data2 = it_node2->SolutionStepData().Data(step);
		
                // Copying this data in the position of the vector we are interested in
                for (unsigned int j = 0; j < step_data_size; j++)
                {
                    new_step_data[j] = 0.5 * (step_data1[j] + step_data2[j]);
                }
            }

            this_model_part.Nodes().push_back(pnode);
        }
        
        this_model_part.Nodes().Sort();
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * Computes the coordinate of the baricenter node of the element (mean of the faces's baricenter)
    * Insert the news nodes in the center of elements and interopolate the variables.
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void Calculate_Coordinate_Center_Node_And_Insert_New_Nodes(ModelPart& this_model_part)
    {
        // Lower face
        array_1d<double, 3 > Coord_Node_1;
        array_1d<double, 3 > Coord_Node_2;
        array_1d<double, 3 > Coord_Node_3;
	
        // Upper face
        array_1d<double, 3 > Coord_Node_4;
        array_1d<double, 3 > Coord_Node_5;
        array_1d<double, 3 > Coord_Node_6;
	
        // Center
        array_1d<double, 3 > Coordinate_center_node;
	
        std::vector<int> node_center;
        NodesArrayType& pNodes = this_model_part.Nodes();
        int Id_Center = pNodes.size() + 1;
        ElementsArrayType& rElements = this_model_part.Elements();
        ElementsArrayType::iterator it_begin = rElements.ptr_begin(); 
        ElementsArrayType::iterator it_end = rElements.ptr_end(); 
        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
        {
            Element::GeometryType& geom = it->GetGeometry();
            noalias(Coord_Node_1) = geom[0].Coordinates();
            noalias(Coord_Node_2) = geom[1].Coordinates();
            noalias(Coord_Node_3) = geom[2].Coordinates();
            noalias(Coord_Node_4) = geom[3].Coordinates();
            noalias(Coord_Node_5) = geom[4].Coordinates();
            noalias(Coord_Node_6) = geom[5].Coordinates();

            unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
            Node < 3 > ::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();
            noalias(Coordinate_center_node) = 0.16666666666666666 * (Coord_Node_1 + Coord_Node_2 + Coord_Node_3 +
                                                                     Coord_Node_4 + Coord_Node_5 + Coord_Node_6);

            /* Inserting the new node in the model part */
            Node < 3 > ::Pointer pnode = this_model_part.CreateNewNode(Id_Center, Coordinate_center_node[0], Coordinate_center_node[1], Coordinate_center_node[2]);
            pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());

            pnode->X0() = 0.16666666666666666 * (geom[0].X0() + geom[1].X0() + geom[2].X0() + geom[3].X0() + geom[4].X0() + geom[5].X0());
            pnode->Y0() = 0.16666666666666666 * (geom[0].Y0() + geom[1].Y0() + geom[2].Y0() + geom[3].Y0() + geom[4].Y0() + geom[5].Y0());
            pnode->Z0() = 0.16666666666666666 * (geom[0].Z0() + geom[1].Z0() + geom[2].Z0() + geom[3].Z0() + geom[4].Z0() + geom[5].Z0());
            
            for (Node < 3 > ::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
            {
                Node < 3 > ::DofType& rDof = *iii;
                Node < 3 > ::DofType::Pointer p_new_dof = pnode->pAddDof(rDof);
                if (geom[0].IsFixed(iii->GetVariable()) == true && geom[1].IsFixed(iii->GetVariable()) == true && geom[2].IsFixed(iii->GetVariable()) == true && geom[3].IsFixed(iii->GetVariable()) == true 
                 && geom[4].IsFixed(iii->GetVariable()) == true && geom[5].IsFixed(iii->GetVariable()) == true)
                {
                    (p_new_dof)->FixDof();
                }
                else
                {
                    (p_new_dof)->FreeDof();
                }
            }

            /* Intepolating the data */
            unsigned int buffer_size = pnode->GetBufferSize();
            for (unsigned int step = 0; step < buffer_size; step++)
            {
                double* new_step_data = pnode->SolutionStepData().Data(step);
		
                // Lower face
                double* step_data1 = geom[0].SolutionStepData().Data(step);
                double* step_data2 = geom[1].SolutionStepData().Data(step);
                double* step_data3 = geom[2].SolutionStepData().Data(step);
		
                // Upper face
                double* step_data4 = geom[3].SolutionStepData().Data(step);
                double* step_data5 = geom[4].SolutionStepData().Data(step);
                double* step_data6 = geom[5].SolutionStepData().Data(step);
		
                // Copying this data in the position of the vector we are interested in
                for (unsigned int j = 0; j < step_data_size; j++)
                {
                    new_step_data[j] = 0.16666666666666666 * (step_data1[j] + step_data2[j] + step_data3[j] +
                                                              step_data4[j] + step_data5[j] + step_data6[j]);
                }
            }
            node_center.push_back(Id_Center);
            Id_Center++;
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * It erases the old elements and it creates the new ones
    * @param Coord: The coordinates of the element
    * @param New_Elements: The new elements created
    * @param interpolate_internal_variables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return this_model_part: The model part of the model (it is the input too)
    */
    
    void Erase_Old_Element_And_Create_New_Prism_Element(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& New_Elements,
            bool interpolate_internal_variables
    )
    {
        ElementsArrayType& rElements = this_model_part.Elements();
        ElementsArrayType::iterator it_begin = rElements.ptr_begin();
        ElementsArrayType::iterator it_end = rElements.ptr_end();
        unsigned int to_be_deleted = 0;
        unsigned int large_id = (rElements.end() - 1)->Id() * 10;
        bool create_element = false;
        int edge_ids[6];
        int t[24];
        int number_elem = 0;
        int splitted_edges = 0;
        int nint = 0;
        array_1d<int, 12 > aux;

        ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();

        std::cout << "****************** REFINING MESH ******************" << std::endl;
        std::cout << "OLD NUMBER ELEMENTS: " << rElements.size() << std::endl;

        unsigned int current_id = (rElements.end() - 1)->Id() + 1;
        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
        {
            for (unsigned int i = 0; i < 24; i++)
            {
                t[i] = -1;
            }
            
            Element::GeometryType& geom = it->GetGeometry();
            Calculate_Edges(geom, Coord, edge_ids, aux);

            // It creates the new conectivities
            create_element = Split_Prism(edge_ids, t, &number_elem, &splitted_edges, &nint);

            // It creates the new elements
            if (create_element == true)
            {
                to_be_deleted++;
                for (int i = 0; i < number_elem; i++)
                {
                    unsigned int base = i * 6;
                    unsigned int i0 = aux[t[base]];
                    unsigned int i1 = aux[t[base + 1]];
                    unsigned int i2 = aux[t[base + 2]];
                    unsigned int i3 = aux[t[base + 3]];
                    unsigned int i4 = aux[t[base + 4]];
                    unsigned int i5 = aux[t[base + 5]];

                    Prism3D6<Node < 3 > > geom(
                        this_model_part.Nodes()(i0),
                        this_model_part.Nodes()(i1),
                        this_model_part.Nodes()(i2),
                        this_model_part.Nodes()(i3),
                        this_model_part.Nodes()(i4),
                        this_model_part.Nodes()(i5)
                    );

                    Element::Pointer p_element;
                    p_element = it->Create(current_id, geom, it->pGetProperties());
                    p_element->Initialize();
                    p_element->InitializeSolutionStep(rCurrentProcessInfo);
                    p_element->FinalizeSolutionStep(rCurrentProcessInfo);

                    // Setting the internal variables in the child elem
                    if (interpolate_internal_variables == true)
                    {
                        InterpolateInteralVariables(number_elem, *it.base(), p_element, rCurrentProcessInfo);
                    }

                    // Transfer elemental variables
                    p_element->Data() = it->Data();
                    p_element->GetValue(SPLIT_ELEMENT) = false;
                    New_Elements.push_back(p_element);

                    current_id++;
                }
                it->SetId(large_id);
                large_id++;
            }
        }

        /* Adding news elements to the model part */
        for (PointerVector< Element >::iterator it_new = New_Elements.begin(); it_new != New_Elements.end(); it_new++)
        {
            rElements.push_back(*(it_new.base()));
        }

        /* All of the elements to be erased are at the end */
        rElements.Sort();

        /* Now remove all of the "old" elements */
        rElements.erase(this_model_part.Elements().end() - to_be_deleted, this_model_part.Elements().end());

        std::cout << "NEW NUMBER ELEMENTS: " << rElements.size() << std::endl;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * Remove the old conditions and creates new ones
    * @param Coord: The coordinates of the nodes of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */
    
    void Erase_Old_Conditions_And_Create_New(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord
    )
    {
        PointerVector< Condition > New_Conditions;

        ConditionsArrayType& rConditions = this_model_part.Conditions();

        if (rConditions.size() > 0)
        {
            ConditionsArrayType::iterator it_begin = rConditions.ptr_begin();
            ConditionsArrayType::iterator it_end   = rConditions.ptr_end();
            unsigned int to_be_deleted = 0;
            unsigned int large_id = (rConditions.end() - 1)->Id() * 7;

            ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();

            unsigned int current_id = (rConditions.end() - 1)->Id() + 1;
            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it)
            {

                Condition::GeometryType& geom = it->GetGeometry();

                if (geom.size() == 2)
                {
                    int index_0 = geom[0].Id() - 1;
                    int index_1 = geom[1].Id() - 1;
                    int new_id;

                    if (index_0 > index_1)
                    {
                        new_id = Coord(index_1, index_0);
                    }
                    else
                    {
                        new_id = Coord(index_0, index_1);
                    }

                    if (new_id > 0) // We need to create a new condition
                    {
                        to_be_deleted++;
                        Line3D2<Node < 3 > > newgeom1(
                            this_model_part.Nodes()(geom[0].Id()),
                            this_model_part.Nodes()(new_id)
                        );

                        Line3D2<Node < 3 > > newgeom2(
                            this_model_part.Nodes()(new_id),
                            this_model_part.Nodes()(geom[1].Id())
                        );

                        Condition::Pointer pcond1 = it->Create(current_id++, newgeom1, it->pGetProperties());
                        Condition::Pointer pcond2 = it->Create(current_id++, newgeom2, it->pGetProperties());

                        pcond1->Data() = it->Data();
                        pcond2->Data() = it->Data();

                        New_Conditions.push_back(pcond1);
                        New_Conditions.push_back(pcond2);

                        it->SetId(large_id);
                        large_id++;
                    }
                }
            }

            /* All of the elements to be erased are at the end */
            this_model_part.Conditions().Sort();

            /* Remove all of the "old" elements */
            this_model_part.Conditions().erase(this_model_part.Conditions().end() - to_be_deleted, this_model_part.Conditions().end());

            unsigned int total_size = this_model_part.Conditions().size() + New_Conditions.size();
            this_model_part.Conditions().reserve(total_size);

            /* Adding news elements to the model part */
            for (PointerVector< Condition >::iterator it_new = New_Conditions.begin(); it_new != New_Conditions.end(); it_new++)
            {
                it_new->Initialize();
                it_new->InitializeSolutionStep(rCurrentProcessInfo);
                it_new->FinalizeSolutionStep(rCurrentProcessInfo);
                this_model_part.Conditions().push_back(*(it_new.base()));
            }

            /* Renumber */
            unsigned int my_index = 1;
            for (ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); it++)
            {
                it->SetId(my_index++);
            }

        }

    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * It calculates the new edges of the new prisms, 
    * first it calculates the new edges correspondign to the lower face (as a triangle),
    * later it added to the upper face
    * @param geom: The prism element geometry
    * @param edge_ids: The ids of the edges
    * @return aux: The vector that includes the index of the new edges
    */
    
    void Calculate_Edges(
            Element::GeometryType& geom,
            const compressed_matrix<int>& Coord,
            int* edge_ids,
            array_1d<int, 12 > & aux
            )
    {
        // Lower face
        int index_0 = geom[0].Id() - 1;
        int index_1 = geom[1].Id() - 1;
        int index_2 = geom[2].Id() - 1;
	
        aux[0] = geom[0].Id();
        aux[1] = geom[1].Id();
        aux[2] = geom[2].Id();
	
        // Upper face
        int index_3 = geom[3].Id() - 1;
        int index_4 = geom[4].Id() - 1;
        int index_5 = geom[5].Id() - 1;

        aux[3] = geom[3].Id();
        aux[4] = geom[4].Id();
        aux[5] = geom[5].Id();
	
        //-------------------------------------------------------------------------
	
        // First node of the triangle face
        if (index_0 > index_1)
        {
            aux[6] = Coord(index_1, index_0);
            aux[9] = Coord(index_4, index_3);
        }
        else
        {
            aux[6] = Coord(index_0, index_1);
            aux[9] = Coord(index_3, index_4);
        }

        // Second node of the triangle face
        if (index_1 > index_2)
        {
            aux[7]  = Coord(index_2, index_1);
            aux[10] = Coord(index_5, index_4);
        }
        else
        {
            aux[7]  = Coord(index_1, index_2);
            aux[10] = Coord(index_4, index_5);
        }

	// Third node of the triangle face
        if (index_2 > index_0)
        {
            aux[8]  = Coord(index_0, index_2);
            aux[11] = Coord(index_3, index_5);
        }
        else
        {
            aux[8]  = Coord(index_2, index_0);
            aux[11] = Coord(index_5, index_3);
        }
	
        //-------------------------------------------------------------------------

        // Edge 01
        if (aux[6] < 0)
        {
            if (index_0 > index_1)
            {
                edge_ids[0] = 0;
                edge_ids[3] = 6;
            }
            else
            {
                edge_ids[0] = 1;
                edge_ids[3] = 7;
            }
        }
        else
        {
            edge_ids[0] = 3;
            edge_ids[3] = 9;
        }

        // Edge 12
        if (aux[7] < 0)
        {
            if (index_1 > index_2)
            {
                edge_ids[1] = 1;
                edge_ids[4] = 7;
            }
            else
            {
                edge_ids[1] = 2;
                edge_ids[4] = 8;
            }
        }
        else
        {
            edge_ids[1] = 4;
            edge_ids[4] = 10;
        }

        // Edge 20
        if (aux[8] < 0)
        {
            if (index_2 > index_0)
            {
                edge_ids[2] = 2;
                edge_ids[5] = 8;
            }
            else
            {
                edge_ids[2] = 0;
                edge_ids[5] = 6;
            }
        }
        else
        {
            edge_ids[2] = 5;
            edge_ids[5] = 11;
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * This process renumerates the elements and nodes  
    * @param New_Elements: Pointers to the new elements created
    * @return this_model_part: The model part of the model (it is the input too)
    */
    
    void Renumering_Elements_And_Nodes(
            ModelPart& this_model_part,
            PointerVector< Element >& New_Elements
    )
    {
        unsigned int id_node = 1;
        unsigned int id_elem = 1;
        NodesArrayType& pNodes = this_model_part.Nodes();
        NodesArrayType::iterator i_begin = pNodes.ptr_begin();
        NodesArrayType::iterator i_end   = pNodes.ptr_end();

        for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i)
        {
            if (i->Id() != id_node)
            {
                i->SetId(id_node);
            }
            id_node++;
        }

        ElementsArrayType& rElements = this_model_part.Elements();
        ElementsArrayType::iterator it_begin = rElements.ptr_begin();
        ElementsArrayType::iterator it_end   = rElements.ptr_end();

        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
        {
            if (it->Id() != id_elem)
            {
                it->SetId(id_elem);
            }
            id_elem++;
        }

    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * It creates a partition of the process between the different threads
    * @param number_of_threads: Number the threads considered in the computation
    * @param number_of_rows: 
    * @return partitions: The vector that contents the partitions corresponding to each thread
    */
    
    inline void CreatePartition(
      unsigned int number_of_threads, 
      const int number_of_rows, 
      vector<unsigned int>& partitions
    )
    {
        partitions.resize(number_of_threads + 1, false);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
        {
            partitions[i] = partitions[i - 1] + partition_size;
        }
    }


    /***********************************************************************************/
    /***********************************************************************************/

protected:
    ModelPart& mr_model_part;

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
    * Interpolates the internal variables
    * @param number_elem: Number of elements
    * @param father_elem: Father element (the original one)
    * @param child_elem: Child element (the new ones created)
    * @param rCurrentProcessInfo: The model part process info
    */
    
    void InterpolateInteralVariables(
            const int& number_elem,
            const Element::Pointer father_elem,
            Element::Pointer child_elem,
            ProcessInfo& rCurrentProcessInfo
            )
    {
        // NOTE: Right now there is not an interpolation at all, it just copying the values
        std::vector<Vector> values;
        father_elem->GetValueOnIntegrationPoints(INTERNAL_VARIABLES, values, rCurrentProcessInfo);
        child_elem->SetValueOnIntegrationPoints(INTERNAL_VARIABLES, values, rCurrentProcessInfo);
    }
};

} // namespace Kratos.

#endif // KRATOS_LOCAL_REFINE_PRISM_MESH  defined 


