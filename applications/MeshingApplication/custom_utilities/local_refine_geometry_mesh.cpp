// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//    Co-authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "local_refine_geometry_mesh.hpp"
#include <unordered_map>

namespace Kratos
{
  void LocalRefineGeometryMesh::LocalRefineMesh(
            bool refine_on_reference,
            bool interpolate_internal_variables
    )
    {
        KRATOS_TRY;

        if (refine_on_reference == true)
        {
            if (!(mModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT)))
            {
                KRATOS_THROW_ERROR(std::logic_error, "DISPLACEMENT Variable is not in the model part -- needed if refine_on_reference = true", "");
            }
        }

        compressed_matrix<int> Coord;                                              // The matrix that stores all the index of the geometry
        boost::numeric::ublas::vector<int> List_New_Nodes;                         // The news nodes
        boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;           // Edges where are the news nodes
        boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node; // The coordinate of the new nodes

        PointerVector< Element > New_Elements;
	New_Elements.reserve(20);

	// Initial renumber of nodes and elemetns
	unsigned int id = 1;
        for (ModelPart::NodesContainerType::iterator it = mModelPart.NodesBegin(); it != mModelPart.NodesEnd(); it++)
	{
	  it->SetId(id++);
	}

	id = 1;
        for (ModelPart::ElementsContainerType::iterator it = mModelPart.ElementsBegin(); it != mModelPart.ElementsEnd(); it++)
	{
	  it->SetId(id++);
	}

	id = 1;
        for (ModelPart::ConditionsContainerType::iterator it = mModelPart.ConditionsBegin(); it != mModelPart.ConditionsEnd(); it++)
	{
	  it->SetId(id++);
	}

        if (refine_on_reference == true)
        {
            for (ModelPart::NodesContainerType::iterator it = mModelPart.NodesBegin(); it != mModelPart.NodesEnd(); it++)
            {
                it->X() = it->X0();
                it->Y() = it->Y0();
                it->Z() = it->Z0();
            }
        }

        this->ResetFatherNodes(mModelPart);

        /* Calling all the functions necessaries to refine the mesh */
        CSRRowMatrix(mModelPart, Coord);

        SearchEdgeToBeRefined(mModelPart, Coord);

        CreateListOfNewNodes(mModelPart, Coord, List_New_Nodes, Position_Node);

        CalculateCoordinateAndInsertNewNodes(mModelPart, Position_Node, List_New_Nodes);

        EraseOldElementAndCreateNewElement(mModelPart, Coord, New_Elements, interpolate_internal_variables);

        EraseOldConditionsAndCreateNew(mModelPart, Coord);

        RenumeringElementsAndNodes(mModelPart, New_Elements);


        if (refine_on_reference == true)
        {
            for (ModelPart::NodesContainerType::iterator it = mModelPart.NodesBegin(); it != mModelPart.NodesEnd(); it++)
            {
                const array_1d<double, 3 > & disp = it->FastGetSolutionStepValue(DISPLACEMENT);
                it->X() = it->X0() + disp[0];
                it->Y() = it->Y0() + disp[1];
                it->Z() = it->Z0() + disp[2];
            }
        }

        UpdateSubModelPartNodes(mModelPart);

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void LocalRefineGeometryMesh::CSRRowMatrix(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord
    )
    {
        KRATOS_TRY;

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

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void LocalRefineGeometryMesh::SearchEdgeToBeRefined(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord
    )
    {
        KRATOS_TRY;

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

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void LocalRefineGeometryMesh::CreateListOfNewNodes(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord,
            boost::numeric::ublas::vector<int> &List_New_Nodes,
            boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node
    )
    {
        KRATOS_TRY;

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

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void LocalRefineGeometryMesh::CalculateCoordinateAndInsertNewNodes(
            ModelPart& this_model_part,
            const boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node,
            const boost::numeric::ublas::vector<int> &List_New_Nodes
    )
    {
	KRATOS_TRY;

        array_1d<double, 3 > Coord_Node_1;
        array_1d<double, 3 > Coord_Node_2;
        boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
        Coordinate_New_Node.resize(Position_Node.size(), false);
        unsigned int step_data_size = this_model_part.GetNodalSolutionStepDataSize();
        Node < 3 > ::DofsContainerType& reference_dofs = (this_model_part.NodesBegin())->GetDofs();

        //fill an auxiliary array with a pointer to the existing nodes
        std::unordered_map< unsigned int, Node<3>::Pointer > aux_node_list;
        for(auto it = this_model_part.NodesBegin(); it!=this_model_part.NodesEnd(); it++)
            aux_node_list[it->Id()] = ( *(it.base()) );


        for (unsigned int i = 0; i < Position_Node.size(); i++)
        {
            /* Calculating the coordinate of the news nodes */
            const int& node_i = Position_Node[i][0];
            const int& node_j = Position_Node[i][1];

            auto& it_node1 = aux_node_list[node_i];
            auto& it_node2 = aux_node_list[node_j];

// //             ModelPart::NodesContainerType::iterator it_node1 = this_model_part.Nodes().find(node_i);
//             std::size_t pos1 = it_node1 - this_model_part.NodesBegin();
            noalias(Coord_Node_1) = it_node1->Coordinates();

// //             ModelPart::NodesContainerType::iterator it_node2 = this_model_part.Nodes().find(node_j);
//             std::size_t pos2 = it_node2 - this_model_part.NodesBegin();
            noalias(Coord_Node_2) = it_node2->Coordinates();

            noalias(Coordinate_New_Node[i]) = 0.50 * (Coord_Node_1 + Coord_Node_2);

            /* Inserting the news node in the model part */
            Node < 3 >::Pointer pnode = Node < 3 >::Pointer(new Node < 3 >(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]));
            pnode->SetSolutionStepVariablesList( this_model_part.NodesBegin()->pGetVariablesList() );
            pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());

//             it_node1 = this_model_part.NodesBegin() + pos1;
//             it_node2 = this_model_part.NodesBegin() + pos2;

            pnode->GetValue(FATHER_NODES).resize(0);
            pnode->GetValue(FATHER_NODES).push_back(Node < 3 > ::WeakPointer( it_node1 ));
            pnode->GetValue(FATHER_NODES).push_back(Node < 3 > ::WeakPointer( it_node2 ));

            pnode->X0() = 0.5 * (it_node1->X0() + it_node2->X0());
            pnode->Y0() = 0.5 * (it_node1->Y0() + it_node2->Y0());
            pnode->Z0() = 0.5 * (it_node1->Z0() + it_node2->Z0());

//             KRATOS_WATCH(__LINE__)

            for (auto iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
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
//             KRATOS_WATCH(__LINE__)
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
//             KRATOS_WATCH(__LINE__)

            aux_node_list[pnode->Id()] = pnode;
            //this_model_part.Nodes().push_back(pnode);
        }

        this_model_part.Nodes().clear();
        this_model_part.Nodes().reserve(aux_node_list.size());
        for(auto & it : aux_node_list)
            this_model_part.Nodes().push_back(it.second);


        this_model_part.Nodes().Sort();

	KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void LocalRefineGeometryMesh::EraseOldElementAndCreateNewElement(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& New_Elements,
            bool interpolate_internal_variables
									     )
    {
      KRATOS_TRY;

      KRATOS_THROW_ERROR( std::logic_error, "Called the virtual function of LocalRefineGeometryMesh for EraseOldElementAndCreateNewElement", "" );

      KRATOS_CATCH( "" );
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void LocalRefineGeometryMesh::EraseOldConditionsAndCreateNew(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord
	 )
    {
      KRATOS_TRY;

      KRATOS_THROW_ERROR( std::logic_error, "Called the virtual function of LocalRefineGeometryMesh for EraseOldConditionsAndCreateNew", "" );

      KRATOS_CATCH( "" );
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void LocalRefineGeometryMesh::CalculateEdges(
            Element::GeometryType& geom,
            const compressed_matrix<int>& Coord,
            int* edge_ids,
            std::vector<int> & aux
            )
    {
      KRATOS_TRY;

      KRATOS_THROW_ERROR( std::logic_error, "Called the virtual function of LocalRefineGeometryMesh for CalculateEdges", "" );

      KRATOS_CATCH( "" );

    }

    /***********************************************************************************/
    /***********************************************************************************/

    void LocalRefineGeometryMesh::RenumeringElementsAndNodes(
            ModelPart& this_model_part,
            PointerVector< Element >& New_Elements
    )
    {
        KRATOS_TRY;

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

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    inline void LocalRefineGeometryMesh::CreatePartition(
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

    void LocalRefineGeometryMesh::InterpolateInteralVariables(
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


    void LocalRefineGeometryMesh::UpdateSubModelPartNodes(ModelPart &rModelPart)
    {
        for (ModelPart::SubModelPartIterator iSubModelPart = rModelPart.SubModelPartsBegin();
                iSubModelPart != rModelPart.SubModelPartsEnd(); iSubModelPart++)
        {
            for (auto iNode = rModelPart.Nodes().ptr_begin();
                    iNode != rModelPart.Nodes().ptr_end(); iNode++)
            {
                WeakPointerVector< Node<3> > &rFatherNodes = (*iNode)->GetValue(FATHER_NODES);
                unsigned int ParentCount = rFatherNodes.size();

                if (ParentCount > 0)
                {
                    unsigned int ParentsInSubModelPart = 0;

                    for ( WeakPointerVector< Node<3> >::iterator iParent = rFatherNodes.begin();
                            iParent != rFatherNodes.end(); iParent++)
                    {
                        unsigned int ParentId = iParent->Id();
                        ModelPart::NodeIterator iFound = iSubModelPart->Nodes().find( ParentId );
                        if ( iFound != iSubModelPart->NodesEnd() )
                            ParentsInSubModelPart++;
                    }

                    if ( ParentCount == ParentsInSubModelPart )
                        iSubModelPart->AddNode( *iNode );
                }
            }
        }
    }


    void LocalRefineGeometryMesh::ResetFatherNodes(ModelPart &rModelPart)
    {
        for (ModelPart::NodeIterator iNode = rModelPart.NodesBegin();
                iNode != rModelPart.NodesEnd(); iNode++)
        {
            ( iNode->GetValue(FATHER_NODES) ).clear();
        }
    }

} // Namespace Kratos.
