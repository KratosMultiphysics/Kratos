// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//

//Project includes
#include "custom_utilities/linear_to_quadratic_tetrahedra_mesh_converter_utility.h"

namespace Kratos {

    void LinearToQuadraticTetrahedraMeshConverter::LocalConvertLinearToQuadraticTetrahedraMesh(
        bool RefineOnReference,
        bool InterpolateInternalVariables)
    {
        //checking all elements to be refined are linear tets
        block_for_each(mModelPart.Elements(), [](const Element& r_element) {
            if(r_element.Has(SPLIT_ELEMENT) && r_element.GetValue(SPLIT_ELEMENT)){
                KRATOS_ERROR_IF_NOT(r_element.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) << "Element #" << r_element.Id() << " is not a linear tetrahedron" << std::endl;
            }
        });
        block_for_each(mModelPart.Conditions(), [](const Condition& r_condition) {
            if(r_condition.Has(SPLIT_ELEMENT) && r_condition.GetValue(SPLIT_ELEMENT)){
                KRATOS_ERROR_IF_NOT(r_condition.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) << "Condition #" << r_condition.Id() << " is not a linear tetrahedron" << std::endl;
            }
        });
        LocalRefineMesh(RefineOnReference, InterpolateInternalVariables);
    }

    Tetrahedra3D10<Node> LinearToQuadraticTetrahedraMeshConverter::GenerateTetrahedra(
        ModelPart& rThisModelPart,
        const std::vector<int>& rNodeIds)
    {
        unsigned int i0 = rNodeIds[0];
        unsigned int i1 = rNodeIds[1];
        unsigned int i2 = rNodeIds[2];
        unsigned int i3 = rNodeIds[3];
        unsigned int i4 = rNodeIds[4];
        unsigned int i6 = rNodeIds[5];
        unsigned int i7 = rNodeIds[6];
        unsigned int i5 = rNodeIds[7];
        unsigned int i8 = rNodeIds[8];
        unsigned int i9 = rNodeIds[9];

        Tetrahedra3D10<Node > geom(
            rThisModelPart.pGetNode(i0),
            rThisModelPart.pGetNode(i1),
            rThisModelPart.pGetNode(i2),
            rThisModelPart.pGetNode(i3),
            rThisModelPart.pGetNode(i4),
            rThisModelPart.pGetNode(i5),
            rThisModelPart.pGetNode(i6),
            rThisModelPart.pGetNode(i7),
            rThisModelPart.pGetNode(i8),
            rThisModelPart.pGetNode(i9)
        );
        return geom;
    }

    Triangle3D6<Node> LinearToQuadraticTetrahedraMeshConverter::GenerateTriangle3D6(
        ModelPart& rThisModelPart,
        const array_1d<int, 6>& rNodeIds)
    {
        unsigned int i0   = rNodeIds[0];
        unsigned int i1   = rNodeIds[1];
        unsigned int i2   = rNodeIds[2];
        unsigned int i3   = rNodeIds[3];
        unsigned int i4   = rNodeIds[4];
        unsigned int i5   = rNodeIds[5];

        Triangle3D6<Node > geom(
            rThisModelPart.pGetNode(i0),
            rThisModelPart.pGetNode(i1),
            rThisModelPart.pGetNode(i2),
            rThisModelPart.pGetNode(i3),
            rThisModelPart.pGetNode(i4),
            rThisModelPart.pGetNode(i5)
        );

        return geom;
    }

    void LinearToQuadraticTetrahedraMeshConverter::EraseOldElementAndCreateNewElement(
            ModelPart& rThisModelPart,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& NewElements,
            bool InterpolateInternalVariables
    )
    {
        auto& r_elements = rThisModelPart.Elements();
        ElementsArrayType::iterator it_begin = r_elements.ptr_begin();
        ElementsArrayType::iterator it_end = r_elements.ptr_end();

        const auto& r_current_process_info = rThisModelPart.GetProcessInfo();
        int edge_ids[6];
        std::vector<int> node_ids;

        const Element& r_elem = KratosComponents<Element>::Get("Element3D10N");
        for (ElementsArrayType::iterator& it = it_begin; it != it_end; ++it)
        {
            if(it->GetValue(SPLIT_ELEMENT)){
                // GlobalPointersVector< Element >& r_child_elements = it->GetValue(NEIGHBOUR_ELEMENTS);
                auto& r_child_elements = it->GetValue(NEIGHBOUR_ELEMENTS);
                r_child_elements.resize(0);

                CalculateEdges(it->GetGeometry(), Coord, edge_ids, node_ids);

                // Generate the new Tetrahedra3D10 element
                Tetrahedra3D10<Node> geom = GenerateTetrahedra(rThisModelPart, node_ids);
                Element::Pointer p_element;
                p_element = r_elem.Create(it->Id(), geom, it->pGetProperties());
                p_element->Initialize(r_current_process_info);
                p_element->InitializeSolutionStep(r_current_process_info);
                p_element->FinalizeSolutionStep(r_current_process_info);

                // Setting the internal variables in the "child" elem (the element replacing the old one)
                if (InterpolateInternalVariables == true)
                {
                    //This method only copies the current information to the new element
                    InterpolateInteralVariables(0, *it.base(), p_element, r_current_process_info);
                }

                // Transfer elemental variables to new element
                p_element->GetData() = it->GetData();
                p_element->GetValue(SPLIT_ELEMENT) = false;
                NewElements.push_back(p_element);

                r_child_elements.push_back( Element::WeakPointer(p_element) );
            } else { //element must not be replaced
                //checking that the element does not have split edges.
                //which would lead to quad tets in contact with linear tets / other elems
                if(it->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4){
                    //linear tet. Using member function to find edges
                    CalculateEdges(it->GetGeometry(), Coord, edge_ids, node_ids);
                    for(unsigned int j=4; j<10; j++){
                        KRATOS_ERROR_IF(node_ids[j]== -2 || node_ids[j]>0)  << "Element #" << it->Id() << " is in contact with a refined edge" << std::endl;
                    }
                } else {
                    //using generic(slow) implementation with edges:
                    const auto& edges = it->GetGeometry().GenerateEdges();
                    for (unsigned int j = 0; j < edges.size(); j++){
                        const GeometryType& rEdge = edges[j];
                        int index_0 = mMapNodeIdToPos[rEdge[0].Id()];
                        int index_1 = mMapNodeIdToPos[rEdge[1].Id()];
                        if(index_0>index_1){
                            const int new_node = Coord(index_1,index_0);
                            KRATOS_ERROR_IF(new_node== -2 || new_node>0)  << "Element #" << it->Id() << " is in contact with a refined edge" << std::endl;
                        } else {
                            const int new_node = Coord(index_0,index_1);
                            KRATOS_ERROR_IF(new_node== -2 || new_node>0)  << "Element #" << it->Id() << " is in contact with a refined edge" << std::endl;
                        }
                    }
                }
            }
        }

        // Now replace the elements in SubModelParts
        if ( NewElements.size() > 0 ) {
            ReplaceElementsInSubModelPart(rThisModelPart.GetRootModelPart());
        }
    }

    void LinearToQuadraticTetrahedraMeshConverter::ReplaceElementsInSubModelPart(ModelPart& rThisModelPart) {
        for(auto& p_element : rThisModelPart.ElementsArray()){
            if( p_element->GetValue(SPLIT_ELEMENT) )
            {
                // GlobalPointersVector< Element >& children = p_element->GetValue(NEIGHBOUR_ELEMENTS);
                auto& children = p_element->GetValue(NEIGHBOUR_ELEMENTS);
                p_element = children[0].shared_from_this();
            }
        }

        //Recursively for all subModelParts
        for (ModelPart::SubModelPartIterator i_submodelpart = rThisModelPart.SubModelPartsBegin();
                i_submodelpart != rThisModelPart.SubModelPartsEnd(); i_submodelpart++)
        {
            ReplaceElementsInSubModelPart(*i_submodelpart);
        }
    }

    void LinearToQuadraticTetrahedraMeshConverter::EraseOldConditionsAndCreateNew(
            ModelPart& rThisModelPart,
            const compressed_matrix<int>& Coord)
    {
        KRATOS_TRY;

        PointerVector< Condition > NewConditions;
        auto& r_conditions = rThisModelPart.Conditions();

        if(r_conditions.size() > 0)
        {
            ConditionsArrayType::iterator it_begin = r_conditions.ptr_begin();
            ConditionsArrayType::iterator it_end = r_conditions.ptr_end();
            int  edge_ids[3];
            array_1d<int, 6> node_ids;

            const auto& r_current_process_info = rThisModelPart.GetProcessInfo();
            const Condition& r_cond = KratosComponents<Condition>::Get("SurfaceCondition3D6N");

            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                if(it->GetValue(SPLIT_ELEMENT)){
                    CalculateEdgesFaces(it->GetGeometry(), Coord, edge_ids, node_ids);

                    // GlobalPointersVector< Condition >& r_child_conditions = it->GetValue(NEIGHBOUR_CONDITIONS);
                    auto& r_child_conditions = it->GetValue(NEIGHBOUR_CONDITIONS);
                    r_child_conditions.resize(0);

                    //Generate the new condition
                    Triangle3D6<Node > geom = GenerateTriangle3D6 (rThisModelPart, node_ids);
                    Condition::Pointer p_cond;
                    p_cond = r_cond.Create(it->Id(), geom, it->pGetProperties());
                    p_cond ->Initialize(r_current_process_info);
                    p_cond ->InitializeSolutionStep(r_current_process_info);
                    p_cond ->FinalizeSolutionStep(r_current_process_info);

                    // Transfer condition variables
                    p_cond->GetData() = it->GetData();
                    p_cond->GetValue(SPLIT_ELEMENT) = false;
                    NewConditions.push_back(p_cond);

                    r_child_conditions.push_back( Condition::WeakPointer( p_cond ) );
                }
            }

            // Replace the conditions in SubModelParts
            if (NewConditions.size() > 0) {
                ReplaceConditionsInSubModelPart(rThisModelPart.GetRootModelPart());
            }
        }
        KRATOS_CATCH("");
    }

    void LinearToQuadraticTetrahedraMeshConverter::ReplaceConditionsInSubModelPart(ModelPart& rThisModelPart) {
        for(auto& p_cond : rThisModelPart.ConditionsArray()){
            if( p_cond->GetValue(SPLIT_ELEMENT) )
            {
                // GlobalPointersVector< Condition >& children = p_cond->GetValue(NEIGHBOUR_CONDITIONS);
                auto& children = p_cond->GetValue(NEIGHBOUR_CONDITIONS);
                p_cond = children[0].shared_from_this();
            }
        }
        for (ModelPart::SubModelPartIterator i_submodelpart = rThisModelPart.SubModelPartsBegin();
                i_submodelpart != rThisModelPart.SubModelPartsEnd(); i_submodelpart++)
        {
            ReplaceConditionsInSubModelPart(*i_submodelpart);
        }
    }
}