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

#if !defined(KRATOS_TETRAHEDRA10_MESH_CONVERTER_UTILITY)
#define  KRATOS_TETRAHEDRA10_MESH_CONVERTER_UTILITY

// System includes

/* Project includes */
#include "custom_utilities/local_refine_tetrahedra_mesh.hpp"
#include "geometries/tetrahedra_3d_10.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{
class Tetrahedra10MeshConverter : public LocalRefineTetrahedraMesh
{
public:

    ///@name Type Definitions
    ///@{
    ///@}

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;

    ///@name Life Cycle
    ///@{

    /// Default constructors
    Tetrahedra10MeshConverter(ModelPart& ModelPart) : LocalRefineTetrahedraMesh(ModelPart)
    {

    }

    /// Destructor
    ~Tetrahedra10MeshConverter() 
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
        
    /**
    * Replaces the geometry of the mesh of Tetrahedra3D4 elements. Resulting is a mesh of Tetrahedra3D10 created
    * by adding intermediate nodes to each Tetrahedra3D4. Does the same for conditions, replacing the Triangle3D3 
    * by Triangle3D6
    * @param refine_on_reference: Boolean that defines if refine or not the mesh according to the reference
    * @param interpolate_internal_variables: Boolean that defines if to interpolate or not the internal variables
    */
    void LocalConvertTetrahedra10Mesh(bool refine_on_reference, bool interpolate_internal_variables) {
        block_for_each(mModelPart.Elements(), [&](Element element) {
            element.SetValue(SPLIT_ELEMENT,true);
        });
        block_for_each(mModelPart.Conditions(), [&](Condition condition) {
            condition.SetValue(SPLIT_ELEMENT,true);
        });
        LocalRefineMesh(refine_on_reference, interpolate_internal_variables);
    } 
  
private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
    * Creates a new tetrahedra3D10
    * @return The new tetrahedra
    */
    Tetrahedra3D10<Node<3>> GenerateTetrahedra(
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

        Tetrahedra3D10<Node < 3 > > Geom(
            rThisModelPart.Nodes()(i0),
            rThisModelPart.Nodes()(i1),
            rThisModelPart.Nodes()(i2),
            rThisModelPart.Nodes()(i3),
            rThisModelPart.Nodes()(i4),
            rThisModelPart.Nodes()(i5),
            rThisModelPart.Nodes()(i6),
            rThisModelPart.Nodes()(i7),
            rThisModelPart.Nodes()(i8),
            rThisModelPart.Nodes()(i9)
        );
        return Geom;
    }

    /**
    * Creates a new triangle3D6
    * @return The new triangle
    */
    Triangle3D6<Node<3>> GenerateTriangle3D6(
        ModelPart& rThisModelPart, 
        const array_1d<int, 6>& rNodeIds) 
    {
        unsigned int i0   = rNodeIds[0];
        unsigned int i1   = rNodeIds[1];
        unsigned int i2   = rNodeIds[2];
        unsigned int i3   = rNodeIds[3];
        unsigned int i4   = rNodeIds[4];
        unsigned int i5   = rNodeIds[5];

        Triangle3D6<Node<3> > Geom(
                rThisModelPart.Nodes()(i0),
                rThisModelPart.Nodes()(i1),
                rThisModelPart.Nodes()(i2),
                rThisModelPart.Nodes()(i3),
                rThisModelPart.Nodes()(i4),
                rThisModelPart.Nodes()(i5)
                );
        return Geom;
    }

        /**
    * It erases the old Tetrahedra3D4 elements and it creates the new Tetrahedra3D10 ones
    * @param Coord: The compressed matrix containing at (i,j) the id of the node created between nodes i,j
    * @param New_Elements: The new elements created
    * @param interpolate_internal_variables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return rThisModelPart: The model part of the model (it is the input too)
    */

    void EraseOldElementAndCreateNewElement(
            ModelPart& rThisModelPart,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& NewElements,
            bool interpolate_internal_variables
    ) override
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
            // GlobalPointersVector< Element >& r_child_elements = it->GetValue(NEIGHBOUR_ELEMENTS);
            auto& r_child_elements = it->GetValue(NEIGHBOUR_ELEMENTS);
            r_child_elements.resize(0);

            CalculateEdges(it->GetGeometry(), Coord, edge_ids, node_ids);

            // Generate the new Tetrahedra3D10 element
            Tetrahedra3D10<Node<3>> geom = GenerateTetrahedra(rThisModelPart, node_ids);
            Element::Pointer p_element;
            p_element = r_elem.Create(it->Id(), geom, it->pGetProperties());
            p_element->Initialize(r_current_process_info);
            p_element->InitializeSolutionStep(r_current_process_info);
            p_element->FinalizeSolutionStep(r_current_process_info);

            // Setting the internal variables in the "child" elem (the element replacing the old one)
            if (interpolate_internal_variables == true)
            {
                //This method only copies the current information to the new element
                InterpolateInteralVariables(0, *it.base(), p_element, r_current_process_info);
            }
            
            // Transfer elemental variables to new element
            p_element->Data() = it->Data();
            p_element->GetValue(SPLIT_ELEMENT) = false;
            NewElements.push_back(p_element);

            r_child_elements.push_back( Element::WeakPointer(p_element) );
        }

        // Now replace the elements in SubModelParts
        if ( NewElements.size() > 0 ) {
            ReplaceElementsInSubModelPart(rThisModelPart);
        }
    }

    /**
    * It replaces the old Tetrahedra3D4 elements in its corresponding submodelpart by its new Tetrahedra3D10
    * @param rThisModelPart: The modelpart or submodelpart to replace the elements
    */
    void ReplaceElementsInSubModelPart(ModelPart& rThisModelPart) {
        for(auto& p_element : rThisModelPart.ElementsArray()){
            if( p_element->GetValue(SPLIT_ELEMENT) )
            {
                // GlobalPointersVector< Element >& children = p_element->GetValue(NEIGHBOUR_ELEMENTS);
                auto& children = p_element->GetValue(NEIGHBOUR_ELEMENTS);
                p_element = children[0].shared_from_this();
            } 
        }

        //Recursively for all subModelParts
        for (ModelPart::SubModelPartIterator iSubModelPart = rThisModelPart.SubModelPartsBegin();
                iSubModelPart != rThisModelPart.SubModelPartsEnd(); iSubModelPart++)
        {
            ReplaceElementsInSubModelPart(*iSubModelPart);
        }
    }

    /**
    * Remove the old Trangle3D3 conditions and replace them by the new Trangle3D6 ones
    * @param Coord: The coordinates of the nodes of the geometry
    * @return rThisModelPart: The model part of the model (it is the input too)
    */

    void EraseOldConditionsAndCreateNew(
            ModelPart& rThisModelPart,
            const compressed_matrix<int>& Coord
            ) override
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
                CalculateEdgesFaces(it->GetGeometry(), Coord, edge_ids, node_ids);

                // GlobalPointersVector< Condition >& r_child_conditions = it->GetValue(NEIGHBOUR_CONDITIONS);
                auto& r_child_conditions = it->GetValue(NEIGHBOUR_CONDITIONS);
                r_child_conditions.resize(0);

                //Generate the new condition
                Triangle3D6<Node<3> > geom = GenerateTriangle3D6 (rThisModelPart, node_ids);
                Condition::Pointer pcond;
                pcond = r_cond.Create(it->Id(), geom, it->pGetProperties());
                pcond ->Initialize(r_current_process_info);
                pcond ->InitializeSolutionStep(r_current_process_info);
                pcond ->FinalizeSolutionStep(r_current_process_info);

                // Transfer condition variables
                pcond->Data() = it->Data();
                pcond->GetValue(SPLIT_ELEMENT) = false;
                NewConditions.push_back(pcond);

                r_child_conditions.push_back( Condition::WeakPointer( pcond ) );
            }

            // Replace the conditions in SubModelParts
            if (NewConditions.size() > 0) {
                ReplaceConditionsInSubModelPart(rThisModelPart);
            }
        }
        KRATOS_CATCH("");
    }


    /**
    * It replaces the old Triangle3D3 elements in its corresponding submodelpart by its new Triangle3D6
    * @param rThisModelPart: The modelpart or submodelpart to replace the elements
    */
    void ReplaceConditionsInSubModelPart(ModelPart& rThisModelPart) {
        for(auto& pCond : rThisModelPart.ConditionsArray()){
            if( pCond->GetValue(SPLIT_ELEMENT) )
            {
                // GlobalPointersVector< Condition >& children = pCond->GetValue(NEIGHBOUR_CONDITIONS);
                auto& children = pCond->GetValue(NEIGHBOUR_CONDITIONS);
                pCond = children[0].shared_from_this();
            } 
        }
        for (ModelPart::SubModelPartIterator iSubModelPart = rThisModelPart.SubModelPartsBegin();
                iSubModelPart != rThisModelPart.SubModelPartsEnd(); iSubModelPart++)
        {
            ReplaceConditionsInSubModelPart(*iSubModelPart);
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{
    ///@}

}; // Class Tetrahedra10MeshConverter

} // namespace Kratos.

#endif // KRATOS_TET10_REFINEMENT_UTILITY  defined
