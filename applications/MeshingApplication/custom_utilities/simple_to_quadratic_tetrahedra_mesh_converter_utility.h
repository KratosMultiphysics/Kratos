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

#if !defined(KRATOS_SIMPLE_TO_QUADRATIC_TETRAHEDRA_MESH_CONVERTER_UTILITY)
#define  KRATOS_SIMPLE_TO_QUADRATIC_TETRAHEDRA_MESH_CONVERTER_UTILITY

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
class SimpleToQuadraticTetrahedraMeshConverter : public LocalRefineTetrahedraMesh
{
public:

    ///@name Type Definitions
    ///@{
    ///@}

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;

     /// Pointer definition of VoxelInsideVolume
    KRATOS_CLASS_POINTER_DEFINITION( SimpleToQuadraticTetrahedraMeshConverter );

    ///@name Life Cycle
    ///@{

    /// Default constructors
    SimpleToQuadraticTetrahedraMeshConverter(ModelPart& ModelPart) : LocalRefineTetrahedraMesh(ModelPart)
    {

    }

    /// Destructor
    ~SimpleToQuadraticTetrahedraMeshConverter() 
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
    * @param RefineOnReference: Boolean that defines if refine or not the mesh according to the reference
    * @param InterpolateInternalVariables: Boolean that defines if to interpolate or not the internal variables
    */
    void LocalConvertSimpleToQuadraticTetrahedraMesh(
        bool RefineOnReference, 
        bool InterpolateInternalVariables) 
    {
        block_for_each(mModelPart.Elements(), [&](Element element) {
            element.SetValue(SPLIT_ELEMENT,true);
        });
        block_for_each(mModelPart.Conditions(), [&](Condition condition) {
            condition.SetValue(SPLIT_ELEMENT,true);
        });
        LocalRefineMesh(RefineOnReference, InterpolateInternalVariables);
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

        Tetrahedra3D10<Node < 3 > > geom(
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

        Triangle3D6<Node<3> > geom(
            rThisModelPart.pGetNode(i0),
            rThisModelPart.pGetNode(i1),
            rThisModelPart.pGetNode(i2),
            rThisModelPart.pGetNode(i3),
            rThisModelPart.pGetNode(i4),
            rThisModelPart.pGetNode(i5)
        );

        return geom;
    }

        /**
    * It erases the old Tetrahedra3D4 elements and it creates the new Tetrahedra3D10 ones
    * @param Coord: The compressed matrix containing at (i,j) the id of the node created between nodes i,j
    * @param New_Elements: The new elements created
    * @param InterpolateInternalVariables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return rThisModelPart: The model part of the model (it is the input too)
    */

    void EraseOldElementAndCreateNewElement(
            ModelPart& rThisModelPart,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& NewElements,
            bool InterpolateInternalVariables
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
        for (ModelPart::SubModelPartIterator i_submodelpart = rThisModelPart.SubModelPartsBegin();
                i_submodelpart != rThisModelPart.SubModelPartsEnd(); i_submodelpart++)
        {
            ReplaceElementsInSubModelPart(*i_submodelpart);
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
                Condition::Pointer p_cond;
                p_cond = r_cond.Create(it->Id(), geom, it->pGetProperties());
                p_cond ->Initialize(r_current_process_info);
                p_cond ->InitializeSolutionStep(r_current_process_info);
                p_cond ->FinalizeSolutionStep(r_current_process_info);

                // Transfer condition variables
                p_cond->Data() = it->Data();
                p_cond->GetValue(SPLIT_ELEMENT) = false;
                NewConditions.push_back(p_cond);

                r_child_conditions.push_back( Condition::WeakPointer( p_cond ) );
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

}; // Class SimpleToQuadraticTetrahedraMeshConverter

} // namespace Kratos.

#endif // KRATOS_SIMPLE_TO_QUADRATIC_TETRAHEDRA_MESH_CONVERTER_UTILITY  defined
