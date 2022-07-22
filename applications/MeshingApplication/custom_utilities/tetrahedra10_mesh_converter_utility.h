// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Ariadna Cortés
//

#if !defined(KRATOS_TETRAHEDRA10_MESH_CONVERTER_UTILITY)
#define  KRATOS_TETRAHEDRA10_MESH_CONVERTER_UTILITY

// System includes

/* Project includes */
#include "includes/node.h"
#include "custom_utilities/local_refine_tetrahedra_mesh.hpp"
#include "containers/model.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/triangle_3d_6.h"

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
    Tetrahedra10MeshConverter(ModelPart& model_part) : LocalRefineTetrahedraMesh(model_part)
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
    * by adding intermediate nodes to each Tetrahedra3D4, with same ids as the originals. Does the same for conditions,
    * replacing the Triangle3D3 by Triangle3D6
    * @param refine_on_reference: Boolean that defines if refine or not the mesh according to the reference
    * @param interpolate_internal_variables: Boolean that defines if to interpolate or not the internal variables
    */
    void LocalConvertTetrahedra10Mesh(bool refine_on_reference, bool interpolate_internal_variables) {
            for (auto element : mModelPart.Elements()) element.SetValue(SPLIT_ELEMENT,true);
            for (auto condition : mModelPart.Conditions()) condition.SetValue(SPLIT_ELEMENT,true);
            LocalRefineMesh(refine_on_reference, interpolate_internal_variables);
        } 
  
    
protected:
    ///@name Protected static Member Variables
    ///@{
    //int mPreviousRefinementLevel;
    //int mCurrentRefinementLevel;
    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

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
    * @return the reference to the new tetrahedra
    */
    Tetrahedra3D10<Node<3>> GenerateTetrahedra(ModelPart& this_model_part, std::vector<int>& aux) {
        unsigned int i0 = aux[0];
        unsigned int i1 = aux[1];
        unsigned int i2 = aux[2];
        unsigned int i3 = aux[3];
        unsigned int i4 = aux[4];
        unsigned int i6 = aux[5];
        unsigned int i7 = aux[6];
        unsigned int i5 = aux[7];
        unsigned int i8 = aux[8];
        unsigned int i9 = aux[9];

        Tetrahedra3D10<Node < 3 > > geom(
            this_model_part.Nodes()(i0),
            this_model_part.Nodes()(i1),
            this_model_part.Nodes()(i2),
            this_model_part.Nodes()(i3),
            this_model_part.Nodes()(i4),
            this_model_part.Nodes()(i5),
            this_model_part.Nodes()(i6),
            this_model_part.Nodes()(i7),
            this_model_part.Nodes()(i8),
            this_model_part.Nodes()(i9)
        );
        return geom;
    }

    /**
    * Creates a new triangle3D6
    * @return the reference to the new triangle
    */
    Triangle3D6<Node<3>> GenerateTriangle3D6(ModelPart& this_model_part, array_1d<int, 6>& aux) {
        unsigned int i0   = aux[0];
        unsigned int i1   = aux[1];
        unsigned int i2   = aux[2];
        unsigned int i3   = aux[3];
        unsigned int i4   = aux[4];
        unsigned int i5   = aux[5];

        Triangle3D6<Node<3> > geom(
                this_model_part.Nodes()(i0),
                this_model_part.Nodes()(i1),
                this_model_part.Nodes()(i2),
                this_model_part.Nodes()(i3),
                this_model_part.Nodes()(i4),
                this_model_part.Nodes()(i5)
                );
        return geom;
    }

        /**
    * It erases the old Tetrahedra3D4 elements and it creates the new Tetrahedra3D10 ones
    * @param Coord: The compressed matrix containing at (i,j) the id of the node created between nodes i,j
    * @param New_Elements: The new elements created
    * @param interpolate_internal_variables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void EraseOldElementAndCreateNewElement(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& NewElements,
            bool interpolate_internal_variables
    ) override
    {
	ElementsArrayType& rElements = this_model_part.Elements();
	ElementsArrayType::iterator it_begin = rElements.ptr_begin();
	ElementsArrayType::iterator it_end = rElements.ptr_end();

	const ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();
	int edge_ids[6];
	std::vector<int> aux;

	for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
	{
        GlobalPointersVector< Element >& rChildElements = it->GetValue(NEIGHBOUR_ELEMENTS);
        rChildElements.resize(0);

        Element::GeometryType& geometry = it->GetGeometry();
        CalculateEdges(geometry, Coord, edge_ids, aux);

        // Generate the new Tetrahedra3D10 element
        Tetrahedra3D10<Node<3>> geom = GenerateTetrahedra(this_model_part, aux);
        Element::Pointer p_element;
        const Element& rElem = KratosComponents<Element>::Get("Element3D10N");
        p_element = rElem.Create(it->Id(), geom, it->pGetProperties());
        p_element->Initialize(rCurrentProcessInfo);
        p_element->InitializeSolutionStep(rCurrentProcessInfo);
        p_element->FinalizeSolutionStep(rCurrentProcessInfo);

        // Setting the internal variables in the child elem
        if (interpolate_internal_variables == true)
        {
            //This method only copies the current information to the new element
            InterpolateInteralVariables(0, *it.base(), p_element, rCurrentProcessInfo);
        }

        // Transfer elemental variables
        p_element->Data() = it->Data();
        p_element->GetValue(SPLIT_ELEMENT) = false;
        NewElements.push_back(p_element);

        rChildElements.push_back( Element::WeakPointer(p_element) );
    }

    // Now update the elements in SubModelParts
    if ( NewElements.size() > 0 ) ReplaceElementsInSubModelPart(this_model_part);
}

    /**
    * It replaces the old Tetrahedra3D4 elements in its corresponding submodelpart by its new Tetrahedra3D10
    * @param this_model_part: The modelpart or submodelpart to replace the elements
    */
    void ReplaceElementsInSubModelPart(ModelPart& this_model_part) {
        for(Element::Pointer& p_element : this_model_part.ElementsArray()){
                if( p_element->GetValue(SPLIT_ELEMENT) )
                {
                    GlobalPointersVector< Element >& children = p_element->GetValue(NEIGHBOUR_ELEMENTS);
                    p_element = children[0].shared_from_this();
                } 
        }

        //Recursively for all subModelParts
        for (ModelPart::SubModelPartIterator iSubModelPart = this_model_part.SubModelPartsBegin();
                iSubModelPart != this_model_part.SubModelPartsEnd(); iSubModelPart++)
        {
            ReplaceElementsInSubModelPart(*iSubModelPart);
        }
    }

    /**
    * Remove the old Trangle3D3 conditions and replace them by the new Trangle3D6 ones
    * @param Coord: The coordinates of the nodes of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void EraseOldConditionsAndCreateNew(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord
            ) override
    {
        KRATOS_TRY;

        PointerVector< Condition > NewConditions;
        ConditionsArrayType& rConditions = this_model_part.Conditions();

        if(rConditions.size() > 0)
        {
            ConditionsArrayType::iterator it_begin = rConditions.ptr_begin();
            ConditionsArrayType::iterator it_end = rConditions.ptr_end();
            int  edge_ids[3];
            array_1d<int, 6> aux;

            const ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();

            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                Condition::GeometryType& geom = it->GetGeometry();
                CalculateEdgesFaces(geom, Coord, edge_ids, aux);

                GlobalPointersVector< Condition >& rChildConditions = it->GetValue(NEIGHBOUR_CONDITIONS);
                rChildConditions.resize(0);

                it->Set(TO_ERASE,true); //Mark them as the "old" conditions for later remove

                Triangle3D6<Node<3> > newgeom = GenerateTriangle3D6 (this_model_part, aux);
            
                // Generate new condition by cloning the base one
                Condition::Pointer pcond;
                const Condition& rCond = KratosComponents<Condition>::Get("SurfaceCondition3D6N");
                pcond = rCond.Create(it->Id(), newgeom, it->pGetProperties());
                pcond ->Initialize(rCurrentProcessInfo);
                pcond ->InitializeSolutionStep(rCurrentProcessInfo);
                pcond ->FinalizeSolutionStep(rCurrentProcessInfo);

                // Transfer condition variables
                pcond->Data() = it->Data();
                pcond->GetValue(SPLIT_ELEMENT) = false;
                NewConditions.push_back(pcond);

                rChildConditions.push_back( Condition::WeakPointer( pcond ) );
            }

            this_model_part.Conditions().reserve(this_model_part.Conditions().size()+ NewConditions.size());

            /// Add the new Conditions to the ModelPart
            for (auto iCond = NewConditions.ptr_begin();
                    iCond != NewConditions.ptr_end(); iCond++)
            {
                this_model_part.Conditions().push_back( *iCond );
            }

            // Now update the conditions in SubModelParts
            if (NewConditions.size() > 0) ReplaceConditionsInSubModelPart(this_model_part);

            //Finally remove the old conditions
            this_model_part.RemoveConditions(TO_ERASE);
        }
        KRATOS_CATCH("");
    }


    /**
    * It replaces the old Tetrahedra3D4 elements in its corresponding submodelpart by its new Tetrahedra3D10
    * @param this_model_part: The modelpart or submodelpart to replace the elements
    */
    void ReplaceConditionsInSubModelPart(ModelPart& this_model_part) {
        for(Condition::Pointer& p_cond : this_model_part.ConditionsArray()){
                if( p_cond->GetValue(SPLIT_ELEMENT) )
                {
                    GlobalPointersVector< Condition >& children = p_cond->GetValue(NEIGHBOUR_CONDITIONS);
                    p_cond = children[0].shared_from_this();
                } 
        }
        for (ModelPart::SubModelPartIterator iSubModelPart = this_model_part.SubModelPartsBegin();
                iSubModelPart != this_model_part.SubModelPartsEnd(); iSubModelPart++)
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

};

} // namespace Kratos.

#endif // KRATOS_TET10_REFINEMENT_UTILITY  defined
