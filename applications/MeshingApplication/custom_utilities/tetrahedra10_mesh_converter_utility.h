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

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

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
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;

    ///@name Life Cycle
    ///@{

    /// Default constructors
    Tetrahedra10MeshConverter(ModelPart& model_part) : LocalRefineTetrahedraMesh(model_part)
    {

    }

    /// Destructor
    ~Tetrahedra10MeshConverter() //TODO maybe {}
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
        
    /**
    * Changes the geometry of the mesh of Tetrahedra3D4 elements. Resulting is a mesh of Tetrahedra3D10 created
    * by adding intermediate nodes to each Tetrahedra3D4
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
	unsigned int to_be_deleted = 0;
	unsigned int current_id = (rElements.end() - 1)->Id() + 1;
	int internal_node = 0;

	const ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();
	int edge_ids[6];
	std::vector<int> aux;

    KRATOS_INFO("") << "************* CONVERTING TO TETRAHEDRA3D8 MESH **************\n" 
                    << "OLD NUMBER ELEMENTS: " << rElements.size() << std::endl;


	for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
	{
        Element::GeometryType& geometry = it->GetGeometry();
        CalculateEdges(geometry, Coord, edge_ids, aux);

        it->Set(TO_ERASE,true);

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

        // Generate the new Tetrahedra3D10 element
        Element::Pointer p_element;
        const Element& rElem = KratosComponents<Element>::Get("Element3D10N");
        p_element = rElem.Create(current_id, geom, it->pGetProperties());
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

        current_id++;
    }

    /* Adding news elements to the model part */
    for (PointerVector< Element >::iterator it_new = NewElements.begin(); it_new != NewElements.end(); it_new++)
    {
    rElements.push_back(*(it_new.base()));
    }

    /* Now remove all of the "old" elements */
    this_model_part.RemoveElements(TO_ERASE);

    KRATOS_INFO("") <<  "NEW NUMBER ELEMENTS: " << rElements.size() << std::endl;

    // Now update the elements in SubModelParts
    if (NewElements.size() > 0)
    {
        for (ModelPart::SubModelPartIterator iSubModelPart = this_model_part.SubModelPartsBegin();
                iSubModelPart != this_model_part.SubModelPartsEnd(); iSubModelPart++)
        {
            to_be_deleted = 0;
            NewElements.clear();

            // Create list of new elements in SubModelPart
            // Count how many elements will be removed
            for (ModelPart::ElementIterator iElem = iSubModelPart->ElementsBegin();
                    iElem != iSubModelPart->ElementsEnd(); iElem++)
            {
                if( iElem->GetValue(SPLIT_ELEMENT) )
                {
                    to_be_deleted++;
                    GlobalPointersVector< Element >& rChildElements = iElem->GetValue(NEIGHBOUR_ELEMENTS);

                    for ( auto iChild = rChildElements.ptr_begin();
                            iChild != rChildElements.ptr_end(); iChild++ )
                    {
                        NewElements.push_back((*iChild)->shared_from_this());
                    }
                }
            }

            // Add new elements to SubModelPart
            iSubModelPart->Elements().reserve( iSubModelPart->Elements().size() + NewElements.size() );
            for (PointerVector< Element >::iterator it_new = NewElements.begin();
                    it_new != NewElements.end(); it_new++)
            {
                iSubModelPart->Elements().push_back(*(it_new.base()));
            }

            // Delete old elements
            iSubModelPart->Elements().Sort();
            iSubModelPart->Elements().erase(iSubModelPart->Elements().end() - to_be_deleted, iSubModelPart->Elements().end());
        }
    }
}

    /**
    * Remove the old Trangle3D3 conditions and creates new Trangle3D6 ones
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
            unsigned int to_be_deleted = 0;
            int  edge_ids[3];
            array_1d<int, 6> aux;

            const ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();

            unsigned int current_id = (rConditions.end() - 1)->Id() + 1;
            for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it)
            {
                Condition::GeometryType& geom = it->GetGeometry();

                CalculateEdgesFaces(geom, Coord, edge_ids, aux);

                // Create the new conditions
                it->Set(TO_ERASE,true);

                unsigned int i0   = aux[0];
                unsigned int i1   = aux[1];
                unsigned int i2   = aux[2];
                unsigned int i3   = aux[3];
                unsigned int i4   = aux[4];
                unsigned int i5   = aux[5];

                Triangle3D6<Node<3> > newgeom(
                        this_model_part.Nodes()(i0),
                        this_model_part.Nodes()(i1),
                        this_model_part.Nodes()(i2),
                        this_model_part.Nodes()(i3),
                        this_model_part.Nodes()(i4),
                        this_model_part.Nodes()(i5)
                        );

                // Generate new condition by cloning the base one
                Condition::Pointer pcond;
                const Condition& rCond = KratosComponents<Condition>::Get("SurfaceCondition3D6N");
                pcond = rCond.Create(current_id, newgeom, it->pGetProperties());
                pcond ->Initialize(rCurrentProcessInfo);
                pcond ->InitializeSolutionStep(rCurrentProcessInfo);
                pcond ->FinalizeSolutionStep(rCurrentProcessInfo);

                // Transfer condition variables
                pcond->Data() = it->Data();
                pcond->GetValue(SPLIT_ELEMENT) = false;
                NewConditions.push_back(pcond);

                current_id++;
            }


            /* Now remove all of the "old" conditions*/
            this_model_part.RemoveConditions(TO_ERASE);

            unsigned int total_size = this_model_part.Conditions().size()+ NewConditions.size();
            this_model_part.Conditions().reserve(total_size);

            /// Add the new Conditions to the ModelPart
            for (auto iCond = NewConditions.ptr_begin();
                    iCond != NewConditions.ptr_end(); iCond++)
            {
                this_model_part.Conditions().push_back( *iCond );
            }

            /* Renumber id */
            unsigned int my_index = 1;
            for(ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); it++)
            {
                it->SetId(my_index++);
            }


            // Now update the conditions in SubModelParts
            if (NewConditions.size() > 0)
            {
                for (ModelPart::SubModelPartIterator iSubModelPart = this_model_part.SubModelPartsBegin();
                        iSubModelPart != this_model_part.SubModelPartsEnd(); iSubModelPart++)
                {
                    to_be_deleted = 0;
                    NewConditions.clear();

                    // Create list of new conditions in SubModelPart
                    // Count how many conditions will be removed
                    for (ModelPart::ConditionIterator iCond = iSubModelPart->ConditionsBegin();
                            iCond != iSubModelPart->ConditionsEnd(); iCond++)
                    {
                        if( iCond->GetValue(SPLIT_ELEMENT) )
                        {
                            to_be_deleted++;
                            GlobalPointersVector< Condition >& rChildConditions = iCond->GetValue(NEIGHBOUR_CONDITIONS);

                            for ( auto iChild = rChildConditions.ptr_begin();
                                    iChild != rChildConditions.ptr_end(); iChild++ )
                            {
                                NewConditions.push_back((*iChild)->shared_from_this());
                            }
                        }
                    }

                    // Add new conditions to SubModelPart
                    iSubModelPart->Conditions().reserve( iSubModelPart->Conditions().size() + NewConditions.size() );
                    for (PointerVector< Condition >::iterator it_new = NewConditions.begin();
                            it_new != NewConditions.end(); it_new++)
                    {
                        iSubModelPart->Conditions().push_back(*(it_new.base()));
                    }

                    // Delete old conditions
                    iSubModelPart->Conditions().Sort();
                    iSubModelPart->Conditions().erase(iSubModelPart->Conditions().end() - to_be_deleted, iSubModelPart->Conditions().end());
                }
            }
        }
        KRATOS_CATCH("");
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
