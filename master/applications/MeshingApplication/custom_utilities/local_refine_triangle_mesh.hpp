// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//    Co-authors:    Vicente Mataix Ferrandiz
//

#pragma once

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

// External includes

// Project includes
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "custom_utilities/local_refine_triangle_mesh_generic.hpp"

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

class LocalRefineTriangleMesh : public LocalRefineTriangleMeshGeneric<Element, ModelPart::ElementsContainerType>
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LocalRefineTriangleMesh(ModelPart& model_part) : LocalRefineTriangleMeshGeneric(model_part)
    {

    }

    /// Destructor
    ~LocalRefineTriangleMesh() override
    = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * It erases the old elements and it creates the new ones
    * @param Coord: The coordinates of the element
    * @param New_Elements: The new elements created
    * @param interpolate_internal_variables: A boolean that defines if it is necessary to interpolate the internal variables
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void EraseOldElementAndCreateNewElement(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& New_Elements,
            bool interpolate_internal_variables
    ) override
    {
        EraseOldObjetcsAndCreateNewObjects(
            this_model_part,
            this_model_part.Elements(),            
            Coord,
            New_Elements,
            interpolate_internal_variables
        );

        // Now update the elements in SubModelParts
        if (New_Elements.size() > 0)
        {
            UpdateSubModelPartElements(this_model_part.GetRootModelPart(), New_Elements);
        }

    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Remove the old conditions and creates new ones
    * @param Coord: The coordinates of the nodes of the geometry
    * @return this_model_part: The model part of the model (it is the input too)
    */

    void EraseOldConditionsAndCreateNew(
	ModelPart& this_model_part,
	const compressed_matrix<int>& Coord
	 ) override
    {
        KRATOS_TRY;

        PointerVector< Condition > New_Conditions;

        ConditionsArrayType& rConditions = this_model_part.Conditions();

        if (rConditions.size() > 0)
        {
            ConditionsArrayType::iterator it_begin = rConditions.ptr_begin();
            ConditionsArrayType::iterator it_end   = rConditions.ptr_end();
            unsigned int to_be_deleted = 0;
            unsigned int large_id = (rConditions.end() - 1)->Id() * 7;

            const ProcessInfo& rCurrentProcessInfo = this_model_part.GetProcessInfo();

	    const unsigned int dimension = it_begin->GetGeometry().WorkingSpaceDimension();

            unsigned int current_id = (rConditions.end() - 1)->Id() + 1;

	    for (ConditionsArrayType::iterator it = it_begin; it != it_end; ++it)
	    {
                Condition::GeometryType& geom = it->GetGeometry();

                if (geom.size() == 2)
                {
                    int index_0 = mMapNodeIdToPos[geom[0].Id()];
                    int index_1 = mMapNodeIdToPos[geom[1].Id()];
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
                        if (dimension == 2)
                        {
                            Line2D2<Node > newgeom1(
                                this_model_part.Nodes()(geom[0].Id()),
                                this_model_part.Nodes()(new_id)
                            );

                            Line2D2<Node > newgeom2(
                                this_model_part.Nodes()(new_id),
                                this_model_part.Nodes()(geom[1].Id())
                            );

                            Condition::Pointer pcond1 = it->Create(current_id++, newgeom1, it->pGetProperties());
                            Condition::Pointer pcond2 = it->Create(current_id++, newgeom2, it->pGetProperties());

                            pcond1->GetData() = it->GetData();
                            pcond2->GetData() = it->GetData();

                            New_Conditions.push_back(pcond1);
                            New_Conditions.push_back(pcond2);
                        }
                        else
                        {
                            Line3D2<Node > newgeom1(
                                this_model_part.Nodes()(geom[0].Id()),
                                this_model_part.Nodes()(new_id)
                            );

                            Line3D2<Node > newgeom2(
                                this_model_part.Nodes()(new_id),
                                this_model_part.Nodes()(geom[1].Id())
                            );

                            Condition::Pointer pcond1 = it->Create(current_id++, newgeom1, it->pGetProperties());
                            Condition::Pointer pcond2 = it->Create(current_id++, newgeom2, it->pGetProperties());

                            pcond1->GetData() = it->GetData();
                            pcond2->GetData() = it->GetData();

                            New_Conditions.push_back(pcond1);
                            New_Conditions.push_back(pcond2);
                        }

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
                it_new->Initialize(rCurrentProcessInfo);
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

        KRATOS_CATCH("");
    }

    /**
    * Updates recursively the elements in the submodelpars
    * @param NewElements: list of elems
    * @return this_model_part: The model part of the model (it is the input too)
    */
    void UpdateSubModelPartElements(ModelPart& this_model_part, PointerVector< Element >& NewElements)
      {
          for (ModelPart::SubModelPartIterator iSubModelPart = this_model_part.SubModelPartsBegin();
                  iSubModelPart != this_model_part.SubModelPartsEnd(); iSubModelPart++)
          {
              unsigned int to_be_deleted = 0;
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
              /*
              KRATOS_WATCH(iSubModelPart->Info());
              KRATOS_WATCH(to_be_deleted);
              KRATOS_WATCH(iSubModelPart->Elements().size());
              KRATOS_WATCH(this_model_part.Elements().size());
              */

            //NEXT LEVEL
            if (NewElements.size() > 0)
            {
               ModelPart &rSubModelPart = *iSubModelPart;
               UpdateSubModelPartElements(rSubModelPart,NewElements);
            }


          }
      }
    /***********************************************************************************/
    /***********************************************************************************/


};

} // namespace Kratos.
