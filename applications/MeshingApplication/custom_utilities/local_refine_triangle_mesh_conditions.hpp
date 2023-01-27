#if !defined(KRATOS_LOCAL_REFINE_TRIANGLE_MESH_CONDITIONS)
#define  KRATOS_LOCAL_REFINE_TRIANGLE_MESH_CONDITIONS

#ifdef _OPENMP
#include <omp.h>
#endif

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

/* Project includes */
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

class LocalRefineTriangleMeshConditions : public LocalRefineTriangleMeshGeneric<Condition, ModelPart::ConditionsContainerType>
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LocalRefineTriangleMeshConditions(ModelPart& model_part) : LocalRefineTriangleMeshGeneric(model_part)
    {

    }

    /// Destructor
    ~LocalRefineTriangleMeshConditions() override
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
    * Nothing done here since only conditions are used 
    */

    void EraseOldElementAndCreateNewElement(
            ModelPart& this_model_part,
            const compressed_matrix<int>& Coord,
            PointerVector< Element >& New_Elements,
            bool interpolate_internal_variables
    ) override
    {
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
        PointerVector< Condition > New_Conditions;

        EraseOldObjetcsAndCreateNewObjects(
            this_model_part,
            this_model_part.Conditions(),            
            Coord,
            New_Conditions,
            false
        );

        // Now update the elements in SubModelParts
        if (New_Conditions.size() > 0)
        {
            KRATOS_WATCH(New_Conditions.size())
            UpdateSubModelPartConditions(this_model_part.GetRootModelPart(), New_Conditions);
        }
    }

    
    void SearchEdgeToBeRefined(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord
    ) override
    {
        KRATOS_TRY;

        ModelPart::ConditionsContainerType& rConditions = this_model_part.Conditions();
        ModelPart::ConditionsContainerType::iterator it_begin = rConditions.ptr_begin();
        ModelPart::ConditionsContainerType::iterator it_end   = rConditions.ptr_end();

        this->SearchEdgeToBeRefinedGeneric<ModelPart::ConditionsContainerType::iterator>(it_begin,it_end,Coord);

        KRATOS_CATCH("");
    }

    void ResetFatherNodes(ModelPart &rModelPart) override
    {
        block_for_each(mModelPart.Nodes(), [&](Node<3>& rNode)
        {
            if(!rNode.Has(FATHER_NODES)){
                GlobalPointersVector<Node<3>> empty_father_vector;
                rNode.SetValue(FATHER_NODES, empty_father_vector);
            }
        });
    }

    /**
    * Updates recursively the elements in the submodelpars
    * @param NewElements: list of elems
    * @return this_model_part: The model part of the model (it is the input too)
    */
    void UpdateSubModelPartConditions(ModelPart& this_model_part, PointerVector< Condition >& NewConditions)
      {
          for (ModelPart::SubModelPartIterator iSubModelPart = this_model_part.SubModelPartsBegin();
                  iSubModelPart != this_model_part.SubModelPartsEnd(); iSubModelPart++)
          {
              unsigned int to_be_deleted = 0;
              NewConditions.clear();

              // Create list of new elements in SubModelPart
              // Count how many elements will be removed
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

              // Add new elements to SubModelPart
              iSubModelPart->Conditions().reserve( iSubModelPart->Conditions().size() + NewConditions.size() );
              for (PointerVector< Condition >::iterator it_new = NewConditions.begin();
                      it_new != NewConditions.end(); it_new++)
              {
                  iSubModelPart->Conditions().push_back(*(it_new.base()));
              }

              // Delete old elements
              iSubModelPart->Conditions().Sort();
              iSubModelPart->Conditions().erase(iSubModelPart->Conditions().end() - to_be_deleted, iSubModelPart->Conditions().end());

            //NEXT LEVEL
            if (NewConditions.size() > 0)
            {
               ModelPart &rSubModelPart = *iSubModelPart;
               UpdateSubModelPartConditions(rSubModelPart,NewConditions);
            }


          }
      }

    void CSRRowMatrix(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord
    ) override
    {
        KRATOS_TRY;

        NodesArrayType& pNodes = this_model_part.Nodes();
        NodesArrayType::iterator it_begin = pNodes.ptr_begin();
        NodesArrayType::iterator it_end   = pNodes.ptr_end();

        Coord.resize(pNodes.size(), pNodes.size(), false);

        for(NodesArrayType::iterator i = it_begin; i!=it_end; i++)
        {
            int index_i = i->Id() - 1; // WARNING: MESH MUST BE IN ORDER
            GlobalPointersVector< Node < 3 > >& neighb_nodes = i->GetValue(NEIGHBOUR_CONDITION_NODES);

            std::vector<unsigned int> aux(neighb_nodes.size());
            unsigned int active = 0;
            for (GlobalPointersVector< Node < 3 > >::iterator inode = neighb_nodes.begin(); inode != neighb_nodes.end(); inode++)
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

};

} // namespace Kratos.

#endif // KRATOS_LOCAL_REFINE_TRIANGLE_MESH_CONDITIONS  defined
