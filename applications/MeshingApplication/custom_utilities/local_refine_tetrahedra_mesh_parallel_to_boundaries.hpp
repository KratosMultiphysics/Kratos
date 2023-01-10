// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Pablo Becker
//

#if !defined(KRATOS_LOCAL_REFINE_TETRAHEDRA_MESH_PARALLEL_TO_BOUNDARIES)
#define  KRATOS_LOCAL_REFINE_TETRAHEDRA_MESH_PARALLEL_TO_BOUNDARIES

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

/* Project includes */
#include "utilities/parallel_utilities.h"
#include "custom_utilities/local_refine_tetrahedra_mesh.hpp"


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
class LocalRefineTetrahedraMeshParallelToBoundaries : public LocalRefineTetrahedraMesh
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LocalRefineTetrahedraMeshParallelToBoundaries(ModelPart& model_part) : LocalRefineTetrahedraMesh(model_part)
    {

    }

    /// Destructor
    ~LocalRefineTetrahedraMeshParallelToBoundaries() //TODO maybe {}
    = default;

    ///@}
    ///@name Operators
    ///@{
    void LocalRefineMesh(
            bool RefineOnReference,
            bool InterpolateInternalVariables
    ) override
    {
        KRATOS_TRY;

        if (RefineOnReference)
        {
            KRATOS_ERROR_IF_NOT(mModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT)) << "Missing DISPLACEMENT variable on solution step data." ;
        }

        compressed_matrix<int> Coord;                                              // The matrix that stores all the index of the geometry
        boost::numeric::ublas::vector<int> List_New_Nodes;                         // The news nodes
        boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;           // Edges where are the news nodes
        boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node; // The coordinate of the new nodes

        PointerVector< Element > New_Elements;
        New_Elements.reserve(20);
        
        // Initial renumber of nodes and elemetns
        //TODO : DIFFERENCE WITH BASE CLASS 1
        mPreviousRefinementLevel=0;
        unsigned int id = 1;
        for (ModelPart::NodesContainerType::iterator it = mModelPart.NodesBegin(); it != mModelPart.NodesEnd(); it++)
        {
            it->SetId(id++);
            int node_refinement_level = it->GetValue(REFINEMENT_LEVEL);
            if(node_refinement_level>mPreviousRefinementLevel)
                    mPreviousRefinementLevel=node_refinement_level;
        }
        mCurrentRefinementLevel = mPreviousRefinementLevel+1;
	
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
	
        if (RefineOnReference)
        {
            block_for_each(mModelPart.Nodes(), [&](Node<3>& rNode)
            {
                rNode.X() = rNode.X0();
                rNode.Y() = rNode.Y0();
                rNode.Z() = rNode.Z0();
            });
        }

        this->ResetFatherNodes(mModelPart); 

        //TODO : DIFFERENCE WITH BASE CLASS 2
        //using the conditions to mark the boundary with the flag boundary
        //note that we DO NOT add the conditions to the model part
        // we also temporarily substract -100 to be able to spot the new ones:
        block_for_each(mModelPart.Nodes(), [&](Node<3>& rNode)
        {
            rNode.GetValue(REFINEMENT_LEVEL)-=100;
            rNode.Set(BOUNDARY,false);
        });
        block_for_each(mModelPart.Conditions(), [&](ModelPart::ConditionType& rCondition) {
            Geometry< Node<3> >& geom = rCondition.GetGeometry();
            for(unsigned int i=0; i<geom.size(); i++){
                geom[i].SetLock();
                geom[i].Set(BOUNDARY,true);
                geom[i].UnSetLock();
            } 
        });

        /* Calling all the functions necessaries to refine the mesh */
        CSRRowMatrix(mModelPart, Coord);

        SearchEdgeToBeRefined(mModelPart, Coord);

        CreateListOfNewNodes(mModelPart, Coord, List_New_Nodes, Position_Node);

        CalculateCoordinateAndInsertNewNodes(mModelPart, Position_Node, List_New_Nodes);

        EraseOldElementAndCreateNewElement(mModelPart, Coord, New_Elements, InterpolateInternalVariables);

        EraseOldConditionsAndCreateNew(mModelPart, Coord);

        RenumeringElementsAndNodes(mModelPart, New_Elements);

        //TODO : DIFFERENCE WITH BASE CLASS 3. fixing refinement level on new and old nodes:
        block_for_each(mModelPart.Nodes(), [&](Node<3>& rNode)
        {
            if(!rNode.Has(REFINEMENT_LEVEL)){
                rNode.SetValue(REFINEMENT_LEVEL,mCurrentRefinementLevel);
            }
            else{
                rNode.GetValue(REFINEMENT_LEVEL)+=100;
            }
        });
        
        
        if (RefineOnReference)
        {
            block_for_each(mModelPart.Nodes(), [&](Node<3>& rNode)
            {
                const array_1d<double, 3 > & disp = rNode.FastGetSolutionStepValue(DISPLACEMENT);
                rNode.X() = rNode.X0() + disp[0];
                rNode.Y() = rNode.Y0() + disp[1];
                rNode.Z() = rNode.Z0() + disp[2];
            });
        }

        UpdateSubModelPartNodes(mModelPart);

        KRATOS_CATCH("");
    }
    ///@}
    ///@name Operations
    ///@{

    /***********************************************************************************/
    /***********************************************************************************/

    void ResetFatherNodes(ModelPart &rModelPart) override
    {
        block_for_each(mModelPart.Nodes(), [&](Node<3>& rNode)
        {
            if(rNode.GetValue(REFINEMENT_LEVEL)==0){
                GlobalPointersVector<Node<3>> empty_father_vector;
                rNode.SetValue(FATHER_NODES, empty_father_vector);
            }
        });
    }

    void SearchEdgeToBeRefined(
            ModelPart& this_model_part,
            compressed_matrix<int>& Coord
    ) override
    {
        KRATOS_TRY;
	
        ElementsArrayType& rElements = this_model_part.Elements();
        ElementsArrayType::iterator it_begin = rElements.ptr_begin();
        ElementsArrayType::iterator it_end   = rElements.ptr_end();

        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it)
        {
            if (it->GetValue(SPLIT_ELEMENT))
            {
                Element::GeometryType& geom = it->GetGeometry(); // Nodes of the element
                for (unsigned int i = 0; i < geom.size(); i++)
                {
                    int index_i = geom[i].Id() - 1;
                    bool is_boundary_i = geom[i].Is(BOUNDARY);
                    for (unsigned int j = 0; j < geom.size(); j++)
                    {
                        int index_j = geom[j].Id() - 1;
                        bool is_boundary_j = geom[j].Is(BOUNDARY);
                        //if (index_j > index_i && is_boundary_j!=is_boundary_i) //old version, only edges that join internal and external nodes
                        if (index_j > index_i && (is_boundary_j||is_boundary_i)) // new version, for single elem in thickness meshes
                        {
                            Coord(index_i, index_j) = -2;
                        }
                    }
                }
            }
        }

        //unmarking edges belonging to the edges of conditions (skin) to avoid refining edges
        ConditionsArrayType& rConditions = this_model_part.Conditions();
        ConditionsArrayType::iterator it_begin_cond = rConditions.ptr_begin();
        ConditionsArrayType::iterator it_end_cond   = rConditions.ptr_end();

        for (ConditionsArrayType::iterator it = it_begin_cond; it != it_end_cond; ++it)
        {
            Condition::GeometryType& geom = it->GetGeometry(); // Nodes of the condition
            for (unsigned int i = 0; i < geom.size(); i++)
            {
                    int index_i = geom[i].Id() - 1;
                    for (unsigned int j = 0; j < geom.size(); j++)
                    {
                        int index_j = geom[j].Id() - 1;
                        if (index_j > index_i) 
                        {
                            Coord(index_i, index_j) = -1;
                        }
                    }
            }
        }

        KRATOS_CATCH("");
    }

    //have to override because we need the refinement level in the ifs ( (*iNode)->GetValue(REFINEMENT_LEVEL)==mcurrent_refinement_level )
    void UpdateSubModelPartNodes(ModelPart &rModelPart) override
    {
        KRATOS_WATCH(rModelPart.Name())

        bool added_nodes=false;

        for (ModelPart::SubModelPartIterator iSubModelPart = rModelPart.SubModelPartsBegin();
                iSubModelPart != rModelPart.SubModelPartsEnd(); iSubModelPart++)
        {
            
            added_nodes=false;
            for (auto iNode = rModelPart.Nodes().ptr_begin();
                    iNode != rModelPart.Nodes().ptr_end(); iNode++)
            {
                GlobalPointersVector< Node<3> > &rFatherNodes = (*iNode)->GetValue(FATHER_NODES);
                unsigned int ParentCount = rFatherNodes.size();

                if (ParentCount > 0 && (*iNode)->GetValue(REFINEMENT_LEVEL)==mCurrentRefinementLevel)
                {
                    unsigned int ParentsInSubModelPart = 0;

                    for ( GlobalPointersVector< Node<3> >::iterator iParent = rFatherNodes.begin();
                            iParent != rFatherNodes.end(); iParent++)
                    {
                        unsigned int ParentId = iParent->Id();
                        ModelPart::NodeIterator iFound = iSubModelPart->Nodes().find( ParentId );
                        if ( iFound != iSubModelPart->NodesEnd() )
                            ParentsInSubModelPart++;
                    }

                    if ( ParentCount == ParentsInSubModelPart )
                    {
                        iSubModelPart->AddNode( *iNode );
                        added_nodes=true;
                    }
                }
            }
            if(added_nodes)
            {
                 ModelPart &rSubModelPart = *iSubModelPart;
                 UpdateSubModelPartNodes(rSubModelPart);
            }    
        }
    }


protected:
    ///@name Protected static Member Variables
    ///@{
    int mPreviousRefinementLevel;
    int mCurrentRefinementLevel;
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

#endif // KRATOS_LOCAL_REFINE_TETRAHEDRA_MESH_PARALLEL_TO_BOUNDARIES  defined
