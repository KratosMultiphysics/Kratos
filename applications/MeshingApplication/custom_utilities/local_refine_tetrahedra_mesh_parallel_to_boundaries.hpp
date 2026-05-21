// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pablo Becker
//

#pragma once

// NOTE: Before compute the remeshing it is necessary to compute the neighbours

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Applicaion includes
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
    explicit LocalRefineTetrahedraMeshParallelToBoundaries(ModelPart& rModelPart)
     : LocalRefineTetrahedraMesh(rModelPart)
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

        //using the conditions to mark the boundary with the flag boundary
        //note that we DO NOT add the conditions to the model part
        block_for_each(mModelPart.Nodes(), [&](Node& rNode) {
            rNode.Set(BOUNDARY,false);
        });

        block_for_each(mModelPart.Conditions(), [&](ModelPart::ConditionType& rCondition) {
            Geometry< Node >& geom = rCondition.GetGeometry();
            for(unsigned int i=0; i<geom.size(); i++){
                geom[i].SetLock();
                geom[i].Set(BOUNDARY,true);
                geom[i].UnSetLock();
            } 
        });

        LocalRefineGeometryMesh::LocalRefineMesh(RefineOnReference,InterpolateInternalVariables);


        KRATOS_CATCH("");
    }
    ///@}
    ///@name Operations
    ///@{

    /***********************************************************************************/
    /***********************************************************************************/

    void ResetFatherNodes(ModelPart &rModelPart) override
    {
        block_for_each(mModelPart.Nodes(), [&](Node& rNode) {
            if(rNode.GetValue(REFINEMENT_LEVEL)==0){
                GlobalPointersVector<Node> empty_father_vector;
                rNode.SetValue(FATHER_NODES, empty_father_vector);
            }
        });
    }

    void SearchEdgeToBeRefined(
            ModelPart& rThisModelPart,
            compressed_matrix<int>& rCoord
    ) override
    {
        KRATOS_TRY;
        for (auto& r_elem: rThisModelPart.Elements()) {
            if (r_elem.GetValue(SPLIT_ELEMENT)) {
                Element::GeometryType& r_geom = r_elem.GetGeometry(); // Nodes of the element
                for (unsigned int i = 0; i < r_geom.size(); i++) {
                    int index_i = mMapNodeIdToPos[r_geom[i].Id()];
                    bool is_boundary_i = r_geom[i].Is(BOUNDARY);
                    for (unsigned int j = 0; j < r_geom.size(); j++) {
                        int index_j = mMapNodeIdToPos[r_geom[j].Id()];
                        bool is_boundary_j = r_geom[j].Is(BOUNDARY);
                        if (index_j > index_i && (is_boundary_j||is_boundary_i)) {
                            rCoord(index_i, index_j) = -2;
                        }
                    }
                }
            }
        }

        //unmarking edges belonging to the edges of conditions (skin) to avoid refining edges
        for (auto& r_cond : rThisModelPart.Conditions()) {
            auto& r_geom = r_cond.GetGeometry(); // Nodes of the condition
            for (unsigned int i = 0; i < r_geom.size(); i++) {
                int index_i = mMapNodeIdToPos[r_geom[i].Id()];
                for (unsigned int j = 0; j < r_geom.size(); j++) {
                    int index_j = mMapNodeIdToPos[r_geom[j].Id()];
                    if (index_j > index_i)  {
                        rCoord(index_i, index_j) = -1;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    //have to override because we need the refinement level in the ifs ( (*iNode)->GetValue(REFINEMENT_LEVEL)==mcurrent_refinement_level )
    void UpdateSubModelPartNodes(ModelPart &rModelPart) override
    {
        bool added_nodes=false;

        for (auto it_submodel_part = rModelPart.SubModelPartsBegin();
                it_submodel_part != rModelPart.SubModelPartsEnd(); it_submodel_part++) {
            added_nodes=false;
            for (auto it_node = rModelPart.Nodes().ptr_begin();
                    it_node != rModelPart.Nodes().ptr_end(); it_node++) {
                auto &r_father_nodes = (*it_node)->GetValue(FATHER_NODES);
                unsigned int parent_count = r_father_nodes.size();
                if (parent_count > 0 && (*it_node)->GetValue(REFINEMENT_LEVEL)==mCurrentRefinementLevel) {
                    unsigned int parents_in_submodel_part = 0;

                    for (auto it_parent = r_father_nodes.begin(); it_parent != r_father_nodes.end(); it_parent++) {
                        unsigned int parent_id = it_parent->Id();
                        ModelPart::NodeIterator iFound = it_submodel_part->Nodes().find( parent_id );
                        if ( iFound != it_submodel_part->NodesEnd() ) {
                            parents_in_submodel_part++;
                        }
                    }

                    if ( parent_count == parents_in_submodel_part ) {
                        it_submodel_part->AddNode( *it_node );
                        added_nodes=true;
                    }
                }
            }
            if(added_nodes) {
                 ModelPart &rSubModelPart = *it_submodel_part;
                 UpdateSubModelPartNodes(rSubModelPart);
            }
        }
    }


protected:
    ///@name Protected static Member Variables
    ///@{
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
