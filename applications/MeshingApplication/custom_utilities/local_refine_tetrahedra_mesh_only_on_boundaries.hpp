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

// Application includes
#include "custom_utilities/local_refine_tetrahedra_mesh_parallel_to_boundaries.hpp"

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
class LocalRefineTetrahedraMeshOnlyOnBoundaries : public LocalRefineTetrahedraMeshParallelToBoundaries
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    explicit LocalRefineTetrahedraMeshOnlyOnBoundaries(ModelPart& rModelPart)
     : LocalRefineTetrahedraMeshParallelToBoundaries(rModelPart)
    {
    }

    /// Destructor
    ~LocalRefineTetrahedraMeshOnlyOnBoundaries() //TODO maybe {}
    = default;

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
                        if (index_j > index_i && (is_boundary_j&&is_boundary_i)) {
                            rCoord(index_i, index_j) = -2;
                        }
                    }
                }
            }
        }

        // check for conditions with SPLIT_ELEMENT flag to refine the corresponding edges
        for (auto& r_cond: rThisModelPart.Conditions()) {
            if (r_cond.GetValue(SPLIT_ELEMENT)) {
                Condition::GeometryType& r_geom = r_cond.GetGeometry(); // Nodes of the condition
                for (unsigned int i = 0; i < r_geom.size(); i++) {
                    int index_i = mMapNodeIdToPos[r_geom[i].Id()];  
                    for (unsigned int j = 0; j < r_geom.size(); j++) {
                        int index_j = mMapNodeIdToPos[r_geom[j].Id()];
                        if (index_j > index_i) {
                            rCoord(index_i, index_j) = -2;
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("");
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
