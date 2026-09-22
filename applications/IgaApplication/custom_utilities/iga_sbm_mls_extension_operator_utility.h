//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "containers/array_1d.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Utilities for the IGA SBM-MLS coupling extension operator E_h.
/**
*   This class computes the Moving-Least-Squares extension operator used to
*   reconstruct a field (and its gradient) at true boundary
*   from a cloud of nearby ACTIVE quadrature points.
*/
class KRATOS_API(IGA_APPLICATION) IgaSbmMlsExtensionOperatorUtility
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    /// Assembled generalized extension operator
    struct ExtensionOperatorResult
    {
        std::vector<IndexType> NodeIds;
        Vector Weights;
        Matrix GradientWeights; 
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Finds the indices of the N closest points to rEvalPoint.
    * @param rCandidatePoints Candidate cloud points, in parametric space 
    * @param rEvalPoint Point to search around, in parametric space
    * @param NClosest Number of closest points to return 
    * @return Indices into rCandidatePoints of the NClosest closest points, sorted by increasing distance
    */
    static std::vector<IndexType> FindNClosestPoints(
        const std::vector<array_1d<double, 3>>& rCandidatePoints,
        const array_1d<double, 3>& rEvalPoint,
        const SizeType NClosest);

    /**
    * @brief Computes the MLS extension operator weights (value and gradient) at one point.
    * @param rCloudPointsParametric Cloud points, in parametric space
    * @param rEvalPointParametric Evaluation point, in parametric space
    * @param MLSOrder Polynomial order of the MLS basis (1 or 2)
    * @param rN Output: MLS weight for each cloud point 
    * @param rDN_DLocal Output: gradient of the MLS weights w.r.t. (xi, eta) for each cloud point
    */
    static void ComputeExtensionOperator(
        const Matrix& rCloudPointsParametric,
        const array_1d<double, 3>& rEvalPointParametric,
        const SizeType MLSOrder,
        Vector& rN,
        Matrix& rDN_DLocal);

    /**
    * @brief Classifies rDomainElements into active (safe for the MLS cloud) and excluded (buffer/intersected).
    * @param rDomainElements Candidate elements
    * @param rEvalPointParametric Point to exclude a buffer band around
    * @param KBuffer Number of nearest elements to exclude
    * @return The remaining (active) elements
    */
    static std::vector<Element::Pointer> ClassifyActiveElements(
        const std::vector<Element::Pointer>& rDomainElements,
        const array_1d<double, 3>& rEvalPointParametric,
        const SizeType KBuffer);

    /**
    * @brief Computes the assembled (per-control-point) extension operator at one Gamma_D point.
    * @param rActiveElements Elements to draw the cloud from
    * @param rEvalPointParametric true boundary point, in this patch's own PARAMETRIC space
    * @param NClosest Number of cloud quadrature points to start from 
    * @param MLSOrder Polynomial order of the MLS basis (1 or 2)
    * @return Assembled extension operator
    */
    static ExtensionOperatorResult ComputeAssembledExtensionOperator(
        const std::vector<Element::Pointer>& rActiveElements,
        const array_1d<double, 3>& rEvalPointParametric,
        const SizeType NClosest,
        const SizeType MLSOrder);

    /**
    * @brief Precomputes and stores the MLS extension operator on every condition
    * @param rCouplingModelPart Model part holding the coupling conditions
    * @param rMasterDomainElements Master patch's own domain elements 
    * @param rSlaveDomainElements Slave patch's own domain elements
    * @param KBuffer Number of nearest elements to exclude around each evaluation point
    * @param NClosest Number of cloud quadrature points to use per evaluation point
    * @param MLSOrder Polynomial order of the MLS basis (1 or 2)
    * @param rExcludedNodeIds Node ids to keep OUT of the MLS cloud entirely 
    */
    static void PrecomputeAndStoreCouplingExtensionOperators(
        ModelPart& rCouplingModelPart,
        const std::vector<Element::Pointer>& rMasterDomainElements,
        const std::vector<Element::Pointer>& rSlaveDomainElements,
        const SizeType KBuffer,
        const SizeType NClosest,
        const SizeType MLSOrder,
        const std::vector<IndexType>& rExcludedNodeIds);

    /**
    * @brief Debug/validation helper: returns the ids of SBM_MLS_ALL_DOF_NODES.
    * @param rCondition Condition previously processed by PrecomputeAndStoreCouplingExtensionOperators
    * @return Node ids, in the exact column order of the stored weight matrices
    */
    static std::vector<IndexType> GetStoredDofNodeIds(const Condition& rCondition);

    /**
    * @brief Wires each CouplingSbmExtensionOperator6pCondition's true boundary adjoint term to its
    *        own patch's surrogate boundary flux source.
    * @param rCouplingModelPart Model part holding the instances
    * @param rMasterSurrogateConditions Master patch's own surrogate conditions
    * @param rSlaveSurrogateConditions Slave patch's own surrogate conditions
    */
    static void PrecomputeAndStoreCouplingShiftSources(
        ModelPart& rCouplingModelPart,
        const std::vector<Condition::Pointer>& rMasterSurrogateConditions,
        const std::vector<Condition::Pointer>& rSlaveSurrogateConditions);

    /**
    * @param rElement Element with exactly one integration point 
    */
    static double GetDifferentialArea(const Element& rElement);

    ///@}

}; // Class IgaSbmMlsExtensionOperatorUtility

///@}

}  // namespace Kratos.
