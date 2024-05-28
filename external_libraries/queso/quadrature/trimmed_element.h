// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MOMENT_FITTING_UTILITIES_INCLUDE_H
#define MOMENT_FITTING_UTILITIES_INCLUDE_H

//// STL includes
#include <vector>
#include <array>
#include <variant>
//// Project includes
#include "embedding/octree.h"
#include "containers/element.hpp"
#include "containers/boundary_integration_point.hpp"
#include "includes/parameters.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  QuadratureTrimmedElement.
 * @author Manuel Messmer
 * @brief  Provides functions to create integration rules for trimmed elements.
 * @tparam TElementType
**/
template<typename TElementType>
class QuadratureTrimmedElement{
public:
    ///@name Type Definition
    ///@{
    typedef TElementType ElementType;
    typedef typename ElementType::IntegrationPointVectorType IntegrationPointVectorType;
    typedef Unique<IntegrationPointVectorType> IntegrationPointVectorPtrType;
    typedef typename ElementType::BoundaryIntegrationPointType BoundaryIntegrationPointType;
    typedef std::vector<BoundaryIntegrationPointType> BoundaryIPsVectorType;
    typedef Unique<BoundaryIPsVectorType> BoundaryIPsVectorPtrType;
    typedef std::vector<double> VectorType;

    ///@}
    ///@name Operations
    ///@{

    ///@brief Creates integration points for trimmed domain.
    ///@details 1. Distributes initial integration points uniformly in trimmed domain.
    ///         2. Computes constant terms of moment fitting equation.
    ///         3. Solves moment fitting equation in iterative point elimination algorithm.
    /// See: M. Meßmer et. al: Efficient CAD-integrated isogeometric analysis of trimmed solids,
    ///      Comput. Methods Appl. Mech. Engrg. 400 (2022) 115584, https://doi.org/10.1016/j.cma.2022.115584.
    ///@param rElement
    ///@param rIntegrationOrder
    ///@param Residual Targeted residual
    ///@param EchoLevel Default: 0
    static double AssembleIPs(ElementType& rElement, const Vector3i& rIntegrationOrder, double Residual, IndexType EchoLevel=0);

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    /// @brief Distributes point within trimmed domain using an octree. In each leaf node, Gauss points according to rIntegrationOrder are generated.
    ///        Only points inside the trimmed domain are considered.
    ///        Every time this function is called the otree is refined and more points are distributed.
    /// @param[out] rIntegrationPoint
    /// @param rOctree
    /// @param MinNumPoints Minimum Number of Points
    /// @param rIntegrationOrder Order of Gauss quadrature.
    static void DistributeIntegrationPoints(IntegrationPointVectorType& rIntegrationPoint, Octree<TrimmedDomain>& rOctree, SizeType MinNumPoints, const Vector3i& rIntegrationOrder);

    /// @brief Computes constant terms of moment fitting equation via volume integration points.
    /// @param[out] rConstantTerms
    /// @param pIntegrationPoints (Unique<T>)
    /// @param rElement
    /// @param rIntegrationOrder
    static void ComputeConstantTerms(VectorType& rConstantTerms, const IntegrationPointVectorPtrType& pIntegrationPoints, const ElementType& rElement, const Vector3i& rIntegrationOrder);

    /// @brief Computes constant terms of moment fitting equation via boundary integration points. This functions uses the divergence theorem
    //         to transform the respective volume integrals to countour/surface integrals.
    /// @param[out] rConstantTerms
    /// @param pBoundaryIPs (Unique<T>)
    /// @param rElement
    /// @param rIntegrationOrder
    static void ComputeConstantTerms(VectorType& rConstantTerms, const BoundaryIPsVectorPtrType& pBoundaryIPs, const ElementType& rElement, const Vector3i& rIntegrationOrder);

    /// @brief Set-Up and solve moment fitting equation. Solve the moment fitting equation for the weights of the integration points.
    ///        Computed weights are directly assigned to rIntegrationPoint.
    /// @param rConstantTerms
    /// @param[out] rIntegrationPoint
    /// @param rElement
    /// @param rIntegrationOrder
    /// @return double Relative residual ||ax -b||_L2 / ||b||_L2
    static double MomentFitting(const VectorType& rConstantTerms, IntegrationPointVectorType& rIntegrationPoint, const ElementType& rElement, const Vector3i& rIntegrationOrder);

    /// @brief Start point elimination algorihtm. Final quadrature rule is stored in rElement.
    /// @param rConstantTerms
    /// @param rIntegrationPoint
    /// @param rElement
    /// @param rIntegrationOrder
    /// @param Residual targeted residual
    /// @return double achieved residual
    static double PointElimination(const VectorType& rConstantTerms, IntegrationPointVectorType& rIntegrationPoint, ElementType& rElement, const Vector3i& rIntegrationOrder, double Residual);
}; // End Class


} // End namespace queso

#endif // MOMENT_FITTING_UTILITIES_INCLUDE_H