//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined(KRATOS_GEOMETRY_METRIC_CALCULATOR_INCLUDED)
#define  KRATOS_GEOMETRY_METRIC_CALCULATOR_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

template<std::size_t TDim, std::size_t TNumNodes>
class KRATOS_API(KRATOS_CORE) GeometryMetricCalculator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryMetricCalculator
    KRATOS_CLASS_POINTER_DEFINITION(GeometryMetricCalculator);

    typedef Geometry<Node<3>> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Deleted default constructor
    GeometryMetricCalculator() = delete;

    /// Deleted copy constructor.
    GeometryMetricCalculator(GeometryMetricCalculator const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    GeometryMetricCalculator& operator=(GeometryMetricCalculator const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Calculate the metric tensor of a given geometry
     * This function calculates the metric tensor and its data for a given geometry.
     * The metric tensor M is computed by solving the problem trans(e)*M*e = 1, that means find
     * the coefficients of the matrix M such that all the geometry edges (e) have unit lenght.
     * @tparam TDim Geometry working dimension
     * @tparam TNumNodes Geometry number of nodes
     * @param rGeometry Reference to the geometry of interest
     * @param rMetricTensor Reference to the metric tensor
     */
    static void CalculateMetricTensor(
        const GeometryType& rGeometry,
        BoundedMatrix<double, TDim, TDim>& rMetricTensor);

    /**
     * @brief Calculate the metric tensor of a given geometry
     * This function calculates the metric tensor and its data for a given geometry.
     * The metric tensor M is computed by solving the problem trans(e)*M*e = 1, that means find
     * the coefficients of the matrix M such that all the geometry edges (e) have unit lenght.
     * The eigenvalues of the metric tensor are also computed in order to obtain the lenghts
     * of the Steiner inertia ellipsis semiaxes. This allows computing the reference element
     * size (as the average of the semiaxes lengths) and the infimum and supremum norms of M.
     * @tparam TDim Geometry working dimension
     * @tparam TNumNodes Geometry number of nodes
     * @param rGeometry Reference to the geometry of interest
     * @param rMetricTensor Reference to the metric tensor
     * @param rReferenceElementSize Reference to the reference element size
     * @param rMetricInfimum Reference to the metric infimum
     * @param rMetricSupremum Reference to the metric supremum
     */
    static void CalculateMetricTensor(
        const GeometryType& rGeometry,
        BoundedMatrix<double, TDim, TDim>& rMetricTensor,
        double& rReferenceElementSize,
        double& rMetricInfimum,
        double& rMetricSupremum);

    /**
     * @brief Calculate the metric tensor of a given geometry
     * This function calculates the metric tensor and its data for a given geometry.
     * The metric tensor M is computed by solving the problem trans(e)*M*e = 1, that means find
     * the coefficients of the matrix M such that all the geometry edges (e) have unit lenght.
     * The eigenvalues of the metric tensor are also computed in order to obtain the lenghts
     * of the Steiner inertia ellipsis semiaxes. This allows computing the reference element
     * size (as the average of the semiaxes lengths) and the infimum and supremum norms of M.
     * Note that this methods returns the dimensionless metric that is the metric tensor
     * multiplied by the square of the reference element size.
     * @tparam TDim Geometry working dimension
     * @tparam TNumNodes Geometry number of nodes
     * @param rGeometry Reference to the geometry of interest
     * @param rMetricTensor Reference to the metric tensor
     * @param rReferenceElementSize Reference to the reference element size
     * @param rMetricInfimum Reference to the metric infimum
     * @param rMetricSupremum Reference to the metric supremum
     */
    static void CalculateMetricTensorDimensionless(
        const GeometryType& rGeometry,
        BoundedMatrix<double, TDim, TDim>& rMetricTensorDimensionless,
        double& rReferenceElementSize,
        double& rMetricInfimum,
        double& rMetricSupremum);

    ///@}
}; // Class GeometryMetricCalculator

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_METRIC_CALCULATOR_INCLUDED  defined

