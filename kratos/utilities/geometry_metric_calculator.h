//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_GEOMETRY_UTILITIES_INCLUDED )
#define  KRATOS_GEOMETRY_UTILITIES_INCLUDED

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

class GeometryMetricCalculator
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

    /// Destructor.
    ~GeometryMetricCalculator();

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    GeometryMetricCalculator& operator=(GeometryMetricCalculator const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    template<unsigned int TDim, unsigned int TNumNodes>
    static void KRATOS_API(KRATOS_CORE) CalculateMetricTensorData(
        const GeometryType& rGeometry,
        BoundedMatrix<double, TDim, TDim>& rMetricTensor,
        double& rReferenceElementSize,
        double& rMetricInfimum,
        double& rMetricSupremum);

    ///@}

private:
    ///@name Private Operations
    ///@{


    ///@}
}; // Class GeometryMetricCalculator

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined


