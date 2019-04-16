//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED

#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/model_part.h"

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

class RansCalculationUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// We create the Pointer related to VariableUtils
    KRATOS_CLASS_POINTER_DEFINITION(RansCalculationUtilities);

    /// Node type
    typedef ModelPart::NodeType NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    RansCalculationUtilities()
    {
    }

    /**
     * Destructor
     */
    ~RansCalculationUtilities()
    {
    }

    ///@}

    ///@name Operations
    ///@{

    void CalculateGeometryData(const GeometryType& rGeometry,
                               const GeometryData::IntegrationMethod& rIntegrationMethod,
                               Vector& rGaussWeights,
                               Matrix& rNContainer,
                               GeometryType::ShapeFunctionsGradientsType& rDN_DX);

    GeometryType::ShapeFunctionsGradientsType CalculateGeometryParameterDerivatives(
        const GeometryType& rGeometry, const GeometryData::IntegrationMethod& rIntegrationMethod);

    double EvaluateInPoint(const GeometryType& rGeometry,
                           const Variable<double>& rVariable,
                           const Vector& rShapeFunction,
                           const int Step = 0);

    array_1d<double, 3> EvaluateInPoint(const GeometryType& rGeometry,
                                        const Variable<array_1d<double, 3>>& rVariable,
                                        const Vector& rShapeFunction,
                                        const int Step = 0);

    template <unsigned int TDim>
    double CalculateMatrixTrace(const BoundedMatrix<double, TDim, TDim>& rMatrix);

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}

    ///@name Private Operations
    ///@{

    ///@}

}; // Class CalculationUtilities

///@}

} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CALCULATION_UTILITIES_H_INCLUDED defined