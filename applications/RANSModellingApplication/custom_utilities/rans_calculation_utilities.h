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
#include "utilities/geometrical_sensitivity_utility.h"

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
    using NodeType = ModelPart::NodeType;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

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

    void CalculateGeometryParameterDerivativesShapeSensitivity(Matrix& rOutput,
                                                               const ShapeParameter& rShapeDerivative,
                                                               const Matrix& rDnDe,
                                                               const Matrix& rDeDx);

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

    template <unsigned int TDim>
    void CalculateGradient(BoundedMatrix<double, TDim, TDim>& rOutput,
                           const Geometry<ModelPart::NodeType>& rGeometry,
                           const Variable<array_1d<double, 3>>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const;

    void CalculateGradient(array_1d<double, 3>& rOutput,
                           const Geometry<ModelPart::NodeType>& rGeometry,
                           const Variable<double>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const;

    template <unsigned int TDim>
    Vector GetVector(const array_1d<double, 3>& rVector) const;

    Vector GetVector(const array_1d<double, 3>& rVector, const unsigned int Dim) const;

    double CalculateLogarithmicYPlusLimit(const double Kappa,
                                          const double Beta,
                                          const int MaxIterations = 20,
                                          const double Tolerance = 1e-6);

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