// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#pragma once

// Project includes

// Application includes
#include "custom_elements/stress_state_policy.h"
#include "custom_retention/retention_law.h"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "geometries/geometry.h"
#include "includes/element.h"
#include "includes/smart_pointers.h"
#include "includes/variables.h"

namespace Kratos
{

class GeoTransportEquationUtilities
{
public:
    template <unsigned int TDim, unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix(
        const Matrix&                            rGradNpT,
        double                                   DynamicViscosityInverse,
        const BoundedMatrix<double, TDim, TDim>& rMaterialPermeabilityMatrix,
        double                                   RelativePermeability,
        double                                   PermeabilityUpdateFactor,
        double                                   IntegrationCoefficient)
    {
        return CalculatePermeabilityMatrix(rGradNpT, DynamicViscosityInverse,
                                           rMaterialPermeabilityMatrix, RelativePermeability,
                                           PermeabilityUpdateFactor, IntegrationCoefficient);
    }

    static inline Matrix CalculatePermeabilityMatrix(const Matrix& rGradNpT,
                                                     double        DynamicViscosityInverse,
                                                     const Matrix& rMaterialPermeabilityMatrix,
                                                     double        RelativePermeability,
                                                     double        PermeabilityUpdateFactor,
                                                     double        IntegrationCoefficient)
    {
        return -PORE_PRESSURE_SIGN_FACTOR * DynamicViscosityInverse *
               prod(rGradNpT, Matrix(prod(rMaterialPermeabilityMatrix, trans(rGradNpT)))) *
               RelativePermeability * PermeabilityUpdateFactor * IntegrationCoefficient;
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes * TDim, TNumNodes> CalculateCouplingMatrix(
        const Matrix& rB, const Vector& rVoigtVector, const Vector& rNp, double BiotCoefficient, double BishopCoefficient, double IntegrationCoefficient)
    {
        return CalculateCouplingMatrix(rB, rVoigtVector, rNp, BiotCoefficient, BishopCoefficient,
                                       IntegrationCoefficient);
    }

    static inline Matrix CalculateCouplingMatrix(const Matrix& rB,
                                                 const Vector& rVoigtVector,
                                                 const Vector& rNp,
                                                 double        BiotCoefficient,
                                                 double        BishopCoefficient,
                                                 double        IntegrationCoefficient)
    {
        return PORE_PRESSURE_SIGN_FACTOR * BiotCoefficient * BishopCoefficient *
               outer_prod(Vector(prod(trans(rB), rVoigtVector)), rNp) * IntegrationCoefficient;
    }

    template <unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TNumNodes> CalculateCompressibilityMatrix(
        const Vector& rNp, double BiotModulusInverse, double IntegrationCoefficient)
    {
        return CalculateCompressibilityMatrix(rNp, BiotModulusInverse, IntegrationCoefficient);
    }

    static inline Matrix CalculateCompressibilityMatrix(const Vector& rNp, double BiotModulusInverse, double IntegrationCoefficient)
    {
        return -PORE_PRESSURE_SIGN_FACTOR * BiotModulusInverse * outer_prod(rNp, rNp) * IntegrationCoefficient;
    }

    static double CalculateSoilDensity(double DegreeOfSaturation, const Properties& rProp)
    {
        return (DegreeOfSaturation * rProp[POROSITY] * rProp[DENSITY_WATER]) +
               (1.0 - rProp[POROSITY]) * rProp[DENSITY_SOLID];
    }

    static Vector CalculateSoilDensityVector(const Geometry<Node>& rGeom,
                                             unsigned int          NumberOfIntegrationPoints,
                                             const Matrix&         rNContainer,
                                             const std::vector<RetentionLaw::Pointer>& rRetentionLawVector,
                                             const Properties&  rProp,
                                             const ProcessInfo& rCurrentProcessInfo)
    {
        const Vector pressure_vector = GeoTransportEquationUtilities::GetSolutionVector(rGeom, WATER_PRESSURE);
        RetentionLaw::Parameters retention_parameters(rProp, rCurrentProcessInfo);
        Vector                   density(NumberOfIntegrationPoints);
        for (unsigned int g_point = 0; g_point < NumberOfIntegrationPoints; ++g_point) {
            const double degree_of_saturation =
                CalculateDegreeOfSaturation(row(rNContainer, g_point), pressure_vector,
                                            retention_parameters, rRetentionLawVector[g_point]);

            density(g_point) = CalculateSoilDensity(degree_of_saturation, rProp);
        }
        return density;
    }

private:
    static double CalculateIntegrationCoefficientInitialConfiguration(
        unsigned int                                      GPoint,
        const Geometry<Node>&                             rGeom,
        const Geometry<Node>::IntegrationPointsArrayType& rIntegrationPoints,
        const GeometryData::IntegrationMethod             IntegrationMethod,
        const StressStatePolicy&                          rStressStatePolicy)
    {
        Matrix J0;
        Matrix inv_J0;
        GeometryUtils::JacobianOnInitialConfiguration(rGeom, rIntegrationPoints[GPoint], J0);
        const Matrix& dN_De = rGeom.ShapeFunctionsLocalGradients(IntegrationMethod)[GPoint];
        double        det_J_initial_configuration;
        MathUtils<double>::InvertMatrix(J0, inv_J0, det_J_initial_configuration);
        Matrix dNu_dX_initial_configuration;
        GeometryUtils::ShapeFunctionsGradients(dN_De, inv_J0, dNu_dX_initial_configuration);

        return rStressStatePolicy.CalculateIntegrationCoefficient(
            rIntegrationPoints[GPoint], det_J_initial_configuration, rGeom);
    }

    static double CalculateDegreeOfSaturation(const Vector&                rNVector,
                                              const Vector&                rPressureVector,
                                              RetentionLaw::Parameters&    rRetentionParameters,
                                              const RetentionLaw::Pointer& rRetentionLaw)
    {
        const double fluid_pressure =
            GeoTransportEquationUtilities::CalculateFluidPressure(rNVector, rPressureVector);
        rRetentionParameters.SetFluidPressure(fluid_pressure);

        return rRetentionLaw->CalculateSaturation(rRetentionParameters);
    }

    static double CalculateFluidPressure(const Vector& rNp, const Vector& rPressureVector)
    {
        return inner_prod(rNp, rPressureVector);
    }

    static Vector GetSolutionVector(const Geometry<Node>& rGeom, const Variable<double>& rSolutionVariable)
    {
        Vector solution_vector(rGeom.size());
        std::transform(rGeom.begin(), rGeom.end(), solution_vector.begin(),
                       [&rSolutionVariable](const auto& node) {
            return node.FastGetSolutionStepValue(rSolutionVariable);
        });
        return solution_vector;
    }
}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
