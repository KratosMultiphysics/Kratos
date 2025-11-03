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
#include "custom_retention/retention_law.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/variables_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/kratos_export_api.h"
#include "includes/variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoTransportEquationUtilities
{
public:
    template <unsigned int TDim, unsigned int TNumNodes>
    static inline BoundedMatrix<double, TNumNodes, TNumNodes> CalculatePermeabilityMatrix(
        const Matrix&                            rGradNpT,
        double                                   DynamicViscosityInverse,
        const BoundedMatrix<double, TDim, TDim>& rMaterialPermeabilityMatrix,
        double                                   RelativePermeability,
        double                                   IntegrationCoefficient)
    {
        return -PORE_PRESSURE_SIGN_FACTOR * DynamicViscosityInverse *
               prod(rGradNpT, BoundedMatrix<double, TDim, TNumNodes>(
                                  prod(rMaterialPermeabilityMatrix, trans(rGradNpT)))) *
               RelativePermeability * IntegrationCoefficient;
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static std::vector<array_1d<double, TDim>> CalculateFluidFluxes(
        const Geometry<Node>&                     rGeometry,
        GeometryData::IntegrationMethod           IntegrationMethod,
        Properties const&                         rProperties,
        const std::vector<RetentionLaw::Pointer>& rRetentionLawVector,
        const std::vector<double>&                rPermeabilityUpdateFactors)
    {
        const IndexType number_of_integration_points = rGeometry.IntegrationPointsNumber(IntegrationMethod);

        std::vector<array_1d<double, TDim>> fluid_fluxes;
        fluid_fluxes.reserve(number_of_integration_points);
        array_1d<double, TNumNodes> pressure_vector;
        VariablesUtilities::GetNodalValues(rGeometry, WATER_PRESSURE, pressure_vector.begin());
        Matrix N_container(number_of_integration_points, TNumNodes);
        N_container = rGeometry.ShapeFunctionsValues(IntegrationMethod);
        BoundedMatrix<double, TDim, TDim> permeability_matrix;
        GeoElementUtilities::FillPermeabilityMatrix(permeability_matrix, rProperties);

        auto relative_permeability_values = RetentionLaw::CalculateRelativePermeabilityValues(
            rRetentionLawVector, rProperties,
            GeoTransportEquationUtilities::CalculateFluidPressures(N_container, pressure_vector));
        std::ranges::transform(relative_permeability_values, rPermeabilityUpdateFactors,
                               relative_permeability_values.begin(), std::multiplies<>{});
        const auto dynamic_viscosity_inverse = 1.0 / rProperties[DYNAMIC_VISCOSITY];
        array_1d<double, TNumNodes * TDim> volume_acceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(volume_acceleration, rGeometry,
                                                                     VOLUME_ACCELERATION);
        array_1d<double, TDim>                      body_acceleration;
        Matrix                                      grad_Np_T(TNumNodes, TDim);
        Vector                                      det_J_Container(number_of_integration_points);
        Geometry<Node>::ShapeFunctionsGradientsType dN_dx_container;
        rGeometry.ShapeFunctionsIntegrationPointsGradients(dN_dx_container, det_J_Container, IntegrationMethod);
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            noalias(grad_Np_T) = dN_dx_container[integration_point];

            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                body_acceleration, N_container, volume_acceleration, integration_point);

            array_1d<double, TDim> GradPressureTerm = prod(trans(grad_Np_T), pressure_vector);
            GradPressureTerm += PORE_PRESSURE_SIGN_FACTOR * rProperties[DENSITY_WATER] * body_acceleration;

            fluid_fluxes.push_back(PORE_PRESSURE_SIGN_FACTOR * dynamic_viscosity_inverse *
                                   relative_permeability_values[integration_point] *
                                   prod(permeability_matrix, GradPressureTerm));
        }

        return fluid_fluxes;
    }

    static Matrix CalculatePermeabilityMatrix(const Matrix& rGradNpT,
                                              double        DynamicViscosityInverse,
                                              const Matrix& rMaterialPermeabilityMatrix,
                                              double        RelativePermeability,
                                              double        IntegrationCoefficient);

    template <typename TMatrixType>
    static inline void CalculateCouplingMatrix(TMatrixType&  rCouplingMatrix,
                                               const Matrix& rB,
                                               const Vector& rVoigtVector,
                                               const Vector& rNp,
                                               double        BiotCoefficient,
                                               double        BishopCoefficient,
                                               double        IntegrationCoefficient)
    {
        const double multiplier =
            PORE_PRESSURE_SIGN_FACTOR * BiotCoefficient * BishopCoefficient * IntegrationCoefficient;

        Vector temp_vector(rB.size2());
        noalias(temp_vector) = prod(trans(rB), rVoigtVector);
        KRATOS_ERROR_IF(rCouplingMatrix.size1() != rB.size2())
            << " Inconsistent sizes: rCouplingMatrix.size1(): " << rCouplingMatrix.size1()
            << " rB.size2(): " << rB.size2() << std::endl;
        noalias(rCouplingMatrix) = multiplier * outer_prod(temp_vector, rNp);
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static BoundedMatrix<double, TNumNodes * TDim, TNumNodes> CalculateCouplingMatrix(const Matrix& rB,
                                                                                      const Vector& rVoigtVector,
                                                                                      const Vector& rNp,
                                                                                      double BiotCoefficient,
                                                                                      double BishopCoefficient,
                                                                                      double IntegrationCoefficient)
    {
        return PORE_PRESSURE_SIGN_FACTOR * BiotCoefficient * BishopCoefficient *
               IntegrationCoefficient * outer_prod(Vector(prod(trans(rB), rVoigtVector)), rNp);
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
        return DegreeOfSaturation * rProp[POROSITY] * rProp[DENSITY_WATER] +
               (1.0 - rProp[POROSITY]) * rProp[DENSITY_SOLID];
    }

    static std::vector<double> CalculateSoilDensities(const std::vector<double>& rDegreesSaturation,
                                                      const Properties&          rProp)
    {
        std::vector<double> result(rDegreesSaturation.size());
        std::transform(rDegreesSaturation.cbegin(), rDegreesSaturation.cend(), result.begin(),
                       [&rProp](const auto& degree_saturation) {
            return CalculateSoilDensity(degree_saturation, rProp);
        });
        return result;
    }

    template <typename VectorType>
    [[nodiscard]] static double CalculateFluidPressure(const VectorType& rN, const Vector& rPressureVector)
    {
        return inner_prod(rN, rPressureVector);
    }

    [[nodiscard]] static std::vector<double> CalculateFluidPressures(const Matrix& rNContainer,
                                                                     const Vector& rPressureVector);

    [[nodiscard]] static std::vector<double> CalculateInverseBiotModuli(const std::vector<double>& rBiotCoefficients,
                                                                        const std::vector<double>& rDegreesOfSaturation,
                                                                        const std::vector<double>& DerivativesOfSaturation,
                                                                        const Properties& rProperties)
    {
        std::vector<double> result;
        result.reserve(rBiotCoefficients.size());
        for (std::size_t i = 0; i < rBiotCoefficients.size(); ++i) {
            result.push_back(CalculateInverseBiotModulus(rBiotCoefficients[i], rDegreesOfSaturation[i],
                                                         DerivativesOfSaturation[i], rProperties));
        }
        return result;
    }

    [[nodiscard]] static double CalculateBulkModulus(const Matrix& rConstitutiveMatrix)
    {
        KRATOS_ERROR_IF(rConstitutiveMatrix.size1() == 0)
            << "Constitutive matrix is empty, aborting bulk modulus calculation.\n";
        const SizeType index_G = rConstitutiveMatrix.size1() - 1;
        return rConstitutiveMatrix(0, 0) - (4.0 / 3.0) * rConstitutiveMatrix(index_G, index_G);
    }

    [[nodiscard]] static std::vector<double> CalculateBiotCoefficients(const std::vector<Matrix>& rConstitutiveMatrices,
                                                                       const Properties& rProperties)
    {
        std::vector<double> result;
        result.reserve(rConstitutiveMatrices.size());
        std::transform(rConstitutiveMatrices.begin(), rConstitutiveMatrices.end(),
                       std::back_inserter(result), [&rProperties](const Matrix& rConstitutiveMatrix) {
            return CalculateBiotCoefficient(rConstitutiveMatrix, rProperties);
        });

        return result;
    }

    [[nodiscard]] static std::vector<double> CalculatePermeabilityUpdateFactors(const std::vector<Vector>& rStrainVectors,
                                                                                const Properties& rProperties)
    {
        auto result = std::vector<double>{};
        result.reserve(rStrainVectors.size());
        std::transform(rStrainVectors.cbegin(), rStrainVectors.cend(), std::back_inserter(result),
                       [&rProperties](const auto& rStrainVector) {
            return CalculatePermeabilityUpdateFactor(rStrainVector, rProperties);
        });
        return result;
    }

    [[nodiscard]] static double CalculateParticleDiameter(const Properties& rProperties);

private:
    [[nodiscard]] static double CalculateBiotCoefficient(const Matrix&     rConstitutiveMatrix,
                                                         const Properties& rProperties)
    {
        return rProperties.Has(BIOT_COEFFICIENT)
                   ? rProperties[BIOT_COEFFICIENT]
                   : 1.0 - CalculateBulkModulus(rConstitutiveMatrix) / rProperties[BULK_MODULUS_SOLID];
    }

    [[nodiscard]] static double CalculateInverseBiotModulus(double BiotCoefficient,
                                                            double DegreeOfSaturation,
                                                            double DerivativeOfSaturation,
                                                            const Properties& rProperties)
    {
        const auto bulk_modulus_fluid = rProperties[IGNORE_UNDRAINED] ? TINY : rProperties[BULK_MODULUS_FLUID];
        double result = (BiotCoefficient - rProperties[POROSITY]) / rProperties[BULK_MODULUS_SOLID] +
                        rProperties[POROSITY] / bulk_modulus_fluid;

        result *= DegreeOfSaturation;
        result -= DerivativeOfSaturation * rProperties[POROSITY];

        return result;
    }

    [[nodiscard]] static double CalculatePermeabilityUpdateFactor(const Vector&     rStrainVector,
                                                                  const Properties& rProperties)
    {
        if (rProperties[PERMEABILITY_CHANGE_INVERSE_FACTOR] > 0.0) {
            const double InverseCK = rProperties[PERMEABILITY_CHANGE_INVERSE_FACTOR];
            const double epsV      = StressStrainUtilities::CalculateTrace(rStrainVector);
            const double ePrevious = rProperties[POROSITY] / (1.0 - rProperties[POROSITY]);
            const double eCurrent  = (1.0 + ePrevious) * std::exp(epsV) - 1.0;
            const double permLog10 = (eCurrent - ePrevious) * InverseCK;
            return std::pow(10.0, permLog10);
        }

        return 1.0;
    }
}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
