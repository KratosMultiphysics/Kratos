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

    static Matrix CalculateMassMatrixDiffOrder(const Geometry<Node>&          rGeom,
                                               const Geometry<Node>::Pointer& rpPressureGeometry,
                                               const GeometryData::IntegrationMethod IntegrationMethod,
                                               const StressStatePolicy& rStressStatePolicy,
                                               const std::vector<RetentionLaw::Pointer>& rRetentionLawVector,
                                               const Properties&  rProp,
                                               const ProcessInfo& rCurrentProcessInfo)
    {
        const SizeType dimension          = rGeom.WorkingSpaceDimension();
        const SizeType number_U_nodes     = rGeom.PointsNumber();
        const SizeType block_element_size = number_U_nodes * dimension;
        const SizeType number_P_nodes     = rpPressureGeometry->PointsNumber();
        const Geometry<Node>::IntegrationPointsArrayType& integration_points =
            rGeom.IntegrationPoints(IntegrationMethod);

        const Matrix& Nu_container = rGeom.ShapeFunctionsValues(IntegrationMethod);
        const Matrix  Np_container = rpPressureGeometry->ShapeFunctionsValues(IntegrationMethod);
        const Vector  pressure_vector =
            GeoTransportEquationUtilities::GetSolutionVector(number_P_nodes, rGeom, WATER_PRESSURE);

        // create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(rProp, rCurrentProcessInfo);

        Matrix         Nu                 = ZeroMatrix(dimension, number_U_nodes * dimension);
        Matrix         aux_density_matrix = ZeroMatrix(dimension, number_U_nodes * dimension);
        Matrix         density_matrix     = ZeroMatrix(dimension, dimension);
        const SizeType element_size       = block_element_size + number_P_nodes;
        Matrix         mass_matrix        = ZeroMatrix(element_size, element_size);
        
        for (unsigned int g_point = 0; g_point < integration_points.size(); ++g_point) {
            const double degree_of_saturation = CalculateDegreeOfSaturation(
                g_point, Np_container, pressure_vector, RetentionParameters, rRetentionLawVector);

            const double density = CalculateSoilDensity(degree_of_saturation, rProp);

            GeoElementUtilities::AssembleDensityMatrix(density_matrix, density);
            GeoElementUtilities::CalculateNuMatrix(dimension, number_U_nodes, Nu, Nu_container, g_point);
            noalias(aux_density_matrix) = prod(density_matrix, Nu);

            // Adding contribution to Mass matrix
            const double integration_coefficient_initial_configuration =
                CalculateIntegrationCoefficientInitialConfiguration(
                    g_point, rGeom, integration_points, IntegrationMethod, rStressStatePolicy);

            GeoElementUtilities::AssembleUUBlockMatrix(
                mass_matrix, prod(trans(Nu), aux_density_matrix) * integration_coefficient_initial_configuration);
        }
        return mass_matrix;
    }

    static double CalculateSoilDensity(double DegreeOfSaturation, const Properties& rProp)
    {
        return (DegreeOfSaturation * rProp[POROSITY] * rProp[DENSITY_WATER]) +
               (1.0 - rProp[POROSITY]) * rProp[DENSITY_SOLID];
    }

    template <unsigned int TDim, unsigned int TNumNodes>
    static Matrix CalculateMassMatrix(unsigned int                              N_DOF,
                                      const Geometry<Node>&                     rGeom,
                                      const GeometryData::IntegrationMethod     IntegrationMethod,
                                      const StressStatePolicy&                  rStressStatePolicy,
                                      const std::vector<RetentionLaw::Pointer>& rRetentionLawVector,
                                      const Properties&                         rProp,
                                      const ProcessInfo&                        rCurrentProcessInfo)
    {
        const Geometry<Node>::IntegrationPointsArrayType& integration_points =
            rGeom.IntegrationPoints(IntegrationMethod);
        const unsigned int number_G_points = integration_points.size();

        const Matrix& N_container = rGeom.ShapeFunctionsValues(IntegrationMethod);
        const Vector  pressure_vector =
            GeoTransportEquationUtilities::GetSolutionVector(TNumNodes, rGeom, WATER_PRESSURE);

        // Create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(rProp, rCurrentProcessInfo);

        BoundedMatrix<double, TDim, TNumNodes * TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
        BoundedMatrix<double, TDim, TNumNodes * TDim> aux_density_matrix = ZeroMatrix(TDim, TNumNodes * TDim);
        BoundedMatrix<double, TDim, TDim> density_matrix = ZeroMatrix(TDim, TDim);
        Matrix                            mass_matrix    = ZeroMatrix(N_DOF, N_DOF);

        for (unsigned int g_point = 0; g_point < number_G_points; ++g_point) {
            const double degree_of_saturation = CalculateDegreeOfSaturation(
                g_point, N_container, pressure_vector, RetentionParameters, rRetentionLawVector);

            const double density = CalculateSoilDensity(degree_of_saturation, rProp);

            GeoElementUtilities::AssembleDensityMatrix(density_matrix, density);
            GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nu, N_container, g_point);
            noalias(aux_density_matrix) = prod(density_matrix, Nu);

            // Adding contribution to Mass matrix
            const double integration_coefficient_initial_configuration =
                CalculateIntegrationCoefficientInitialConfiguration(
                    g_point, rGeom, integration_points, IntegrationMethod, rStressStatePolicy);
            GeoElementUtilities::AssembleUUBlockMatrix(
                mass_matrix, prod(trans(Nu), aux_density_matrix) * integration_coefficient_initial_configuration);
        }
        return mass_matrix;
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

    static double CalculateDegreeOfSaturation(unsigned int              GPoint,
                                              const Matrix&             rNContainer,
                                              const Vector&             rPressureVector,
                                              RetentionLaw::Parameters& rRetentionParameters,
                                              const std::vector<RetentionLaw::Pointer>& rRetentionLawVector)
    {
        const double fluid_pressure = GeoTransportEquationUtilities::CalculateFluidPressure(
            row(rNContainer, GPoint), rPressureVector);
        rRetentionParameters.SetFluidPressure(fluid_pressure);

        return rRetentionLawVector[GPoint]->CalculateSaturation(rRetentionParameters);
    }

    static double CalculateFluidPressure(const Vector& rNp, const Vector& rPressureVector)
    {
        return inner_prod(rNp, rPressureVector);
    }

    static Vector GetSolutionVector(SizeType NumberPNodes, const Geometry<Node>& rGeom, const Variable<double>& rSolutionVariable)
    {
        Vector solution_vector(NumberPNodes);
        std::transform(rGeom.begin(), rGeom.end(), solution_vector.begin(),
                       [&rSolutionVariable](const auto& node) {
            return node.FastGetSolutionStepValue(rSolutionVariable);
        });
        return solution_vector;
    }
}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
