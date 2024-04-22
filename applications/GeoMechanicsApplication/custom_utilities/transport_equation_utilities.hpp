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

    static Matrix CalculateMassMatrix(const Geometry<Node>&                     rGeom,
                                      const Geometry<Node>::Pointer&            rpPressureGeometry,
                                      const GeometryData::IntegrationMethod     IntegrationMethod,
                                      std::unique_ptr<StressStatePolicy> const& mpStressStatePolicy,
                                      std::vector<RetentionLaw::Pointer>&       rRetentionLawVector,
                                      const Properties&                         rProp,
                                      const ProcessInfo&                        rCurrentProcessInfo)
    {
        const SizeType dimension          = rGeom.WorkingSpaceDimension();
        const SizeType number_U_nodes     = rGeom.PointsNumber();
        const SizeType block_element_size = number_U_nodes * dimension;
        const SizeType number_P_nodes     = rpPressureGeometry->PointsNumber();
        const Geometry<Node>::IntegrationPointsArrayType& integration_points =
            rGeom.IntegrationPoints(IntegrationMethod);

        Matrix nu_container = rGeom.ShapeFunctionsValues(IntegrationMethod);
        Matrix np_container = rpPressureGeometry->ShapeFunctionsValues(IntegrationMethod);
        Vector pressure_vector =
            GeoTransportEquationUtilities::GetSolutionVector(number_P_nodes, rGeom, WATER_PRESSURE);

        // create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(rProp, rCurrentProcessInfo);

        Matrix mass_matrix_contribution = ZeroMatrix(block_element_size, block_element_size);

        // Defining shape functions and the determinant of the jacobian at all integration points
        // Loop over integration points
        Matrix Nu                 = ZeroMatrix(dimension, number_U_nodes * dimension);
        Matrix aux_density_matrix = ZeroMatrix(dimension, number_U_nodes * dimension);
        Matrix density_matrix     = ZeroMatrix(dimension, dimension);

        for (unsigned int g_point = 0; g_point < integration_points.size(); ++g_point) {
            Matrix J0;
            Matrix inv_J0;
            GeometryUtils::JacobianOnInitialConfiguration(rGeom, integration_points[g_point], J0);
            const Matrix& dN_De = rGeom.ShapeFunctionsLocalGradients(IntegrationMethod)[g_point];
            double        det_J_initial_configuration;
            MathUtils<double>::InvertMatrix(J0, inv_J0, det_J_initial_configuration);
            Matrix dNu_dX_initial_configuration;
            GeometryUtils::ShapeFunctionsGradients(dN_De, inv_J0, dNu_dX_initial_configuration);

            // calculating weighting coefficient for integration
            const double integration_coefficient_initial_configuration =
                mpStressStatePolicy->CalculateIntegrationCoefficient(
                    integration_points[g_point], det_J_initial_configuration, rGeom);

            const double fluid_pressure = GeoTransportEquationUtilities::CalculateFluidPressure(
                row(np_container, g_point), pressure_vector);
            RetentionParameters.SetFluidPressure(fluid_pressure);

            const double degree_of_saturation =
                rRetentionLawVector[g_point]->CalculateSaturation(RetentionParameters);

            const double density = CalculateSoilDensity(degree_of_saturation, rProp);

            // Setting the shape function matrix
            SizeType index = 0;
            for (SizeType i = 0; i < number_U_nodes; ++i) {
                for (SizeType i_dim = 0; i_dim < dimension; ++i_dim) {
                    Nu(i_dim, index) = nu_container(g_point, i);
                    index++;
                }
            }

            GeoElementUtilities::AssembleDensityMatrix(density_matrix, density);

            noalias(aux_density_matrix) = prod(density_matrix, Nu);

            // Adding contribution to Mass matrix
            noalias(mass_matrix_contribution) +=
                prod(trans(Nu), aux_density_matrix) * integration_coefficient_initial_configuration;
        }

        // Distribute mass block matrix into the elemental matrix
        const SizeType element_size = block_element_size + number_P_nodes;
        Matrix         mass_matrix  = ZeroMatrix(element_size, element_size);

        for (SizeType i = 0; i < number_U_nodes * dimension; ++i) {
            for (SizeType j = 0; j < number_U_nodes * dimension; ++j) {
                mass_matrix(i, j) += mass_matrix_contribution(i, j);
            }
        }
        return mass_matrix;
    }

    static double CalculateFluidPressure(const Vector& rNp, const Vector& rPressureVector)
    {
        return inner_prod(rNp, rPressureVector);
    }

    static Vector GetSolutionVector(SizeType NumberPNodes, const Geometry<Node>& rGeom, const Variable<double>& SolutionVariable)
    {
        Vector solution_vector(NumberPNodes);
        for (SizeType i = 0; i < NumberPNodes; ++i) {
            solution_vector[i] = rGeom[i].FastGetSolutionStepValue(SolutionVariable);
        }
        return solution_vector;
    }

    static double CalculateSoilDensity(double DegreeOfSaturation, const Properties& rProp)
    {
        return (DegreeOfSaturation * rProp[POROSITY] * rProp[DENSITY_WATER]) +
               (1.0 - rProp[POROSITY]) * rProp[DENSITY_SOLID];
    }

}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
