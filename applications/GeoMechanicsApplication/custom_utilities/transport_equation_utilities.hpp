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
                                      const SizeType                            NumPNodes,
                                      const GeometryData::IntegrationMethod     IntegrationMethod,
                                      const std::unique_ptr<StressStatePolicy>& mpStressStatePolicy,
                                      std::vector<RetentionLaw::Pointer>&       rRetentionLawVector,
                                      const Properties&                         rProp,
                                      const ProcessInfo&                        rCurrentProcessInfo)
    {
        Matrix                                            rMassMatrix;
        const SizeType                                    Dim       = rGeom.WorkingSpaceDimension();
        const SizeType                                    NumUNodes = rGeom.PointsNumber();
        const SizeType                                    BlockElementSize = NumUNodes * Dim;
        const Geometry<Node>::IntegrationPointsArrayType& IntegrationPoints =
            rGeom.IntegrationPoints(IntegrationMethod);

        // ElementVariables Variables;
        // this->InitializeElementVariables(Variables, rCurrentProcessInfo);
        //  from InitializeElementVariables
        Matrix NuContainer = rGeom.ShapeFunctionsValues(IntegrationMethod);
        Vector vector_Nu(NumUNodes);
        Vector Np;
        Vector PressureVector =
            GeoTransportEquationUtilities::GetSolutionVector(NumPNodes, rGeom, WATER_PRESSURE);

        // create general parameters of retention law
        RetentionLaw::Parameters RetentionParameters(rProp, rCurrentProcessInfo); // do I really need it?

        Matrix MassMatrixContribution = ZeroMatrix(BlockElementSize, BlockElementSize);

        // Defining shape functions and the determinant of the jacobian at all integration points

        // Loop over integration points
        Matrix Nu               = ZeroMatrix(Dim, NumUNodes * Dim);
        Matrix AuxDensityMatrix = ZeroMatrix(Dim, NumUNodes * Dim);
        Matrix DensityMatrix    = ZeroMatrix(Dim, Dim);

        for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
            double detJInitialConfiguration;
            Matrix DNu_DXInitialConfiguration;
            // compute element kinematics (Np, gradNpT, |J|, B)
            // this->CalculateKinematics(Variables, GPoint);
            {
                // Setting the vector of shape functions and the matrix of the shape functions global gradients
                noalias(vector_Nu) = row(NuContainer, GPoint); // it was Variables.Nu
                // noalias(Np) = row(NpContainer, GPoint);

                // noalias(DNu_DX) = DNu_DXContainer[GPoint];
                // noalias(DNp_DX) = DNp_DXContainer[GPoint];

                // Compute the deformation matrix B
                // this->CalculateBMatrix(B, DNu_DX, Nu);

                // detJ = rVariables.detJuContainer[GPoint];

                Matrix J0;
                Matrix InvJ0;

                // this->CalculateDerivativesOnInitialConfiguration(
                //     detJInitialConfiguration, J0, InvJ0, DNu_DXInitialConfiguration, GPoint);
                //  double& detJ, Matrix& J0, Matrix& InvJ0, Matrix& DNu_DX0, unsigned int GPoint) const
                {
                    GeometryUtils::JacobianOnInitialConfiguration(rGeom, IntegrationPoints[GPoint], J0);
                    const Matrix& DN_De = rGeom.ShapeFunctionsLocalGradients(IntegrationMethod)[GPoint];
                    MathUtils<double>::InvertMatrix(J0, InvJ0, detJInitialConfiguration);
                    GeometryUtils::ShapeFunctionsGradients(DN_De, InvJ0, DNu_DXInitialConfiguration);
                }
            }

            // calculating weighting coefficient for integration
            double IntegrationCoefficientInitialConfiguration =
                //    this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, detJInitialConfiguration);
                mpStressStatePolicy->CalculateIntegrationCoefficient(
                    IntegrationPoints[GPoint], detJInitialConfiguration, rGeom);

            double DegreeOfSaturation;
            double Density;
            // CalculateRetentionResponse(Variables, RetentionParameters, GPoint);
            //  void SmallStrainUPwDiffOrderElement::CalculateRetentionResponse(ElementVariables& rVariables,
            //                                                                  RetentionLaw::Parameters& rRetentionParameters,
            //                                                                  unsigned int GPoint)
            {
                double FluidPressure = GeoTransportEquationUtilities::CalculateFluidPressure(
                    Np, PressureVector); // it needs rVariables.Np, rVariables.PressureVector
                // SetRetentionParameters(rVariables, RetentionParameters);
                {
                    RetentionParameters.SetFluidPressure(FluidPressure);
                }

                DegreeOfSaturation = rRetentionLawVector[GPoint]->CalculateSaturation(RetentionParameters); // it needs FluidPressure
                // DerivativeOfSaturation = mRetentionLawVector[GPoint]->CalculateDerivativeOfSaturation(RetentionParameters);
                // RelativePermeability = mRetentionLawVector[GPoint]->CalculateRelativePermeability(RetentionParameters);
                // BishopCoefficient = mRetentionLawVector[GPoint]->CalculateBishopCoefficient(RetentionParameters);
            }

            // this->CalculateSoilDensity(Variables);
            Density = CalculateSoilDensity(DegreeOfSaturation, rProp);

            // Setting the shape function matrix
            SizeType Index = 0;
            for (SizeType i = 0; i < NumUNodes; ++i) {
                for (SizeType iDim = 0; iDim < Dim; ++iDim) {
                    Nu(iDim, Index++) = vector_Nu(i);
                }
            }

            GeoElementUtilities::AssembleDensityMatrix(DensityMatrix, Density);

            noalias(AuxDensityMatrix) = prod(DensityMatrix, Nu);

            // Adding contribution to Mass matrix
            noalias(MassMatrixContribution) +=
                prod(trans(Nu), AuxDensityMatrix) * IntegrationCoefficientInitialConfiguration;
        }

        // Distribute mass block matrix into the elemental matrix
        const SizeType ElementSize = BlockElementSize + NumPNodes;

        if (rMassMatrix.size1() != ElementSize || rMassMatrix.size2() != ElementSize)
            rMassMatrix.resize(ElementSize, ElementSize, false);
        noalias(rMassMatrix) = ZeroMatrix(ElementSize, ElementSize);

        for (SizeType i = 0; i < NumUNodes; ++i) {
            SizeType Index_i = i * Dim;

            for (SizeType j = 0; j < NumUNodes; ++j) {
                SizeType Index_j = j * Dim;
                for (SizeType idim = 0; idim < Dim; ++idim) {
                    for (SizeType jdim = 0; jdim < Dim; ++jdim) {
                        rMassMatrix(Index_i + idim, Index_j + jdim) +=
                            MassMatrixContribution(Index_i + idim, Index_j + jdim);
                    }
                }
            }
        }
        return rMassMatrix;
    }

    static double CalculateFluidPressure(const Vector& rNp, const Vector& rPressureVector)
    {
        return inner_prod(rNp, rPressureVector);
    }

    static Vector GetSolutionVector(SizeType NumPNodes, const Geometry<Node>& rGeom, const Variable<double>& Solution)
    {
        Vector SolutionVector(NumPNodes);
        for (SizeType i = 0; i < NumPNodes; ++i) {
            SolutionVector[i] = rGeom[i].FastGetSolutionStepValue(Solution);
        }
        return SolutionVector;
    }

    static double CalculateSoilDensity(double DegreeOfSaturation, const Properties& rProp)
    {
        return (DegreeOfSaturation * rProp[POROSITY] * rProp[DENSITY_WATER]) +
               (1.0 - rProp[POROSITY]) * rProp[DENSITY_SOLID];
    }

}; /* Class GeoTransportEquationUtilities*/
} /* namespace Kratos.*/
