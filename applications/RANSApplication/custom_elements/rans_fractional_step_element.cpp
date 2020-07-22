//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//  Extended by :    Suneth Warnakulasuriya
//

// Include base h
#include "rans_fractional_step_element.h"

namespace Kratos
{
template <unsigned int TDim>
void RansFractionalStepElement<TDim>::CalculateLocalFractionalVelocitySystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = TDim * NumNodes;

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize);

    rRightHandSideVector = ZeroVector(LocalSize);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    MatrixType MassMatrix = ZeroMatrix(LocalSize, LocalSize);

    const double eta = rCurrentProcessInfo[PRESSURE_COEFFICIENT];

    // Stabilization parameters
    double ElemSize = this->ElementSize();
    double TauOne;
    double TauTwo;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& N = row(NContainer, g);
        const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

        // Evaluate required variables at the integration point
        double Density;
        double MassProjection;
        array_1d<double, 3> BodyForce = ZeroVector(3);
        array_1d<double, 3> MomentumProjection = ZeroVector(3);

        this->EvaluateInPoint(Density, DENSITY, N);
        this->EvaluateInPoint(MassProjection, DIVPROJ, N);
        this->EvaluateInPoint(BodyForce, BODY_FORCE, N);
        //        this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
        this->EvaluateInPoint(MomentumProjection, CONV_PROJ, N);

        // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
        double OldPressure;
        this->EvaluateInPoint(OldPressure, PRESSURE, N, 0);

        // For ALE: convective velocity
        array_1d<double, 3> ConvVel = ZeroVector(3);
        this->EvaluateConvVelocity(ConvVel, N);

        double Viscosity = this->EffectiveViscosity(
            Density, N, rDN_DX, ElemSize, rCurrentProcessInfo);
        this->CalculateTau(TauOne, TauTwo, ElemSize, ConvVel, Density,
                           Viscosity, rCurrentProcessInfo);

        // Evaluate convection operator Velocity * Grad(N)
        Vector UGradN(NumNodes);
        this->ConvectionOperator(UGradN, ConvVel, rDN_DX);

        // Add integration point contribution to the local mass matrix
        this->AddMomentumMassTerm(MassMatrix, N, GaussWeight * Density);

        // Add convection, stabilization and RHS contributions to the local system equation
        this->AddMomentumSystemTerms(rLeftHandSideMatrix, rRightHandSideVector,
                                     Density, UGradN, BodyForce, OldPressure * eta,
                                     TauOne, TauTwo, MomentumProjection,
                                     MassProjection, N, rDN_DX, GaussWeight);

        // Add viscous term
        const double ViscousCoeff = Viscosity * GaussWeight;
        this->AddViscousTerm(rLeftHandSideMatrix, rDN_DX, ViscousCoeff);
    }

    // Add residual of previous iteration to RHS
    VectorType LastValues = ZeroVector(LocalSize);
    this->GetVelocityValues(LastValues, 0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, LastValues);

    // Add dynamic term
    const Vector& rBDFCoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
    noalias(rLeftHandSideMatrix) += rBDFCoeffs[0] * MassMatrix;

    VectorType TimeTerm = rBDFCoeffs[0] * LastValues;
    for (SizeType i = 1; i < rBDFCoeffs.size(); i++)
    {
        this->GetVelocityValues(LastValues, i);
        noalias(TimeTerm) += rBDFCoeffs[i] * LastValues;
    }

    noalias(rRightHandSideVector) -= prod(MassMatrix, TimeTerm);
}

template <unsigned int TDim>
void RansFractionalStepElement<TDim>::CalculateLocalPressureSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    // Stabilization parameters
    double ElemSize = this->ElementSize();
    double TauOne;
    double TauTwo;

    const double eta = rCurrentProcessInfo[PRESSURE_COEFFICIENT];

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& N = row(NContainer, g);
        const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

        // Evaluate required variables at the integration point
        double Density;
        //        array_1d<double,3> Velocity = ZeroVector(3);
        //        array_1d<double,3> MeshVelocity = ZeroVector(3);
        array_1d<double, 3> BodyForce = ZeroVector(3);
        array_1d<double, 3> MomentumProjection = ZeroVector(3);

        this->EvaluateInPoint(Density, DENSITY, N);
        //        this->EvaluateInPoint(Velocity,VELOCITY,N);
        //        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
        this->EvaluateInPoint(BodyForce, BODY_FORCE, N);
        //        this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
        this->EvaluateInPoint(MomentumProjection, PRESS_PROJ, N);

        //        // Evaluate the pressure and pressure gradient at this point
        //        (for the G * P_n term) double OldPressure;
        //        this->EvaluateInPoint(OldPressure,PRESSURE,N,0);

        array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
        this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

        //        // For ALE: convective velocity
        //        array_1d<double,3> ConvVel = Velocity - MeshVelocity;

        // Stabilization parameters
        array_1d<double, 3> ConvVel = ZeroVector(3);
        this->EvaluateConvVelocity(ConvVel, N);
        double Viscosity = this->EffectiveViscosity(
            Density, N, rDN_DX, ElemSize, rCurrentProcessInfo);
        this->CalculateTau(TauOne, TauTwo, ElemSize, ConvVel, Density,
                           Viscosity, rCurrentProcessInfo);

        //        // Evaluate convection operator Velocity * Grad(N)
        //        Vector UGradN(NumNodes);
        //        this->EvaluateConvection(UGradN,ConvVel,mDN_DX);

        double DivU;
        this->EvaluateDivergenceInPoint(DivU, VELOCITY, rDN_DX);

        // constant coefficient multiplying the pressure Laplacian (See Codina, Badia 2006 paper for details in case of a BDF2 time scheme)
        const double LaplacianCoeff =
            1.0 / (Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]);

        // Add convection, stabilization and RHS contributions to the local system equation
        for (SizeType i = 0; i < NumNodes; ++i)
        {
            // LHS contribution
            for (SizeType j = 0; j < NumNodes; ++j)
            {
                double Lij = 0.0;
                for (SizeType d = 0; d < TDim; ++d)
                    Lij += rDN_DX(i, d) * rDN_DX(j, d);
                Lij *= (LaplacianCoeff + TauOne);

                rLeftHandSideMatrix(i, j) += GaussWeight * Lij;
            }

            // RHS contribution

            // Velocity divergence
            double RHSi = -N[i] * DivU;

            for (SizeType d = 0; d < TDim; ++d)
            {
                //                double Conv = UGradN[0] * rGeom[0].FastGetSolutionStepValue(VELOCITY)[d];
                //                for (SizeType j = 1; j < NumNodes; ++j) Conv += UGradN[i] * rGeom[i].FastGetSolutionStepValue(VELOCITY)[d];
                // Momentum stabilization
                RHSi += rDN_DX(i, d) * TauOne *
                        (Density * (BodyForce[d] /* - Conv*/) -
                         OldPressureGradient[d] - MomentumProjection[d]);
                RHSi += (eta - 1) * LaplacianCoeff * rDN_DX(i, d) * OldPressureGradient[d];
            }

            rRightHandSideVector[i] += GaussWeight * RHSi;
        }
    }
}

/*
 * Template class definition (this should allow us to compile the desired template instantiations)
 */

template class RansFractionalStepElement<2>;
template class RansFractionalStepElement<3>;

} // namespace Kratos
