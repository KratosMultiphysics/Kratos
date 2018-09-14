//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
#include "vms.h"

namespace Kratos
{
///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see VMS::EquationIdVector
 */
template <>
void VMS<2>::EquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(3),LocalSize(9);
    unsigned int LocalIndex = 0;

    unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X,vpos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y,vpos+1).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE,ppos).EquationId();
    }
}

/**
 * @see VMS::EquationIdVector
 */
template <>
void VMS<3>::EquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(4),LocalSize(16);
    unsigned int LocalIndex = 0;
    unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X,vpos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y,vpos+1).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z,vpos+2).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE,ppos).EquationId();
    }
}

/**
 * @see VMS::GetDofList
 */
template <>
void VMS<2>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(3),LocalSize(9);
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see VMS::GetDofList
 */
template <>
void VMS<3>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(4),LocalSize(16);
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see VMS::GetFirstDerivativesVector
 */
template <>
void VMS<2>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    const unsigned int NumNodes(3),LocalSize(9);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        Values[LocalIndex++] = rVelocity[0];
        Values[LocalIndex++] = rVelocity[1];
        Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

/**
 * @see VMS::GetFirstDerivativesVector
 */
template <>
void VMS<3>::GetFirstDerivativesVector(Vector& Values, int Step)
{
    const unsigned int NumNodes(4),LocalSize(16);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        Values[LocalIndex++] = rVelocity[0];
        Values[LocalIndex++] = rVelocity[1];
        Values[LocalIndex++] = rVelocity[2];
        Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

/**
 * @see VMS::GetSecondDerivativesVector
 */
template <>
void VMS<2>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    const unsigned int NumNodes(3),LocalSize(9);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        array_1d<double,3>& rAcceleration = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        Values[LocalIndex++] = rAcceleration[0];
        Values[LocalIndex++] = rAcceleration[1];
        Values[LocalIndex++] = 0.0; // Pressure Dof
    }
}

/**
 * @see VMS::GetSecondDerivativesVector
 */
template <>
void VMS<3>::GetSecondDerivativesVector(Vector& Values, int Step)
{
    const unsigned int NumNodes(4),LocalSize(16);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        array_1d<double,3>& rAcceleration = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        Values[LocalIndex++] = rAcceleration[0];
        Values[LocalIndex++] = rAcceleration[1];
        Values[LocalIndex++] = rAcceleration[2];
        Values[LocalIndex++] = 0.0; // Pressure Dof
    }
}

/**
 * The size of the 2D element is estimated as the diameter of a circle of the same area.
 * Area = Pi * (h/2)^2
 * @see VMS::ElementSize
 */
template <>
double VMS<2,3>::ElementSize(const double Area)
{
    return 1.128379167 * sqrt(Area); //Diameter of circumference of given Area
}

/**
 * The size of the 3D element is estimated as the diameter of the sphere
 * circumscribed to a regular tetrahedron with the same volume.
 * @see VMS::ElementSize
 */
template <>
double VMS<3,4>::ElementSize(const double Volume)
{
    return 0.60046878 * pow(Volume,0.333333333333333333333);
}


/**
 * @see VMS::GetValueOnIntegrationPoints
 */
template <>
void VMS<2>::GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
        std::vector<array_1d<double,3> >& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int Dim(2),NumNodes(3);
    if (rVariable == VORTICITY)
    {
        // Set output vector (for a single integration point)
        rOutput.resize(1);
        array_1d<double, 3 > & rVorticity = rOutput[0];
        rVorticity[0] = 0.0;
        rVorticity[1] = 0.0;
        rVorticity[2] = 0.0;

        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            const array_1d<double, 3 > & rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
            rVorticity[2] += DN_DX(iNode,0)*rVelocity[1] - DN_DX(iNode,1)*rVelocity[0];
        }
    }
    else if (rVariable == SUBSCALE_VELOCITY)
    {
        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        array_1d<double,3> AdvVel;
        this->GetAdvectiveVel(AdvVel,N);

        double Density;
        this->EvaluateInPoint(Density,DENSITY,N);

        double ElemSize = this->ElementSize(Area);
        double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

        // stabilization parameters
        double TauOne, TauTwo;
        this->CalculateTau(TauOne,TauTwo,AdvVel,ElemSize,Density,Viscosity,rCurrentProcessInfo);

        // Set output vector (for a single integration point)
        rOutput.resize(1);
        array_1d<double,3> MomError = ZeroVector(3);
        if (rCurrentProcessInfo[OSS_SWITCH]==1)
        {
            this->OSSMomResidual(AdvVel,Density,MomError,N,DN_DX,1.0);
        }
        else
        {
            this->ASGSMomResidual(AdvVel,Density,MomError,N,DN_DX,1.0);
        }
        MomError *= TauOne;
        array_1d<double,3>& rSubscale = rOutput[0];
        rSubscale[0] = MomError[0];
        rSubscale[1] = MomError[1];
        rSubscale[2] = 0.0;
    }
    else // Default behaviour (returns elemental data)
    {
        rOutput.resize(1);
        /*
         The cast is done to avoid modification of the element's data. Data modification
         would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         with associated value of 0.0). This is catastrophic if the variable referenced
         goes out of scope.
         */
        const VMS<Dim,NumNodes>* const_this = static_cast< const VMS<Dim,NumNodes>* >(this);
        rOutput[0] = const_this->GetValue(rVariable);
    }
}

/**
 * @see VMS::GetValueOnIntegrationPoints
 */
template <>
void VMS<3>::GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
        std::vector<array_1d<double,3> >& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int Dim(3),NumNodes(4);
    if (rVariable == VORTICITY)
    {
        // Set output vector (for a single integration point)
        rOutput.resize(1);
        array_1d<double, 3 > & rVorticity = rOutput[0];
        rVorticity[0] = 0.0;
        rVorticity[1] = 0.0;
        rVorticity[2] = 0.0;

        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            const array_1d<double, 3 > & rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
            rVorticity[0] += DN_DX(iNode,1)*rVelocity[2] - DN_DX(iNode,2)*rVelocity[1];
            rVorticity[1] += DN_DX(iNode,2)*rVelocity[0] - DN_DX(iNode,0)*rVelocity[2];
            rVorticity[2] += DN_DX(iNode,0)*rVelocity[1] - DN_DX(iNode,1)*rVelocity[0];
        }
    }
    else if(rVariable == SUBSCALE_VELOCITY)
    {
        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        array_1d<double,3> AdvVel;
        this->GetAdvectiveVel(AdvVel,N);

        double Density;
        this->EvaluateInPoint(Density,DENSITY,N);

        double ElemSize = this->ElementSize(Area);
        double Viscosity = this->EffectiveViscosity(Density,N,DN_DX,ElemSize,rCurrentProcessInfo);

        // stabilization parameters
        double TauOne, TauTwo;
        this->CalculateTau(TauOne,TauTwo,AdvVel,ElemSize,Density,Viscosity,rCurrentProcessInfo);

        // Set output vector (for a single integration point)
        rOutput.resize(1);
        array_1d<double,3> MomError = ZeroVector(3);
        if (rCurrentProcessInfo[OSS_SWITCH]==1)
        {
            this->OSSMomResidual(AdvVel,Density,MomError,N,DN_DX,1.0);
        }
        else
        {
            this->ASGSMomResidual(AdvVel,Density,MomError,N,DN_DX,1.0);
        }
        MomError *= TauOne;
        array_1d<double,3>& rSubscale = rOutput[0];
        rSubscale[0] = MomError[0];
        rSubscale[1] = MomError[1];
        rSubscale[2] = MomError[2];
    }
    else // Default behaviour (returns elemental data)
    {
        rOutput.resize(1);
        /*
         The cast is done to avoid modification of the element's data. Data modification
         would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         with associated value of 0.0). This is catastrophic if the variable referenced
         goes out of scope.
         */
        const VMS<Dim,NumNodes>* const_this = static_cast< const VMS<Dim,NumNodes>* >(this);
        rOutput[0] = const_this->GetValue(rVariable);
    }
}


template <>
double VMS<2,3>::EquivalentStrainRate(const BoundedMatrix<double,3,2> &rDN_DX) const
{
    const GeometryType& rGeom = this->GetGeometry();

    // Calculate Symetric gradient (Voigt notation)
    array_1d<double,3> S = ZeroVector(3);
    for (unsigned int n = 0; n < 3; ++n)
    {
        const array_1d<double,3>& rVel = rGeom[n].FastGetSolutionStepValue(VELOCITY);
        S[0] += rDN_DX(n,0)*rVel[0];
        S[1] += rDN_DX(n,1)*rVel[1];
        S[2] += rDN_DX(n,0)*rVel[1]+rDN_DX(n,1)*rVel[0];
    }

    // Norm of symetric gradient (cross terms don't get the 2)
    return std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + S[2]*S[2]);
}


template <>
double VMS<3,4>::EquivalentStrainRate(const BoundedMatrix<double,4,3> &rDN_DX) const
{
    const GeometryType& rGeom = this->GetGeometry();

    // Calculate Symetric gradient (Voigt notation)
    array_1d<double,6> S = ZeroVector(6);
    for (unsigned int n = 0; n < 4; ++n)
    {
        const array_1d<double,3>& rVel = rGeom[n].FastGetSolutionStepValue(VELOCITY);
        S[0] += rDN_DX(n,0)*rVel[0];
        S[1] += rDN_DX(n,1)*rVel[1];
        S[2] += rDN_DX(n,2)*rVel[2];
        S[3] += rDN_DX(n,2)*rVel[1]+rDN_DX(n,1)*rVel[2];
        S[4] += rDN_DX(n,2)*rVel[0]+rDN_DX(n,0)*rVel[2];
        S[5] += rDN_DX(n,1)*rVel[0]+rDN_DX(n,0)*rVel[1];
    }

    // Norm of symetric gradient (cross terms don't get the 2)
    return std::sqrt(2.*S[0]*S[0] + 2.*S[1]*S[1] + 2.*S[2]*S[2] + S[3]*S[3] + S[4]*S[4] + S[5]*S[5]);
}


/**
 * See VMS::CalculateB
 */
template <>
void VMS<2,3>::CalculateB( BoundedMatrix<double, 3, 6 >& rB,
                           const BoundedMatrix<double, 3, 2 >& rShapeDeriv)
{
    for (unsigned int i = 0; i < 3; i++)
    {
        unsigned int index = 2 * i;

        rB(0, index) = rShapeDeriv(i, 0);
        rB(0, index + 1) = 0.0;
        rB(1, index) = 0.0;
        rB(1, index + 1) = rShapeDeriv(i, 1);
        rB(2, index) = rShapeDeriv(i, 1);
        rB(2, index + 1) = rShapeDeriv(i, 0);
    }
}

/**
 * See VMS::CalculateB
 */
template <>
void VMS<3,4>::CalculateB( BoundedMatrix<double, 6, 12 >& rB,
                           const BoundedMatrix<double, 4, 3 >& rShapeDeriv)
{
    const unsigned int Dim = 3;
    const unsigned int NumNodes = 4;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        unsigned int index = Dim*i;

        rB(0, index) = rShapeDeriv(i, 0);
        rB(0, index + 1) = 0.0;
        rB(0, index + 2) = 0.0;
        rB(1, index) = 0.0;
        rB(1, index + 1) = rShapeDeriv(i, 1);
        rB(1, index + 2) =0.0;
        rB(2, index) = 0.0;
        rB(2, index + 1) = 0.0;
        rB(2, index + 2) = rShapeDeriv(i, 2);
        rB(3, index) = rShapeDeriv(i, 1);
        rB(3, index + 1) = rShapeDeriv(i, 0);
        rB(3, index + 2) = 0.0;
        rB(4, index) = 0.0;
        rB(4, index + 1) = rShapeDeriv(i, 2);
        rB(4, index + 2) = rShapeDeriv(i, 1);
        rB(5, index) = rShapeDeriv(i, 2);
        rB(5, index + 1) = 0.0;
        rB(5, index + 2) = rShapeDeriv(i, 0);
    }
}

/**
 * See VMS::CalculateC
 */
template <>
void VMS<2,3>::CalculateC(BoundedMatrix<double, 3, 3 > & rC,
                          const double Viscosity)
{
    rC(0, 0) =  Viscosity*(1.3333333333333333333333333333333);
    rC(0, 1) = -Viscosity*(0.666666666666666666666666666667);
    rC(0, 2) = 0.0;
    rC(1, 0) = -Viscosity*(0.666666666666666666666666666667);
    rC(1, 1) =  Viscosity*(1.3333333333333333333333333333);
    rC(1, 2) = 0.0;
    rC(2, 0) = 0.0;
    rC(2, 1) = 0.0;
    rC(2, 2) = Viscosity;
}

/**
 * See VMS::CalculateC
 */
template <>
void VMS<3,4>::CalculateC(BoundedMatrix<double, 6,6 > & rC,
                          const double Viscosity)
{
    noalias(rC) = ZeroMatrix(6,6);
    // First row
    rC(0, 0) =  Viscosity*(1.3333333333333333333333333333333);
    rC(0, 1) = -Viscosity*(0.666666666666666666666666666667);
    rC(0, 2) = -Viscosity*(0.666666666666666666666666666667);

    // Second row
    rC(1, 0) = -Viscosity*(0.666666666666666666666666666667);
    rC(1, 1) =  Viscosity*(1.3333333333333333333333333333);
    rC(1, 2) = -Viscosity*(0.666666666666666666666666666667);

    // Third row
    rC(2, 0) = -Viscosity*(0.666666666666666666666666666667);
    rC(2, 1) = -Viscosity*(0.666666666666666666666666666667);
    rC(2, 2) =  Viscosity*(1.3333333333333333333333333333);

    // Fourth row
    rC(3, 3) = Viscosity;

    // Fifth row
    rC(4, 4) = Viscosity;

    // Sixth row
    rC(5, 5) = Viscosity;

}

/**
 * @see VMS::AddViscousTerm
 */
template <>
void VMS<2,3>::AddViscousTerm(MatrixType& rDampingMatrix,
                              const BoundedMatrix<double,3,2>& rShapeDeriv,
                              const double Weight)
{
    const unsigned int NumNodes = 3;

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0),FirstCol(0);

    for (unsigned int j = 0; j < NumNodes; ++j)
    {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            // First Row
            rDampingMatrix(FirstRow,FirstCol) += Weight * ( FourThirds * rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) );
            rDampingMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );

            // Second Row
            rDampingMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
            rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( FourThirds * rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,0) * rShapeDeriv(j,0) );

            // Update Counter
            FirstRow += 3;
        }
        FirstRow = 0;
        FirstCol += 3;
    }
}

/**
 * @see VMS::AddViscousTerm
 */
template <>
void VMS<3,4>::AddViscousTerm(MatrixType& rDampingMatrix,
                              const BoundedMatrix<double,4,3>& rShapeDeriv,
                              const double Weight)
{
    const unsigned int NumNodes = 4;

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0),FirstCol(0);

    for (unsigned int j = 0; j < NumNodes; ++j)
    {
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            // (dN_i/dx_k dN_j/dx_k)
            const double Diag =  rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,2) * rShapeDeriv(j,2);

            // First Row
            rDampingMatrix(FirstRow,FirstCol) += Weight * ( OneThird * rShapeDeriv(i,0) * rShapeDeriv(j,0) + Diag );
            rDampingMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );
            rDampingMatrix(FirstRow,FirstCol+2) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,2) + rShapeDeriv(i,2) * rShapeDeriv(j,0) );

            // Second Row
            rDampingMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
            rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( OneThird * rShapeDeriv(i,1) * rShapeDeriv(j,1) + Diag );
            rDampingMatrix(FirstRow+1,FirstCol+2) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,2) + rShapeDeriv(i,2) * rShapeDeriv(j,1) );

            // Third Row
            rDampingMatrix(FirstRow+2,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,2) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,2) );
            rDampingMatrix(FirstRow+2,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,2) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,2) );
            rDampingMatrix(FirstRow+2,FirstCol+2) += Weight * ( OneThird * rShapeDeriv(i,2) * rShapeDeriv(j,2) + Diag );

            // Update Counter
            FirstRow += 4;
        }
        FirstRow = 0;
        FirstCol += 4;
    }
}

template<>
double VMS<2,3>::ConsistentMassCoef(const double Area)
{
    const double Coef = 1.0/12.0;
    return Area * Coef;
}

template<>
double VMS<3,4>::ConsistentMassCoef(const double Volume)
{
    const double Coef = 1.0/20.0;
    return Volume * Coef;
}

///@} // Specialized implementations
}
