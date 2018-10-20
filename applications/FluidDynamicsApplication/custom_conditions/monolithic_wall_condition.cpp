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

#include "monolithic_wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see MonolithicWallCondition::EquationIdVector
 */
template <>
void MonolithicWallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 6;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see MonolithicWallCondition::EquationIdVector
 */
template <>
void MonolithicWallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see MonolithicWallCondition::GetDofList
 */
template <>
void MonolithicWallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 6;

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
 * @see MonolithicWallCondition::GetDofList
 */
template <>
void MonolithicWallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;

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




template<unsigned int TDim, unsigned int TNumNodes>
void MonolithicWallCondition<TDim,TNumNodes>::CalculateLocalVelocityContribution(MatrixType &rDampMatrix,
                                                                                 VectorType &rRightHandSideVector,
                                                                                 ProcessInfo &rCurrentProcessInfo)
{
    // Initialize local contributions
    const SizeType LocalSize = (TDim + 1) * TNumNodes;

    if (rDampMatrix.size1() != LocalSize)
        rDampMatrix.resize(LocalSize,LocalSize);
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize);

    noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

//    // Pressure Neumann boundary contribution
//    const double NodeFactor = 1.0/double(TDim);
//    GeometryType& rGeom = this->GetGeometry();
//    double ExternalPressure = rGeom[0].FastGetSolutionStepValue(EXTERNAL_PRESSURE);
//    for (unsigned int d = 1; d < TDim; d++)
//        ExternalPressure += rGeom[d].FastGetSolutionStepValue(EXTERNAL_PRESSURE);
//    ExternalPressure *= NodeFactor;

//    array_1d<double,3> Normal = ZeroVector(3);
//    this->CalculateNormal(Normal);

//    unsigned int Row = 0;
//    for (unsigned int i = 0; i < TNumNodes; i++)
//    {
//        Row = i*(TDim+1);
//        for (unsigned int d = 0; d < TDim; d++)
//            rRightHandSideVector[Row+d] -= NodeFactor*Normal[d]*ExternalPressure;
//    }

    this->ApplyNeumannCondition(rDampMatrix,rRightHandSideVector);

    this->ApplyWallLaw(rDampMatrix,rRightHandSideVector,rCurrentProcessInfo);
}

template<unsigned int TDim, unsigned int TNumNodes>
void MonolithicWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> > &rVariable,
                                                                          std::vector<array_1d<double,3> > &rValues,
                                                                          const ProcessInfo &rCurrentProcessInfo)
{
    rValues.resize(1);
    if (rVariable == NORMAL)
    {
        this->CalculateNormal(rValues[0]);
    }
    else
    {
        /* The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const MonolithicWallCondition* const_this = static_cast< const MonolithicWallCondition* >(this);
        rValues[0] = const_this->GetValue(rVariable);
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void MonolithicWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                                          std::vector<double>& rValues,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    /*
     The cast is done to avoid modification of the element's data. Data modification
     would happen if rVariable is not stored now (would initialize a pointer to &rVariable
     with associated value of 0.0). This is catastrophic if the variable referenced
     goes out of scope.
     */
    const MonolithicWallCondition* const_this = static_cast< const MonolithicWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void MonolithicWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                                                          std::vector<array_1d<double, 6 > >& rValues,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const MonolithicWallCondition* const_this = static_cast< const MonolithicWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void MonolithicWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                          std::vector<Vector>& rValues,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const MonolithicWallCondition* const_this = static_cast< const MonolithicWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void MonolithicWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                          std::vector<Matrix>& rValues,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const MonolithicWallCondition* const_this = static_cast< const MonolithicWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


// protected funcions


template <>
void MonolithicWallCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;

}

template <>
void MonolithicWallCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}

template<unsigned int TDim, unsigned int TNumNodes>
void MonolithicWallCondition<TDim,TNumNodes>::ApplyNeumannCondition(MatrixType &rLocalMatrix, VectorType &rLocalVector)
{
    const MonolithicWallCondition<TDim,TNumNodes>& rConstThis = *this;
    if (rConstThis.GetValue(IS_STRUCTURE) == 0.0)
    {
        const unsigned int LocalSize = TDim+1;
        const GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();

        MatrixType NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        array_1d<double,3> Normal;
        this->CalculateNormal(Normal); //this already contains the area
        double A = std::sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
        Normal /= A;

        // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
        double J = (TDim == 2) ? 0.5*A : 2.0*A;

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            Vector N = row(NContainer,g);
            double Weight = J * IntegrationPoints[g].Weight();

            // Neumann boundary condition
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                //unsigned int row = i*LocalSize;
                const NodeType& rConstNode = this->GetGeometry()[i];
                const double pext = rConstNode.FastGetSolutionStepValue(EXTERNAL_PRESSURE);

                for (unsigned int j = 0; j < TNumNodes; j++)
                {
                    unsigned int row = j*LocalSize;
                    for (unsigned int d = 0; d < TDim;d++)
                        rLocalVector[row+d] -= Weight*N[j]*N[i]*pext*Normal[d];
                }
            }

            // Velocity inflow correction
            array_1d<double,3> Vel = ZeroVector(3);
            double Density = 0.0;

            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                const NodeType& rConstNode = this->GetGeometry()[i];
                Vel += N[i]*rConstNode.FastGetSolutionStepValue(VELOCITY);
                Density += N[i]*rConstNode.FastGetSolutionStepValue(DENSITY);
            }

            double Proj = 0.0; //Vel[0]*Normal[0] + Vel[1]*Normal[1] + Vel[2]*Normal[2];

            if (Proj < 0)
            {
                const double W = Weight*Density*Proj;
                for (unsigned int i = 0; i < TNumNodes; i++)
                {
                    double row = i*LocalSize;
                    for (unsigned int j = 0; j < TNumNodes; j++)
                    {
                        double col = j*LocalSize;
                        const array_1d<double,3>& rVel = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
                        for (unsigned int d = 0; d < TDim; d++)
                        {
                            double Tij = W*N[i]*N[j];
                            rLocalMatrix(row+d,col+d) -= Tij;
                            rLocalVector[row+d] += Tij*rVel[d];
                        }
                    }
                }
            }
        }
    }
}

template class MonolithicWallCondition<2,2>;
template class MonolithicWallCondition<3,3>;

} // namespace Kratos
