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
#include "hypo.h"

namespace Kratos
{
///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see VMS::EquationIdVector
 */
template <>
void HYPO<2>::EquationIdVector(EquationIdVectorType& rResult,
                              const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int NumNodes(3),LocalSize(6);
    unsigned int LocalIndex = 0;

    unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    //unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);
    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X,vpos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y,vpos+1).EquationId();
        //KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl; 
        //rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE,ppos).EquationId();
    }
}

/**
 * @see HYPO::EquationIdVector
 */
template <>
void HYPO<3>::EquationIdVector(EquationIdVectorType& rResult,
                              const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int NumNodes(4),LocalSize(12);
    unsigned int LocalIndex = 0;
    unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    //unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X,vpos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y,vpos+1).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z,vpos+2).EquationId();
        //rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE,ppos).EquationId();
    }
}

/**
 * @see HYPO::GetDofList
 */
template <>
void HYPO<2>::GetDofList(DofsVectorType& rElementalDofList,
                        const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int NumNodes(3),LocalSize(6);
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        //rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see HYPO::GetDofList
 */
template <>
void HYPO<3>::GetDofList(DofsVectorType& rElementalDofList,
                        const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int NumNodes(4),LocalSize(12);
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        //rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see HYPO::GetFirstDerivativesVector
 */
template <>
void HYPO<2>::GetFirstDerivativesVector(Vector& Values, int Step) const
{
    const unsigned int NumNodes(3),LocalSize(6);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        const array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        Values[LocalIndex++] = rVelocity[0];
        Values[LocalIndex++] = rVelocity[1];
        //Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

/**
 * @see HYPO::GetFirstDerivativesVector
 */
template <>
void HYPO<3>::GetFirstDerivativesVector(Vector& Values, int Step) const
{
    const unsigned int NumNodes(4),LocalSize(12);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        const array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        Values[LocalIndex++] = rVelocity[0];
        Values[LocalIndex++] = rVelocity[1];
        Values[LocalIndex++] = rVelocity[2];
        //Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
    }
}

/**
 * @see HYPO::GetSecondDerivativesVector
 */
template <>
void HYPO<2>::GetSecondDerivativesVector(Vector& Values, int Step) const
{
    const unsigned int NumNodes(3),LocalSize(6);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        const array_1d<double,3>& rAcceleration = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        Values[LocalIndex++] = rAcceleration[0];
        Values[LocalIndex++] = rAcceleration[1];
        
        
    }
    //KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha","");
}

/**
 * @see HYPO::GetSecondDerivativesVector
 */
template <>
void HYPO<3>::GetSecondDerivativesVector(Vector& Values, int Step) const
{
    const unsigned int NumNodes(4),LocalSize(12);
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        const array_1d<double,3>& rAcceleration = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
        Values[LocalIndex++] = rAcceleration[0];
        Values[LocalIndex++] = rAcceleration[1];
        Values[LocalIndex++] = rAcceleration[2];
        //Values[LocalIndex++] = 0.0; // Pressure Dof
    }
}

/**
 * The size of the 2D element is estimated as the diameter of a circle of the same area.
 * Area = Pi * (h/2)^2
 * @see HYPO::ElementSize
 */
template <>
double HYPO<2,3>::ElementSize(const double Area)
{
    return 1.128379167 * sqrt(Area); //Diameter of circumference of given Area
}

/**
 * The size of the 3D element is estimated as the diameter of the sphere
 * circumscribed to a regular tetrahedron with the same volume.
 * @see HYPO::ElementSize
 */
template <>
double HYPO<3,4>::ElementSize(const double Volume)
{
    return 0.60046878 * pow(Volume,0.333333333333333333333);
}


/**
 * @see HYPO::CalculateOnIntegrationPoints
 */
template <>
void HYPO<2>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3> >& rVariable,
    std::vector<array_1d<double,3> >& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{}

/**
 * @see VMS::CalculateOnIntegrationPoints
 */
template <>
void HYPO<3>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3> >& rVariable,
    std::vector<array_1d<double,3> >& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{}


template<>
double HYPO<2,3>::ConsistentMassCoef(const double Area)
{
    const double Coef = 1.0/12.0;
    return Area * Coef;
}

template<>
double HYPO<3,4>::ConsistentMassCoef(const double Volume)
{
    const double Coef = 1.0/20.0;
    return Volume * Coef;
}

///@} // Specialized implementations
}
