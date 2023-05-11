//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//
#include "HypoElasticSolid.h"

namespace Kratos
{
///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see VMS::EquationIdVector
 */
template <>
void HYPOELASTICSOLID<2>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) 
{
    const unsigned int NumNodes(3),LocalSize(9);
    unsigned int LocalIndex = 0;

KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X,vpos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y,vpos+1).EquationId();
       // rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE,ppos).EquationId();
    }
}


template <>
void HYPOELASTICSOLID<3>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(4),LocalSize(12);
    unsigned int LocalIndex = 0;
    unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
 //   unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
//KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
//KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
//KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X,vpos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y,vpos+1).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z,vpos+2).EquationId();
        //rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE,ppos).EquationId();
    }
    
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha",""); 	

}

/*
 * @see HYPOELASTICSOLID::GetDofList
 */
/*template <>
void HYPOELASTICSOLID<2>::GetDofList(DofsVectorType& rElementalDofList,
                        const ProcessInfo& rCurrentProcessInfo) 
{
    KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha",""); 	
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;
 

    if(rElementalDofList.size() != number_of_nodes*dim)
        rElementalDofList.resize(number_of_nodes*dim);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rElementalDofList[i*dim] = GetGeometry()[i].pGetDof(VELOCITY_X);
        rElementalDofList[i*dim+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
    }
}
*/
/**
 * @see HYPOELASTICSOLID::GetDofList
 */
/*template <>
void HYPOELASTICSOLID<3>::GetDofList(DofsVectorType& rElementalDofList,
                        const ProcessInfo& rCurrentProcessInfo) 
{
    const unsigned int NumNodes(4),LocalSize(12);//LocalSize(16)
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

   KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha",""); 	
    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        //rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}*/
*/
/**
 * @see HYPOELASTICSOLID::GetFirstDerivativesVector
 */
/*template <>
void HYPOELASTICSOLID<2>::GetFirstDerivativesVector(Vector& Values, int Step) 
{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;
    if(Values.size() != MatSize)   Values.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i*dim;
        const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,Step);
        Values[index] = vel[0];
        Values[index + 1] = vel[1];
    }
}*/

/**
 * @see HYPOELASTICSOLID::GetFirstDerivativesVector
 */
/*template <>
void HYPOELASTICSOLID<3>::GetFirstDerivativesVector(Vector& Values, int Step) 
{

   KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha",""); 	

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
    KRATOS_WATCH("fluiddddddddddddddddddo")
    KRATOS_WATCH("fluiddddddddddddddddddo")
    KRATOS_WATCH("fluiddddddddddddddddddo")
    KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha",""); 	

}
*/
/**
 * @see HYPOELASTICSOLID::GetSecondDerivativesVector
 */
/*template <>
void HYPOELASTICSOLID<2>::GetSecondDerivativesVector(Vector& Values, int Step) 
{
KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;
    if(Values.size() != MatSize) Values.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i*dim;
        const array_1d<double,3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION,Step);
        Values[index] = acc[0];
        Values[index + 1] = acc[1];
    }
}
*/

/**
 * @see HYPOELASTICSOLID::GetSecondDerivativesVector
 */
/*template <>
void HYPOELASTICSOLID<3>::GetSecondDerivativesVector(Vector& Values, int Step) 
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes*dim;
    if(Values.size() != MatSize) Values.resize(MatSize,false);

    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
    KRATOS_WATCH("SOLIDOOOOOOOOOOOOOOOOOOOOOOOOOO")
KRATOS_THROW_ERROR(std::logic_error,"not dereeeeeeeeeeeeeeeeeeeeeeecha",""); 	
    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {

        unsigned int index = i*dim;
        const array_1d<double,3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION,Step);
        Values[index] = acc[0];
        Values[index + 1] = acc[1];
        Values[index + 2] = acc[2];
    }
}
*/
/**
 * The size of the 2D element is estimated as the diameter of a circle of the same area.
 * Area = Pi * (h/2)^2
 * @see HYPOELASTICSOLID::ElementSize
 */

/**
 * The size of the 3D element is estimated as the diameter of the sphere
 * circumscribed to a regular tetrahedron with the same volume.
 * @see HYPOELASTICSOLID::ElementSize
 */


/**
 * @see HYPOELASTICSOLID::CalculateOnIntegrationPoints
 */
/*template <>
void HYPOELASTICSOLID<2>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3> >& rVariable,
    std::vector<array_1d<double,3> >& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int Dim(2),NumNodes(3);
    
}*/

/**
 * @see HYPOELASTICSOLID::CalculateOnIntegrationPoints
 */
/*template <>
void HYPOELASTICSOLID<3>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3> >& rVariable,
    std::vector<array_1d<double,3> >& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int Dim(3),NumNodes(4);
    
}
*/



/**
 * See HYPOELASTICSOLID::CalculateB
 */


/**
 * See HYPOELASTICSOLID::CalculateB
 */


/**
 * See HYPOELASTICSOLID::CalculateC
 */

/**
 * See HYPOELASTICSOLID::CalculateC
 */


/**
 * @see HYPOELASTICSOLID::AddViscousTerm
 */


/**
 * @see HYPOELASTICSOLID::AddViscousTerm
 */

///@} // Specialized implementations
}
