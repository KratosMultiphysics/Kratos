//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/generate_dem_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

template <SizeType TDim>
GenerateDemProcess<TDim>::GenerateDemProcess(
    ModelPart& rModelPart,
    ModelPart& rDemModelPart)
    : mrModelPart(rModelPart),
      mrDEMModelPart(rDemModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void GenerateDemProcess<2>::Execute() 
{
    ParticleCreatorDestructor particle_creator = ParticleCreatorDestructor();



}

/***********************************************************************************/
/***********************************************************************************/

template <>
void GenerateDemProcess<3>::Execute() 
{
    ParticleCreatorDestructor particle_creator = ParticleCreatorDestructor();



}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
double GenerateDemProcess<TDim>::CalculateDistanceBetweenNodes(
    const NodeType Node1, 
    const NodeType Node2
    )
{
    const double X1 = Node1.X();
    const double X2 = Node2.X();
    const double Y1 = Node1.Y();
    const double Y2 = Node2.Y();
    const double Z1 = Node1.Z();
    const double Z2 = Node2.Z();
    return std::sqrt(std::pow(X1-X2, 2.0) + std::pow(Y1-Y2, 2.0) + std::pow(Z1-Z2, 2.0));
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
array_1d<double,3> GenerateDemProcess<TDim>::GetNodeCoordinates(
    const NodeType Node
    )
{
    array_1d<double,3> coordinates;
    coordinates[0] = Node.X();
    coordinates[1] = Node.Y();
    coordinates[2] = Node.Z();
    return coordinates;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
double GenerateDemProcess<TDim>::GetMinimumValue2(
    const double a, 
    const double b
    )
{
    if (a < b) return a;
    else return b;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
double GenerateDemProcess<TDim>::GetMinimumValue3(
    const double a, 
    const double b,
    const double c
    )
{
    double aux = a;
    if (aux > b) aux = b;
    if (aux > c) aux = c;
    return aux;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenerateDemProcess<2>;
template class GenerateDemProcess<3>;

}  // namespace Kratos