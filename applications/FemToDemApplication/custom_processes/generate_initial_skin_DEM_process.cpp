//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "custom_processes/generate_initial_skin_DEM_process.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

GenerateInitialSkinDEMProcess::GenerateInitialSkinDEMProcess(
    ModelPart& rModelPart,
    ModelPart& rDemModelPart)
    : mrModelPart(rModelPart),
      mrDEMModelPart(rDemModelPart)
{
}

/***********************************************************************************/
/***********************************************************************************/

void GenerateInitialSkinDEMProcess::Execute() 
{
    auto nodal_neigh_process = FindNodalNeighboursProcess(mrModelPart, 5, 5);
    nodal_neigh_process.Execute();

    const auto it_node_begin = mrModelPart.NodesBegin();
    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); i++) {
        auto it_node = it_node_begin + i;

    }
}


/***********************************************************************************/
/***********************************************************************************/

void GenerateInitialSkinDEMProcess::CreateDEMParticle(
    const int Id,
    const array_1d<double, 3> Coordinates,
    const Properties::Pointer pProperties,
    const double Radius,
    NodeType& rNode
)
{
    auto spheric_particle = mParticleCreator.CreateSphericParticleRaw(mrDEMModelPart, Id, Coordinates, pProperties, Radius, "SphericParticle3D");
    rNode.SetValue(IS_DEM, true);
    rNode.SetValue(RADIUS, Radius);
    rNode.SetValue(DEM_PARTICLE_POINTER, spheric_particle);
}


/***********************************************************************************/
/***********************************************************************************/

double GenerateInitialSkinDEMProcess::CalculateDistanceBetweenNodes(
    const NodeType& rNode1, 
    const NodeType& rNode2
    )
{
    const double X1 = rNode1.X();
    const double X2 = rNode2.X();
    const double Y1 = rNode1.Y();
    const double Y2 = rNode2.Y();
    const double Z1 = rNode1.Z();
    const double Z2 = rNode2.Z();
    return std::sqrt(std::pow(X1-X2, 2) + std::pow(Y1-Y2, 2) + std::pow(Z1-Z2, 2));
}

/***********************************************************************************/
/***********************************************************************************/

double GenerateInitialSkinDEMProcess::GetMinimumValue(
    const Vector& rValues
    )
{ // this method assumes that the Vector is NOT full of 0.0's
    double aux = 1.0e10;
    for (int i = 0; i < rValues.size(); i++) 
        if (aux > rValues[i] && rValues[i] != 0.0) 
            aux = rValues[i];
    return aux;
}

/***********************************************************************************/
/***********************************************************************************/

int GenerateInitialSkinDEMProcess::GetMaximumDEMId()
{
    int max_id = 0;
    const auto it_DEM_begin = mrDEMModelPart.ElementsBegin();
    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrDEMModelPart.Elements().size()); i++) {
        auto it_DEM = it_DEM_begin + i;
        auto& r_geometry = it_DEM->GetGeometry();
        const int DEM_id = r_geometry[0].Id();
        max_id = (max_id < DEM_id) ? DEM_id : max_id;
    }
    return max_id;
}

}  // namespace Kratos
