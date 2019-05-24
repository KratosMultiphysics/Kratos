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
    const auto it_element_begin = mrModelPart.ElementsBegin();
    //#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++) {
        auto it_elem = it_element_begin + i;
        auto& r_geom = it_elem->GetGeometry();

        bool is_active = true;
        if (it_elem->IsDefined(ACTIVE))
            is_active = it_elem->Is(ACTIVE);

        bool dem_generated = it_elem->GetValue(DEM_GENERATED);
        const int number_of_dem = this->GetNumberOfDemOnElement(it_elem);

        if (!is_active && !dem_generated) {
            // std::cout << "Elemento eliminado -- >" << it_elem->Id() << std::endl;
            auto p_DEM_properties = mrDEMModelPart.pGetProperties(1);
			auto& node0 = r_geom[0];
			auto& node1 = r_geom[1];
			auto& node2 = r_geom[2];
            const double dist01 = this->CalculateDistanceBetweenNodes(node0, node1);
            const double dist02 = this->CalculateDistanceBetweenNodes(node0, node2);
            const double dist12 = this->CalculateDistanceBetweenNodes(node1, node2);

            if (number_of_dem == 0) { // We must create 3 DEM
                // Node 0
                const double r0 = this->GetMinimumValue2(dist01, dist02)*0.5;
                const array_1d<double, 3> coordinates0 = this->GetNodeCoordinates(node0);
                this->CreateDEMParticle(node0.Id(), coordinates0, p_DEM_properties, r0, node0);
                // Node 1
                const double r1 = dist01 - r0;
                const array_1d<double, 3> coordinates1 = this->GetNodeCoordinates(node1);
                this->CreateDEMParticle(node1.Id(), coordinates1, p_DEM_properties, r1, node1);
                // Node 2
                const double r2 = this->GetMinimumValue2(dist02 - r0, dist12 - r1);
                const array_1d<double, 3> coordinates2 = this->GetNodeCoordinates(node2);
                this->CreateDEMParticle(node2.Id(), coordinates2, p_DEM_properties, r2, node2);

            } else if (number_of_dem == 2) { // We must create 1 DEM
                const int local_id_no_DEM = this->GetLocalIdWithoutDEM(it_elem);
                auto& node = r_geom[local_id_no_DEM];
                Matrix distances(3,3);
                this->CreateDistancesMatrix(distances, dist01, dist02, dist12);
                int local_id_dem_1, local_id_dem_2;
                this->Get2LocalIdFrom1(local_id_no_DEM, local_id_dem_1, local_id_dem_2);
                const double r1 = r_geom[local_id_dem_1].GetValue(RADIUS);
                const double r2 = r_geom[local_id_dem_2].GetValue(RADIUS);
                const double radius = this->GetMinimumValue2(distances(local_id_dem_1, local_id_no_DEM)-r1, distances(local_id_dem_2, local_id_no_DEM)-r2);
                const array_1d<double, 3> coordinates = this->GetNodeCoordinates(node);
                this->CreateDEMParticle(node.Id(), coordinates, p_DEM_properties, radius, node);

            } else if (number_of_dem == 1) { // We must create 2 DEM
                const int local_id_DEM = this->GetLocalIdWithDEM(it_elem);

            }









        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void GenerateDemProcess<3>::Execute() 
{

}

/***********************************************************************************/
/***********************************************************************************/

template <>
void GenerateDemProcess<2>::CreateDistancesMatrix(
    Matrix& rDistancesMatrix,
    const double d01,
    const double d02,
    const double d12
    )
{
    rDistancesMatrix(0,0) = 0.0;
    rDistancesMatrix(1,1) = 0.0;
    rDistancesMatrix(2,2) = 0.0;
    rDistancesMatrix(0,1) = d01;
    rDistancesMatrix(0,2) = d02;
    rDistancesMatrix(1,2) = d12;

    // Symmetric
    rDistancesMatrix(1,0) = rDistancesMatrix(0,1);
    rDistancesMatrix(2,0) = rDistancesMatrix(0,2);
    rDistancesMatrix(2,1) = rDistancesMatrix(1,2);
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void GenerateDemProcess<2>::Get2LocalIdFrom1(
    const int FixedLocalId,
    int& LocalIDWithDEM1,
    int& LocalIDWithDEM2
    )
{
    LocalIDWithDEM1 = (FixedLocalId == 0) ? 1 : (FixedLocalId == 1) ? 0 : 1;
    LocalIDWithDEM2 = (FixedLocalId == 0) ? 2 : (FixedLocalId == 1) ? 2 : 0;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
int GenerateDemProcess<TDim>::GetLocalIdWithoutDEM(
    ElementIteratorType ItElem
    )
{
    // We use this method when 1 DEM is missing
    auto& r_geom = ItElem->GetGeometry();
    for (int i = 0; i < r_geom.size(); i++) {
        if (!r_geom[i].GetValue(IS_DEM)) {
            return i;
        }
    }
    return -1;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
int GenerateDemProcess<TDim>::GetLocalIdWithDEM(
    ElementIteratorType ItElem
    )
{
    // We use this method when there is only 1 DEM
    auto& r_geom = ItElem->GetGeometry();
    for (int i = 0; i < r_geom.size(); i++) {
        if (r_geom[i].GetValue(IS_DEM)) {
            return i;
        }
    }
    return -1;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
void GenerateDemProcess<TDim>::CreateDEMParticle(
    const int Id,
    const array_1d<double, 3> Coordinates,
    const Properties::Pointer pProperties,
    const double Radius,
    NodeType& rNode
)
{
    mParticleCreator.CreateSphericParticle(mrDEMModelPart, Id, Coordinates, pProperties, Radius, "SphericParticle3D");
    rNode.SetValue(IS_DEM, true);
    rNode.SetValue(RADIUS, Radius);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
int GenerateDemProcess<TDim>::GetNumberOfDemOnElement(
    ElementIteratorType ItElem
    )
{
    int number_dem = 0;
    auto& r_geom = ItElem->GetGeometry();
    for (int i = 0; i < r_geom.size(); i++) {
        auto& r_node = r_geom[i];
        if (r_node.GetValue(IS_DEM)) 
            number_dem++;
    }
    return number_dem;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
double GenerateDemProcess<TDim>::CalculateDistanceBetweenNodes(
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
    return std::sqrt(std::pow(X1-X2, 2.0) + std::pow(Y1-Y2, 2.0) + std::pow(Z1-Z2, 2.0));
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDim>
array_1d<double,3> GenerateDemProcess<TDim>::GetNodeCoordinates(
    const NodeType& rNode
    )
{
    array_1d<double,3> coordinates;
    coordinates[0] = rNode.X();
    coordinates[1] = rNode.Y();
    coordinates[2] = rNode.Z();
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