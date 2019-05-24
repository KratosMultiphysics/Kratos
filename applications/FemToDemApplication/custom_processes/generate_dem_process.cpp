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
    auto nodal_neigh_process = FindNodalNeighboursProcess(mrModelPart, 5, 5);
    nodal_neigh_process.Execute();

    const auto it_element_begin = mrModelPart.ElementsBegin();
    // #pragma omp parallel for
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

            // We perform a loop over the nodes to create its DEM
            for (int i = 0; i < r_geom.size(); i++) {
                auto& r_node = r_geom[i];
                if (!r_node.GetValue(IS_DEM)) {
                    auto& r_neigh_nodes = r_node.GetValue(NEIGHBOUR_NODES);
                    Vector potential_radii(r_neigh_nodes.size());
                    Vector distances(r_neigh_nodes.size());

                    // Loop over the neighbours of that node
                    bool has_dem_neigh = false;
                    for (int neigh = 0; neigh < r_neigh_nodes.size(); neigh++) {
                        auto& r_neighbour = r_neigh_nodes[neigh];
                        const double dist_between_nodes = this->CalculateDistanceBetweenNodes(r_node, r_neighbour);
                        distances(neigh) = dist_between_nodes;
                        if (r_neighbour.GetValue(IS_DEM)) {
                            has_dem_neigh = true;
                            potential_radii(neigh) = dist_between_nodes - r_neighbour.GetValue(RADIUS);
                            if (potential_radii(neigh) < 0.0) { // Houston-> We have a problem
                                const double new_radius = dist_between_nodes*0.5;
                                auto& r_radius_neigh_old = mrDEMModelPart.GetNode(r_neighbour.Id()).GetSolutionStepValue(RADIUS);
                                r_radius_neigh_old = new_radius;
                                r_neighbour.SetValue(RADIUS, new_radius);
                                potential_radii(neigh) = new_radius;
                            }
                        } else {
                            potential_radii(neigh) = 0.0;
                        }
                    }
                    // Let's compute the Radius of the new DEM
                    double radius;
                    if (has_dem_neigh) {
                        radius = this->GetMinimumValue(potential_radii);
                    } else {
                        radius = this->GetMinimumValue(distances)*0.5;
                    }
                    const array_1d<double,3>& r_coordinates = this->GetNodeCoordinates(r_node);
                    this->CreateDEMParticle(r_node.Id(), r_coordinates, p_DEM_properties, radius, r_node);
                }
            }

            it_elem->SetValue(DEM_GENERATED, true);
            it_elem->Set(TO_ERASE, true);
        } else if (!is_active && dem_generated)
            it_elem->Set(TO_ERASE, true);
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

template <SizeType TDim>
void GenerateDemProcess<TDim>::ModifyRadiusToNodes(
    NodeType& rNode1,
    NodeType& rNode2,
    const double NewR1,
    const double NewR2
    )
{
    auto& r_radius_1_old = mrDEMModelPart.GetNode(rNode1.Id()).GetSolutionStepValue(RADIUS);
    r_radius_1_old = NewR1;
    auto& r_radius_2_old = mrDEMModelPart.GetNode(rNode2.Id()).GetSolutionStepValue(RADIUS);
    r_radius_2_old = NewR2;
    rNode1.SetValue(RADIUS, NewR1);
    rNode2.SetValue(RADIUS, NewR2);
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
double GenerateDemProcess<TDim>::GetMinimumValue(
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