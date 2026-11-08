#pragma once 

#include "custom_utilities/compute_volume_utilities.h"

namespace Kratos
{

void ComputeVolumeUtilities::CalculateBoundaryNetNormal(Geometry<Node>& rGeom)
{

    std::size_t domain_size = rGeom.WorkingSpaceDimension();

    if (domain_size == 2){
        
        for (auto& edge : rGeom.GenerateEdges()){

            double area = edge.DomainSize();
            Geometry<Node>::CoordinatesArrayType local_coords;
        
            array_1d<double, 3> normal = edge.Normal(local_coords);
            int number_of_points_edge = edge.PointsNumber();
            auto& node_point = edge.Points();
            array_1d<double, 3> aux = (area * normal) / number_of_points_edge;
        
            Vector area_aux(domain_size);
        
            for (IndexType d = 0; d < domain_size; ++d) area_aux[d] = aux[d];
        
            for (IndexType i = 0; i < number_of_points_edge; ++i){
                auto& boundary_normal_area = node_point[i].GetValue(BOUNDARY_NORMAL_AREA);
                boundary_normal_area += area_aux;
                node_point[i].SetValue(BOUNDARY_NORMAL_AREA, boundary_normal_area);
            }
        }

    } else if (domain_size == 3){
        
        for (auto& face : rGeom.GenerateFaces()){

            double area = face.DomainSize();
            Geometry<Node>::CoordinatesArrayType local_coords;
        
            array_1d<double, 3> normal = face.Normal(local_coords);
            int number_of_points_face = face.PointsNumber();
            auto& node_point = face.Points();
            array_1d<double, 3> aux = (area * normal) / number_of_points_face;
        
            Vector area_aux(domain_size);
        
            for (IndexType d = 0; d < domain_size; ++d) area_aux[d] = aux[d];
        
            for (IndexType i = 0; i < number_of_points_face; ++i){
                auto& boundary_normal_area = node_point[i].GetValue(BOUNDARY_NORMAL_AREA);
                boundary_normal_area += area_aux;
                node_point[i].SetValue(BOUNDARY_NORMAL_AREA, boundary_normal_area);
            }
        }
    }
}

void ComputeVolumeUtilities::CheckBoundaryNetNormal(ModelPart::NodesContainerType& rNodes)
{
    std::size_t domain_size = rNodes.begin()->GetValue(BOUNDARY_NORMAL_AREA).size();
    Vector sum = ZeroVector(domain_size);

    for (auto& node : rNodes){
        auto& temp = node.GetValue(BOUNDARY_NORMAL_AREA);
        sum += temp;
    }

    if (norm_2(sum) > 1e-10){
        KRATOS_WARNING("ComputeVolumeProcess")<< "The domain boundary does not satisfy the zero net normal condition" <<std::endl;
        return;
    }
    KRATOS_INFO("ComputeVolumeProcess")<< "The domain boundary does satisfy the zero net normal condition" <<std::endl;


}

}