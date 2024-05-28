// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

// STL includes
#include <stdexcept>
// Project includes
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"
#include "quadrature/single_element.h"

namespace queso {

template<typename TElementType>
void QuadratureSingleElement<TElementType>::AssembleIPs(ElementType& rElement, const Vector3i& rOrder, IntegrationMethodType Method)
{
    auto& integration_points = rElement.GetIntegrationPoints();

    PointType lower_bound_param = rElement.GetBoundsUVW().first;
    PointType upper_bound_param = rElement.GetBoundsUVW().second;

    AssembleIPs(integration_points, lower_bound_param, upper_bound_param, rOrder, Method);
}

template<typename TElementType>
void QuadratureSingleElement<TElementType>::AssembleIPs(IntegrationPointVectorType& rIntegrationPoints, const PointType& rLowerBoundParam, const PointType& rUpperBoundParam,
        const Vector3i& rOrder, IntegrationMethodType Method ){

    const auto p_ip_list_u = IntegrationPointFactory1D::GetGauss(rOrder[0], Method);
    const auto p_ip_list_v = IntegrationPointFactory1D::GetGauss(rOrder[1], Method);
    const auto p_ip_list_w = IntegrationPointFactory1D::GetGauss(rOrder[2], Method);

    const SizeType n_point_u = p_ip_list_u->size();
    const SizeType n_point_v = p_ip_list_v->size();
    const SizeType n_point_w = p_ip_list_w->size();

    const SizeType n_point = (n_point_u)*(n_point_v)*(n_point_w);
    rIntegrationPoints.clear();
    rIntegrationPoints.reserve(n_point);

    const double length_u = std::abs(rUpperBoundParam[0] - rLowerBoundParam[0]);
    const double length_v = std::abs(rUpperBoundParam[1] - rLowerBoundParam[1]);
    const double length_w = std::abs(rUpperBoundParam[2] - rLowerBoundParam[2]);

    for (SizeType u = 0; u < p_ip_list_u->size(); ++u) {
        for (SizeType v = 0; v < p_ip_list_v->size(); ++v) {
            for( SizeType w = 0; w < p_ip_list_w->size(); ++w) {
                rIntegrationPoints.push_back( IntegrationPoint( rLowerBoundParam[0] + length_u * (*p_ip_list_u)[u][0],
                                                                rLowerBoundParam[1] + length_v * (*p_ip_list_v)[v][0],
                                                                rLowerBoundParam[2] + length_w * (*p_ip_list_w)[w][0],
                                                                (*p_ip_list_u)[u][1] * length_u *
                                                                (*p_ip_list_v)[v][1] * length_v *
                                                                (*p_ip_list_w)[w][1] * length_w ) );
            }
        }
    }
}

/// Explicit class instantiation
template class QuadratureSingleElement<Element<IntegrationPoint, BoundaryIntegrationPoint>>;

} // End namespace queso

// void QuadratureSingleElement::Assemble(
//     IntersectionTest& rInsideTest,
//     IntegrationPointType& rIntegrationPoints,
//     std::array<double,3> LocalPointA,
//     std::array<double,3> LocalPointB,
//     SizeType OrderU,
//     SizeType OrderV,
//     SizeType OrderW)
// {

//     const SizeType number_of_integration_points = (OrderU + 1)*(OrderV + 1)*(OrderW + 1);

//     IntegrationPointType tmp_integration_points;
//     tmp_integration_points.reserve(number_of_integration_points);
//     IntegrationPoints3D(
//         rInsideTest,
//         tmp_integration_points,
//         OrderU+1, OrderV+1, OrderW+1,
//         LocalPointA[0], LocalPointB[0],
//         LocalPointA[1], LocalPointB[1],
//         LocalPointA[2], LocalPointB[2]);

//     //#pragma omp critical
//     rIntegrationPoints.insert(rIntegrationPoints.end(), tmp_integration_points.begin(), tmp_integration_points.end());
// }


// void QuadratureSingleElement::IntegrationPoints3D(
//         IntersectionTest& rInsideTest,
//         IntegrationPointType& rIntegrationPoints,
//         SizeType PointsInU, SizeType PointsInV, SizeType PointsInW,
//         double U0, double U1, double V0, double V1, double W0, double W1)
// {

//    if(PointsInU < 1 || PointsInV < 1 || PointsInW < 1){
//         throw  std::invalid_argument("PointsInU, -V and -W need to be bigger than 0 - PointsInU:" + std::to_string(PointsInU)
//             + ", PointsInV:" + std::to_string(PointsInV) + " and PointsInW:" + std::to_string(PointsInW) + "\n");
//    }

//     const double distance_u = U1 - U0;
//     const double length_u = std::abs(U1 - U0);
//     const double distance_v = V1 - V0;
//     const double length_v = std::abs(V1 - V0);
//     const double distance_w = W1 - W0;
//     const double length_w = std::abs(W1 - W0);

//     const auto p_ip_list_u = IntegrationPointFactory1D::GetGauss(PointsInU, IntegrationPointFactory1D::IntegrationMethod::Gauss);
//     const auto integration_point_list_v = IntegrationPointFactory1D::GetGauss(PointsInV, IntegrationPointFactory1D::IntegrationMethod::Gauss);
//     const auto integration_point_list_w = IntegrationPointFactory1D::GetGauss(PointsInW, IntegrationPointFactory1D::IntegrationMethod::Gauss);

//     for (SizeType u = 0; u < PointsInU; ++u) {
//         for (SizeType v = 0; v < PointsInV; ++v) {
//             for( SizeType w = 0; w < PointsInW; ++w) {
//                 std::array<double,3> tmp_point = {  U0 + distance_u * (*p_ip_list_u)[u][0],
//                                                     V0 + distance_v * (*integration_point_list_v)[v][0],
//                                                     W0 + distance_w * (*integration_point_list_w)[w][0] };

//                 if( !rInsideTest.IsInsideLocalCoordinates(tmp_point) ){
//                     rIntegrationPoints.push_back( IntegrationPoint( tmp_point[0],
//                                                                     tmp_point[1],
//                                                                     tmp_point[2],
//                                                                     (*p_ip_list_u)[u][1] * length_u *
//                                                                     (*integration_point_list_v)[v][1] * length_v *
//                                                                     (*integration_point_list_w)[w][1] * length_w ) );
//                 }
//             }
//         }
//     }
// }