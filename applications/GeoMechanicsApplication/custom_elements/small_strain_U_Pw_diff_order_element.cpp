// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// Project includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/node_utilities.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

namespace Kratos
{

Vector GetPressures(const Geometry<Node>& rGeometry, size_t n_nodes)
{
    Vector pressure(n_nodes);
    std::transform(rGeometry.begin(), rGeometry.begin() + n_nodes, pressure.begin(),
                   [](const auto& node) { return node.FastGetSolutionStepValue(WATER_PRESSURE); });
    return pressure;
}

template <>
void SmallStrainUPwDiffOrderElement<2, 6>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // 2D T6P3
    mpPressureGeometry = make_shared<Triangle2D3<Node>>(r_geometry(0), r_geometry(1), r_geometry(2));
}

template <>
void SmallStrainUPwDiffOrderElement<2, 8>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // 2D Q8P4
    mpPressureGeometry =
        make_shared<Quadrilateral2D4<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3));
}

template <>
void SmallStrainUPwDiffOrderElement<2, 9>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // 2D Q9P4
    mpPressureGeometry =
        make_shared<Quadrilateral2D4<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3));
}

template <>
void SmallStrainUPwDiffOrderElement<3, 10>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    //  3D T10P4
    mpPressureGeometry =
        make_shared<Tetrahedra3D4<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3));
}

template <>
void SmallStrainUPwDiffOrderElement<2, 10>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // 2D T10P6
    mpPressureGeometry = make_shared<Triangle2D6<Node>>(
        r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3), r_geometry(4), r_geometry(5));
}

template <>
void SmallStrainUPwDiffOrderElement<2, 15>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // 2D T15P10
    mpPressureGeometry = make_shared<Triangle2D10<Node>>(
        r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3), r_geometry(4), r_geometry(5),
        r_geometry(6), r_geometry(7), r_geometry(8), r_geometry(9));
}

template <>
void SmallStrainUPwDiffOrderElement<3, 20>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // 3D H20P8
    mpPressureGeometry =
        make_shared<Hexahedra3D8<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3),
                                        r_geometry(4), r_geometry(5), r_geometry(6), r_geometry(7));
}

template <>
void SmallStrainUPwDiffOrderElement<3, 27>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // 3D H27P8
    mpPressureGeometry =
        make_shared<Hexahedra3D8<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3),
                                        r_geometry(4), r_geometry(5), r_geometry(6), r_geometry(7));
}

void set_arithmetic_average_pressure(Geometry<Node>&                               rGeometry,
                                     const Vector&                                 rPressure,
                                     const std::vector<std::pair<size_t, size_t>>& rIndexPpairs,
                                     size_t DestinationOffset = 0)
{
    for (size_t i = 0; const auto& [first_index, second_index] : rIndexPpairs) {
        NodeUtilities::ThreadSafeNodeWrite(rGeometry[DestinationOffset + i], WATER_PRESSURE,
                                           0.5 * (rPressure[first_index] + rPressure[second_index]));
        ++i;
    }
}

void set_arithmetic_average_pressure(Geometry<Node>& rGeometry,
                                     const Vector&   rPressure,
                                     const std::vector<std::tuple<size_t, size_t, size_t, size_t>>& rIndices,
                                     size_t DestinationOffset = 0)
{
    for (size_t i = 0; const auto& [first_index, second_index, third_index, fourth_index] : rIndices) {
        NodeUtilities::ThreadSafeNodeWrite(rGeometry[DestinationOffset + i], WATER_PRESSURE,
                                           0.25 * (rPressure[first_index] + rPressure[second_index] +
                                                   rPressure[third_index] + rPressure[fourth_index]));
        ++i;
    }
}

template <>
void SmallStrainUPwDiffOrderElement<2, 6>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();

    // 2D T6P3
    const Vector                                 pressure = GetPressures(this->GetGeometry(), 3);
    const std::vector<std::pair<size_t, size_t>> pairs    = {{0, 1}, {1, 2}, {2, 0}};
    set_arithmetic_average_pressure(r_geom, pressure, pairs, 3);
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<2, 8>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 2D Q8P4
    const Vector                                 pressure = GetPressures(this->GetGeometry(), 4);
    const std::vector<std::pair<size_t, size_t>> pairs    = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    set_arithmetic_average_pressure(r_geom, pressure, pairs, 4);
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<2, 9>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 2D Q9P4
    const Vector                                 pressure = GetPressures(this->GetGeometry(), 4);
    const std::vector<std::pair<size_t, size_t>> pairs    = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    set_arithmetic_average_pressure(r_geom, pressure, pairs, 4);
    const std::vector<std::tuple<size_t, size_t, size_t, size_t>>& indices = {{0, 1, 2, 3}};
    set_arithmetic_average_pressure(r_geom, pressure, indices, 8);
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<3, 10>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 3D T10P4
    const Vector                                 pressure = GetPressures(this->GetGeometry(), 4);
    const std::vector<std::pair<size_t, size_t>> pairs    = {{0, 1}, {1, 2}, {2, 0},
                                                             {0, 3}, {1, 3}, {2, 3}};
    set_arithmetic_average_pressure(r_geom, pressure, pairs, 4);
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<2, 10>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 2D T10P6
    constexpr double c1 = 1.0 / 9.0;
    const Vector     p  = GetPressures(this->GetGeometry(), 6);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[3], WATER_PRESSURE, (2.0 * p[0] - p[1] + 8.0 * p[3]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE, (2.0 * p[1] - p[0] + 8.0 * p[3]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE, (2.0 * p[1] - p[2] + 8.0 * p[4]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE, (2.0 * p[2] - p[1] + 8.0 * p[4]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE, (2.0 * p[2] - p[0] + 8.0 * p[5]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE, (2.0 * p[0] - p[2] + 8.0 * p[5]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[9], WATER_PRESSURE,
                                       (4.0 * (p[3] + p[4] + p[5]) - (p[0] + p[1] + p[2])) * c1);
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<2, 15>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 2D T15P10
    constexpr double c1 = 0.0390625;
    const Vector     p  = GetPressures(this->GetGeometry(), 10);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[3], WATER_PRESSURE,
                                       (3.0 * p[0] + p[1] + 27.0 * p[3] - 5.4 * p[4]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE,
                                       (14.4 * (p[3] + p[4]) - 1.6 * (p[0] + p[1])) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE,
                                       (3.0 * p[1] + p[0] + 27.0 * p[4] - 5.4 * p[3]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE,
                                       (3.0 * p[1] + p[2] + 27.0 * p[5] - 5.4 * p[6]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE,
                                       (14.4 * (p[5] + p[6]) - 1.6 * (p[1] + p[2])) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE,
                                       (3.0 * p[2] + p[1] + 27.0 * p[6] - 5.4 * p[5]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[9], WATER_PRESSURE,
                                       (3.0 * p[2] + p[0] + 27.0 * p[7] - 5.4 * p[8]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[10], WATER_PRESSURE,
                                       (14.4 * (p[7] + p[8]) - 1.6 * (p[0] + p[2])) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[11], WATER_PRESSURE,
                                       (3.0 * p[0] + p[2] + 27.0 * p[8] - 5.4 * p[7]) * c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[12], WATER_PRESSURE,
                                       (p[1] + p[2] + 7.2 * (p[3] + p[8]) - 3.6 * (p[4] + p[7]) -
                                        1.8 * (p[5] + p[6]) + 21.6 * p[9] - 1.6 * p[0]) *
                                           c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[13], WATER_PRESSURE,
                                       (p[0] + p[2] + 7.2 * (p[4] + p[5]) - 3.6 * (p[3] + p[6]) -
                                        1.8 * (p[7] + p[8]) + 21.6 * p[9] - 1.6 * p[1]) *
                                           c1);
    NodeUtilities::ThreadSafeNodeWrite(r_geom[14], WATER_PRESSURE,
                                       (p[0] + p[1] + 7.2 * (p[6] + p[7]) - 3.6 * (p[5] + p[8]) -
                                        1.8 * (p[3] + p[4]) + 21.6 * p[9] - 1.6 * p[2]) *
                                           c1);
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<3, 20>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 3D H20P8
    const Vector                                 pressure = GetPressures(this->GetGeometry(), 8);
    const std::vector<std::pair<size_t, size_t>> pairs =
        // edges -- bottom
        {{0, 1},
         {1, 2},
         {2, 3},
         {3, 0},
         // edges -- middle
         {4, 0},
         {5, 1},
         {6, 2},
         {7, 3},
         // edges -- top
         {4, 5},
         {5, 6},
         {6, 7},
         {7, 4}};
    set_arithmetic_average_pressure(r_geom, pressure, pairs, 8);
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<3, 27>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 3D H27P8
    const Vector                                 pressure = GetPressures(this->GetGeometry(), 8);
    const std::vector<std::pair<size_t, size_t>> pairs =
        // edges -- bottom
        {{0, 1},
         {1, 2},
         {2, 3},
         {3, 0},
         // edges -- middle
         {4, 0},
         {5, 1},
         {6, 2},
         {7, 3},
         // edges -- top
         {4, 5},
         {5, 6},
         {6, 7},
         {7, 0}};
    set_arithmetic_average_pressure(r_geom, pressure, pairs, 8);
    // face centers
    const std::vector<std::tuple<size_t, size_t, size_t, size_t>>& indices = {
        {0, 1, 2, 3}, {0, 1, 4, 5}, {1, 2, 5, 6}, {2, 3, 6, 7}, {3, 0, 7, 4}, {4, 5, 6, 7}};
    set_arithmetic_average_pressure(r_geom, pressure, indices, 20);
    // element center
    NodeUtilities::ThreadSafeNodeWrite(r_geom[26], WATER_PRESSURE,
                                       0.125 * (pressure[0] + pressure[1] + pressure[2] + pressure[3] +
                                                pressure[4] + pressure[5] + pressure[6] + pressure[7]));
    KRATOS_CATCH("")
}

template class SmallStrainUPwDiffOrderElement<2, 6>;
template class SmallStrainUPwDiffOrderElement<2, 8>;
template class SmallStrainUPwDiffOrderElement<2, 9>;
template class SmallStrainUPwDiffOrderElement<2, 10>;
template class SmallStrainUPwDiffOrderElement<2, 15>;
template class SmallStrainUPwDiffOrderElement<3, 10>;
template class SmallStrainUPwDiffOrderElement<3, 20>;
template class SmallStrainUPwDiffOrderElement<3, 27>;
} // Namespace Kratos
