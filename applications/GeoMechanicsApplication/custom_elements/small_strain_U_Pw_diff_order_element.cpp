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
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

namespace Kratos
{
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

template <>
void SmallStrainUPwDiffOrderElement<2, 6>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();

    // 2D T6P3
    const Vector p = GetPressures(3);
    ThreadSafeNodeWrite(r_geom[3], WATER_PRESSURE, 0.5 * (p[0] + p[1]));
    ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE, 0.5 * (p[1] + p[2]));
    ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE, 0.5 * (p[2] + p[0]));
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<2, 8>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 2D Q8P4
    const Vector p = GetPressures(4);
    ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE, 0.5 * (p[0] + p[1]));
    ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE, 0.5 * (p[1] + p[2]));
    ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE, 0.5 * (p[2] + p[3]));
    ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE, 0.5 * (p[3] + p[0]));
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<2, 9>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 2D Q9P4
    const Vector p = GetPressures(4);
    ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE, 0.5 * (p[0] + p[1]));
    ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE, 0.5 * (p[1] + p[2]));
    ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE, 0.5 * (p[2] + p[3]));
    ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE, 0.5 * (p[3] + p[0]));
    ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE, 0.25 * (p[0] + p[1] + p[2] + p[3]));
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<3, 10>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 3D T10P4  //2D T10P6
    const Vector p = GetPressures(4);
    ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE, 0.5 * (p[0] + p[1]));
    ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE, 0.5 * (p[1] + p[2]));
    ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE, 0.5 * (p[2] + p[0]));
    ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE, 0.5 * (p[0] + p[3]));
    ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE, 0.5 * (p[1] + p[3]));
    ThreadSafeNodeWrite(r_geom[9], WATER_PRESSURE, 0.5 * (p[2] + p[3]));
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<2, 10>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType&    r_geom = GetGeometry();
    constexpr double c1     = 1.0 / 9.0;
    const Vector     p      = GetPressures(6);
    ThreadSafeNodeWrite(r_geom[3], WATER_PRESSURE, (2.0 * p[0] - p[1] + 8.0 * p[3]) * c1);
    ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE, (2.0 * p[1] - p[0] + 8.0 * p[3]) * c1);
    ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE, (2.0 * p[1] - p[2] + 8.0 * p[4]) * c1);
    ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE, (2.0 * p[2] - p[1] + 8.0 * p[4]) * c1);
    ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE, (2.0 * p[2] - p[0] + 8.0 * p[5]) * c1);
    ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE, (2.0 * p[0] - p[2] + 8.0 * p[5]) * c1);
    ThreadSafeNodeWrite(r_geom[9], WATER_PRESSURE, (4.0 * (p[3] + p[4] + p[5]) - (p[0] + p[1] + p[2])) * c1);
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
    const Vector     p  = GetPressures(10);
    ThreadSafeNodeWrite(r_geom[3], WATER_PRESSURE, (3.0 * p[0] + p[1] + 27.0 * p[3] - 5.4 * p[4]) * c1);
    ThreadSafeNodeWrite(r_geom[4], WATER_PRESSURE, (14.4 * (p[3] + p[4]) - 1.6 * (p[0] + p[1])) * c1);
    ThreadSafeNodeWrite(r_geom[5], WATER_PRESSURE, (3.0 * p[1] + p[0] + 27.0 * p[4] - 5.4 * p[3]) * c1);
    ThreadSafeNodeWrite(r_geom[6], WATER_PRESSURE, (3.0 * p[1] + p[2] + 27.0 * p[5] - 5.4 * p[6]) * c1);
    ThreadSafeNodeWrite(r_geom[7], WATER_PRESSURE, (14.4 * (p[5] + p[6]) - 1.6 * (p[1] + p[2])) * c1);
    ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE, (3.0 * p[2] + p[1] + 27.0 * p[6] - 5.4 * p[5]) * c1);
    ThreadSafeNodeWrite(r_geom[9], WATER_PRESSURE, (3.0 * p[2] + p[0] + 27.0 * p[7] - 5.4 * p[8]) * c1);
    ThreadSafeNodeWrite(r_geom[10], WATER_PRESSURE, (14.4 * (p[7] + p[8]) - 1.6 * (p[0] + p[2])) * c1);
    ThreadSafeNodeWrite(r_geom[11], WATER_PRESSURE, (3.0 * p[0] + p[2] + 27.0 * p[8] - 5.4 * p[7]) * c1);
    ThreadSafeNodeWrite(r_geom[12], WATER_PRESSURE,
                        (p[1] + p[2] + 7.2 * (p[3] + p[8]) - 3.6 * (p[4] + p[7]) -
                         1.8 * (p[5] + p[6]) + 21.6 * p[9] - 1.6 * p[0]) *
                            c1);
    ThreadSafeNodeWrite(r_geom[13], WATER_PRESSURE,
                        (p[0] + p[2] + 7.2 * (p[4] + p[5]) - 3.6 * (p[3] + p[6]) -
                         1.8 * (p[7] + p[8]) + 21.6 * p[9] - 1.6 * p[1]) *
                            c1);
    ThreadSafeNodeWrite(r_geom[14], WATER_PRESSURE,
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
    const Vector p = GetPressures(8);
    // edges -- bottom
    ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE, 0.5 * (p[0] + p[1]));
    ThreadSafeNodeWrite(r_geom[9], WATER_PRESSURE, 0.5 * (p[1] + p[2]));
    ThreadSafeNodeWrite(r_geom[10], WATER_PRESSURE, 0.5 * (p[2] + p[3]));
    ThreadSafeNodeWrite(r_geom[11], WATER_PRESSURE, 0.5 * (p[3] + p[0]));
    // edges -- middle
    ThreadSafeNodeWrite(r_geom[12], WATER_PRESSURE, 0.5 * (p[4] + p[0]));
    ThreadSafeNodeWrite(r_geom[13], WATER_PRESSURE, 0.5 * (p[5] + p[1]));
    ThreadSafeNodeWrite(r_geom[14], WATER_PRESSURE, 0.5 * (p[6] + p[2]));
    ThreadSafeNodeWrite(r_geom[15], WATER_PRESSURE, 0.5 * (p[7] + p[3]));
    // edges -- top
    ThreadSafeNodeWrite(r_geom[16], WATER_PRESSURE, 0.5 * (p[4] + p[5]));
    ThreadSafeNodeWrite(r_geom[17], WATER_PRESSURE, 0.5 * (p[5] + p[6]));
    ThreadSafeNodeWrite(r_geom[18], WATER_PRESSURE, 0.5 * (p[6] + p[7]));
    ThreadSafeNodeWrite(r_geom[19], WATER_PRESSURE, 0.5 * (p[7] + p[4]));
    KRATOS_CATCH("")
}

template <>
void SmallStrainUPwDiffOrderElement<3, 27>::AssignPressureToIntermediateNodes()
{
    // Assign pressure values to the intermediate nodes for post-processing
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();
    // 3D H27P8
    const Vector p = GetPressures(8);
    // edges -- bottom
    ThreadSafeNodeWrite(r_geom[8], WATER_PRESSURE, 0.5 * (p[0] + p[1]));
    ThreadSafeNodeWrite(r_geom[9], WATER_PRESSURE, 0.5 * (p[1] + p[2]));
    ThreadSafeNodeWrite(r_geom[10], WATER_PRESSURE, 0.5 * (p[2] + p[3]));
    ThreadSafeNodeWrite(r_geom[11], WATER_PRESSURE, 0.5 * (p[3] + p[0]));
    // edges -- middle
    ThreadSafeNodeWrite(r_geom[12], WATER_PRESSURE, 0.5 * (p[4] + p[0]));
    ThreadSafeNodeWrite(r_geom[13], WATER_PRESSURE, 0.5 * (p[5] + p[1]));
    ThreadSafeNodeWrite(r_geom[14], WATER_PRESSURE, 0.5 * (p[6] + p[2]));
    ThreadSafeNodeWrite(r_geom[15], WATER_PRESSURE, 0.5 * (p[7] + p[3]));
    // edges -- top
    ThreadSafeNodeWrite(r_geom[16], WATER_PRESSURE, 0.5 * (p[4] + p[5]));
    ThreadSafeNodeWrite(r_geom[17], WATER_PRESSURE, 0.5 * (p[5] + p[6]));
    ThreadSafeNodeWrite(r_geom[18], WATER_PRESSURE, 0.5 * (p[6] + p[7]));
    ThreadSafeNodeWrite(r_geom[19], WATER_PRESSURE, 0.5 * (p[7] + p[0]));
    // face centers
    ThreadSafeNodeWrite(r_geom[20], WATER_PRESSURE, 0.25 * (p[0] + p[1] + p[2] + p[3]));
    ThreadSafeNodeWrite(r_geom[21], WATER_PRESSURE, 0.25 * (p[0] + p[1] + p[4] + p[5]));
    ThreadSafeNodeWrite(r_geom[22], WATER_PRESSURE, 0.25 * (p[1] + p[2] + p[5] + p[6]));
    ThreadSafeNodeWrite(r_geom[23], WATER_PRESSURE, 0.25 * (p[2] + p[3] + p[6] + p[7]));
    ThreadSafeNodeWrite(r_geom[24], WATER_PRESSURE, 0.25 * (p[3] + p[0] + p[7] + p[4]));
    ThreadSafeNodeWrite(r_geom[25], WATER_PRESSURE, 0.25 * (p[4] + p[5] + p[6] + p[7]));
    // element center
    ThreadSafeNodeWrite(r_geom[26], WATER_PRESSURE,
                        0.125 * (p[0] + p[1] + p[2] + p[3] + p[4] + p[5] + p[6] + p[7]));
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
