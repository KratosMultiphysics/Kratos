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
void SmallStrainUPwDiffOrderElement<2, 6>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    //    case 6: // 2D T6P3
    mpPressureGeometry = make_shared<Triangle2D3<Node>>(r_geometry(0), r_geometry(1), r_geometry(2));
}

void SmallStrainUPwDiffOrderElement<2, 8>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // case 8: // 2D Q8P4
    mpPressureGeometry =
        make_shared<Quadrilateral2D4<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3));
}

void SmallStrainUPwDiffOrderElement<2, 9>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // case 9: // 2D Q9P4
    mpPressureGeometry =
        make_shared<Quadrilateral2D4<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3));
}

void SmallStrainUPwDiffOrderElement<3, 10>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // case 10:
    //     if (dimension == 3) // 3D T10P4
    mpPressureGeometry =
        make_shared<Tetrahedra3D4<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3));
}

void SmallStrainUPwDiffOrderElement<2, 10>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    //       else if (dimension == 2) // 2D T10P6
    mpPressureGeometry = make_shared<Triangle2D6<Node>>(
        r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3), r_geometry(4), r_geometry(5));
}

void SmallStrainUPwDiffOrderElement<2, 15>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // case 15: // 2D T15P10
    mpPressureGeometry = make_shared<Triangle2D10<Node>>(
        r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3), r_geometry(4), r_geometry(5),
        r_geometry(6), r_geometry(7), r_geometry(8), r_geometry(9));
}

void SmallStrainUPwDiffOrderElement<3, 20>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // case 20: // 3D H20P8
    mpPressureGeometry =
        make_shared<Hexahedra3D8<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3),
                                        r_geometry(4), r_geometry(5), r_geometry(6), r_geometry(7));
}

void SmallStrainUPwDiffOrderElement<3, 27>::SetUpPressureGeometryPointer()
{
    const auto& r_geometry = GetGeometry();
    // case 27: // 3D H27P8
    mpPressureGeometry =
        make_shared<Hexahedra3D8<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3),
                                        r_geometry(4), r_geometry(5), r_geometry(6), r_geometry(7));
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
