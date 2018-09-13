// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes
//

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/line_2d.h"
#include "convection_diffusion_application.h"
#include "includes/variables.h"

namespace Kratos {



KratosConvectionDiffusionApplication::KratosConvectionDiffusionApplication()
    : KratosApplication("ConvectionDiffusionApplication"),
      mEulerianConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mEulerianConvDiff2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConvDiff3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
      mEulerianDiffusion2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mEulerianDiffusion3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      mConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      mLaplacian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mLaplacian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
      mLaplacian3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
      mLaplacian3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27)))),
      mThermalFace2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
      mThermalFace3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mFluxCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
      mFluxCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mFluxCondition3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))) {}

void KratosConvectionDiffusionApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << 
    " KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ " << std::endl <<
    "       / __/ _ || || | | / /__|   |_ _| __| __|" << std::endl <<
    "      | (_| (_) | .` || V /___| |) | || _|| _| " << std::endl <<
    "       |___|___/|_||_| |_/    |___/___|_| |_|  APPLICATION" << std::endl;

    // Registering variables
    KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_1)
    KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_2)
    KRATOS_REGISTER_VARIABLE(BFECC_ERROR)
    KRATOS_REGISTER_VARIABLE(BFECC_ERROR_1)

    KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
    KRATOS_REGISTER_VARIABLE(PROJECTED_SCALAR1)
    KRATOS_REGISTER_VARIABLE(DELTA_SCALAR1)
    KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)

    KRATOS_REGISTER_VARIABLE(THETA)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

    // Registering elements and conditions here
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff2D", mEulerianConvDiff2D);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff2D4N", mEulerianConvDiff2D4N);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff3D", mEulerianConvDiff3D);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff3D8N", mEulerianConvDiff3D8N);
    KRATOS_REGISTER_ELEMENT("EulerianDiffusion2D", mEulerianDiffusion2D);
    KRATOS_REGISTER_ELEMENT("EulerianDiffusion3D", mEulerianDiffusion3D);
    KRATOS_REGISTER_ELEMENT("ConvDiff2D", mConvDiff2D);
    KRATOS_REGISTER_ELEMENT("ConvDiff3D", mConvDiff3D);
    KRATOS_REGISTER_ELEMENT("LaplacianElement2D3N", mLaplacian2D3N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D4N", mLaplacian3D4N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D8N", mLaplacian3D8N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D27N", mLaplacian3D27N);

    KRATOS_REGISTER_CONDITION("ThermalFace2D", mThermalFace2D);
    KRATOS_REGISTER_CONDITION("ThermalFace3D", mThermalFace3D);
    KRATOS_REGISTER_CONDITION("FluxCondition2D2N", mFluxCondition2D2N);
    KRATOS_REGISTER_CONDITION("FluxCondition3D3N", mFluxCondition3D3N);
    KRATOS_REGISTER_CONDITION("FluxCondition3D4N", mFluxCondition3D4N);
}

}  // namespace Kratos.
