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
#include "convection_diffusion_application.h"
#include "includes/variables.h"

namespace Kratos {



KratosConvectionDiffusionApplication::KratosConvectionDiffusionApplication()
    : KratosApplication("ConvectionDiffusionApplication"),


      mAxisymmetricEulerianConvectionDiffusion2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mAxisymmetricEulerianConvectionDiffusion2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConvDiff2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mEulerianConvDiff2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConvDiff3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConvDiff3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node >(Element::GeometryType::PointsArrayType(8)))),
      mEulerianDiffusion2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mEulerianDiffusion3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mLaplacian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mLaplacian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mLaplacian3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node >(Element::GeometryType::PointsArrayType(8)))),
      mLaplacian3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node >(Element::GeometryType::PointsArrayType(27)))),
      mLaplacianShiftedBoundary2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
      mLaplacianShiftedBoundary3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
      mEmbeddedLaplacian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mEmbeddedLaplacian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mMixedLaplacianElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mMixedLaplacianElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mMixedLaplacianShiftedBoundary2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
      mMixedLaplacianShiftedBoundary3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
      mAdjointDiffusionElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mAdjointDiffusionElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mAxisymmetricThermalFace2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mThermalFace2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mThermalFace3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mThermalFace3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mFluxCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mFluxCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mFluxCondition3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mLaplacianShiftedBoundaryCondition(0, Element::GeometryType::Pointer(new Geometry<Node>())),
      mMixedLaplacianShiftedBoundaryCondition(0, Element::GeometryType::Pointer(new Geometry<Node>())),
      mAdjointThermalFace2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
      mAdjointThermalFace3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mQSConvectionDiffusionExplicit2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mQSConvectionDiffusionExplicit3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mDConvectionDiffusionExplicit2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mDConvectionDiffusionExplicit3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))) {}

void KratosConvectionDiffusionApplication::Register() {
    KRATOS_INFO("") <<
    " KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ " << std::endl <<
    "       / __/ _ || || | | / /__|   |_ _| __| __|" << std::endl <<
    "      | (_| (_) | .` || V /___| |) | || _|| _| " << std::endl <<
    "       |___|___/|_||_| |_/    |___/___|_| |_|  APPLICATION" << std::endl;

    // Registering variables
    KRATOS_REGISTER_VARIABLE(AUX_FLUX)
    KRATOS_REGISTER_VARIABLE(AUX_TEMPERATURE)
    KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_1)
    KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_2)
    KRATOS_REGISTER_VARIABLE(BFECC_ERROR)
    KRATOS_REGISTER_VARIABLE(BFECC_ERROR_1)

    KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
    KRATOS_REGISTER_VARIABLE(PROJECTED_SCALAR1)
    KRATOS_REGISTER_VARIABLE(DELTA_SCALAR1)
    KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)

    KRATOS_REGISTER_VARIABLE(TRANSFER_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(ADJOINT_HEAT_TRANSFER)
    KRATOS_REGISTER_VARIABLE(SCALAR_PROJECTION)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

    // Registering elements and conditions here
    KRATOS_REGISTER_ELEMENT("AxisymmetricEulerianConvectionDiffusion2D3N", mAxisymmetricEulerianConvectionDiffusion2D3N);
    KRATOS_REGISTER_ELEMENT("AxisymmetricEulerianConvectionDiffusion2D4N", mAxisymmetricEulerianConvectionDiffusion2D4N);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff2D", mEulerianConvDiff2D3N); //TODO: To be removed as it does not follow the naming convention
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff2D3N", mEulerianConvDiff2D3N);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff2D4N", mEulerianConvDiff2D4N);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff3D", mEulerianConvDiff3D4N); //TODO: To be removed as it does not follow the naming convention
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff3D4N", mEulerianConvDiff3D4N);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff3D8N", mEulerianConvDiff3D8N);
    KRATOS_REGISTER_ELEMENT("EulerianDiffusion2D3N", mEulerianDiffusion2D3N);
    KRATOS_REGISTER_ELEMENT("EulerianDiffusion3D4N", mEulerianDiffusion3D4N);
    KRATOS_REGISTER_ELEMENT("ConvDiff2D", mConvDiff2D);
    KRATOS_REGISTER_ELEMENT("ConvDiff3D", mConvDiff3D);
    KRATOS_REGISTER_ELEMENT("LaplacianElement2D3N", mLaplacian2D3N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D4N", mLaplacian3D4N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D8N", mLaplacian3D8N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D27N", mLaplacian3D27N);
    KRATOS_REGISTER_ELEMENT("LaplacianShiftedBoundaryElement2D3N", mLaplacianShiftedBoundary2D3N);
    KRATOS_REGISTER_ELEMENT("LaplacianShiftedBoundaryElement3D4N", mLaplacianShiftedBoundary3D4N);
    KRATOS_REGISTER_ELEMENT("MixedLaplacianElement2D3N", mMixedLaplacianElement2D3N);
    KRATOS_REGISTER_ELEMENT("MixedLaplacianElement3D4N", mMixedLaplacianElement3D4N);
    KRATOS_REGISTER_ELEMENT("MixedLaplacianShiftedBoundaryElement2D3N", mMixedLaplacianShiftedBoundary2D3N);
    KRATOS_REGISTER_ELEMENT("MixedLaplacianShiftedBoundaryElement3D4N", mMixedLaplacianShiftedBoundary3D4N);
    KRATOS_REGISTER_ELEMENT("EmbeddedLaplacianElement2D3N", mEmbeddedLaplacian2D3N);
    KRATOS_REGISTER_ELEMENT("EmbeddedLaplacianElement3D4N", mEmbeddedLaplacian3D4N);
    KRATOS_REGISTER_ELEMENT("QSConvectionDiffusionExplicit2D3N", mQSConvectionDiffusionExplicit2D3N);
    KRATOS_REGISTER_ELEMENT("QSConvectionDiffusionExplicit3D4N", mQSConvectionDiffusionExplicit3D4N);
    KRATOS_REGISTER_ELEMENT("DConvectionDiffusionExplicit2D3N", mDConvectionDiffusionExplicit2D3N);
    KRATOS_REGISTER_ELEMENT("DConvectionDiffusionExplicit3D4N", mDConvectionDiffusionExplicit3D4N);

    KRATOS_REGISTER_ELEMENT("AdjointDiffusionElement2D3N", mAdjointDiffusionElement2D3N);
    KRATOS_REGISTER_ELEMENT("AdjointDiffusionElement3D4N", mAdjointDiffusionElement3D4N);

    KRATOS_REGISTER_CONDITION("AxisymmetricThermalFace2D2N", mAxisymmetricThermalFace2D2N);
    KRATOS_REGISTER_CONDITION("ThermalFace2D2N", mThermalFace2D2N);
    KRATOS_REGISTER_CONDITION("ThermalFace3D3N", mThermalFace3D3N);
    KRATOS_REGISTER_CONDITION("ThermalFace3D4N", mThermalFace3D4N);
    KRATOS_REGISTER_CONDITION("FluxCondition2D2N", mFluxCondition2D2N);
    KRATOS_REGISTER_CONDITION("FluxCondition3D3N", mFluxCondition3D3N);
    KRATOS_REGISTER_CONDITION("FluxCondition3D4N", mFluxCondition3D4N);
    KRATOS_REGISTER_CONDITION("LaplacianShiftedBoundaryCondition", mLaplacianShiftedBoundaryCondition);
    KRATOS_REGISTER_CONDITION("MixedLaplacianShiftedBoundaryCondition", mMixedLaplacianShiftedBoundaryCondition);

    KRATOS_REGISTER_CONDITION("AdjointThermalFace2D2N", mAdjointThermalFace2D2N);
    KRATOS_REGISTER_CONDITION("AdjointThermalFace3D3N", mAdjointThermalFace3D3N);
}

}  // namespace Kratos.
