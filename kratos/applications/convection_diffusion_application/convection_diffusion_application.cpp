//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-12-15 15:41:36 $
//   Revision:            $Revision: 1.6 $
//
//



// System includes


// External includes
//

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/line_2d.h"
#include "convection_diffusion_application.h"
#include "includes/variables.h"


namespace Kratos
{

KRATOS_CREATE_VARIABLE(double,  MELT_TEMPERATURE_1)
KRATOS_CREATE_VARIABLE(double,  MELT_TEMPERATURE_2)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

KratosConvectionDiffusionApplication::KratosConvectionDiffusionApplication():
    mEulerianConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mEulerianConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),   
    mConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),    
    mLaplacian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mLaplacian3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
    mLaplacian3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27, Node<3>())))),
    mThermalFace2D(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
    mThermalFace3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
{}



void KratosConvectionDiffusionApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosConvectionDiffusionApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_1)
    KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_2)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

    // Registering elements and conditions here
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff2D", mEulerianConvDiff2D);
    KRATOS_REGISTER_ELEMENT("EulerianConvDiff3D", mEulerianConvDiff3D);
    KRATOS_REGISTER_ELEMENT("ConvDiff2D", mConvDiff2D);
    KRATOS_REGISTER_ELEMENT("ConvDiff3D", mConvDiff3D);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D4N", mLaplacian3D4N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D8N", mLaplacian3D8N);
    KRATOS_REGISTER_ELEMENT("LaplacianElement3D27N", mLaplacian3D27N);

    KRATOS_REGISTER_CONDITION("ThermalFace2D", mThermalFace2D);
    KRATOS_REGISTER_CONDITION("ThermalFace3D", mThermalFace3D);

}

}  // namespace Kratos.


