//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
//



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d_2.h"
#include "thermo_mechanical_application.h"
#include "includes/variables.h"


namespace Kratos
{

KratosThermoMechanicalApplication::KratosThermoMechanicalApplication():
    KratosApplication("ThermoMechanicalApplication"),
// 		mElem2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
/*		mMonolithic2DNeumann(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),*/
    mHeatContact2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mHeatContact3D(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mThermalFace2D(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mThermalFace3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mEnvironmentContact(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mSUPGConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mSUPGConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSUPGConv3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSUPGConv2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
	mSUPGConvLevelSet(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
   // mPoisson3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
{}


void KratosThermoMechanicalApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosThermoMechanicalApplication... " << std::endl;

    // KRATOS_REGISTER_ELEMENT("Elem2D", mElem2D);
    // KRATOS_REGISTER_ELEMENT("Elemt3D", mElem3D);
    KRATOS_REGISTER_CONDITION("HeatContact2D", mHeatContact2D);
    KRATOS_REGISTER_CONDITION("HeatContact3D", mHeatContact3D);
    KRATOS_REGISTER_CONDITION("ThermalFace2D", mThermalFace2D);
    KRATOS_REGISTER_CONDITION("ThermalFace3D", mThermalFace3D);
    KRATOS_REGISTER_CONDITION("EnvironmentContact", mEnvironmentContact);

    KRATOS_REGISTER_ELEMENT("SUPGConvDiff2D", mSUPGConvDiff2D);
    KRATOS_REGISTER_ELEMENT("SUPGConvDiff3D", mSUPGConvDiff3D);
    KRATOS_REGISTER_ELEMENT("SUPGConv3D", mSUPGConv3D);
    KRATOS_REGISTER_ELEMENT("SUPGConv2D", mSUPGConv2D);
    KRATOS_REGISTER_ELEMENT("SUPGConv3D", mSUPGConv3D);
    // KRATOS_REGISTER_ELEMENT("Poisson3D", mPoisson3D);
	KRATOS_REGISTER_ELEMENT("SUPGConvLevelSet", mSUPGConvLevelSet);
}

}  // namespace Kratos.


