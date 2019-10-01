//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/variables.h"

// Application includes
#include "fluid_transport_application.h"

namespace Kratos
{

KratosFluidTransportApplication::KratosFluidTransportApplication():
    KratosApplication("FluidTransportApplication"),
    mSteadyConvectionDiffusionFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType (3) ) ) ),
    mSteadyConvectionDiffusionFICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3>>( Element::GeometryType::PointsArrayType (4) ) ) ),
    mSteadyConvectionDiffusionFICElement3D4N( 0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSteadyConvectionDiffusionFICElement3D8N( 0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),

    mTransientConvectionDiffusionFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType (3) ) ) ),
    mTransientConvectionDiffusionFICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3>>( Element::GeometryType::PointsArrayType (4) ) ) ),
    mTransientConvectionDiffusionFICElement3D4N( 0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTransientConvectionDiffusionFICElement3D8N( 0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),

    mTransientConvectionDiffusionFICExplicitElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType (3) ) ) ),
    mTransientConvectionDiffusionFICExplicitElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3>>( Element::GeometryType::PointsArrayType (4) ) ) ),
    mTransientConvectionDiffusionFICExplicitElement3D4N( 0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTransientConvectionDiffusionFICExplicitElement3D8N( 0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),

    mTransientConvectionDiffusionPFEM2FICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>>( Element::GeometryType::PointsArrayType (3) ) ) ),
    mTransientConvectionDiffusionPFEM2FICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3>>( Element::GeometryType::PointsArrayType (4) ) ) ),
    mTransientConvectionDiffusionPFEM2FICElement3D4N( 0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTransientConvectionDiffusionPFEM2FICElement3D8N( 0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8))))

    {}

void KratosFluidTransportApplication::Register()
{
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosFluidTransportApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE(PECLET);
    KRATOS_REGISTER_VARIABLE(THETA);
    KRATOS_REGISTER_VARIABLE(PHI_THETA);
    KRATOS_REGISTER_VARIABLE(NODAL_ANALYTIC_SOLUTION);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_PHI_GRADIENT);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PHI_GRADIENT);

    KRATOS_REGISTER_ELEMENT( "SteadyConvectionDiffusionFICElement2D3N", mSteadyConvectionDiffusionFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SteadyConvectionDiffusionFICElement2D4N", mSteadyConvectionDiffusionFICElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SteadyConvectionDiffusionFICElement3D4N", mSteadyConvectionDiffusionFICElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SteadyConvectionDiffusionFICElement3D8N", mSteadyConvectionDiffusionFICElement3D8N )

    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionFICElement2D3N", mTransientConvectionDiffusionFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionFICElement2D4N", mTransientConvectionDiffusionFICElement2D4N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionFICElement3D4N", mTransientConvectionDiffusionFICElement3D4N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionFICElement3D8N", mTransientConvectionDiffusionFICElement3D8N )

    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionFICExplicitElement2D3N", mTransientConvectionDiffusionFICExplicitElement2D3N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionFICExplicitElement2D4N", mTransientConvectionDiffusionFICExplicitElement2D4N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionFICExplicitElement3D4N", mTransientConvectionDiffusionFICExplicitElement3D4N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionFICExplicitElement3D8N", mTransientConvectionDiffusionFICExplicitElement3D8N )

    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionPFEM2FICElement2D3N", mTransientConvectionDiffusionPFEM2FICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionPFEM2FICElement2D4N", mTransientConvectionDiffusionPFEM2FICElement2D4N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionPFEM2FICElement3D4N", mTransientConvectionDiffusionPFEM2FICElement3D4N )
    KRATOS_REGISTER_ELEMENT( "TransientConvectionDiffusionPFEM2FICElement3D8N", mTransientConvectionDiffusionPFEM2FICElement3D8N )
}

}  // namespace Kratos.
