//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


// System includes


// External includes


// Project includes
#include "stabilized_cfd_application.h"
#include "stabilized_cfd_application_variables.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_2.h"

namespace Kratos {

KratosStabilizedCFDApplication::KratosStabilizedCFDApplication():
    KratosApplication("StabilizedCFDApplication"),
    mDSS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDSS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mDSS2D_FIC(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDSS3D_FIC(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS2D4N_FIC(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS3D8N_FIC(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mDSS2D_FIC_LIMITED(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDSS3D_FIC_LIMITED(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS2D4N_FIC_LIMITED(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS3D8N_FIC_LIMITED(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mDSS2D_GLS(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDSS3D_GLS(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS2D4N_GLS(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS3D8N_GLS(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mDynSS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDynSS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDynSS2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDynSS3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mDSS2D_notau2(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDSS3D_notau2(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS2D4N_notau2(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS3D8N_notau2(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mDYNSS2D_NOTAU2(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDYNSS3D_NOTAU2(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDYNSS2D4N_NOTAU2(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDYNSS3D8N_NOTAU2(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mDSS2D_PS(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDSS3D_PS(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS2D4N_PS(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDSS3D8N_PS(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mDSSFace2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mDSSFace3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mDSSFace3D4N(0, Element::GeometryType::Pointer( new Quadrilateral3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
    {}

void KratosStabilizedCFDApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosStabilizedCFDApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE( FIC_BETA )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DIRECTIONAL_BETA )
    KRATOS_REGISTER_VARIABLE( RECORDED_STEPS )
    KRATOS_REGISTER_VARIABLE( MEAN_KINETIC_ENERGY )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MEAN_VELOCITY )
    KRATOS_REGISTER_VARIABLE( MEAN_PRESSURE )
    KRATOS_REGISTER_VARIABLE( VELOCITY_COVARIANCES )
    KRATOS_REGISTER_VARIABLE( TURBULENCE_STATISTICS )
    KRATOS_REGISTER_VARIABLE( TRACE_XI )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DIV_XI )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_PROJECTION )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MOMENTUM_PROJECTION_RHS )
    KRATOS_REGISTER_VARIABLE( MASS_PROJECTION )
    KRATOS_REGISTER_VARIABLE( MASS_PROJECTION_RHS )

    KRATOS_REGISTER_ELEMENT("DSS2D",mDSS2D);
    KRATOS_REGISTER_ELEMENT("DSS3D",mDSS3D);
    KRATOS_REGISTER_ELEMENT("DSS2D4N",mDSS2D4N);
    KRATOS_REGISTER_ELEMENT("DSS3D8N",mDSS3D8N);

    KRATOS_REGISTER_ELEMENT("DSS2D_FIC",mDSS2D_FIC);
    KRATOS_REGISTER_ELEMENT("DSS3D_FIC",mDSS3D_FIC);
    KRATOS_REGISTER_ELEMENT("DSS2D4N_FIC",mDSS2D4N_FIC);
    KRATOS_REGISTER_ELEMENT("DSS3D8N_FIC",mDSS3D8N_FIC);

    KRATOS_REGISTER_ELEMENT("DSS2D_FIC_LIMITED",mDSS2D_FIC_LIMITED);
    KRATOS_REGISTER_ELEMENT("DSS3D_FIC_LIMITED",mDSS3D_FIC_LIMITED);
    KRATOS_REGISTER_ELEMENT("DSS2D4N_FIC_LIMITED",mDSS2D4N_FIC_LIMITED);
    KRATOS_REGISTER_ELEMENT("DSS3D8N_FIC_LIMITED",mDSS3D8N_FIC_LIMITED);

    KRATOS_REGISTER_ELEMENT("DSS2D_GLS",mDSS2D_GLS);
    KRATOS_REGISTER_ELEMENT("DSS3D_GLS",mDSS3D_GLS);
    KRATOS_REGISTER_ELEMENT("DSS2D4N_GLS",mDSS2D4N_GLS);
    KRATOS_REGISTER_ELEMENT("DSS3D8N_GLS",mDSS3D8N_GLS);

    KRATOS_REGISTER_ELEMENT("DynSS2D",mDynSS2D);
    KRATOS_REGISTER_ELEMENT("DynSS3D",mDynSS3D);
    KRATOS_REGISTER_ELEMENT("DynSS2D4N",mDynSS2D4N);
    KRATOS_REGISTER_ELEMENT("DynSS3D8N",mDynSS3D8N);

    KRATOS_REGISTER_ELEMENT("DSS2D_notau2",mDSS2D_notau2);
    KRATOS_REGISTER_ELEMENT("DSS3D_notau2",mDSS3D_notau2);
    KRATOS_REGISTER_ELEMENT("DSS2D4N_notau2",mDSS2D4N_notau2);
    KRATOS_REGISTER_ELEMENT("DSS3D8N_notau2",mDSS3D8N_notau2);

    KRATOS_REGISTER_ELEMENT("DYNSS2D_NOTAU2",mDYNSS2D_NOTAU2);
    KRATOS_REGISTER_ELEMENT("DYNSS3D_NOTAU2",mDYNSS3D_NOTAU2);
    KRATOS_REGISTER_ELEMENT("DYNSS2D4N_NOTAU2",mDYNSS2D4N_NOTAU2);
    KRATOS_REGISTER_ELEMENT("DYNSS3D8N_NOTAU2",mDYNSS3D8N_NOTAU2);

    KRATOS_REGISTER_ELEMENT("DSS_PS2D",mDSS2D_PS);
    KRATOS_REGISTER_ELEMENT("DSS_PS3D",mDSS3D_PS);
    KRATOS_REGISTER_ELEMENT("DSS_PS2D4N",mDSS2D4N_PS);
    KRATOS_REGISTER_ELEMENT("DSS_PS3D8N",mDSS3D8N_PS);

    KRATOS_REGISTER_CONDITION("DSSFace2D",mDSSFace2D);
    KRATOS_REGISTER_CONDITION("DSSFace3D",mDSSFace3D);
    KRATOS_REGISTER_CONDITION("DSSFace3D4N",mDSSFace3D4N);

}
}  // namespace Kratos.
