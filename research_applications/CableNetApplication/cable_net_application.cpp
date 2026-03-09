//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//


// System includes


// External includes


// Project includes
#include "cable_net_application.h"
#include "cable_net_application_variables.h"


#include "geometries/triangle_3d_3.h"
#include "geometries/line_3d_2.h"
#include "custom_geometries/line_3d_n.h"

namespace Kratos {

    typedef Node NodeType;

KratosCableNetApplication::KratosCableNetApplication():
    KratosApplication("CableNetApplication"),
    mWeakSlidingElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
    mSlidingCableElement3D3N(0, Element::GeometryType::Pointer(new Line3DN<NodeType >(Element::GeometryType::PointsArrayType(3)))),
    mRingElement3D4N(0, Element::GeometryType::Pointer(new Line3DN<NodeType >(Element::GeometryType::PointsArrayType(4)))),
    mRingElement3D3N(0, Element::GeometryType::Pointer(new Line3DN<NodeType >(Element::GeometryType::PointsArrayType(3)))),
    mEmpiricalSpringElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2))))
    {}

void KratosCableNetApplication::Register()
{
    KRATOS_INFO("")  <<  "    KRATOS    ___      _     _          __     _\n"
                     <<  "             / __\\__ _| |__ | | ___  /\\ \\ \\___| |_\n"
                     <<  "            / /  / _` | '_ \\| |/ _ \\/  \\/ / _ \\ __|\n"
                     <<  "           / /__| (_| | |_) | |  __/ /\\  /  __/ |_\n"
                     <<  "           \\____/\\__,_|_.__/|_|\\___\\_\\ \\/ \\___|\\__| Application\n"
                     <<    "Initializing KratosCableNetApplication..." << std::endl;



    KRATOS_REGISTER_ELEMENT("WeakSlidingElement3D3N", mWeakSlidingElement3D3N)
    KRATOS_REGISTER_ELEMENT("SlidingCableElement3D3N", mSlidingCableElement3D3N)
    KRATOS_REGISTER_ELEMENT("RingElement3D4N", mRingElement3D4N)
    KRATOS_REGISTER_ELEMENT("RingElement3D3N", mRingElement3D3N)
    KRATOS_REGISTER_ELEMENT("EmpiricalSpringElement3D2N", mEmpiricalSpringElement3D2N)

    KRATOS_REGISTER_VARIABLE(SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL)
}
}  // namespace Kratos.
