//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "cable_net_application.h"
#include "cable_net_application_variables.h"


#include "geometries/triangle_3d_3.h"
#include "custom_geometries/line_3d_n.h"

namespace Kratos {

    typedef Node<3> NodeType;

KratosCableNetApplication::KratosCableNetApplication():
    KratosApplication("CableNetApplication"),
    mWeakSlidingElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
    mSlidingCableElement3D3N(0, Element::GeometryType::Pointer(new Line3DN<NodeType >(Element::GeometryType::PointsArrayType(3)))),
    mRingElement3D4N(0, Element::GeometryType::Pointer(new Line3DN<NodeType >(Element::GeometryType::PointsArrayType(4)))),
    mRingElement3D3N(0, Element::GeometryType::Pointer(new Line3DN<NodeType >(Element::GeometryType::PointsArrayType(3))))
    {}

void KratosCableNetApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();

    KRATOS_INFO("") <<    "KRATOS  ___/\\/\\/\\/\\/\\______________/\\/\\________/\\/\\________________/\\/\\____/\\/\\________________/\\/\\_____\n"
                    <<    "       _/\\/\\__________/\\/\\/\\______/\\/\\________/\\/\\______/\\/\\/\\____/\\/\\/\\__/\\/\\____/\\/\\/\\____/\\/\\/\\/\\/\\_\n"
                    <<    "      _/\\/\\______________/\\/\\____/\\/\\/\\/\\____/\\/\\____/\\/\\/\\/\\/\\__/\\/\\/\\/\\/\\/\\__/\\/\\/\\/\\/\\____/\\/\\_____\n"
                    <<    "     _/\\/\\__________/\\/\\/\\/\\____/\\/\\__/\\/\\__/\\/\\____/\\/\\________/\\/\\__/\\/\\/\\__/\\/\\__________/\\/\\_____\n"
                    <<    "    ___/\\/\\/\\/\\/\\__/\\/\\/\\/\\/\\__/\\/\\/\\/\\____/\\/\\/\\____/\\/\\/\\/\\__/\\/\\____/\\/\\____/\\/\\/\\/\\____/\\/\\/\\___\n"
                    <<    "   ________________________________________________________________________________________________ Application\n"
                    <<    "Initializing KratosCableNetApplication..." << std::endl;

    KRATOS_REGISTER_ELEMENT("WeakSlidingElement3D3N", mWeakSlidingElement3D3N)
    KRATOS_REGISTER_ELEMENT("SlidingCableElement3D3N", mSlidingCableElement3D3N)
    KRATOS_REGISTER_ELEMENT("RingElement3D4N", mRingElement3D4N)
    KRATOS_REGISTER_ELEMENT("RingElement3D3N", mRingElement3D3N)
}
}  // namespace Kratos.
