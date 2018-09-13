/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

// System includes

// External includes

// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"

namespace Kratos {

KratosIgaApplication::KratosIgaApplication()
    : KratosApplication("IgaApplication")
    , mTrussDiscreteElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
    , mShellKLDiscreteElement(0, Element::GeometryType::Pointer(
        new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1))))
{
}

void KratosIgaApplication::Register() {
    KratosApplication::Register();
    std::cout << "Initializing KratosIgaApplication... " << std::endl;

    KRATOS_REGISTER_ELEMENT("TrussDiscreteElement", mTrussDiscreteElement)
    KRATOS_REGISTER_ELEMENT("ShellKLDiscreteElement", mShellKLDiscreteElement)
}

}  // namespace Kratos
