// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"

#include "meshing_application.h"
#include "includes/variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d_2.h"

namespace Kratos {

typedef array_1d<double, 3> Vector3;

KRATOS_CREATE_VARIABLE(double, AVERAGE_NODAL_ERROR);                          // The average nodal error
KRATOS_CREATE_VARIABLE(double, ANISOTROPIC_RATIO);                            // The anisotropic aspect ratio
KRATOS_CREATE_VARIABLE(Vector3, AUXILIAR_GRADIENT);                           // An auxiliar gradient needed to compute the metric
KRATOS_CREATE_VARIABLE(Vector, AUXILIAR_HESSIAN);                             // An auxiliar hessian needed to compute the metric
KRATOS_CREATE_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(METRIC_TENSOR_2D); // A 2D metric vector
KRATOS_CREATE_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(METRIC_TENSOR_3D); // A 3D metric vector

KRATOS_CREATE_VARIABLE(int, NUMBER_OF_DIVISIONS)
KRATOS_CREATE_VARIABLE(int, SUBSCALE_INDEX)
KRATOS_CREATE_VARIABLE(Node<3>::Pointer, SLAVE_NODE)
KRATOS_CREATE_VARIABLE(Element::Pointer, FATHER_ELEMENT)
KRATOS_CREATE_VARIABLE(Condition::Pointer, FATHER_CONDITION)
KRATOS_CREATE_VARIABLE(std::vector<double>, FATHER_NODES_WEIGHTS)

//for ULF (surface_tension) application:
KRATOS_CREATE_VARIABLE(double, TRIPLE_POINT)
KRATOS_CREATE_VARIABLE(double, CONTACT_ANGLE)

KratosMeshingApplication::KratosMeshingApplication()
    : KratosApplication("MeshingApplication"),
      mTestElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mTestElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))) {}

void KratosMeshingApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing Kratos MeshingApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE(AVERAGE_NODAL_ERROR);                                  // The average nodal error
    KRATOS_REGISTER_VARIABLE(ANISOTROPIC_RATIO);                                    // The anisotropic aspect ratio
    KRATOS_REGISTER_VARIABLE(AUXILIAR_GRADIENT);                                    // An auxiliar gradient needed to compute the metric
    KRATOS_REGISTER_VARIABLE(AUXILIAR_HESSIAN);                                     // An auxiliar hessian needed to compute the metric
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(METRIC_TENSOR_2D); // A 2D metric vector
    KRATOS_REGISTER_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(METRIC_TENSOR_3D); // A 3D metric vector

    KRATOS_REGISTER_VARIABLE(NUMBER_OF_DIVISIONS)
    KRATOS_REGISTER_VARIABLE(SUBSCALE_INDEX)
    KRATOS_REGISTER_VARIABLE(SLAVE_NODE)
    KRATOS_REGISTER_VARIABLE(FATHER_ELEMENT)
    KRATOS_REGISTER_VARIABLE(FATHER_CONDITION)
    KRATOS_REGISTER_VARIABLE(FATHER_NODES_WEIGHTS)

    //--------------- ULF Application (surface_tension) -------------------//
    KRATOS_REGISTER_VARIABLE(TRIPLE_POINT)
    KRATOS_REGISTER_VARIABLE(CONTACT_ANGLE)


    KRATOS_REGISTER_ELEMENT("TestElement2D", mTestElement2D);
    KRATOS_REGISTER_ELEMENT("TestElement3D", mTestElement3D);
}

}  // namespace Kratos.
