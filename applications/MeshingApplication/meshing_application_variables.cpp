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
#include "meshing_application_variables.h"

namespace Kratos
{
typedef array_1d<double, 3> Vector3;

KRATOS_CREATE_VARIABLE(double, AVERAGE_NODAL_ERROR);                          // The average nodal error
KRATOS_CREATE_VARIABLE(double, ANISOTROPIC_RATIO);                            // The anisotropic aspect ratio
KRATOS_CREATE_VARIABLE(Vector3, AUXILIAR_GRADIENT);                           // An auxiliar gradient needed to compute the metric
KRATOS_CREATE_VARIABLE(Vector, AUXILIAR_HESSIAN);                             // An auxiliar hessian needed to compute the metric
KRATOS_CREATE_VARIABLE(double, METRIC_SCALAR);                                // A single scalar metric
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
}
