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

#if !defined(KRATOS_MESHING_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_MESHING_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "includes/condition.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef array_1d<double, 3> Vector3;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

    // Variables definition
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, double, AVERAGE_NODAL_ERROR);                          // The average nodal error
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, double, ANISOTROPIC_RATIO);                            // The anisotropic aspect ratio
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, Vector3, AUXILIAR_GRADIENT);                           // An auxiliar gradient needed to compute the metric
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, Vector,  AUXILIAR_HESSIAN);                            // An auxiliar hessian needed to compute the metric
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, double, METRIC_SCALAR);                                // A single scalar metric
    KRATOS_DEFINE_SYMMETRIC_2D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(MESHING_APPLICATION, METRIC_TENSOR_2D); // A 2D metric vector
    KRATOS_DEFINE_SYMMETRIC_3D_TENSOR_APPLICATION_VARIABLE_WITH_COMPONENTS(MESHING_APPLICATION, METRIC_TENSOR_3D); // A 3D metric vector

    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, int, NUMBER_OF_DIVISIONS);      // The number of divisions for the multi scale refining
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, int, SUBSCALE_INDEX)
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, Node<3>::Pointer, SLAVE_NODE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, Element::Pointer, FATHER_ELEMENT)
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, Condition::Pointer, FATHER_CONDITION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, std::vector<double>, FATHER_NODES_WEIGHTS)

    //for ULF (surface_tension) application:
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, double, TRIPLE_POINT)
    KRATOS_DEFINE_APPLICATION_VARIABLE(MESHING_APPLICATION, double, CONTACT_ANGLE)
}

#endif /* KRATOS_MESHING_APPLICATION_VARIABLES_H_INCLUDED */
