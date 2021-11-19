//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/kernel.h"

#include "includes/kratos_flags.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE(GlobalPointersVector<Node<3> >, NEIGHBOUR_NODES)
    KRATOS_CREATE_VARIABLE(GlobalPointersVector<Node<3> >, NEIGHBOUR_CONDITION_NODES)
    KRATOS_CREATE_VARIABLE(GlobalPointersVector<Node<3> >, FATHER_NODES)

  void KratosApplication::RegisterGlobalPointerVariables()
  {
        KRATOS_REGISTER_VARIABLE(  NEIGHBOUR_NODES )
        KRATOS_REGISTER_VARIABLE(  NEIGHBOUR_CONDITION_NODES )
        KRATOS_REGISTER_VARIABLE(  FATHER_NODES )
  }


}  // namespace Kratos.
