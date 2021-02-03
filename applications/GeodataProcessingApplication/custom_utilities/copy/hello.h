//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//                   Simon Wenczowski
//

#if !defined(KRATOS_HELLO_H )
#define  KRATOS_HELLO_H

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"

#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// A set of functions to compute quantities of interest for turbulent flows.
class Hello {

public:

Hello(){};

void Greet();

void test_nicola();

ModelPart& CheckIfInternal(ModelPart& VolumeModelPart, ModelPart& GeometryModelPart);

};  // Class Hello

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HELLO_H  defined
