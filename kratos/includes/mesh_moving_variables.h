//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    msandre
//

#if !defined(KRATOS_ALE_VARIABLES_H_INCLUDED)
#define KRATOS_ALE_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "includes/kratos_components.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
#include "containers/weak_pointer_vector.h"
#include "containers/periodic_variables_container.h"

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

namespace Kratos
{

   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT);
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_ACCELERATION);
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION);
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS);

   KRATOS_DEFINE_VARIABLE(int, LAPLACIAN_DIRECTION);

} // namespace Kratos

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_ALE_VARIABLES_H_INCLUDED defined
