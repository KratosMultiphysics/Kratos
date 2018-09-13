//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Jordi Cotela
//                  Ruben Zorrilla
//

#if !defined(KRATOS_FSI_VARIABLES_H_INCLUDED )
#define  KRATOS_FSI_VARIABLES_H_INCLUDED



// System includes
#include <string>
#include <iostream>

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

//TODO: move to the FSI application
namespace Kratos
{
  // Variables definition
  KRATOS_DEFINE_VARIABLE(int, CONVERGENCE_ACCELERATOR_ITERATION)
  KRATOS_DEFINE_VARIABLE(double, MAPPER_SCALAR_PROJECTION_RHS)
  KRATOS_DEFINE_VARIABLE(double, SCALAR_PROJECTED)
  KRATOS_DEFINE_VARIABLE(double, FICTITIOUS_FLUID_DENSITY)
  KRATOS_DEFINE_VARIABLE(double, FSI_INTERFACE_RESIDUAL_NORM)
  KRATOS_DEFINE_VARIABLE(double, FSI_INTERFACE_MESH_RESIDUAL_NORM)

  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAPPER_VECTOR_PROJECTION_RHS)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VAUX_EQ_TRACTION)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VECTOR_PROJECTED)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(RELAXED_DISP)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(FSI_INTERFACE_RESIDUAL)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(FSI_INTERFACE_MESH_RESIDUAL)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(POSITIVE_MAPPED_VECTOR_VARIABLE)
  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(NEGATIVE_MAPPED_VECTOR_VARIABLE)

} // namespace Kratos

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT

#endif // KRATOS_FSI_VARIABLES_H_INCLUDED  defined
