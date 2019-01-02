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


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/fsi_variables.h"
#include "python/add_fsi_variables_to_python.h"



namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void  AddFSIVariablesToPython(pybind11::module& m)
{

  //registering variables in python
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,CONVERGENCE_ACCELERATOR_ITERATION);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,MAPPER_SCALAR_PROJECTION_RHS);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FICTITIOUS_FLUID_DENSITY);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SCALAR_PROJECTED);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FSI_INTERFACE_RESIDUAL_NORM);
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,FSI_INTERFACE_MESH_RESIDUAL_NORM);

  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,MAPPER_VECTOR_PROJECTION_RHS);
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,VAUX_EQ_TRACTION);
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,VECTOR_PROJECTED);
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,RELAXED_DISP);
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,FSI_INTERFACE_RESIDUAL);
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,FSI_INTERFACE_MESH_RESIDUAL);
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,POSITIVE_MAPPED_VECTOR_VARIABLE);
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,NEGATIVE_MAPPED_VECTOR_VARIABLE);

}
}  // namespace Python.
} // Namespace Kratos

