// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#if defined(KRATOS_PYTHON)

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <pybind11/pybind11.h>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define_python.h"
#include "shape_optimization_application.h"
#include "custom_python/add_custom_utilities_to_python.h"

// ==============================================================================

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosShapeOptimizationApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosShapeOptimizationApplication,
        KratosShapeOptimizationApplication::Pointer,
        KratosApplication >(m, "KratosShapeOptimizationApplication")
        .def(py::init<>())
        ;

    AddCustomUtilitiesToPython(m);

    //registering variables in python

    // Geometry variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, NORMALIZED_SURFACE_NORMAL);

    // Optimization variables
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DF1DX);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DF1DX_MAPPED);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC1DX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC2DX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC3DX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC4DX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC5DX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC6DX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC7DX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC8DX);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC9DX);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC1DX_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC2DX_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC3DX_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC4DX_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC5DX_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC6DX_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC7DX_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC8DX_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DC9DX_MAPPED);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SEARCH_DIRECTION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONTROL_POINT_UPDATE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONTROL_POINT_CHANGE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHAPE_UPDATE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SHAPE_CHANGE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MESH_CHANGE);

    // For edge damping
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DAMPING_FACTOR);

    // For mapping
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MAPPING_ID);

    // Bead optimization
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ALPHA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ALPHA_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DF1DALPHA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DF1DALPHA_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DPDALPHA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DPDALPHA_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DLDALPHA);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, BEAD_DIRECTION);

    // For auxiliary operations
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_VARIABLE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_VARIABLE_MAPPED);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VECTOR_VARIABLE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VECTOR_VARIABLE_MAPPED);
  }

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
