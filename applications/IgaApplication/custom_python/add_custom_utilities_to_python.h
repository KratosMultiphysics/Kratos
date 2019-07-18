/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED

// System includes

// External includes
#include <pybind11/pybind11.h>
// Project includes
#include "includes/define_python.h"

#include "iga_application_variables.h"
#include "custom_utilities/anurbs.h"
#include "custom_utilities/node_curve_geometry_3d.h"
#include "custom_utilities/node_surface_geometry_3d.h"

#include "custom_utilities/brep_json_io.h"
#include "custom_utilities/nurbs_brep_modeler.h"

#include "custom_utilities/iga_flags.h"
#include "includes/define.h"

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m);

} // namespace Python
} // namespace Kratos

#endif // !defined(KRATOS_ADD_UTILITIES_TO_PYTHON_H_INCLUDED)
