//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//

#if !defined(KRATOS_ADD_EDGEBASED_LEVELSET_TO_PYTHON_H_INCLUDED)
#define KRATOS_ADD_EDGEBASED_LEVELSET_TO_PYTHON_H_INCLUDED

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define_python.h"
#include "processes/process.h"

#include "custom_utilities/edge_data.h"
#include "custom_utilities/edge_data_c2c.h"
#include "custom_utilities/edgebased_levelset.h"
#include "custom_utilities/edgebased_levelset_substep.h"
#include "custom_python/add_custom_edgebased_levelset_solver_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

    namespace Python
    {

        void AddCustomEdgeBasedLevelSetToPython(pybind11::module &pymodule);

    } // namespace Python.

} // namespace Kratos.

#endif // KRATOS_ADD_EDGEBASED_LEVELSET_TO_PYTHON_H_INCLUDED  defined
