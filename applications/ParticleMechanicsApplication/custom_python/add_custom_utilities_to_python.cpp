//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/mpm_search_element_utility.h"
#include "custom_utilities/mpm_create_element_condition_utility.h"

namespace Kratos{
namespace Python{

    void  AddCustomUtilitiesToPython(pybind11::module& m)
    {
        m.def("SearchElement2D", &MPMSearchElementUtility::SearchElement< 2 >);
        m.def("SearchElement3D", &MPMSearchElementUtility::SearchElement< 3 >);
        m.def("CreateMaterialPointElement", &MpmCreateElementConditionUtility::CreateMaterialPointElement);
        m.def("CreateMaterialPointCondition", &MpmCreateElementConditionUtility::CreateMaterialPointCondition);
    }

}  // namespace Python.
} // Namespace Kratos

