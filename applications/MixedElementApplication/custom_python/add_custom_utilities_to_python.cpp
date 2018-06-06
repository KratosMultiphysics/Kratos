//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
//#include "linear_solvers/linear_solver.h"
#include "custom_utilities/plot_damage.h"
#include "custom_utilities/check_constitutive_law.h"

#include "includes/model_part.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<PlotUtility > (m,"PlotUtils")
        .def(init<>())
        .def("PlotVariable", &PlotUtility::PlotVariable)
    ;

    class_<CheckConstitutiveLawUtility > (m,"CheckConstitutiveLawUtils")
        .def(init<Properties::Pointer , ProcessInfo&, double >())
        .def("CheckConstitutiveLaw", &CheckConstitutiveLawUtility::CheckConstitutiveLaw)
        .def("GetScalarVar",&CheckConstitutiveLawUtility::GetScalarVar)
        .def("InitializeSolutionStep",&CheckConstitutiveLawUtility::InitializeSolutionStep)
        .def("FinalizeSolutionStep",&CheckConstitutiveLawUtility::FinalizeSolutionStep)
    ;

}





}  // namespace Python.

} // Namespace Kratos

