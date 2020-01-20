//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"

//Utilities
#include "custom_utilities/generalized_eigenvalue_utility.hpp"


namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

	m.def("ComputeGeneralizedEigenvalues", &GeneralizedEigenvalueUtility::Compute<TUblasDenseSpace<double>>);
       
}

}  // namespace Python.
} // Namespace Kratos

