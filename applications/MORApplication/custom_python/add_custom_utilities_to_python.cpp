// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/define.h"
#include <pybind11/pybind11.h>

//Utilities
#include "custom_utilities/eigen_qr_utility.hpp"
#include "custom_utilities/ublas_wrapper.h"
#include "custom_utilities/mor_qr_utility.hpp"


namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

	//typedef EigenQrUtility EigenQrUtilityType;

		py::class_<EigenQrUtility, EigenQrUtility::Pointer>(m, "EigenQrUtility")
		.def(py::init<unsigned int, const std::size_t>());

        /*py::class_< storage_order_mor>(m, "MorQrUtility")
		.def(py::init<>());*/
       
}

}  // namespace Python.
} // Namespace Kratos

