//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// Elastic laws
#include "custom_constitutive/helmholtz_jacobian_stiffened_3d.h"

namespace Kratos::Python {

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< HelmholtzJacobianStiffened3D, typename HelmholtzJacobianStiffened3D::Pointer, ConstitutiveLaw >
    (m, "HelmholtzJacobianStiffened3DLaw").def(py::init<>() )
    ;

}

}  // namespace Kratos::Python.
