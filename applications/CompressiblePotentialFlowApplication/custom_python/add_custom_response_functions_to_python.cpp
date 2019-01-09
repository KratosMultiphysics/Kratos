// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:        BSD License
//	                license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "custom_python/add_custom_response_functions_to_python.h"

#include "boost/numeric/ublas/vector.hpp"
#include "spaces/ublas_space.h"

// Response Functions
// #include "custom_response_functions/adjoint_structural_response_function.h"
#include "custom_response_functions/adjoint_lift_response_function.h"
#include "custom_response_functions/adjoint_potential_static_scheme.h"
#include "custom_response_functions/adjoint_postprocess.h"


namespace Kratos {
namespace Python {

void  AddCustomResponseFunctionUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    typedef AdjointPotentialStaticScheme< SparseSpaceType, LocalSpaceType > AdjointPotentialStaticSchemeType;


    // Response Functions
    py::class_<AdjointLiftResponseFunction, AdjointLiftResponseFunction::Pointer, AdjointResponseFunction>
        (m, "AdjointLiftResponseFunction")
        .def(py::init<ModelPart&, Parameters>());
    py::class_<AdjointPotentialStaticSchemeType, AdjointPotentialStaticSchemeType::Pointer, BaseSchemeType>
    (m, "AdjointPotentialStaticScheme")
        .def(py::init<Parameters, AdjointResponseFunction::Pointer>());

    py::class_<AdjointPostprocess, AdjointPostprocess::Pointer>
    (m, "AdjointPostprocess")
        .def(py::init<ModelPart&, AdjointResponseFunction&, Parameters>())
        .def("Initialize", &AdjointPostprocess::Initialize)
        .def("InitializeSolutionStep", &AdjointPostprocess::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &AdjointPostprocess::FinalizeSolutionStep)
        .def("UpdateSensitivities", &AdjointPostprocess::UpdateSensitivities);

}

}  // namespace Python.
} // Namespace Kratos

