//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes

// External includes

//Application includes
#include "custom_python/add_custom_utilities_to_python.h"

// Project includes
#include "includes/node.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

#include "custom_utilities/two_step_v_p_settings.h"
#include "custom_utilities/postprocess_utilities.h"
#include "custom_utilities/pfem_fluid_gid_io.h"

namespace Kratos {
namespace Python {
    namespace py = pybind11;

    void  AddCustomUtilitiesToPython(pybind11::module& m) {

        py::class_<PostProcessUtilities, PostProcessUtilities::Pointer>(m, "PostProcessUtilities")
            .def(py::init<>())
            .def("RebuildPostProcessModelPart", &PostProcessUtilities::RebuildPostProcessModelPart)
        ;

        py::class_<PfemFluidGidIO<>, PfemFluidGidIO<>::Pointer, GidIO<> >(m,
        "PfemFluidGidIO")
            .def(py::init<std::string const&, GiD_PostMode,
                MultiFileFlag,
                WriteDeformedMeshFlag,
                WriteConditionsFlag>())
            ;

    }

}  // namespace Python.
} // Namespace Kratos

