/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/utility_functions.h"


namespace Kratos {

    namespace Python {

        using namespace pybind11;

        void  AddCustomUtilitiesToPython(pybind11::module& m) {

            class_<DEMFEMVolumeCouplingUtilities> (m, "DEMFEMVolumeCouplingUtilities")
                .def(init<>())
                .def("SetNodalCouplingWeightsOnFEMLinearly", &DEMFEMVolumeCouplingUtilities::SetNodalCouplingWeightsOnFEMLinearly)
                .def("CalculateDisplacementDifference", &DEMFEMVolumeCouplingUtilities::CalculateDisplacementDifference)
                .def("AssignPointLoads", &DEMFEMVolumeCouplingUtilities::AssignPointLoads)
                .def("CalculateNodalCouplingForces", &DEMFEMVolumeCouplingUtilities::CalculateNodalCouplingForces)
                .def("CalculateNodalDEMCouplingForces", &DEMFEMVolumeCouplingUtilities::CalculateNodalDEMCouplingForces)
                .def("CalculateMomentum", &DEMFEMVolumeCouplingUtilities::CalculateMomentum)
                .def("CalculateDEMForces", &DEMFEMVolumeCouplingUtilities::CalculateDEMForces)
            ;



        }
    }  // namespace Python
} // Namespace Kratos
