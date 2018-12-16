/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/dem_structures_coupling_utilities.h"
#include "custom_utilities/compute_dem_face_load_utility.h"
#include "custom_utilities/interpolate_structural_solution_for_dem_utility.h"

namespace Kratos {

    namespace Python {

        using namespace pybind11;

        void  AddCustomUtilitiesToPython(pybind11::module& m) {

            class_<DemStructuresCouplingUtilities> (m, "DemStructuresCouplingUtilities")
                .def(init<>())
                .def("TransferStructuresSkinToDem", &DemStructuresCouplingUtilities::TransferStructuresSkinToDem)
                .def("CheckProvidedProperties", &DemStructuresCouplingUtilities::CheckProvidedProperties)
            ;
            
            class_<ComputeDEMFaceLoadUtility> (m, "ComputeDEMFaceLoadUtility")
                .def(init<>())
                .def("ClearDEMFaceLoads", &ComputeDEMFaceLoadUtility::ClearDEMFaceLoads)
                .def("CalculateDEMFaceLoads", &ComputeDEMFaceLoadUtility::CalculateDEMFaceLoads)
            ;

            class_<InterpolateStructuralSolutionForDEM> (m, "InterpolateStructuralSolutionForDEM")
                .def(init<>())
                .def("SaveStructuralSolution", &InterpolateStructuralSolutionForDEM::SaveStructuralSolution)
                .def("InterpolateStructuralSolution", &InterpolateStructuralSolutionForDEM::InterpolateStructuralSolution)
                .def("RestoreStructuralSolution", &InterpolateStructuralSolutionForDEM::RestoreStructuralSolution)
            ;
        }
    }  // namespace Python
} // Namespace Kratos
