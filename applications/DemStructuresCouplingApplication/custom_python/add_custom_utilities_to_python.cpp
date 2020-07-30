/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/dem_structures_coupling_utilities.h"
#include "custom_utilities/compute_dem_face_load_utility.h"
#include "custom_utilities/interpolate_structural_solution_for_dem_utility.h"
#include "custom_utilities/control_module_fem_dem_utilities.hpp"
#include "custom_utilities/stress_failure_check_utilities.hpp"
#include "custom_utilities/post_process_utilities.hpp"
#include "custom_utilities/sand_production_utilities.hpp"

namespace Kratos {

    namespace Python {

        using namespace pybind11;

        void  AddCustomUtilitiesToPython(pybind11::module& m) {

            class_<DemStructuresCouplingUtilities> (m, "DemStructuresCouplingUtilities")
                .def(init<>())
                .def("TransferStructuresSkinToDem", &DemStructuresCouplingUtilities::TransferStructuresSkinToDem)
                .def("CheckProvidedProperties", &DemStructuresCouplingUtilities::CheckProvidedProperties)
                .def("SmoothLoadTrasferredToFem", &DemStructuresCouplingUtilities::SmoothLoadTrasferredToFem)
                .def("ComputeSandProduction", &DemStructuresCouplingUtilities::ComputeSandProduction)
                .def("ComputeSandProductionWithDepthFirstSearch", &DemStructuresCouplingUtilities::ComputeSandProductionWithDepthFirstSearch)
                .def("ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation", &DemStructuresCouplingUtilities::ComputeSandProductionWithDepthFirstSearchNonRecursiveImplementation)
                .def("ComputeTriaxialSandProduction", &DemStructuresCouplingUtilities::ComputeTriaxialSandProduction)
                .def("MarkBrokenSpheres", &DemStructuresCouplingUtilities::MarkBrokenSpheres)
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

            class_<ControlModuleFemDemUtilities> (m, "ControlModuleFemDemUtilities")
                .def(init<ModelPart&,ModelPart&,Parameters&>())
                .def("ExecuteInitialize", &ControlModuleFemDemUtilities::ExecuteInitialize)
                .def("ExecuteInitializeSolutionStep", &ControlModuleFemDemUtilities::ExecuteInitializeSolutionStep)
                .def("ExecuteFinalizeSolutionStep", &ControlModuleFemDemUtilities::ExecuteFinalizeSolutionStep)
            ;

            class_<StressFailureCheckUtilities> (m, "StressFailureCheckUtilities")
                .def(init<ModelPart&,Parameters&>())
                .def("ExecuteFinalizeSolutionStep", &StressFailureCheckUtilities::ExecuteFinalizeSolutionStep)
            ;

            class_<PostProcessUtilities,PostProcessUtilities::Pointer>(m, "PostProcessUtilities", module_local())
                .def(init<ModelPart&>())
                .def("GetStickyStatus", &PostProcessUtilities::GetStickyStatus)
                .def("GetInitialContinuumBonds", &PostProcessUtilities::GetInitialContinuumBonds)
                .def("GetCurrentContinuumBonds", &PostProcessUtilities::GetCurrentContinuumBonds)
                ;

            class_<SandProductionUtilities, SandProductionUtilities::Pointer>(m, "SandProductionUtilities")
                .def(init<>())
                .def("MarkSandProductionParticlesForErasing", &SandProductionUtilities::MarkSandProductionParticlesForErasing)
                ;



        }
    }  // namespace Python
} // Namespace Kratos
