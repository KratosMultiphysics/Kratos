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
#include "custom_utilities/control_module_fem_dem_2d_utilities.hpp"
#include "custom_utilities/stress_failure_check_utilities.hpp"
#include "custom_utilities/post_process_utilities.hpp"
#include "custom_utilities/sand_production_utilities.hpp"
#include "custom_utilities/multiaxial_control_module_fem_dem_generalized_2d_utilities.hpp"
#include "custom_utilities/effective_stresses_communicator_utility.hpp"
#include "custom_utilities/pore_pressure_communicator_utility.hpp"
#include "custom_utilities/permeability_tensor_communicator_utility.hpp"

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

            class_<EffectiveStressesCommunicatorUtility> (m, "EffectiveStressesCommunicatorUtility")
                .def(init<ModelPart&,ModelPart&>())
                .def("Initialize", &EffectiveStressesCommunicatorUtility::Initialize)
                .def("CopyWallCurrentEffectiveStressesToOldEffectiveStresses", &EffectiveStressesCommunicatorUtility::CopyWallCurrentEffectiveStressesToOldEffectiveStresses)
                .def("CommunicateCurrentRadialEffectiveStressesToDemWalls", &EffectiveStressesCommunicatorUtility::CommunicateCurrentRadialEffectiveStressesToDemWalls)
                .def("CommunicateGivenRadialEffectiveStressesToDemWalls", &EffectiveStressesCommunicatorUtility::CommunicateGivenRadialEffectiveStressesToDemWalls)
            ;

            class_<PorePressureCommunicatorUtility> (m, "PorePressureCommunicatorUtility")
                .def(init<ModelPart&,ModelPart&>())
                .def("Initialize", &PorePressureCommunicatorUtility::Initialize)
                .def("ComputeForceOnParticlesDueToPorePressureGradient", &PorePressureCommunicatorUtility::ComputeForceOnParticlesDueToPorePressureGradient)
            ;

            class_<PermeabilityTensorCommunicatorUtility> (m, "PermeabilityTensorCommunicatorUtility")
                .def(init<ModelPart&,ModelPart&>())
                .def("Initialize", &PermeabilityTensorCommunicatorUtility::Initialize)
                .def("TrasferUpdatedPermeabilityTensor", &PermeabilityTensorCommunicatorUtility::TrasferUpdatedPermeabilityTensor)
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

            class_<ControlModuleFemDem2DUtilities> (m, "ControlModuleFemDem2DUtilities")
                .def(init<ModelPart&,ModelPart&,Parameters&>())
                .def("ExecuteInitialize", &ControlModuleFemDem2DUtilities::ExecuteInitialize)
                .def("ExecuteInitializeSolutionStep", &ControlModuleFemDem2DUtilities::ExecuteInitializeSolutionStep)
                .def("ExecuteFinalizeSolutionStep", &ControlModuleFemDem2DUtilities::ExecuteFinalizeSolutionStep)
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

            class_<MultiaxialControlModuleFEMDEMGeneralized2DUtilities> (m, "MultiaxialControlModuleFEMDEMGeneralized2DUtilities")
                .def(init<ModelPart&,ModelPart&,Parameters&>())
                .def("ExecuteInitialize", &MultiaxialControlModuleFEMDEMGeneralized2DUtilities::ExecuteInitialize)
                .def("ExecuteInitializeSolutionStep", &MultiaxialControlModuleFEMDEMGeneralized2DUtilities::ExecuteInitializeSolutionStep)
                .def("ExecuteFinalizeSolutionStep", &MultiaxialControlModuleFEMDEMGeneralized2DUtilities::ExecuteFinalizeSolutionStep)
            ;


        }
    }  // namespace Python
} // Namespace Kratos
