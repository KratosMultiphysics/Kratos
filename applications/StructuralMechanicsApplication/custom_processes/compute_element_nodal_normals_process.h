// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

/**
 * @class ComputeElementNodalNormalsProcess
 * @ingroup StructuralMechanicsApplication
 * @brief Computes element-based averaged nodal unit normals and stores them in NORMAL.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ComputeElementNodalNormalsProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ComputeElementNodalNormalsProcess); // same pointer as for base Proccess class. Because base API expects Process::Pointer (shared_ptr)

    ComputeElementNodalNormalsProcess() = default;

    ComputeElementNodalNormalsProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters
        ) override;

    ~ComputeElementNodalNormalsProcess() override = default;

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    const Parameters GetDefaultParameters() const override;

    std::string Info() const override
    {
        return "ComputeElementNodalNormalsProcess";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

private:
    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.StructuralMechanicsApplication", Process, ComputeElementNodalNormalsProcess) // The base Process class uses "Processes.KratosMultiphysics" because it's part of the core Kratos library, not an application. My class belongs to StructuralMechanicsApplication. (Processes.StructuralMechanicsApplication is a registry key like a folder path)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.All", Process, ComputeElementNodalNormalsProcess)                            // Core processes → "Processes.KratosMultiphysics", Application processes → "Processes.ApplicationName". Both also register in → "Processes.All" for global discovery

    // process internal state

    // target data (which model part)
    ModelPart* mpModelPart = nullptr;       // which model part the process will operate on. It is a pointer (initialized to nullptr) so the class can be default-constructed for registry prototype creation. it gets set in constructor that receives ModelPart&
    // when to execute the computation
    bool mComputeOnInitialize = true;       // Option flag: run normal computation in ExecuteInitialize(). Default true => compute once at the start unless user disables it in parameters
    bool mRecomputeEachStep = false;        // Option flag: run normal computation in ExecuteInitializeSolutionStep(). Default false => do not recompute each step unless user enables it in parameters

    void ComputeNormals();
};

} // namespace Kratos
