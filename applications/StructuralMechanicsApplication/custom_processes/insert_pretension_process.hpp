// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// --- Core Includes ---
#include "processes/process.h" // Process
#include "containers/model.h" // Model


namespace Kratos {


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertPretensionProcess final : public Process {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertPretensionProcess);

    InsertPretensionProcess(Model& rModel, Parameters Settings);

    InsertPretensionProcess(InsertPretensionProcess&&) noexcept;

    ~InsertPretensionProcess();

    void ExecuteBeforeSolutionLoop() override;

    const Parameters GetDefaultParameters() const override;

private:
    InsertPretensionProcess(const InsertPretensionProcess&) = delete;

    InsertPretensionProcess& operator=(const InsertPretensionProcess&) = delete;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class InsertPretensionProcess


} // namespace Kratos
