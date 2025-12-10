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
#include "operations/operation.h" // Operation
#include "containers/model.h" // Model


namespace Kratos {


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertPretensionOperation final : public Operation {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertPretensionOperation);

    InsertPretensionOperation(Model& rModel, Parameters Settings);

    InsertPretensionOperation(InsertPretensionOperation&&) noexcept;

    ~InsertPretensionOperation();

    void Execute() override;

    const Parameters GetDefaultParameters() const override;

private:
    InsertPretensionOperation(const InsertPretensionOperation&) = delete;

    InsertPretensionOperation& operator=(const InsertPretensionOperation&) = delete;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class InsertPretensionOperation


} // namespace Kratos
