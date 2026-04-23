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

// --- STL Includes ---
#include <unordered_set>
#include <unordered_map>


namespace Kratos {


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertPreTensionOperation : public Operation {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertPreTensionOperation);

    InsertPreTensionOperation(Model& rModel, Parameters Settings);

    InsertPreTensionOperation(InsertPreTensionOperation&&) noexcept;

    ~InsertPreTensionOperation();

    void Execute() override;

    const Parameters GetDefaultParameters() const override;

    virtual std::string Info() const;

protected:
    virtual void InsertControlNodeConstraints(
        ModelPart& rModelPart,
        array_1d<double,3> SurfaceNormal,
        const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
        Node::Pointer pControlNode,
        const std::unordered_set<const Dof<double>*> rPositiveSideDofs) const = 0;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;

private:
    InsertPreTensionOperation(const InsertPreTensionOperation&) = delete;

    InsertPreTensionOperation& operator=(const InsertPreTensionOperation&) = delete;
}; // class InsertPreTensionOperation


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertDirichletPreTensionOperation final : public InsertPreTensionOperation {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertDirichletPreTensionOperation);

    using InsertPreTensionOperation::InsertPreTensionOperation;

    std::string Info() const override;

protected:
    /// @details Inserts a constraint that ties the average out-of-plane
    ///          relative displacement of duplicated nodes to a prescribed value.
    void InsertControlNodeConstraints(
        ModelPart& rModelPart,
        array_1d<double,3> SurfaceNormal,
        const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
        Node::Pointer pControlNode,
        const std::unordered_set<const Dof<double>*> rPositiveSideDofs) const override;
}; // class InsertDirichletPreTensionOperation


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) InsertNeumannPreTensionOperation final : public InsertPreTensionOperation {
public:
    KRATOS_CLASS_POINTER_DEFINITION(InsertNeumannPreTensionOperation);

    using InsertPreTensionOperation::InsertPreTensionOperation;

    std::string Info() const override;

protected:
    /// @details Loads the average out-of-plane displacement with a prescribed
    ///          value on both sides of the pre-tension surface in opposite
    ///          directions.
    void InsertControlNodeConstraints(
        ModelPart& rModelPart,
        array_1d<double,3> SurfaceNormal,
        const std::unordered_map<Node*,Node::Pointer> rDuplicateNodeMap,
        Node::Pointer pControlNode,
        const std::unordered_set<const Dof<double>*> rPositiveSideDofs) const override;
}; // class InsertNeumannPreTensionOperation


} // namespace Kratos
