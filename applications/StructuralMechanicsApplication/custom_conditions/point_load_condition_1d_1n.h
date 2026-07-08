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

// --- Core Includes ---
#include "includes/condition.h"
#include "geometries/point_2d.h"
#include "includes/variables.h"

// --- Structural Includes ---
#include "structural_mechanics_application_variables.h" // POINT_LOAD_X


namespace Kratos {


/// @brief Simple condition loading a single degree-of-freedom.
/// @details This condition is primarily used to load virtual DoFs belonging to @ref MultfreedomConstraint "constraints".
/// @see @ref pre_tensioning "Pre-Tensioning"
/// @ingroup pre_tensioning
class PointLoadCondition1D1N : public Condition {
public:
    KRATOS_CLASS_POINTER_DEFINITION(PointLoadCondition1D1N);

    PointLoadCondition1D1N() = default;

    PointLoadCondition1D1N(
        Condition::IndexType Id,
        Geometry<Node>::Pointer pGeometry,
        Properties::Pointer pProperties)
            : Condition(Id, pGeometry, pProperties)
    {}

    Condition::Pointer Clone(
        Condition::IndexType Id,
        const Condition::NodesArrayType& rNodes) const override {
            KRATOS_ERROR_IF_NOT(rNodes.size() == 1)
                << "PointLoadCondition1D1N::Clone expects 1 node but got " << rNodes.size();
            return Condition::Pointer(new PointLoadCondition1D1N(
                Id,
                Geometry<Node>::Pointer(new Point2D<Node>(rNodes)),
                this->pGetProperties()));
    }

    void EquationIdVector(
        Condition::EquationIdVectorType& rIndices,
        const ProcessInfo&) const override {
            rIndices.resize(1);
            rIndices[0] = this->FindDof().EquationId();
    }

    void GetDofList(
        Condition::DofsVectorType& rDofs,
        const ProcessInfo&) const override {
            rDofs.resize(1);
            rDofs[0] = &this->FindDof();
    }

    void GetValuesVector(
        Vector& rOutput,
        int Step) const override {
            rOutput.resize(1);
            rOutput[0] = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
    }

    void GetFirstDerivativesVector(
        Vector& rOutput,
        int Step) const override {
            rOutput.resize(1);
            rOutput[0] = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X, Step);
    }

    void GetSecondDerivativesVector(
        Vector& rOutput,
        int Step) const override {
            rOutput.resize(1);
            rOutput[0] = this->GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION_X, Step);
    }

    void CalculateRightHandSide(
        Condition::VectorType& rRhs,
        const ProcessInfo&) override {
            rRhs.resize(1);
            rRhs[0] = this->GetLoadFactor() * this->GetValue(POINT_LOAD_X);
    }

    void CalculateLocalSystem(
        Condition::MatrixType& rLhs,
        Condition::VectorType& rRhs,
        const ProcessInfo& rProcessInfo) override {
            rLhs.resize(1, 1, false);
            rLhs(0, 0) = 0.0;
            this->CalculateRightHandSide(rRhs, rProcessInfo);
    }

    void CalculateMassMatrix(
        Condition::MatrixType& rMatrix,
        const ProcessInfo&) override {
            rMatrix.resize(1, 1, false);
            rMatrix(0, 0) = 0.0;
    }

    void CalculateDampingMatrix(
        Condition::MatrixType& rMatrix,
        const ProcessInfo&) override {
            rMatrix.resize(1, 1, false);
            rMatrix(0, 0) = 0.0;
    }

    std::string Info() const override {
        return "PointLoadCondition1D1N";
    }

private:
    Dof<double>& FindDof() const {
        auto& r_dofs = this->GetGeometry()[0].GetDofs();
        const auto it_dof = std::find_if(
            r_dofs.begin(),
            r_dofs.end(),
            [] (const auto& rp_dof) {return rp_dof->GetVariable().Key() == DISPLACEMENT_X.Key();});
        KRATOS_ERROR_IF(it_dof == r_dofs.end())
            << this->Info() << ": " << DISPLACEMENT_X.Name()
            << " is not a DoF in node " << this->GetGeometry()[0].Id();
        return **it_dof;
    }

    double GetLoadFactor() const {
        return this->GetProperties().Has(CONSTRAINT_SCALE_FACTOR)
            ? this->GetProperties()[CONSTRAINT_SCALE_FACTOR]
            : 1.0;
    }
}; // class PointLoadCondition1D1N

} // namespace Kratos
