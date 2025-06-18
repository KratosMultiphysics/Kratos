//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_VELOCITY_BOSSAK_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED)
#define KRATOS_VELOCITY_BOSSAK_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_strategies/schemes/simple_steady_sensitivity_builder_scheme.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class VelocityBossakSensitivityBuilderScheme : public SimpleSteadySensitivityBuilderScheme
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(VelocityBossakSensitivityBuilderScheme);

    using BaseType = SimpleSteadySensitivityBuilderScheme;

    using NodeType = typename BaseType::NodeType;

    using ConditionType = typename BaseType::ConditionType;

    using ElementType = typename BaseType::ElementType;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    VelocityBossakSensitivityBuilderScheme(
        const double NewAlphaBossak,
        const IndexType Dimension,
        const IndexType BlockSize)
        : BaseType(Dimension, BlockSize)
    {
        //default values for the Newmark Scheme
        mAlphaBossak = NewAlphaBossak;
        mBetaNewmark = 0.25 * pow((1.00 - mAlphaBossak), 2);
        mGammaNewmark = 0.5 - mAlphaBossak;

        const int number_of_threads = ParallelUtilities::GetNumThreads();
        mMassMatrices.resize(number_of_threads);
        mDampingMatrices.resize(number_of_threads);
        mSecondDerivativesVectors.resize(number_of_threads);
    }

    /// Destructor.
    ~VelocityBossakSensitivityBuilderScheme() = default;

    ///@}
    ///@name Operations
    ///@{

    void Clear() override
    {
        BaseType::Clear();
        mMassMatrices.clear();
        mDampingMatrices.clear();
        mSecondDerivativesVectors.clear();
    }

    ///@}
    ///@name Operations
    ///@{

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        ModelPart& rSensitivityModelPart,
        AdjointResponseFunction& rResponseFunction) override
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, rSensitivityModelPart, rResponseFunction);

        const double delta_time = -rModelPart.GetProcessInfo()[DELTA_TIME];
        KRATOS_ERROR_IF(delta_time < 1.0e-12)
            << "Detected delta_time > 0 in the adjoint Bossak scheme. Adjoints "
               "are calculated in reverse time, therefore DELTA_TIME should be "
               "negative."
            << std::endl;

        // initializing constants
        ma1 = delta_time * mBetaNewmark / mGammaNewmark;
        mam = (1.0 - mAlphaBossak) / (mGammaNewmark * delta_time);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "VelocityBossakSensitivityBuilderScheme";
    }

    ///@}

private:
    ///@name Private Members
    ///@{

    double mAlphaBossak;
    double mBetaNewmark;
    double mGammaNewmark;

    double ma1;
    double mam;

    std::vector<Matrix> mMassMatrices;
    std::vector<Matrix> mDampingMatrices;
    std::vector<Vector> mSecondDerivativesVectors;

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateLHSAndRHS(
        ElementType& rElement,
        Matrix& rLHS,
        Vector& rRHS,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateEntityLHSAndRHS(rElement, rLHS, rRHS, rProcessInfo);
    }

    void CalculateLHSAndRHS(
        ConditionType& rCondition,
        Matrix& rLHS,
        Vector& rRHS,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateEntityLHSAndRHS(rCondition, rLHS, rRHS, rProcessInfo);
    }

    template<class EntityType>
    void CalculateEntityLHSAndRHS(
        EntityType& rEntity,
        Matrix& rLHS,
        Vector& rRHS,
        const ProcessInfo& rProcessInfo)
    {
        const int k = OpenMPUtils::ThisThread();

        auto& mass_matrix = mMassMatrices[k];
        auto& damping_matrix = mDampingMatrices[k];
        auto& second_derivatives_vector = mSecondDerivativesVectors[k];

        // following calls uses the same method calls as in the primal scheme to be consistent
        rEntity.CalculateLocalSystem(rLHS, rRHS, rProcessInfo);
        rEntity.CalculateMassMatrix(mass_matrix, rProcessInfo);
        rEntity.CalculateLocalVelocityContribution(damping_matrix, rRHS, rProcessInfo);

        AddDynamicsToLHS(rLHS, damping_matrix, mass_matrix);
        AddDynamicsToRHS(rRHS, second_derivatives_vector, mass_matrix, rEntity, rProcessInfo);
    }

    void AddDynamicsToLHS(
        Matrix& rLHS,
        const Matrix& rDampingMatrix,
        const Matrix& rMassMatrix) const
    {
        // multiplying time scheme factor
        rLHS *= ma1;

        // adding mass contribution to the dynamic stiffness
        if (rMassMatrix.size1() != 0) {
            noalias(rLHS) += mam * rMassMatrix;
        }

        // adding  damping contribution
        if (rDampingMatrix.size1() != 0) {
            noalias(rLHS) += rDampingMatrix;
        }
    }

    template <class TEntityType>
    void AddDynamicsToRHS(
        Vector& rRHS,
        Vector& rSecondDerivativesVector,
        const Matrix& rMassMatrix,
        TEntityType& rEntity,
        const ProcessInfo& rProcessInfo)
    {
        // adding inertia contribution
        if (rMassMatrix.size1() != 0) {
            rEntity.Calculate(PRIMAL_RELAXED_SECOND_DERIVATIVE_VALUES, rSecondDerivativesVector, rProcessInfo);
            noalias(rRHS) -= prod(rMassMatrix, rSecondDerivativesVector);
        }
    }

    ///@}

}; /* Class VelocityBossakSensitivityBuilderScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_VELOCITY_BOSSAK_SENSITIVITY_BUILDER_SCHEME_H_INCLUDED defined */
