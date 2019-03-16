//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Jordi Cotela
//                  Suneth Warnakulasuriya
//

#if !defined(KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED

// System includes
#include <unordered_set>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/indirect_scalar.h"
#include "utilities/scheme_extension.h"
#include "utilities/time_discretization.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for steady and dynamic equations, using Bossak time integration.
/**
 * It can be used for either first- or second-order time derivatives. Elements
 * and conditions must provide a specialization of SchemeExtension via
 * their data value container, which allows the scheme to operate independently
 * of the variable arrangements in the element or condition.
 */
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedBossakVelocityScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBossakVelocityScheme);

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TSystemMatrixType SystemMatrixType;

    typedef typename BaseType::TSystemVectorType SystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.

    ResidualBasedBossakVelocityScheme(const double AlphaBossak,
                                      const bool UpdateAcceleration = true,
                                      const bool UpdateDisplacement = true)
        : mAlphaBossak(AlphaBossak),
          mUpdateAcceleration(UpdateAcceleration),
          mUpdateDisplacement(UpdateDisplacement)
    {
        KRATOS_INFO("ResidualBasedBossakVelocityScheme")
            << " Using bossak velocity scheme with alpha_bossak = " << std::scientific
            << mAlphaBossak << " [UpdateAcceleration: " << mUpdateAcceleration
            << ", UpdateDisplacement: " << mUpdateDisplacement << "]\n";

        // Allocate auxiliary memory.
        int num_threads = OpenMPUtils::GetNumThreads();

        mMassMatrix.resize(num_threads);
        mDampingMatrix.resize(num_threads);
        mValuesVector.resize(num_threads);
        mSecondDerivativeValuesVector.resize(num_threads);
        mSecondDerivativeValuesVectorOld.resize(num_threads);
        mIndirectCurrentVelocityVector.resize(num_threads);
        mIndirectCurrentAccelerationVector.resize(num_threads);
        mIndirectCurrentDisplacementVector.resize(num_threads);
        mIndirectOldVelocityVector.resize(num_threads);
        mIndirectOldAccelerationVector.resize(num_threads);
        mIndirectOldDisplacementVector.resize(num_threads);
    }

    /// Destructor.
    ~ResidualBasedBossakVelocityScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void InitializeSolutionStep(ModelPart& rModelPart,
                                SystemMatrixType& rA,
                                SystemVectorType& rDx,
                                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

        KRATOS_ERROR_IF(delta_time < std::numeric_limits<double>::epsilon())
            << "detected delta_time = 0 in the Bossak Scheme ... "
               "check if the time step is created correctly for "
               "the current model part.";

        ResidualBasedBossakVelocityScheme::CalculateBossakConstants(
            mBossak, mAlphaBossak, delta_time);

        KRATOS_CATCH("");
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        this->UpdateTimeSchemeVariables(rModelPart);

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const int k = OpenMPUtils::ThisThread();

        (pCurrentElement)->InitializeNonLinearIteration(rCurrentProcessInfo);
        (pCurrentElement)->CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentElement)->CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);

        if (mUpdateAcceleration)
        {
            (pCurrentElement)->CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(pCurrentElement, rRHS_Contribution, mDampingMatrix[k],
                             mMassMatrix[k], rCurrentProcessInfo);
        }
        AddDynamicsToLHS(rLHS_Contribution, mDampingMatrix[k], mMassMatrix[k],
                         rCurrentProcessInfo);

        (pCurrentElement)->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Calculate_RHS_Contribution(Element::Pointer pCurrentElement,
                                    LocalSystemVectorType& rRHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        int k = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current element
        (pCurrentElement)->InitializeNonLinearIteration(rCurrentProcessInfo);

        // basic operations for the element considered
        (pCurrentElement)->CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentElement)->CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentElement)->EquationIdVector(rEquationId, rCurrentProcessInfo);

        // adding the dynamic contributions (static is already included)
        if (mUpdateAcceleration)
        {
            (pCurrentElement)->CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(pCurrentElement, rRHS_Contribution, mDampingMatrix[k],
                             mMassMatrix[k], rCurrentProcessInfo);
        }
    }

    void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                LocalSystemMatrixType& rLHS_Contribution,
                                                LocalSystemVectorType& rRHS_Contribution,
                                                Condition::EquationIdVectorType& rEquationId,
                                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        int k = OpenMPUtils::ThisThread();

        (pCurrentCondition)->InitializeNonLinearIteration(rCurrentProcessInfo);
        (pCurrentCondition)->CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentCondition)->CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentCondition)->EquationIdVector(rEquationId, rCurrentProcessInfo);

        if (mUpdateAcceleration)
        {
            (pCurrentCondition)->CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(pCurrentCondition, rRHS_Contribution,
                             mDampingMatrix[k], mMassMatrix[k], rCurrentProcessInfo);
        }
        AddDynamicsToLHS(rLHS_Contribution, mDampingMatrix[k], mMassMatrix[k],
                         rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Condition_Calculate_RHS_Contribution(Condition::Pointer pCurrentCondition,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int k = OpenMPUtils::ThisThread();

        (pCurrentCondition)->InitializeNonLinearIteration(rCurrentProcessInfo);
        (pCurrentCondition)->CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentCondition)->CalculateLocalVelocityContribution(mDampingMatrix[k], rRHS_Contribution, rCurrentProcessInfo);
        (pCurrentCondition)->EquationIdVector(rEquationId, rCurrentProcessInfo);

        // adding the dynamic contributions (static is already included)
        if (mUpdateAcceleration)
        {
            (pCurrentCondition)->CalculateMassMatrix(mMassMatrix[k], rCurrentProcessInfo);
            AddDynamicsToRHS(pCurrentCondition, rRHS_Contribution,
                             mDampingMatrix[k], mMassMatrix[k], rCurrentProcessInfo);
        }
        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedBossakVelocityScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }
    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    struct BossakConstants
    {
        double Alpha;
        double Gamma;
        double Beta;
        double C0;
        double C1;
        double C2;
        double C3;
        double C4;
        double C5;
        double C6;
    };

    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<LocalSystemVectorType> mSecondDerivativeValuesVectorOld;
    std::vector<LocalSystemVectorType> mSecondDerivativeValuesVector;
    std::vector<LocalSystemVectorType> mValuesVector;

    std::vector<LocalSystemMatrixType> mMassMatrix;
    std::vector<LocalSystemMatrixType> mDampingMatrix;

    double mAlphaBossak;
    bool mUpdateAcceleration;
    bool mUpdateDisplacement;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    //****************************************************************************

    /**
    Kdyn = am*M + D + a1*K
     */
    void AddDynamicsToLHS(LocalSystemMatrixType& rLHS_Contribution,
                          LocalSystemMatrixType& rDampingMatrix,
                          LocalSystemMatrixType& rMassMatrix,
                          ProcessInfo& CurrentProcessInfo)
    {
        // multipling time scheme factor
        rLHS_Contribution *= mBossak.C1;

        // adding mass contribution to the dynamic stiffness
        if (rMassMatrix.size1() != 0 && mUpdateAcceleration) // if M matrix declared
        {
            noalias(rLHS_Contribution) += mBossak.C0 * rMassMatrix;
        }

        // adding  damping contribution
        if (rDampingMatrix.size1() != 0) // if M matrix declared
        {
            noalias(rLHS_Contribution) += rDampingMatrix;
        }
    }

    //****************************************************************************

    /// Add Bossak contributions from the inertial term to the RHS vector.
    /** This essentially performs bdyn = b - M*acc for the current element.
     *  Note that viscous/pressure contributions to the RHS are expected to be added by the element itself.
     *  @param[in] rCurrentElement The fluid element we are assembling.
     *  @param[in/out] rRHS_Contribution The right hand side term where the contribution will be added.
     *  @param[in] rD The elemental velocity/pressure LHS matrix.
     *  @param[in] rM The elemental acceleration LHS matrix.
     *  @param[in] rCurrentProcessInfo ProcessInfo instance for the containing ModelPart.
     */
    void AddDynamicsToRHS(Element::Pointer rCurrentElement,
                          LocalSystemVectorType& rRHS_Contribution,
                          LocalSystemMatrixType& rDampingMatrix,
                          LocalSystemMatrixType& rMassMatrix,
                          ProcessInfo& rCurrentProcessInfo)
    {
        // adding inertia contribution
        if (rMassMatrix.size1() != 0)
        {
            int k = OpenMPUtils::ThisThread();
            rCurrentElement->GetSecondDerivativesVector(
                mSecondDerivativeValuesVector[k], 0);
            (mSecondDerivativeValuesVector[k]) *= (1.00 - mBossak.Alpha);
            rCurrentElement->GetSecondDerivativesVector(
                mSecondDerivativeValuesVectorOld[k], 1);
            noalias(mSecondDerivativeValuesVector[k]) +=
                mBossak.Alpha * mSecondDerivativeValuesVectorOld[k];
            noalias(rRHS_Contribution) -=
                prod(rMassMatrix, mSecondDerivativeValuesVector[k]);
        }
    }

    /// Add Bossak contributions from the inertial term to the RHS vector.
    /** This essentially performs bdyn = b - M*acc for the current condition.
     *  Note that viscous/pressure contributions to the RHS are expected to be added by the element condition.
     *  @param[in] rCurrentCondition The fluid condition we are assembling.
     *  @param[in/out] rRHS_Contribution The right hand side term where the contribution will be added.
     *  @param[in] rD The elemental velocity/pressure LHS matrix.
     *  @param[in] rM The elemental acceleration LHS matrix.
     *  @param[in] rCurrentProcessInfo ProcessInfo instance for the containing ModelPart.
     */
    void AddDynamicsToRHS(Condition::Pointer rCurrentCondition,
                          LocalSystemVectorType& rRHS_Contribution,
                          LocalSystemMatrixType& rDampingMatrix,
                          LocalSystemMatrixType& rMassMatrix,
                          ProcessInfo& rCurrentProcessInfo)
    {
        // adding inertia contribution
        if (rMassMatrix.size1() != 0)
        {
            int k = OpenMPUtils::ThisThread();
            rCurrentCondition->GetSecondDerivativesVector(
                mSecondDerivativeValuesVector[k], 0);
            (mSecondDerivativeValuesVector[k]) *= (1.00 - mBossak.Alpha);
            rCurrentCondition->GetSecondDerivativesVector(
                mSecondDerivativeValuesVectorOld[k], 1);
            noalias(mSecondDerivativeValuesVector[k]) +=
                mBossak.Alpha * mSecondDerivativeValuesVectorOld[k];

            noalias(rRHS_Contribution) -=
                prod(rMassMatrix, mSecondDerivativeValuesVector[k]);
        }
    }

    void UpdateTimeSchemeVariables(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        UpdateAcceleration(rModelPart);
        UpdateDisplacement(rModelPart);

        KRATOS_CATCH("");
    }

    void UpdateAcceleration(double& rCurrentAcceleration,
                            const double CurrentVelocity,
                            const double OldVelocity,
                            const double OldAcceleration)
    {
        rCurrentAcceleration = mBossak.C2 * (CurrentVelocity - OldVelocity) -
                               mBossak.C3 * OldAcceleration;
    }

    void UpdateDisplacement(double& rCurrentDisplacement,
                            const double OldDisplacement,
                            const double OldVelocity,
                            const double CurrentAcceleration,
                            const double OldAcceleration)
    {
        rCurrentDisplacement = OldDisplacement + mBossak.C6 * OldVelocity +
                               mBossak.C4 * OldAcceleration + mBossak.C5 * CurrentAcceleration;
    }

    static void CalculateBossakConstants(BossakConstants& rBossakConstants,
                                         const double Alpha,
                                         const double DeltaTime)
    {
        TimeDiscretization::Bossak bossak(Alpha, 0.25, 0.5);
        rBossakConstants.Alpha = bossak.GetAlphaM();
        rBossakConstants.Gamma = bossak.GetGamma();
        rBossakConstants.Beta = bossak.GetBeta();

        rBossakConstants.C0 =
            (1.0 - rBossakConstants.Alpha) / (rBossakConstants.Gamma * DeltaTime);
        rBossakConstants.C1 =
            DeltaTime / (rBossakConstants.Beta * rBossakConstants.Gamma);
        rBossakConstants.C2 = 1.0 / (rBossakConstants.Gamma * DeltaTime);
        rBossakConstants.C3 = (1.0 - rBossakConstants.Gamma) / rBossakConstants.Gamma;
        rBossakConstants.C4 =
            std::pow(DeltaTime, 2) * (-2.0 * rBossakConstants.Beta + 1.0) / 2.0;
        rBossakConstants.C5 = std::pow(DeltaTime, 2) * rBossakConstants.Beta;
        rBossakConstants.C6 = DeltaTime;
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater =
        TSparseSpace::CreateDofUpdater();

    BossakConstants mBossak;

    std::vector<std::vector<IndirectScalar<double>>> mIndirectCurrentVelocityVector;
    std::vector<std::vector<IndirectScalar<double>>> mIndirectCurrentAccelerationVector;
    std::vector<std::vector<IndirectScalar<double>>> mIndirectCurrentDisplacementVector;
    std::vector<std::vector<IndirectScalar<double>>> mIndirectOldVelocityVector;
    std::vector<std::vector<IndirectScalar<double>>> mIndirectOldAccelerationVector;
    std::vector<std::vector<IndirectScalar<double>>> mIndirectOldDisplacementVector;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    // class to hold all the derivatives for updated target variable
    template <int TNumberOfVariables>
    struct DerivativesHolder
    {
        double* mpTargetDerivative;
        double mRequiredDerivatives[TNumberOfVariables];

        bool operator==(const DerivativesHolder& value) const
        {
            return (this->mpTargetDerivative == value.mpTargetDerivative);
        }
    };

    // class for hash function
    template <int TNumberOfVariables>
    class DerivativesHasher
    {
    public:
        // id is returned as hash function
        std::size_t operator()(const DerivativesHolder<TNumberOfVariables>& value) const
        {
            return (std::size_t)value.mpTargetDerivative;
        }
    };

    void UpdateAcceleration(ModelPart& rModelPart)
    {
        if (!mUpdateAcceleration)
            return;

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int number_of_elements = rModelPart.NumberOfElements();

        typedef DerivativesHolder<3> derivatives_holder;
        typedef std::unordered_set<derivatives_holder, DerivativesHasher<3>> derivatives_set;

        derivatives_set global_derivatives_list;
        global_derivatives_list.reserve(number_of_elements * 20);

#pragma omp parallel
        {
            derivatives_set derivatives_tmp_set;
            derivatives_tmp_set.reserve(20000);

#pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_elements; ++i)
            {
                Element& r_element = *(rModelPart.ElementsBegin() + i);
                Geometry<Node<3>>& r_geometry = r_element.GetGeometry();
                SchemeExtension& r_extensions = *r_element.GetValue(SCHEME_EXTENSION);

                const int k = OpenMPUtils::ThisThread();

                for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
                {
                    r_extensions.GetFirstDerivativesVector(
                        i_node, mIndirectCurrentVelocityVector[k], 0, r_current_process_info);
                    r_extensions.GetFirstDerivativesVector(
                        i_node, mIndirectOldVelocityVector[k], 1, r_current_process_info);
                    r_extensions.GetSecondDerivativesVector(
                        i_node, mIndirectCurrentAccelerationVector[k], 0, r_current_process_info);
                    r_extensions.GetSecondDerivativesVector(
                        i_node, mIndirectOldAccelerationVector[k], 1, r_current_process_info);

                    for (unsigned int i_var = 0;
                         i_var < mIndirectCurrentAccelerationVector[k].size(); ++i_var)
                    {
                        derivatives_holder current_derivative;
                        current_derivative.mpTargetDerivative =
                            mIndirectCurrentAccelerationVector[k][i_var].pGetValue();
                        current_derivative.mRequiredDerivatives[0] =
                            mIndirectCurrentVelocityVector[k][i_var];
                        current_derivative.mRequiredDerivatives[1] =
                            mIndirectOldVelocityVector[k][i_var];
                        current_derivative.mRequiredDerivatives[2] =
                            mIndirectOldAccelerationVector[k][i_var];

                        derivatives_tmp_set.insert(current_derivative);
                    }
                }
            }

#pragma omp critical
            {
                global_derivatives_list.insert(derivatives_tmp_set.begin(),
                                               derivatives_tmp_set.end());
            }
        }

        const int number_of_update_variables = global_derivatives_list.size();
        std::vector<derivatives_holder> derivatives_holder_list(number_of_update_variables);

        int local_index = 0;
        for (auto r_derivatives : global_derivatives_list)
        {
            derivatives_holder_list[local_index++] = r_derivatives;
        }

#pragma omp parallel for
        for (int i = 0; i < number_of_update_variables; ++i)
        {
            derivatives_holder& r_derivatives = derivatives_holder_list[i];
            UpdateAcceleration(*r_derivatives.mpTargetDerivative,
                               r_derivatives.mRequiredDerivatives[0],
                               r_derivatives.mRequiredDerivatives[1],
                               r_derivatives.mRequiredDerivatives[2]);
        }
    }

    void UpdateDisplacement(ModelPart& rModelPart)
    {
        if (!mUpdateDisplacement)
            return;

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const int number_of_elements = rModelPart.NumberOfElements();

        typedef DerivativesHolder<4> derivatives_holder;
        typedef std::unordered_set<derivatives_holder, DerivativesHasher<4>> derivatives_set;

        derivatives_set global_derivatives_list;
        global_derivatives_list.reserve(number_of_elements * 20);

#pragma omp parallel
        {
            derivatives_set derivatives_tmp_set;
            derivatives_tmp_set.reserve(20000);

#pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_elements; ++i)
            {
                Element& r_element = *(rModelPart.ElementsBegin() + i);
                Geometry<Node<3>>& r_geometry = r_element.GetGeometry();
                SchemeExtension& r_extensions = *r_element.GetValue(SCHEME_EXTENSION);

                const int k = OpenMPUtils::ThisThread();

                for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
                {
                    r_extensions.GetZeroDerivativesVector(
                        i_node, mIndirectCurrentDisplacementVector[k], 0, r_current_process_info);
                    r_extensions.GetZeroDerivativesVector(
                        i_node, mIndirectOldDisplacementVector[k], 1, r_current_process_info);
                    r_extensions.GetFirstDerivativesVector(
                        i_node, mIndirectOldVelocityVector[k], 1, r_current_process_info);
                    r_extensions.GetSecondDerivativesVector(
                        i_node, mIndirectCurrentAccelerationVector[k], 0, r_current_process_info);
                    r_extensions.GetSecondDerivativesVector(
                        i_node, mIndirectOldAccelerationVector[k], 1, r_current_process_info);

                    for (unsigned int i_var = 0;
                         i_var < mIndirectCurrentDisplacementVector[k].size(); ++i_var)
                    {
                        derivatives_holder current_derivative;
                        current_derivative.mpTargetDerivative =
                            mIndirectCurrentDisplacementVector[k][i_var].pGetValue();
                        current_derivative.mRequiredDerivatives[0] =
                            mIndirectOldDisplacementVector[k][i_var];
                        current_derivative.mRequiredDerivatives[1] =
                            mIndirectOldVelocityVector[k][i_var];
                        current_derivative.mRequiredDerivatives[2] =
                            mIndirectCurrentAccelerationVector[k][i_var];
                        current_derivative.mRequiredDerivatives[3] =
                            mIndirectOldAccelerationVector[k][i_var];

                        derivatives_tmp_set.insert(current_derivative);
                    }
                }
            }

#pragma omp critical
            {
                global_derivatives_list.insert(derivatives_tmp_set.begin(),
                                               derivatives_tmp_set.end());
            }
        }

        const int number_of_update_variables = global_derivatives_list.size();
        std::vector<derivatives_holder> derivatives_holder_list(number_of_update_variables);

        int local_index = 0;
        for (auto r_derivatives : global_derivatives_list)
        {
            derivatives_holder_list[local_index++] = r_derivatives;
        }

#pragma omp parallel for
        for (int i = 0; i < number_of_update_variables; ++i)
        {
            derivatives_holder& r_derivatives = derivatives_holder_list[i];
            UpdateDisplacement(*r_derivatives.mpTargetDerivative,
                               r_derivatives.mRequiredDerivatives[0],
                               r_derivatives.mRequiredDerivatives[1],
                               r_derivatives.mRequiredDerivatives[2],
                               r_derivatives.mRequiredDerivatives[3]);
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedBossakVelocityScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME_H_INCLUDED defined */
