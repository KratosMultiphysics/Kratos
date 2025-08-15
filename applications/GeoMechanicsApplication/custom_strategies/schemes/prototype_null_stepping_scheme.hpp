// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "backward_euler_quasistatic_U_Pw_scheme.hpp"

namespace Kratos
{
template <class TSparseSpace, class TDenseSpace>
class PrototypeNullSteppingScheme : public BackwardEulerQuasistaticUPwScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrototypeNullSteppingScheme);

    void CalculateSystemContributions(
        Element& rCurrentElement,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemMatrixType& LHS_Contribution,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_INFO("CalculateSystemContributions for element ") << rCurrentElement.GetId() << std::endl;
        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.Calculate(INTERNAL_FORCES_VECTOR, RHS_Contribution, CurrentProcessInfo);
        auto fraction_of_unbalance = (CurrentProcessInfo[TIME] - CurrentProcessInfo[START_TIME]) /
                                     (CurrentProcessInfo[END_TIME] - CurrentProcessInfo[START_TIME]);

        KRATOS_INFO("fraction_of_unbalance") << fraction_of_unbalance << std::endl;
        KRATOS_INFO("mInternalForcesAtStartbyElementId[rCurrentElement.GetId()]") << mInternalForcesAtStartbyElementId[rCurrentElement.GetId()] << std::endl;
        KRATOS_INFO("mExternalForcesAtStartbyElementId[rCurrentElement.GetId()]") << mExternalForcesAtStartbyElementId[rCurrentElement.GetId()] << std::endl;
        const auto f_ext = -mInternalForcesAtStartbyElementId[rCurrentElement.GetId()] +
                           (mExternalForcesAtStartbyElementId[rCurrentElement.GetId()] +
                            mInternalForcesAtStartbyElementId[rCurrentElement.GetId()]) *
                               fraction_of_unbalance;
        KRATOS_INFO("f_ext") << f_ext << std::endl;

        RHS_Contribution += f_ext;

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    void CalculateRHSContribution(Element& rCurrentElement,
                                  typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_INFO("CalculateRHSContribution for element") << rCurrentElement.GetId() << std::endl;

        rCurrentElement.Calculate(INTERNAL_FORCES_VECTOR, RHS_Contribution, CurrentProcessInfo);
        KRATOS_INFO("RHS_Contribution") << RHS_Contribution << std::endl;
        auto fraction_of_unbalance = (CurrentProcessInfo[TIME] - CurrentProcessInfo[START_TIME]) /
                                     (CurrentProcessInfo[END_TIME] - CurrentProcessInfo[START_TIME]);

        KRATOS_INFO("fraction_of_unbalance") << fraction_of_unbalance << std::endl;
        KRATOS_INFO("mInternalForcesAtStartbyElementId[rCurrentElement.GetId()]") << mInternalForcesAtStartbyElementId[rCurrentElement.GetId()] << std::endl;
        KRATOS_INFO("mExternalForcesAtStartbyElementId[rCurrentElement.GetId()]") << mExternalForcesAtStartbyElementId[rCurrentElement.GetId()] << std::endl;
        const auto f_ext = -mInternalForcesAtStartbyElementId[rCurrentElement.GetId()] +
                           (mExternalForcesAtStartbyElementId[rCurrentElement.GetId()] +
                            mInternalForcesAtStartbyElementId[rCurrentElement.GetId()]) *
                               fraction_of_unbalance;
        KRATOS_INFO("f_ext") << f_ext << std::endl;

        RHS_Contribution += f_ext;
        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);
    }

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemMatrixType& A,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemVectorType& dX,
        typename GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::TSystemVectorType& b) override
    {
        KRATOS_INFO("InitializeSolutionStep of prototype") << std::endl;
        BackwardEulerQuasistaticUPwScheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(
            rModelPart, A, dX, b);
        KRATOS_INFO("Initialized base") << std::endl;

        if (!mIsInitialized) {
            for (auto& rElement : rModelPart.Elements()) {
                KRATOS_INFO("Initializing element") << rElement.GetId() << std::endl;
                mExternalForcesAtStartbyElementId.insert({rElement.GetId(), Vector{}});
                rElement.Calculate(EXTERNAL_FORCES_VECTOR,
                                   mExternalForcesAtStartbyElementId[rElement.GetId()],
                                   rModelPart.GetProcessInfo());
                KRATOS_INFO("Initialized external forces") << std::endl;

                mInternalForcesAtStartbyElementId.insert({rElement.GetId(), Vector{}});
                rElement.Calculate(INTERNAL_FORCES_VECTOR,
                                   mInternalForcesAtStartbyElementId[rElement.GetId()],
                                   rModelPart.GetProcessInfo());
                KRATOS_INFO("Initialized internal forces") << std::endl;

            }
            mIsInitialized = true;
        }
    }

private:
    std::map<std::size_t, Vector> mExternalForcesAtStartbyElementId;
    std::map<std::size_t, Vector> mInternalForcesAtStartbyElementId;

    bool mIsInitialized = false;
};
} // namespace Kratos