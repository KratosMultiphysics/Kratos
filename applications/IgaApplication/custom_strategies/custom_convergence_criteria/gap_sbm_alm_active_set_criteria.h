// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ /
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

#pragma once

// Project includes
#include "custom_conditions/gap_sbm_alm_contact_condition.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
class GapSbmALMActiveSetCriteria
    : public ConvergenceCriteria<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GapSbmALMActiveSetCriteria);

    using BaseType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using ClassType = GapSbmALMActiveSetCriteria<TSparseSpace, TDenseSpace>;
    using IndexType = std::size_t;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    GapSbmALMActiveSetCriteria()
        : BaseType()
        , mParameters(GetDefaultParameters())
    {
        this->SetEchoLevel(mParameters["echo_level"].GetInt());
    }

    explicit GapSbmALMActiveSetCriteria(Parameters ThisParameters)
        : BaseType()
    {
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        mParameters = ThisParameters;
        this->SetEchoLevel(ThisParameters["echo_level"].GetInt());
    }

    ~GapSbmALMActiveSetCriteria() override = default;

    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb) override
    {
        UpdateActiveSet(rModelPart);
        return true;
    }

    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb) override
    {
        const auto [n_active, n_changes] = UpdateActiveSet(rModelPart);
        const bool is_converged = n_changes == 0;

        KRATOS_INFO_IF("GapSbmALMActiveSetCriteria", this->GetEchoLevel() > 0)
            << "active_conditions=" << n_active
            << " active_set_changes=" << n_changes
            << " converged=" << is_converged << std::endl;

        return is_converged;
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "echo_level" : 0,
            "activity_tolerance" : 1.0e-12
        })");

        default_parameters.RecursivelyAddMissingParameters(BaseType::GetDefaultParameters());
        return default_parameters;
    }

    std::string Info() const override
    {
        return "GapSbmALMActiveSetCriteria";
    }

protected:
    std::pair<IndexType, IndexType> UpdateActiveSet(ModelPart& rModelPart)
    {
        if (!rModelPart.HasSubModelPart("ContactInterface")) {
            return {0, 0};
        }

        // const double activity_tolerance = mParameters["activity_tolerance"].GetDouble();

        const double activity_tolerance = 1e-3;
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        IndexType number_of_active_conditions = 0;
        IndexType number_of_changes = 0;

        ModelPart& r_contact_interface_model_part = rModelPart.GetSubModelPart("ContactInterface");
        for (auto& r_contact_model_part : r_contact_interface_model_part.SubModelParts()) {
            if (!r_contact_model_part.HasSubModelPart("contact")) {
                continue;
            }

            ModelPart& r_contact_sub_model_part = r_contact_model_part.GetSubModelPart("contact");
            for (auto& r_condition : r_contact_sub_model_part.Conditions()) {
                const auto* p_alm_condition = dynamic_cast<const GapSbmALMContactCondition*>(&r_condition);
                if (p_alm_condition == nullptr) {
                    continue;
                }

                const double augmented_traction =
                    p_alm_condition->CalculateAugmentedTraction(r_process_info);
                const int new_activation_level = augmented_traction > activity_tolerance ? 3 : 0;
                const int old_activation_level = r_condition.GetValue(ACTIVATION_LEVEL);

                if (new_activation_level != old_activation_level) {
                    r_condition.SetValue(OLD_ACTIVATION_LEVEL, old_activation_level);
                    r_condition.SetValue(ACTIVATION_LEVEL, new_activation_level);
                    ++number_of_changes;
                }

                if (new_activation_level == 3) {
                    ++number_of_active_conditions;
                }
            }
        }

        return {number_of_active_conditions, number_of_changes};
    }

private:
    Parameters mParameters;
};

} // namespace Kratos
