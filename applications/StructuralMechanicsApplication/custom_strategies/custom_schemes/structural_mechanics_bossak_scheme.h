//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "structural_mechanics_application_variables.h"
#include "includes/cfd_variables.h" //TODO: For OSS_SWITCH (think about a better place for this variable)
#include "processes/calculate_nodal_area_process.h"

// Application includes
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @brief Structural mechanics Bossak scheme
 * This scheme extends the implementation in the ResidualBasedBossakDisplacementScheme to do some structural mechanics specifics
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 */
template<class TSparseSpace,  class TDenseSpace >
class StructuralMechanicsBossakScheme
    : public ResidualBasedBossakDisplacementScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(StructuralMechanicsBossakScheme);

    /// Base type for the scheme
    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    /// Implicit base type for the scheme
    using BossakBaseType = ResidualBasedBossakDisplacementScheme<TSparseSpace, TDenseSpace>;

    /// Class type for the scheme
    using ClassType = StructuralMechanicsBossakScheme<TSparseSpace, TDenseSpace>;

    /// Type for the dofs array within BossakBaseType
    using DofsArrayType = typename BossakBaseType::DofsArrayType;

    /// Type for the system matrix within BossakBaseType
    using TSystemMatrixType = typename BossakBaseType::TSystemMatrixType;

    /// Type for the system vector within BossakBaseType
    using TSystemVectorType = typename BossakBaseType::TSystemVectorType;

    /// Pointer type for the BaseType
    using BaseTypePointer = typename BaseType::Pointer;

    /// Component type as 'double'
    using ComponentType = double;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit StructuralMechanicsBossakScheme(Parameters ThisParameters)
        : BossakBaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // For pure Newmark Scheme
        this->mNewmark.gamma = 0.5;

        this->AuxiliarInitializeBossak();
    }

    explicit StructuralMechanicsBossakScheme(StructuralMechanicsBossakScheme& rOther)
        :BossakBaseType(rOther)
    {}

    ~StructuralMechanicsBossakScheme() override
    {}

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    BaseTypePointer Clone() override
    {
        return BaseTypePointer(new StructuralMechanicsBossakScheme(*this));
    }

    void Initialize(ModelPart& rModelPart) override
    {
        // Call the base Initialize method
        BossakBaseType::Initialize(rModelPart);

        // Allocate the OSS projection variables
        const auto &r_process_info = rModelPart.GetProcessInfo();
        const bool oss_switch = r_process_info.Has(OSS_SWITCH) ? r_process_info[OSS_SWITCH] : false;
        if (oss_switch) {
            const array_1d<double,3> aux_zero = ZeroVector(3);
            block_for_each(rModelPart.Nodes(), aux_zero, [this](Node& rNode, array_1d<double,3>& rAuxZero){
                for (auto& rp_var : mProjectionScalarVariablesList) {
                    rNode.SetValue(*rp_var, 0.0);
                }
                for (auto& rp_var : mProjectionVectorVariablesList) {
                    rNode.SetValue(*rp_var, rAuxZero);
                }
            });
        }
    }

    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Call the base Predict method
        BossakBaseType::Predict(rModelPart, rDofSet, rA, rDx, rb);

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Predicting angular time derivatives if they are available (nodally for efficiency)
        if(rModelPart.HasNodalSolutionStepVariable(ROTATION) && rModelPart.Nodes().size() > 0)
        {
            const auto it_node_begin = rModelPart.Nodes().begin();

            // Getting position
            KRATOS_ERROR_IF_NOT(it_node_begin->HasDofFor(ROTATION_X)) << "StructuralMechanicsBossakScheme:: ROTATION is not added" << std::endl;
            const int rot_pos = it_node_begin->GetDofPosition(ROTATION_X);

            // Getting dimension
            const std::size_t dimension = r_current_process_info.Has(DOMAIN_SIZE) ? r_current_process_info.GetValue(DOMAIN_SIZE) : 3;

            //Auxiliar variable
            const std::array<const Variable<ComponentType>*, 3> rot_components = {&ROTATION_X, &ROTATION_Y, &ROTATION_Z};

            block_for_each(rModelPart.Nodes(), array_1d<double, 3>(), [&](Node& rNode, array_1d<double, 3>& rDeltaRotationTLS){
                //Predicting
                const array_1d<double, 3>& r_previous_rotation = rNode.FastGetSolutionStepValue(ROTATION, 1);
                const array_1d<double, 3>& r_previous_angular_velocity = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY, 1);
                const array_1d<double, 3>& r_previous_angular_acceleration = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION, 1);
                array_1d<double, 3>& r_current_rotation = rNode.FastGetSolutionStepValue(ROTATION);
                array_1d<double, 3>& r_current_angular_velocity = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3>& r_current_angular_acceleration = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION);

                for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                    if (!rNode.GetDof(*rot_components[i_dim], rot_pos + i_dim).IsFixed()) {
                        r_current_rotation[i_dim] = r_previous_rotation[i_dim] + delta_time * r_previous_angular_velocity[i_dim] + 0.5 * std::pow(delta_time, 2) * r_previous_angular_acceleration[i_dim];
                    }
                }

                // Updating
                noalias(rDeltaRotationTLS) = r_current_rotation - r_previous_rotation;
                UpdateAngularVelocity(r_current_angular_velocity, rDeltaRotationTLS, r_previous_angular_velocity, r_previous_angular_acceleration);
                UpdateAngularAcceleration(r_current_angular_acceleration, rDeltaRotationTLS, r_previous_angular_velocity, r_previous_angular_acceleration);

            });
        }

        KRATOS_CATCH("")
    }

    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Call the base Update method
        BossakBaseType::Update(rModelPart, rDofSet, rA, rDx, rb);

        // Updating angular time derivatives if they are available (nodally for efficiency)
        if(rModelPart.HasNodalSolutionStepVariable(ROTATION))
        {
            block_for_each(rModelPart.Nodes(), array_1d<double,3>(), [&](Node& rNode, array_1d<double,3>& rDeltaRotationTLS){
                    noalias(rDeltaRotationTLS) = rNode.FastGetSolutionStepValue(ROTATION) - rNode.FastGetSolutionStepValue(ROTATION, 1);

                    array_1d<double, 3>& r_current_angular_velocity = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY);
                    const array_1d<double, 3>& r_previous_angular_velocity = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY, 1);

                    array_1d<double, 3>& r_current_angular_acceleration = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION);
                    const array_1d<double, 3>& r_previous_angular_acceleration = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION, 1);

                    UpdateAngularVelocity(r_current_angular_velocity, rDeltaRotationTLS, r_previous_angular_velocity, r_previous_angular_acceleration);
                    UpdateAngularAcceleration(r_current_angular_acceleration, rDeltaRotationTLS, r_previous_angular_velocity, r_previous_angular_acceleration);
            });
        }

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        // Call the base InitializeSolutionStep method
        BossakBaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        // Update the NODAL_AREA (note that this strictly required only in the updated Lagrangian case)
        const auto& r_current_process_info = rModelPart.GetProcessInfo();
        const bool oss_switch = r_current_process_info.Has(OSS_SWITCH) ? r_current_process_info[OSS_SWITCH] : false;
        if (oss_switch && mUpdateNodalArea) {
            CalculateNodalAreaProcess<false>(rModelPart).Execute();
        }
    }

    void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY

        // Call base class FinalizeNonLinIteration
        BossakBaseType::FinalizeNonLinIteration(rModelPart, rA, rDx, rb);

        // Check if the Orthogonal SubScales (OSS) are active
        const auto& r_process_info = rModelPart.GetProcessInfo();
        const bool oss_switch = r_process_info.Has(OSS_SWITCH) ? r_process_info[OSS_SWITCH] : false;

        // Calculate the OSS projections
        if (oss_switch) {
            // Initialize the projection values
            const array_1d<double,3> aux_zero = ZeroVector(3);
            block_for_each(rModelPart.Nodes(), aux_zero, [this](Node& rNode, array_1d<double,3>& rAuxZero){
                for (auto& rp_var : mProjectionScalarVariablesList) {
                    rNode.FastGetSolutionStepValue(*rp_var) = 0.0;
                }
                for (auto& rp_var : mProjectionVectorVariablesList) {
                    rNode.FastGetSolutionStepValue(*rp_var) = rAuxZero;
                }
            });

            // Calculate the element residuals projection
            std::tuple<double, array_1d<double,3>> oss_proj_tls;
            block_for_each(rModelPart.Elements(), oss_proj_tls, [&, this](Element& rElement, std::tuple<double, array_1d<double,3>>& rOssProjTLS){
                double& r_scalar_proj = std::get<0>(rOssProjTLS);
                array_1d<double,3>& r_vector_proj = std::get<1>(rOssProjTLS);
                for (auto& rp_var : mProjectionScalarVariablesList) {
                    rElement.Calculate(*rp_var, r_scalar_proj, r_process_info);
                }
                for (auto& rp_var : mProjectionVectorVariablesList) {
                    rElement.Calculate(*rp_var, r_vector_proj, r_process_info);
                }
            });

            // Do the nodal weighting
            //TODO: We the weighted L2 projection with the density should be done in the multimaterial case
            //TODO: Note that this would require customizing the weights calculation in here (standard NODAL_AREA cannot be used)
            block_for_each(rModelPart.Nodes(), [this](Node& rNode){
                const double nodal_area = rNode.GetValue(NODAL_AREA);
                for (auto& rp_var : mProjectionScalarVariablesList) {
                    rNode.FastGetSolutionStepValue(*rp_var) /= nodal_area;
                }
                for (auto& rp_var : mProjectionVectorVariablesList) {
                    rNode.FastGetSolutionStepValue(*rp_var) /= nodal_area;
                }
            });
        }

        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"({
            "name" : "structural_mechanics_bossak_scheme",
            "update_nodal_area" : true,
            "projection_variables_list" : []
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BossakBaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /// @brief Return the name of the class as used in the settings (@a snake_case).
    static std::string Name()
    {
        return "structural_mechanics_bossak_scheme";
    }

    ///@}
    ///@name Input and output
    ///@{

    /// @brief Return information as a string.
    std::string Info() const override
    {
        return "StructuralMechanicsBossakScheme";
    }

    /// @brief Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// @brief Print the instance's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        BossakBaseType::AssignSettings(ThisParameters);

        mUpdateNodalArea = ThisParameters["update_nodal_area"].GetBool();
        const auto& r_proj_var_list = ThisParameters["projection_variables_list"].GetStringArray();
        for (const std::string& r_var_name : r_proj_var_list) {
            if (KratosComponents<Variable<double>>::Has(r_var_name)) {
                mProjectionScalarVariablesList.push_back(&KratosComponents<Variable<double>>::Get(r_var_name));
            } else if (KratosComponents<Variable<array_1d<double,3>>>::Has(r_var_name)) {
                mProjectionVectorVariablesList.push_back(&KratosComponents<Variable<array_1d<double, 3>>>::Get(r_var_name));
            } else {
                KRATOS_ERROR << "Wrong projection variable '" << r_var_name << "'." << std::endl;
            }
        }
    }

    inline void UpdateAngularVelocity(
        array_1d<double, 3>& rCurrentAngularVelocity,
        const array_1d<double, 3>& rDeltaRotation,
        const array_1d<double, 3>& rPreviousAngularVelocity,
        const array_1d<double, 3>& rPreviousAngularAcceleration
        )
    {
        noalias(rCurrentAngularVelocity) = (this->mBossak.c1 * rDeltaRotation - this->mBossak.c4 * rPreviousAngularVelocity - this->mBossak.c5 * rPreviousAngularAcceleration);
    }

    inline void UpdateAngularAcceleration(
        array_1d<double, 3>& rCurrentAngularAcceleration,
        const array_1d<double, 3>& rDeltaRotation,
        const array_1d<double, 3>& rPreviousAngularVelocity,
        const array_1d<double, 3>& rPreviousAngularAcceleration
        )
    {
        noalias(rCurrentAngularAcceleration) = (this->mBossak.c0 * rDeltaRotation - this->mBossak.c2 * rPreviousAngularVelocity - this->mBossak.c3 * rPreviousAngularAcceleration);
    }

    ///@}
private:
    ///@name Private Member Variables
    ///@{

    bool mUpdateNodalArea;

    std::vector<const Variable<double>*> mProjectionScalarVariablesList;

    std::vector<const Variable<array_1d<double,3>>*> mProjectionVectorVariablesList;

    ///@}
}; // Class StructuralMechanicsBossakScheme

///@}

} // namespace Kratos
