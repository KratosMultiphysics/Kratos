//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_SHALLOW_WATER_RESIDUAL_BASED_BDF_SCHEME_H_INCLUDED
#define KRATOS_SHALLOW_WATER_RESIDUAL_BASED_BDF_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "shallow_water_application_variables.h"
#include "custom_utilities/flow_rate_slip_utility.h"
#include "solving_strategies/schemes/residual_based_bdf_scheme.h"
#include "processes/calculate_nodal_area_process.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

/**
 * @class ShallowWaterResidualBasedBDFScheme
 * @ingroup KratosShallowWaterApplication
 * @brief BDF integration scheme (for dynamic problems)
 * @details The \f$n\f$ order Backward Differentiation Formula (BDF) method is a two step \f$n\f$ order accurate method.
 * This scheme is designed to solve a system of the type:
 * \f[
 *   \mathbf{M} \frac{du_{n0}}{dt} + \mathbf{K} u_{n0} = \mathbf{f}_{ext}
 * \f]
 * @author Miguel Maso Sotomayor
 */
template<class TSparseSpace,  class TDenseSpace>
class ShallowWaterResidualBasedBDFScheme
    : public ResidualBasedBDFScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ShallowWaterResidualBasedBDFScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                             BaseType;

    typedef ResidualBasedBDFScheme<TSparseSpace,TDenseSpace>          BDFBaseType;

    typedef typename BDFBaseType::DofsArrayType                     DofsArrayType;

    typedef typename BDFBaseType::TSystemMatrixType             TSystemMatrixType;

    typedef typename BDFBaseType::TSystemVectorType             TSystemVectorType;

    typedef typename BDFBaseType::LocalSystemVectorType     LocalSystemVectorType;

    typedef typename BDFBaseType::LocalSystemMatrixType     LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                          NodesArrayType;

    typedef typename ModelPart::NodeType                                 NodeType;

    typedef FlowRateSlipUtility<LocalSystemMatrixType,LocalSystemVectorType,double>FlowRateSlipToolType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor
    explicit ShallowWaterResidualBasedBDFScheme(const std::size_t Order = 2, bool UpdateVelocities = false)
        : BDFBaseType(Order)
        , mRotationTool()
        , mUpdateVelocities(UpdateVelocities)
        , mProjectDispersiveField(false)
    {
        mVariables = {&MOMENTUM_X, &MOMENTUM_Y, &HEIGHT};
        mDerivativeVariables = {&ACCELERATION_X, &ACCELERATION_Y, &VERTICAL_VELOCITY};
    }

    // Constructor with parameters
    explicit ShallowWaterResidualBasedBDFScheme(Parameters ThisParameters)
        : BDFBaseType(ThisParameters["integration_order"].GetDouble())
        , mRotationTool()
    {
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    // Copy Constructor
    explicit ShallowWaterResidualBasedBDFScheme(ShallowWaterResidualBasedBDFScheme& rOther)
        : BDFBaseType(rOther)
        , mRotationTool()
        , mUpdateVelocities(rOther.mUpdateVelocities)
        , mProjectDispersiveField(rOther.mProjectDispersiveField)
        , mVariables(rOther.mVariables)
        , mDerivativeVariables(rOther.mDerivativeVariables)
    {}

    /**
     * Clone
     */
    typename BaseType::Pointer Clone() override
    {
        return typename BaseType::Pointer(new ShallowWaterResidualBasedBDFScheme(*this));
    }

    // Destructor
    ~ShallowWaterResidualBasedBDFScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Performing the update of the solution within newton iteration
     * @param rModelPart The model of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param rA LHS matrix
     * @param rDx incremental update of primary variables
     * @param rb RHS Vector
     */
    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        mRotationTool.RotateVelocities(rModelPart);

        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        mRotationTool.RecoverVelocities(rModelPart);

        BDFBaseType::UpdateDerivatives(rModelPart, rDofSet, rA, rDx, rb);

        if (mUpdateVelocities) UpdateVelocities(rModelPart);

        KRATOS_CATCH("ShallowWaterResidualBasedBDFScheme.Update");
    }

    /**
     * @brief Performing the prediction of the solution
     * @details It predicts the solution for the current step
     * @param rModelPart The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     */
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        const double delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];

        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        const auto it_node_begin = rModelPart.Nodes().begin();

        IndexPartition<std::size_t>(num_nodes).for_each([&](std::size_t i){
            auto it_node = it_node_begin + i;

            for (std::size_t j = 0; j < 3; ++j)
            {
                if (!it_node->IsFixed(*mVariables[j])) {
                    double& un0 = it_node->FastGetSolutionStepValue(*mVariables[j]);
                    double un1 = it_node->FastGetSolutionStepValue(*mVariables[j], 1);
                    double dot_un1 = it_node->FastGetSolutionStepValue(*mDerivativeVariables[j], 1);
                    un0 = un1 + delta_time * dot_un1;
                }
            }

            UpdateFirstDerivative(it_node);
        });

        KRATOS_CATCH("ShallowWaterResidualBasedBDFScheme.Predict");
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rCurrentElement The element to compute
     * @param rLHS_Contribution The LHS matrix contribution
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BDFBaseType::CalculateSystemContributions(
            rCurrentElement,
            rLHS_Contribution,
            rRHS_Contribution,
            rEquationId,
            rCurrentProcessInfo);

        mRotationTool.Rotate(rLHS_Contribution,rRHS_Contribution,rCurrentElement.GetGeometry());
        mRotationTool.ApplySlipCondition(rLHS_Contribution,rRHS_Contribution,rCurrentElement.GetGeometry());
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentElement The element to compute
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BDFBaseType::CalculateRHSContribution(
            rCurrentElement,
            rRHS_Contribution,
            rEquationId,
            rCurrentProcessInfo);

        mRotationTool.Rotate(rRHS_Contribution,rCurrentElement.GetGeometry());
        mRotationTool.ApplySlipCondition(rRHS_Contribution,rCurrentElement.GetGeometry());
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rCurrentCondition The condition to compute
     * @param rLHS_Contribution The LHS matrix contribution
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BDFBaseType::CalculateSystemContributions(
            rCurrentCondition,
            rLHS_Contribution,
            rRHS_Contribution,
            rEquationId,
            rCurrentProcessInfo);

        mRotationTool.Rotate(rLHS_Contribution,rRHS_Contribution,rCurrentCondition.GetGeometry());
        mRotationTool.ApplySlipCondition(rLHS_Contribution,rRHS_Contribution,rCurrentCondition.GetGeometry());
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentCondition The condition to compute
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        BDFBaseType::CalculateRHSContribution(
            rCurrentCondition,
            rRHS_Contribution,
            rEquationId,
            rCurrentProcessInfo);

        mRotationTool.Rotate(rRHS_Contribution,rCurrentCondition.GetGeometry());
        mRotationTool.ApplySlipCondition(rRHS_Contribution,rCurrentCondition.GetGeometry());
    }

    /**
     * @brief Initialize the nodal area and the derivatives recovery.
     * @details The nodal area is used to apply the explicit prediction.
     * @param rModelPart The model part of the problem to solve
     */
    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);
        if (mProjectDispersiveField) {
            CalculateNodalAreaProcess<true>(rModelPart).Execute();
        }
    }

    /**
     * @brief Calculate the global projection of the dispersive field.
     * @param rModelPart The model part of the problem to solve
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     */
    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    ) override
    {
        if (mProjectDispersiveField) {
            const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
            block_for_each(rModelPart.Nodes(), [](NodeType& rNode){
                rNode.FastGetSolutionStepValue(DISPERSION_H) = ZeroVector(3);
                rNode.FastGetSolutionStepValue(DISPERSION_V) = ZeroVector(3);
            });
            block_for_each(rModelPart.Elements(), [&](Element& rElement){
                rElement.InitializeNonLinearIteration(r_process_info);
            });
            block_for_each(rModelPart.Conditions(), [&](Condition& rCondition){
                rCondition.InitializeNonLinearIteration(r_process_info);
            });
            block_for_each(rModelPart.Nodes(), [](NodeType& rNode){
                const double nodal_area = rNode.FastGetSolutionStepValue(NODAL_AREA);
                rNode.FastGetSolutionStepValue(DISPERSION_H) /= nodal_area;
                rNode.FastGetSolutionStepValue(DISPERSION_V) /= nodal_area;
            });
            ApplyLaplacianBoundaryConditions(rModelPart);
        }
    }

    /**
     * @brief Free memory allocated by this class.
     */
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                     : "shallow_water_residual_based_bdf_scheme",
            "integration_order"        : 2,
            "update_velocities"        : false,
            "project_dispersive_field" : false,
            "solution_variables"       : ["MOMENTUM","HEIGHT"]
        })" );

        // Getting base class default parameters
        const Parameters base_default_parameters = BDFBaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
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
        return "ShallowWaterResidualBasedBDFScheme";
    }

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected Operations
    ///@{

    /**
     * @brief This method validate and assign default parameters
     * @param rParameters Parameters to be validated
     * @param DefaultParameters The default parameters
     * @return Returns validated Parameters
     */
    Parameters ValidateAndAssignParameters(
        Parameters ThisParameters,
        const Parameters DefaultParameters
        ) const override
    {
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        return ThisParameters;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        mUpdateVelocities = ThisParameters["update_velocities"].GetBool();
        mProjectDispersiveField = ThisParameters["project_dispersive_field"].GetBool();

        const auto variable_names = ThisParameters["solution_variables"].GetStringArray();
        for (std::string variable_name : variable_names)
        {
            if (KratosComponents<Variable<double>>::Has(variable_name))
            {
                const auto& r_var = KratosComponents<Variable<double>>::Get(variable_name);
                mVariables.push_back(&r_var);
            }
            else if (KratosComponents<Variable<array_1d<double,3>>>::Has(variable_name)) {
                const auto& r_var_x = KratosComponents<Variable<double>>::Get(variable_name+"_X");
                const auto& r_var_y = KratosComponents<Variable<double>>::Get(variable_name+"_Y");
                mVariables.push_back(&r_var_x);
                mVariables.push_back(&r_var_y);
            }
            else
            {
                KRATOS_ERROR << "Only double and component variables are allowed in the solution variables list." ;
            }
        }
        mDerivativeVariables.push_back(&ACCELERATION_X);
        mDerivativeVariables.push_back(&ACCELERATION_Y);
        mDerivativeVariables.push_back(&VERTICAL_VELOCITY);
    }

    ///@}
    ///@name Protected member Variables
    ///@{

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    FlowRateSlipToolType mRotationTool;

    bool mUpdateVelocities;
    bool mProjectDispersiveField;

    std::vector<const Variable<double>*> mVariables;
    std::vector<const Variable<double>*> mDerivativeVariables;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Updating first time derivative
     * @param itNode the node interator
     */
    void UpdateFirstDerivative(NodesArrayType::iterator itNode) override
    {
        for (std::size_t i_var = 0; i_var < 3; ++i_var)
        {
            double& dot_un0 = itNode->FastGetSolutionStepValue(*mDerivativeVariables[i_var]);
            dot_un0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(*mVariables[i_var]);

            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
            {
                dot_un0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(*mVariables[i_var], i_order);
            }
        }
    }

    /**
     * @brief Updating second time derivative
     * @param itNode the node interator
     */
    void UpdateSecondDerivative(NodesArrayType::iterator itNode) override {}

    /**
     * @brief Updating the velocities
     * @param rModelPart The model part to compute
     */
    void UpdateVelocities(ModelPart& rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [&](NodeType& r_node){
            auto& vel = r_node.FastGetSolutionStepValue(VELOCITY);
            const auto& q = r_node.FastGetSolutionStepValue(MOMENTUM);
            const auto& h = r_node.FastGetSolutionStepValue(HEIGHT);
            vel = q / h;
        });
    }

    /**
     * @brief It adds the dynamic LHS contribution of the elements
     * @param rLHS_Contribution The dynamic contribution for the LHS
     * @param rD The damping matrix
     * @param rM The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding mass contribution to the dynamic stiffness
        if (rM.size1() != 0) { // if M matrix declared
            noalias(rLHS_Contribution) += rM * BDFBaseType::mBDF[0];
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements
     * @param rElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (rM.size1() != 0) {
            rElement.GetFirstDerivativesVector(BDFBaseType::mVector.dotun0[this_thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, BDFBaseType::mVector.dotun0[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition
     * @param rCondition The condition to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Condition& rCondition,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (rM.size1() != 0) {
            rCondition.GetFirstDerivativesVector(BDFBaseType::mVector.dotun0[this_thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, BDFBaseType::mVector.dotun0[this_thread]);
        }
    }


    void ApplyLaplacianBoundaryConditions(ModelPart& rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [](NodeType& rNode){
            if (rNode.IsFixed(VELOCITY_X)) {
                rNode.FastGetSolutionStepValue(DISPERSION_H_X) = 0.0;
                rNode.FastGetSolutionStepValue(DISPERSION_V_X) = 0.0;
            }
            if (rNode.IsFixed(VELOCITY_Y)) {
                rNode.FastGetSolutionStepValue(DISPERSION_H_Y) = 0.0;
                rNode.FastGetSolutionStepValue(DISPERSION_V_Y) = 0.0;
            }
        });
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
    ///@{

}; // Class ShallowWaterResidualBasedBDFScheme

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // Namespace Kratos

#endif // KRATOS_SHALLOW_WATER_RESIDUAL_BASED_BDF_SCHEME_H_INCLUDED defined
