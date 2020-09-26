//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


#if !defined(KRATOS_RESIDUAL_BASED_BDF_CUSTOM_SCHEME )
#define  KRATOS_RESIDUAL_BASED_BDF_CUSTOM_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/residual_based_bdf_scheme.h"
#include "includes/variables.h"
#include "includes/kratos_parameters.h"
#include "includes/checks.h"

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
 * @class ResidualBasedBDFCustomScheme
 * @ingroup KratosCore
 * @brief BDF integration scheme (for dynamic problems)
 * @details The second order Backward Differentiation Formula (BDF) method is a two step second order accurate method.
 * This scheme is a generalization of the only displacement scheme, where any list of variables and its derivatives can be considered instead
 * Look at the base class for more details
 * @see ResidualBasedBDFScheme
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,  class TDenseSpace>
class ResidualBasedBDFCustomScheme
    : public ResidualBasedBDFScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBDFCustomScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                  BaseType;

    typedef ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace> ImplicitBaseType;

    typedef ResidualBasedBDFScheme<TSparseSpace,TDenseSpace>               BDFBaseType;

    typedef ResidualBasedBDFCustomScheme<TSparseSpace, TDenseSpace>          ClassType;

    typedef typename ImplicitBaseType::TDataType                             TDataType;

    typedef typename ImplicitBaseType::DofsArrayType                     DofsArrayType;

    typedef typename Element::DofsVectorType                            DofsVectorType;

    typedef typename ImplicitBaseType::TSystemMatrixType             TSystemMatrixType;

    typedef typename ImplicitBaseType::TSystemVectorType             TSystemVectorType;

    typedef typename ImplicitBaseType::LocalSystemVectorType     LocalSystemVectorType;

    typedef typename ImplicitBaseType::LocalSystemMatrixType     LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                               NodesArrayType;

    typedef ModelPart::ElementsContainerType                         ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                     ConditionsArrayType;

    typedef typename BaseType::Pointer                                 BaseTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. The BDF method
     * @param ThisParameters The parameters containing the list of variables to consider
     * @todo The ideal would be to use directly the dof or the variable itself to identify the type of variable and is derivatives
     */
    explicit ResidualBasedBDFCustomScheme(Parameters ThisParameters)
        :BDFBaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // Creating variables list
        CreateVariablesList(ThisParameters);
    }

    /**
     * @brief Constructor. The BDF method
     * @param Order The integration order
     * @param ThisParameters The parameters containing the list of variables to consider
     * @todo The ideal would be to use directly the dof or the variable itself to identify the type of variable and is derivatives
     */
    explicit ResidualBasedBDFCustomScheme(
        const std::size_t Order = 2,
        Parameters ThisParameters =  Parameters(R"({})")
        )
        :BDFBaseType(Order)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // Creating variables list
        CreateVariablesList(ThisParameters);
    }

    /** Copy Constructor.
     */
    explicit ResidualBasedBDFCustomScheme(ResidualBasedBDFCustomScheme& rOther)
        :BDFBaseType(rOther)
        ,mDoubleVariable(rOther.mDoubleVariable)
        ,mFirstDoubleDerivatives(rOther.mFirstDoubleDerivatives)
        ,mSecondDoubleDerivatives(rOther.mSecondDoubleDerivatives)
    {
    }

    /**
     * Clone
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedBDFCustomScheme(*this) );
    }

    /** Destructor.
     */
    ~ResidualBasedBDFCustomScheme
    () override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief This is the place to initialize the Scheme.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        BDFBaseType::Initialize(rModelPart);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Getting dimension
        KRATOS_WARNING_IF("ResidualBasedBDFCustomScheme", !r_current_process_info.Has(DOMAIN_SIZE)) << "DOMAIN_SIZE not defined. Please define DOMAIN_SIZE. 3D case will be assumed" << std::endl;
        const std::size_t domain_size = r_current_process_info.Has(DOMAIN_SIZE) ? r_current_process_info.GetValue(DOMAIN_SIZE) : 3;
        if (domain_size != mDomainSize) {
            const std::size_t total_number_of_variables = mDoubleVariable.size();

            // We remove the third component
            if (domain_size == 2) {
                const std::size_t number_variables_added = total_number_of_variables/3;
                for (std::size_t i = 0; i < number_variables_added; ++i) {
                    mDoubleVariable.erase(mDoubleVariable.begin() + (2 + 2 * i));
                    mFirstDoubleDerivatives.erase(mDoubleVariable.begin() + (2 + 2 * i));
                    mSecondDoubleDerivatives.erase(mDoubleVariable.begin() + (2 + 2 * i));
                }
            } else if (domain_size == 3) { // We need to add the third component
                const std::size_t number_variables_added = total_number_of_variables/2;
                for (std::size_t i = 0; i < number_variables_added; ++i) {
                    const std::string variable_name = ((*(mDoubleVariable.begin() + 2 * i))->GetSourceVariable()).Name();
                    const auto& r_var_z = KratosComponents<Variable<double>>::Get(variable_name + "_Z");
                    mDoubleVariable.push_back(&r_var_z);
                    mFirstDoubleDerivatives.push_back(&(r_var_z.GetTimeDerivative()));
                    mSecondDoubleDerivatives.push_back(&((r_var_z.GetTimeDerivative()).GetTimeDerivative()));
                }
            } else {
                KRATOS_ERROR << "DOMAIN_SIZE can onbly be 2 or 3. It is: " << domain_size << std::endl;
            }

            mDomainSize = domain_size;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart The model part of the problem to solve
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        BDFBaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        const auto it_node_begin = rModelPart.Nodes().begin();

        // Auxiliar fixed value
        bool fixed = false;

        #pragma omp parallel for private(fixed)
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            std::size_t counter = 0;
            for (auto p_var : mDoubleVariable) {

                fixed = false;

                // Derivatives
                const auto& dvar = *mFirstDoubleDerivatives[counter];
                const auto& d2var = *mSecondDoubleDerivatives[counter];

                if (it_node->HasDofFor(d2var)) {
                    if (it_node->IsFixed(d2var)) {
                        it_node->Fix(*p_var);
                        fixed = true;
                    }
                }

                if (it_node->HasDofFor(dvar)) {
                    if (it_node->IsFixed(dvar) && !fixed) {
                        it_node->Fix(*p_var);
                    }
                }

                counter++;
            }
        }

        KRATOS_CATCH("ResidualBasedBDFCustomScheme.InitializeSolutionStep");
    }

    /**
     * @brief Performing the prediction of the solution
     * @details It predicts the solution for the current step x = xold + vold * Dt
     * @param rModelPart The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        // Getting process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Getting delta time
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );

        // Getting first node iterator
        const auto it_node_begin = rModelPart.Nodes().begin();

        #pragma omp parallel for
        for(int i = 0;  i< num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            std::size_t counter = 0;
            for (auto p_var : mDoubleVariable) {
                // Derivatives
                const auto& dvar = *mFirstDoubleDerivatives[counter];
                const auto& d2var = *mSecondDoubleDerivatives[counter];

                ComputePredictComponent(it_node, *p_var, dvar, d2var, delta_time);

                counter++;
            }

            // Updating time derivatives
            UpdateFirstDerivative(it_node);
            UpdateSecondDerivative(it_node);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided.
     * @details Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */

    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY;

        const int err = BDFBaseType::Check(rModelPart);
        if(err!=0) return err;

        // Check for variables keys
        // Verify that the variables are correctly initialized
        for ( auto p_var : mDoubleVariable)
            KRATOS_CHECK_VARIABLE_KEY((*p_var))
        for ( auto p_var : mFirstDoubleDerivatives)
            KRATOS_CHECK_VARIABLE_KEY((*p_var))
        for ( auto p_var : mSecondDoubleDerivatives)
            KRATOS_CHECK_VARIABLE_KEY((*p_var))

        // Check that variables are correctly allocated
        for(auto& r_node : rModelPart.Nodes()) {
            for ( auto p_var : mDoubleVariable)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_var), r_node)
            for ( auto p_var : mFirstDoubleDerivatives)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_var), r_node)
            for ( auto p_var : mSecondDoubleDerivatives)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA((*p_var), r_node)

            for ( auto p_var : mDoubleVariable)
                KRATOS_CHECK_DOF_IN_NODE((*p_var), r_node)
        }

        KRATOS_CATCH( "" );

        return 0;
    }

        /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"               : "bdf_scheme",
            "domain_size"        : 3,
            "integration_order"  : 2,
            "solution_variables" : ["DISPLACEMENT"]
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BDFBaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "bdf_scheme";
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
        return "ResidualBasedBDFCustomScheme";
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

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<const Variable<double>*> mDoubleVariable;          /// The double variables
    std::vector<const Variable<double>*> mFirstDoubleDerivatives;   /// The first derivative double variable to compute
    std::vector<const Variable<double>*> mSecondDoubleDerivatives; /// The second derivative double variable to compute

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Updating first time derivative (velocity)
     * @param itNode the node interator
     */
    inline void UpdateFirstDerivative(NodesArrayType::iterator itNode) override
    {
        // DOUBLES
        std::size_t counter = 0;
        for (auto p_var : mDoubleVariable) {
            double& dotun0 = itNode->FastGetSolutionStepValue(*mFirstDoubleDerivatives[counter]);
            dotun0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(*p_var);
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                dotun0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(*p_var, i_order);
            counter++;
        }
    }

    /**
     * @brief Updating second time derivative (acceleration)
     * @param itNode the node interator
     */
    inline void UpdateSecondDerivative(NodesArrayType::iterator itNode) override
    {
        // DOUBLES
        std::size_t counter = 0;
        for (auto p_var : mFirstDoubleDerivatives) {
            double& dot2un0 = itNode->FastGetSolutionStepValue(*mSecondDoubleDerivatives[counter]);
            dot2un0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(*p_var);
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                dot2un0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(*p_var, i_order);
            counter++;
        }
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

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::size_t mDomainSize = 3; /// This auxiliar variable is used to store the domain size of the problem

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method reduces the code duplication for each components when computing the prediction
     * @param itNode The node iterator of the node currently being computed
     * @param rVariable The variable currently being integrated
     * @param rDerivedVariable The first time derivative of the current variable
     * @param rDerived2Variable The second time derivative of the current variable
     * @param DeltaTime The increment of time for the time integration
     */
    template<class TClassVar>
    void ComputePredictComponent(
        NodesArrayType::iterator itNode,
        const TClassVar& rVariable,
        const TClassVar& rDerivedVariable,
        const TClassVar& rDerived2Variable,
        const double DeltaTime
        )
    {
        // Values
        const double dot2un1 = itNode->FastGetSolutionStepValue(rDerived2Variable, 1);
        const double dotun1 = itNode->FastGetSolutionStepValue(rDerivedVariable, 1);
        const double un1 = itNode->FastGetSolutionStepValue(rVariable, 1);
        const double dot2un0 = itNode->FastGetSolutionStepValue(rDerived2Variable);
        double& dotun0 = itNode->FastGetSolutionStepValue(rDerivedVariable);
        double& un0 = itNode->FastGetSolutionStepValue(rVariable);

        if (itNode->HasDofFor(rDerived2Variable) && itNode->IsFixed(rDerived2Variable)) {
            dotun0 = dot2un0;
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                dotun0 -= BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(rDerivedVariable, i_order);
            dotun0 /= BDFBaseType::mBDF[0];

            un0 = dotun0;
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                un0 -= BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(rVariable, i_order);
            un0 /= BDFBaseType::mBDF[0];
        } else if (itNode->HasDofFor(rDerivedVariable) && itNode->IsFixed(rDerivedVariable)) {
            un0 = dotun0;
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                un0 -= BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(rVariable, i_order);
            un0 /= BDFBaseType::mBDF[0];
        } else if (!itNode->IsFixed(rVariable)) {
            un0 = un1 + DeltaTime * dotun1 + 0.5 * std::pow(DeltaTime, 2) * dot2un1;
        }
    }

    /**
     * @brief This method creates the list of variables
     * @param ThisParameters The configuration parameters
     */
    void CreateVariablesList(Parameters ThisParameters)
    {
        const std::size_t n_variables = ThisParameters["solution_variables"].size();

        // The current dimension
        mDomainSize = ThisParameters["domain_size"].GetInt();

        const auto variable_names = ThisParameters["solution_variables"].GetStringArray();

        for (std::size_t p_var = 0; p_var < n_variables; ++p_var){
            const std::string& variable_name = variable_names[p_var];

            if(KratosComponents<Variable<double>>::Has(variable_name)){
                const auto& r_var = KratosComponents<Variable<double>>::Get(variable_name);
                mDoubleVariable.push_back(&r_var);
                mFirstDoubleDerivatives.push_back(&(r_var.GetTimeDerivative()));
                mSecondDoubleDerivatives.push_back(&((r_var.GetTimeDerivative()).GetTimeDerivative()));
            } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(variable_name)) {
                // Components
                const auto& r_var_x = KratosComponents<Variable<double>>::Get(variable_name+"_X");
                const auto& r_var_y = KratosComponents<Variable<double>>::Get(variable_name+"_Y");
                mDoubleVariable.push_back(&r_var_x);
                mDoubleVariable.push_back(&r_var_y);
                mFirstDoubleDerivatives.push_back(&(r_var_x.GetTimeDerivative()));
                mFirstDoubleDerivatives.push_back(&(r_var_y.GetTimeDerivative()));
                mSecondDoubleDerivatives.push_back(&((r_var_x.GetTimeDerivative()).GetTimeDerivative()));
                mSecondDoubleDerivatives.push_back(&((r_var_y.GetTimeDerivative()).GetTimeDerivative()));
                if (mDomainSize == 3) {
                    const auto& r_var_z = KratosComponents<Variable<double>>::Get(variable_name+"_Z");
                    mDoubleVariable.push_back(&r_var_z);
                    mFirstDoubleDerivatives.push_back(&(r_var_z.GetTimeDerivative()));
                    mSecondDoubleDerivatives.push_back(&((r_var_z.GetTimeDerivative()).GetTimeDerivative()));
                }
            } else {
                KRATOS_ERROR << "Only double and vector variables are allowed in the variables list." ;
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
}; /* Class ResidualBasedBDFCustomScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BDF_CUSTOM_SCHEME defined */
