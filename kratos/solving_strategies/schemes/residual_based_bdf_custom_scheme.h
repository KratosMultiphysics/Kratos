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

    typedef VectorComponentAdaptor< array_1d< double, 3 > >              ComponentType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. The BDF method
     * @param ThisParameters The parameters containing the list of variables to consider
     * @todo The ideal would be to use directly the dof or the variable itself to identify the type of variable and is derivatives
     */
    explicit ResidualBasedBDFCustomScheme(Parameters ThisParameters)
    {
        // Getting default parameters
        Parameters default_parameters = GetDefaultParameters();
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // Now here call the base class constructor
        BDFBaseType( ThisParameters["integration_order"].GetInt());

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
        // Getting default parameters
        Parameters default_parameters = GetDefaultParameters();
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // Creating variables list
        CreateVariablesList(ThisParameters);
    }

    /** Copy Constructor.
     */
    explicit ResidualBasedBDFCustomScheme(ResidualBasedBDFCustomScheme& rOther)
        :BDFBaseType(rOther)
        ,mDoubleVariable(rOther.mDoubleVariable)
        ,mFirtsDoubleDerivatives(rOther.mFirtsDoubleDerivatives)
        ,mSecondDoubleDerivatives(rOther.mSecondDoubleDerivatives)
        ,mArrayVariable(rOther.mArrayVariable)
        ,mFirtsArrayDerivatives(rOther.mFirtsArrayDerivatives)
        ,mSecondArrayDerivatives(rOther.mSecondArrayDerivatives)
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

            //ATTENTION::: the prediction is performed only on free nodes
            std::size_t counter = 0;
            for ( const auto& var : mDoubleVariable) {
                // Derivatives
                const auto& dvar = mFirtsDoubleDerivatives[counter];
                const auto& d2var = mSecondDoubleDerivatives[counter];

                ComputePredictComponent(it_node, var, dvar, d2var, delta_time);

                counter++;
            }
            counter = 0;
            for ( const auto& var : mArrayVariable) {
                // Derivatives
                const auto& dvar = mFirtsArrayDerivatives[counter];
                const auto& d2var = mSecondArrayDerivatives[counter];

                ComputePredictComponent(it_node, var, dvar, d2var, delta_time);
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

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const int err = BDFBaseType::Check(rModelPart);
        if(err!=0) return err;

        // Check for variables keys
        // Verify that the variables are correctly initialized
        for ( const auto& i_var : mDoubleVariable)
            KRATOS_CHECK_VARIABLE_KEY(i_var)
        for ( const auto& i_var : mFirtsDoubleDerivatives)
            KRATOS_CHECK_VARIABLE_KEY(i_var)
        for ( const auto& i_var : mSecondDoubleDerivatives)
            KRATOS_CHECK_VARIABLE_KEY(i_var)
        for ( const auto& i_var : mArrayVariable)
            KRATOS_CHECK_VARIABLE_KEY(i_var)
        for ( const auto& i_var : mFirtsArrayDerivatives)
            KRATOS_CHECK_VARIABLE_KEY(i_var)
        for ( const auto& i_var : mSecondArrayDerivatives)
            KRATOS_CHECK_VARIABLE_KEY(i_var)

        // Check that variables are correctly allocated
        for(auto& r_node : rModelPart.Nodes()) {
            for ( const auto& i_var : mDoubleVariable)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(i_var, r_node)
            for ( const auto& i_var : mFirtsDoubleDerivatives)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(i_var, r_node)
            for ( const auto& i_var : mSecondDoubleDerivatives)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(i_var, r_node)
            for ( const auto& i_var : mArrayVariable)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(i_var, r_node)
            for ( const auto& i_var : mFirtsArrayDerivatives)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(i_var, r_node)
            for ( const auto& i_var : mSecondArrayDerivatives)
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(i_var, r_node)

            for ( const auto& i_var : mDoubleVariable)
                KRATOS_CHECK_DOF_IN_NODE(i_var, r_node)

            for ( const auto& i_var : mArrayVariable) {
                KRATOS_CHECK_DOF_IN_NODE(i_var, r_node)
            }
        }

        KRATOS_CATCH( "" );

        return 0;
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

    std::vector<Variable<double>> mDoubleVariable;                         /// The double variables
    std::vector<Variable<double>> mFirtsDoubleDerivatives;                 /// The first derivative double variable to compute
    std::vector<Variable<double>> mSecondDoubleDerivatives;                /// The second derivative double variable to compute
    std::vector<VariableComponent<ComponentType>> mArrayVariable;          /// The array variables to compute
    std::vector<VariableComponent<ComponentType>> mFirtsArrayDerivatives;  /// The first derivative array variable to compute
    std::vector<VariableComponent<ComponentType>> mSecondArrayDerivatives; /// The second derivative array variable to compute

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
        for ( const auto& i_var : mDoubleVariable) {
            if (!itNode->IsFixed(mFirtsDoubleDerivatives[counter])) {
                double& dotun0 = itNode->FastGetSolutionStepValue(mFirtsDoubleDerivatives[counter]);
                dotun0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(i_var);
                for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                    dotun0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(i_var, i_order);
            }
            counter++;
        }

        // ARRAYS
        counter = 0;
        for ( const auto& i_var : mArrayVariable) {
            if (!itNode->IsFixed(mFirtsArrayDerivatives[counter])) {
                double& dotun0 = itNode->FastGetSolutionStepValue(mFirtsArrayDerivatives[counter]);
                dotun0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(i_var);
                for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                    dotun0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(i_var, i_order);
            }
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
        for ( const auto& i_var : mFirtsDoubleDerivatives) {
            if (!itNode->IsFixed(mSecondDoubleDerivatives[counter])) {
                double& dot2un0 = itNode->FastGetSolutionStepValue(mSecondDoubleDerivatives[counter]);
                dot2un0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(i_var);
                for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                    dot2un0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(i_var, i_order);
            }
            counter++;
        }

        // ARRAYS
        counter = 0;
        for ( const auto& i_var : mFirtsArrayDerivatives) {
            if (!itNode->IsFixed(mSecondArrayDerivatives[counter])) {
                double& dot2un0 = itNode->FastGetSolutionStepValue(mSecondArrayDerivatives[counter]);
                dot2un0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(i_var);
                for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                    dot2un0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(i_var, i_order);
            }
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method reduces the code duplication for each components when computing the prediction
     * @param itNode The node iterator of the node currently being computed
     * @param iVar The variable currently being integrated
     * @param DerivedVariable The first time derivative of the current variable
     * @param Derived2Variable The second time derivative of the current variable
     * @param DeltaTime The increment of time for the time integration
     */
    template<class TClassVar>
    void ComputePredictComponent(
        NodesArrayType::iterator itNode,
        const TClassVar& iVar,
        const TClassVar& DerivedVariable,
        const TClassVar& Derived2Variable,
        const double DeltaTime
        )
    {
        // Values
        const double dot2un1 = itNode->FastGetSolutionStepValue(Derived2Variable, 1);
        const double dotun1 = itNode->FastGetSolutionStepValue(DerivedVariable, 1);
        const double un1 = itNode->FastGetSolutionStepValue(iVar, 1);
        const double dot2un0 = itNode->FastGetSolutionStepValue(Derived2Variable);
        double& dotun0 = itNode->FastGetSolutionStepValue(DerivedVariable);
        double& un0 = itNode->FastGetSolutionStepValue(iVar);

        if (itNode->IsFixed(Derived2Variable)) {
            dotun0 = dot2un0;
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                dotun0 -= BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(DerivedVariable, i_order);
            dotun0 /= BDFBaseType::mBDF[0];

            un0 = dotun0;
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                un0 -= BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(iVar, i_order);
            un0 /= BDFBaseType::mBDF[0];
        } else if (itNode->IsFixed(DerivedVariable)) {
            un0 = dotun0;
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                un0 -= BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(iVar, i_order);
            un0 /= BDFBaseType::mBDF[0];
        } else if (!itNode->IsFixed(iVar)) {
            un0 = un1 + DeltaTime * dotun1 + 0.5 * std::pow(DeltaTime, 2) * dot2un1;
        }
    }

    /**
     * @brief This method creates the list of variables
     * @param ThisParameters The configuration parameters
     */
    void CreateVariablesList(Parameters ThisParameters)
    {
        const std::size_t n_variables = ThisParameters["variable"].size();
        const std::size_t n_first_derivative = ThisParameters["first_derivative"].size();
        const std::size_t n_second_derivative = ThisParameters["second_derivative"].size();

        // Size check
        KRATOS_ERROR_IF(n_variables != n_first_derivative) << "Your list of variables is not the same size as the list of first derivatives variables" << std::endl;
        KRATOS_ERROR_IF(n_variables != n_second_derivative) << "Your list of variables is not the same size as the list of second derivatives variables" << std::endl;

        // The current dimension
        const std::size_t domain_size = ThisParameters["domain_size"].GetInt();

        for (std::size_t i_var = 0; i_var < n_variables; ++i_var){
            const std::string& variable_name = ThisParameters["variable"].GetArrayItem(i_var).GetString();
            const std::string& first_derivative_name = ThisParameters["first_derivative"].GetArrayItem(i_var).GetString();
            const std::string& second_derivative_name = ThisParameters["second_derivative"].GetArrayItem(i_var).GetString();

            if(KratosComponents<Variable<double>>::Has(variable_name)){
                Variable<double> variable = KratosComponents< Variable<double> >::Get(variable_name);
                Variable<double> first_derivative = KratosComponents< Variable<double> >::Get(first_derivative_name);
                Variable<double> second_derivative = KratosComponents< Variable<double> >::Get(second_derivative_name);
                mDoubleVariable.push_back(variable);
                mFirtsDoubleDerivatives.push_back(first_derivative);
                mSecondDoubleDerivatives.push_back(second_derivative);
            } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(variable_name)) {
                Variable<array_1d< double, 3>> variable = KratosComponents< Variable<array_1d< double, 3>> >::Get(variable_name);
                Variable<array_1d< double, 3>> first_derivative = KratosComponents< Variable<array_1d< double, 3>> >::Get(first_derivative_name);
                Variable<array_1d< double, 3>> second_derivative = KratosComponents< Variable<array_1d< double, 3>> >::Get(second_derivative_name);

                // Components
                mArrayVariable.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(variable_name+"_X"));
                mArrayVariable.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(variable_name+"_Y"));
                if (domain_size == 3)
                    mArrayVariable.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(variable_name+"_Z"));

                mFirtsArrayDerivatives.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(first_derivative_name+"_X"));
                mFirtsArrayDerivatives.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(first_derivative_name+"_Y"));
                if (domain_size == 3)
                    mFirtsArrayDerivatives.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(first_derivative_name+"_Z"));

                mSecondArrayDerivatives.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(second_derivative_name+"_X"));
                mSecondArrayDerivatives.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(second_derivative_name+"_Y"));
                if (domain_size == 3)
                    mSecondArrayDerivatives.push_back(KratosComponents< VariableComponent<ComponentType>>::Get(second_derivative_name+"_Z"));
            } else {
                KRATOS_ERROR << "Only double and vector variables are allowed in the variables list." ;
            }
        }
    }

    /**
     * @brief This method returns the defaulr parameters in order to avoid code duplication
     * @return Returns the default parameters
     */
    Parameters GetDefaultParameters()
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                  : "ResidualBasedBDFCustomScheme",
            "domain_size"           : 3,
            "integration_order"     : 2,
            "variable"              : ["DISPLACEMENT"],
            "first_derivative"      : ["VELOCITY"],
            "second_derivative"     : ["ACCELERATION"]
        })" );

        return default_parameters;
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
