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
#include "unordered_map"

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
     * Constructor.
     * The BDF method
     * @param Order The integration order
     * @param rParameters The parameters containing the list of variables to consider
     * @todo The ideal would be to use directly the dof or the variable itself to identify the type of variable and is derivatives
     */
    ResidualBasedBDFCustomScheme(
        const std::size_t Order = 2,
        Parameters rParameters =  Parameters(R"({})")
        )
        :BDFBaseType(Order)
    {
        Parameters default_parameters = GetDefaultParameters(); 
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        const std::size_t n_variables = rParameters["variable"].size();
        const std::size_t n_first_derivative = rParameters["first_derivative"].size();
        const std::size_t n_second_derivative = rParameters["second_derivative"].size();
    
        // Size check
        KRATOS_ERROR_IF(n_variables != n_first_derivative) << "Your list of variables is not the same size as the list of first derivatives variables" << std::endl;
        KRATOS_ERROR_IF(n_variables != n_second_derivative) << "Your list of variables is not the same size as the list of second derivatives variables" << std::endl;
        
        for (unsigned int i_var = 0; i_var < n_variables; ++i_var){
            std::string variable_name = rParameters["variable"].GetArrayItem(i_var).GetString();
            std::string first_derivative_name = rParameters["first_derivative"].GetArrayItem(i_var).GetString();
            std::string second_derivative_name = rParameters["second_derivative"].GetArrayItem(i_var).GetString();
        
            if(KratosComponents<Variable<double>>::Has(variable_name)){
                Variable<double> variable = KratosComponents< Variable<double> >::Get(variable_name);
                Variable<double> first_derivative = KratosComponents< Variable<double> >::Get(first_derivative_name);
                Variable<double> second_derivative = KratosComponents< Variable<double> >::Get(second_derivative_name);
                const std::size_t key = variable.Key();
                mDoubleVariable[key] = variable;
                mFirtsDoubleDerivatives[key] = first_derivative;
                mSecondDoubleDerivatives[key] = second_derivative;
            } else if (KratosComponents< Variable< array_1d< double, 3> > >::Has(variable_name)) {
                Variable<array_1d< double, 3>> variable = KratosComponents< Variable<array_1d< double, 3>> >::Get(variable_name);
                Variable<array_1d< double, 3>> first_derivative = KratosComponents< Variable<array_1d< double, 3>> >::Get(first_derivative_name);
                Variable<array_1d< double, 3>> second_derivative = KratosComponents< Variable<array_1d< double, 3>> >::Get(second_derivative_name);
                const std::size_t key = variable.Key();
                mArrayVariable[key] = variable;
                mFirtsArrayDerivatives[key] = first_derivative;
                mSecondArrayDerivatives[key] = second_derivative;
            } else {
                KRATOS_ERROR << "Only double, component and vector variables are allowed in the visualization variables list." ;
            }
            
        }
    }

    /** Copy Constructor.
     */
    ResidualBasedBDFCustomScheme(ResidualBasedBDFCustomScheme& rOther)
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

        ProcessInfo& current_process_info = rModelPart.GetProcessInfo();
        const double delta_time = current_process_info[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );

        #pragma omp parallel for
        for(int i = 0;  i< num_nodes; ++i) {
            auto it_node = rModelPart.Nodes().begin() + i;

            //ATTENTION::: the prediction is performed only on free nodes
            
            for ( auto var_it = mDoubleVariable.begin(); var_it!= mDoubleVariable.end(); ++var_it ) {
                // Derivatives
                Variable<double> dvar = mFirtsDoubleDerivatives[var_it->first];
                Variable<double> d2var = mSecondDoubleDerivatives[var_it->first];
                
                // Values
                const double dot2un1 = it_node->FastGetSolutionStepValue(d2var, 1);
                const double dotun1 = it_node->FastGetSolutionStepValue(dvar, 1);
                const double un1 = it_node->FastGetSolutionStepValue(var_it->second, 1);
                const double dot2un0 = it_node->FastGetSolutionStepValue(d2var);
                double& dotun0 = it_node->FastGetSolutionStepValue(dvar);
                double& un0 = it_node->FastGetSolutionStepValue(var_it->second);
                
                if (it_node->HasDofFor(d2var)) {
                    if (it_node -> IsFixed(d2var)) {
                        dotun0 = (dot2un0 - BDFBaseType::mBDF[1] * dotun1)/BDFBaseType::mBDF[0];
                        un0 = (dotun0 - BDFBaseType::mBDF[1] * un1)/BDFBaseType::mBDF[0];
                } } else if (it_node->HasDofFor(dvar)) {
                    if (it_node -> IsFixed(dvar)) {
                        un0 = (dotun1 - BDFBaseType::mBDF[1] * un1)/BDFBaseType::mBDF[0];
                } } else if (it_node -> IsFixed(var_it->second) == false) {
                    un0 = un1 + delta_time * dotun1 + 0.5 * std::pow(delta_time, 2) * dot2un1;
                }
                
                for (std::size_t i_order = 2; i_order < BDFBaseType::mOrder + 1; ++i_order) {
                    const double dotun = it_node->FastGetSolutionStepValue(dvar, i_order);
                    const double un = it_node->FastGetSolutionStepValue(var_it->second, i_order);
                    
                    if (it_node->HasDofFor(d2var)) {
                        if (it_node -> IsFixed(d2var)) {
                            dotun0 -= (BDFBaseType::mBDF[i_order] * dotun)/BDFBaseType::mBDF[0];
                            un0 -= (BDFBaseType::mBDF[i_order] * un)/BDFBaseType::mBDF[0];
                    } } else if (it_node->HasDofFor(dvar)) {
                        if (it_node -> IsFixed(dvar)) {
                            un0 -= (BDFBaseType::mBDF[i_order] * un)/BDFBaseType::mBDF[0];
                    } }
                }
            }
            
            for ( auto var_it = mArrayVariable.begin(); var_it!= mArrayVariable.end(); ++var_it ) {
                // Values
                const array_1d<double, 3>& dot2un1 = it_node->FastGetSolutionStepValue(mSecondArrayDerivatives[var_it->first], 1);
                const array_1d<double, 3>& dotun1 = it_node->FastGetSolutionStepValue(mFirtsArrayDerivatives[var_it->first], 1);
                const array_1d<double, 3>& un1 = it_node->FastGetSolutionStepValue(var_it->second, 1);
                const array_1d<double, 3>& dot2un0 = it_node->FastGetSolutionStepValue(mSecondArrayDerivatives[var_it->first]);
                array_1d<double, 3>& dotun0 = it_node->FastGetSolutionStepValue(mFirtsArrayDerivatives[var_it->first]);
                array_1d<double, 3>& un0 = it_node->FastGetSolutionStepValue(var_it->second);
                
                // Components
                std::string variable_name = (var_it->second).Name();
                VariableComponent<ComponentType> var_x = KratosComponents< VariableComponent<ComponentType>>::Get(variable_name.append("_X"));
                VariableComponent<ComponentType> var_y = KratosComponents< VariableComponent<ComponentType>>::Get(variable_name.append("_Y"));
                VariableComponent<ComponentType> var_z = KratosComponents< VariableComponent<ComponentType>>::Get(variable_name.append("_Z"));
                std::string dvariable_name = (mFirtsArrayDerivatives[var_it->first]).Name();
                VariableComponent<ComponentType> dvar_x = KratosComponents< VariableComponent<ComponentType>>::Get(dvariable_name.append("_X"));
                VariableComponent<ComponentType> dvar_y = KratosComponents< VariableComponent<ComponentType>>::Get(dvariable_name.append("_Y"));
                VariableComponent<ComponentType> dvar_z = KratosComponents< VariableComponent<ComponentType>>::Get(dvariable_name.append("_Z"));
                std::string d2variable_name = (mSecondArrayDerivatives[var_it->first]).Name();
                VariableComponent<ComponentType> d2var_x = KratosComponents< VariableComponent<ComponentType>>::Get(d2variable_name.append("_X"));
                VariableComponent<ComponentType> d2var_y = KratosComponents< VariableComponent<ComponentType>>::Get(d2variable_name.append("_Y"));
                VariableComponent<ComponentType> d2var_z = KratosComponents< VariableComponent<ComponentType>>::Get(d2variable_name.append("_Z"));
                
                if (it_node->HasDofFor(d2var_x)) {
                    if (it_node -> IsFixed(d2var_x)) {
                        dotun0[0] = (dot2un0[0] - BDFBaseType::mBDF[1] * dotun1[0])/BDFBaseType::mBDF[0];
                        un0[0] = (dotun0[0] - BDFBaseType::mBDF[1] * un1[0])/BDFBaseType::mBDF[0];
                } } else if (it_node->HasDofFor(dvar_x)) {
                    if (it_node -> IsFixed(dvar_x)) {
                        un0[0] = (dotun1[0] - BDFBaseType::mBDF[1] * un1[0])/BDFBaseType::mBDF[0];
                } } else if (it_node -> IsFixed(var_x) == false) {
                    un0[0] = un1[0] + delta_time * dotun1[0] + 0.5 * std::pow(delta_time, 2) * dot2un1[0];
                }

                if (it_node->HasDofFor(d2var_y)) {
                    if (it_node -> IsFixed(d2var_y)) {
                        dotun0[1] = (dot2un0[1] - BDFBaseType::mBDF[1] * dotun1[1])/BDFBaseType::mBDF[0];
                        un0[1] = (dotun0[1] - BDFBaseType::mBDF[1] * un1[1])/BDFBaseType::mBDF[0];
                } } else if (it_node->HasDofFor(dvar_y)) {
                    if (it_node -> IsFixed(dvar_y)) {
                        un0[1] = (dotun1[1] - BDFBaseType::mBDF[1] * un1[1])/BDFBaseType::mBDF[0];
                } } else if (it_node -> IsFixed(var_y) == false) {
                    un0[1] = un1[1] + delta_time * dotun1[1] + 0.5 * std::pow(delta_time, 2) * dot2un1[1];
                }

                // For 3D cases
                if (it_node -> HasDofFor(var_z)) {
                    if (it_node->HasDofFor(d2var_z)) {
                        if (it_node -> IsFixed(d2var_z)) {
                            dotun0[2] = (dot2un0[2] - BDFBaseType::mBDF[1] * dotun1[2])/BDFBaseType::mBDF[0];
                            un0[2] = (dotun0[2] - BDFBaseType::mBDF[1] * un1[2])/BDFBaseType::mBDF[0];
                    } } else if (it_node->HasDofFor(dvar_y)) {
                        if (it_node -> IsFixed(dvar_y)) {
                            un0[2] = (dotun1[2] - BDFBaseType::mBDF[1] * un1[2])/BDFBaseType::mBDF[0];
                    } } else if (it_node -> IsFixed(var_z) == false) {
                        un0[2] = un1[2] + delta_time * dotun1[2] + 0.5 * std::pow(delta_time, 2) * dot2un1[2];
                    }
                }
                
                for (std::size_t i_order = 2; i_order < BDFBaseType::mOrder + 1; ++i_order) {
                    const array_1d<double, 3>& dotun = it_node->FastGetSolutionStepValue(mFirtsArrayDerivatives[var_it->first], i_order);
                    const array_1d<double, 3>& un = it_node->FastGetSolutionStepValue(var_it->second, i_order);
                    
                    if (it_node->HasDofFor(d2var_x)) {
                        if (it_node -> IsFixed(d2var_x)) {
                            dotun0[0] -= (BDFBaseType::mBDF[i_order] * dotun[0])/BDFBaseType::mBDF[0];
                            un0[0] -= (BDFBaseType::mBDF[i_order] * un[0])/BDFBaseType::mBDF[0];
                    } } else if (it_node->HasDofFor(dvar_x)) {
                        if (it_node -> IsFixed(dvar_x)) {
                            un0[0] -= (BDFBaseType::mBDF[i_order] * un[0])/BDFBaseType::mBDF[0];
                    } }

                    if (it_node->HasDofFor(d2var_y)) {
                        if (it_node -> IsFixed(d2var_y)) {
                            dotun0[1] -= (BDFBaseType::mBDF[i_order] * dotun[1])/BDFBaseType::mBDF[0];
                            un0[1] -= (BDFBaseType::mBDF[i_order] * un[1])/BDFBaseType::mBDF[0];
                    } } else if (it_node->HasDofFor(dvar_y)) {
                        if (it_node -> IsFixed(dvar_x)) {
                            un0[1] -= (BDFBaseType::mBDF[i_order] * un[1])/BDFBaseType::mBDF[0];
                    } }

                    // For 3D cases
                    if (it_node -> HasDofFor(var_z)) {
                        if (it_node->HasDofFor(d2var_z)) {
                            if (it_node -> IsFixed(d2var_z)) {
                                dotun0[1] -= (BDFBaseType::mBDF[i_order] * dotun[2])/BDFBaseType::mBDF[0];
                                un0[1] -= (BDFBaseType::mBDF[i_order] * un[2])/BDFBaseType::mBDF[0];
                        } } else if (it_node->HasDofFor(dvar_y)) {
                            if (it_node -> IsFixed(dvar_x)) {
                                un0[1] -= (BDFBaseType::mBDF[i_order] * un[2])/BDFBaseType::mBDF[0];
                        } }
                    }
                }
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
        for ( auto var_it = mDoubleVariable.begin(); var_it!= mDoubleVariable.end(); ++var_it )
            KRATOS_CHECK_VARIABLE_KEY(var_it->second) 
        for ( auto var_it = mFirtsDoubleDerivatives.begin(); var_it!= mFirtsDoubleDerivatives.end(); ++var_it )
            KRATOS_CHECK_VARIABLE_KEY(var_it->second) 
        for ( auto var_it = mSecondDoubleDerivatives.begin(); var_it!= mSecondDoubleDerivatives.end(); ++var_it )
            KRATOS_CHECK_VARIABLE_KEY(var_it->second) 
        for ( auto var_it = mArrayVariable.begin(); var_it!= mArrayVariable.end(); ++var_it )
            KRATOS_CHECK_VARIABLE_KEY(var_it->second) 
        for ( auto var_it = mFirtsArrayDerivatives.begin(); var_it!= mFirtsArrayDerivatives.end(); ++var_it )
            KRATOS_CHECK_VARIABLE_KEY(var_it->second) 
        for ( auto var_it = mSecondArrayDerivatives.begin(); var_it!= mSecondArrayDerivatives.end(); ++var_it )
            KRATOS_CHECK_VARIABLE_KEY(var_it->second) 

        // Check that variables are correctly allocated
        for(auto& rnode : rModelPart.Nodes()) {
            for ( auto var_it = mDoubleVariable.begin(); var_it!= mDoubleVariable.end(); ++var_it )
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(var_it->second, rnode) 
            for ( auto var_it = mFirtsDoubleDerivatives.begin(); var_it!= mFirtsDoubleDerivatives.end(); ++var_it )
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(var_it->second, rnode)
            for ( auto var_it = mSecondDoubleDerivatives.begin(); var_it!= mSecondDoubleDerivatives.end(); ++var_it )
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(var_it->second, rnode) 
            for ( auto var_it = mArrayVariable.begin(); var_it!= mArrayVariable.end(); ++var_it )
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(var_it->second, rnode) 
            for ( auto var_it = mFirtsArrayDerivatives.begin(); var_it!= mFirtsArrayDerivatives.end(); ++var_it )
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(var_it->second, rnode) 
            for ( auto var_it = mSecondArrayDerivatives.begin(); var_it!= mSecondArrayDerivatives.end(); ++var_it )
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(var_it->second, rnode) 
    
            for ( auto var_it = mDoubleVariable.begin(); var_it!= mDoubleVariable.end(); ++var_it )
                KRATOS_CHECK_DOF_IN_NODE(var_it->second, rnode)
                
            for ( auto var_it = mArrayVariable.begin(); var_it!= mArrayVariable.end(); ++var_it ) {
                std::string variable_name = (var_it->second).Name();
                VariableComponent<ComponentType> var_x = KratosComponents< VariableComponent<ComponentType>>::Get(variable_name.append("_X"));
                VariableComponent<ComponentType> var_y = KratosComponents< VariableComponent<ComponentType>>::Get(variable_name.append("_Y"));
                VariableComponent<ComponentType> var_z = KratosComponents< VariableComponent<ComponentType>>::Get(variable_name.append("_Z"));
                KRATOS_CHECK_DOF_IN_NODE(var_x, rnode)
                KRATOS_CHECK_DOF_IN_NODE(var_y, rnode)
                KRATOS_CHECK_DOF_IN_NODE(var_z, rnode)
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

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    
    std::unordered_map<std::size_t, Variable< double>> mDoubleVariable;                     /// The double variables  
    std::unordered_map<std::size_t, Variable< double>> mFirtsDoubleDerivatives;             /// The first derivative double variable to compute 
    std::unordered_map<std::size_t, Variable< double>> mSecondDoubleDerivatives;            /// The second derivative double variable to compute
    std::unordered_map<std::size_t, Variable<array_1d<double, 3>>> mArrayVariable;          /// The array variables to compute 
    std::unordered_map<std::size_t, Variable<array_1d<double, 3>>> mFirtsArrayDerivatives;  /// The first derivative array variable to compute 
    std::unordered_map<std::size_t, Variable<array_1d<double, 3>>> mSecondArrayDerivatives; /// The second derivative array variable to compute

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
        for ( auto var_it = mDoubleVariable.begin(); var_it!= mDoubleVariable.end(); ++var_it ) {
            double& dotun0 = itNode->FastGetSolutionStepValue(mFirtsDoubleDerivatives[var_it->first]);
            dotun0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(var_it->second);
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                dotun0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(var_it->second, i_order);
        }
        
        // ARRAYS
        for ( auto var_it = mArrayVariable.begin(); var_it!= mArrayVariable.end(); ++var_it ) {
            array_1d<double, 3>& dotun0 = itNode->FastGetSolutionStepValue(mFirtsArrayDerivatives[var_it->first]);
            noalias(dotun0) = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(var_it->second);
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                noalias(dotun0) += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(var_it->second, i_order);
        }
    }

    /**
     * @brief Updating second time derivative (acceleration)
     * @param itNode the node interator
     */
    
    inline void UpdateSecondDerivative(NodesArrayType::iterator itNode) override
    {
        // DOUBLES
        for ( auto var_it = mFirtsDoubleDerivatives.begin(); var_it!= mFirtsDoubleDerivatives.end(); ++var_it ) {
            double& dot2un0 = itNode->FastGetSolutionStepValue(mSecondDoubleDerivatives[var_it->first]);
            dot2un0 = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(var_it->second);
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                dot2un0 += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(var_it->second, i_order);
        }
        
        // ARRAYS
        for ( auto var_it = mFirtsArrayDerivatives.begin(); var_it!= mFirtsArrayDerivatives.end(); ++var_it ) {
            array_1d<double, 3>& dot2un0 = itNode->FastGetSolutionStepValue(mSecondArrayDerivatives[var_it->first]);
            noalias(dot2un0) = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(var_it->second);
            for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
                noalias(dot2un0) += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(var_it->second, i_order);
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
     * @brief This method returns the defaulr parameters in order to avoid code duplication
     * @return Returns the default parameters
     */
    
    Parameters GetDefaultParameters()
    {
        Parameters default_parameters = Parameters(R"(
        {
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
