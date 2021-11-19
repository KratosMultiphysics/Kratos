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


#if !defined(KRATOS_RESIDUAL_BASED_BDF_DISPLACEMENT_SCHEME )
#define  KRATOS_RESIDUAL_BASED_BDF_DISPLACEMENT_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/schemes/residual_based_bdf_scheme.h"
#include "includes/variables.h"
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
 * @class ResidualBasedBDFDisplacementScheme
 * @ingroup KratosCore
 * @brief BDF integration scheme (displacement based)
 * @details The \f$ n \f$ order Backward Differentiation Formula (BDF) method is a two step \f$ n \f$ order accurate method.
 * Look at the base class for more details
 * @see ResidualBasedBDFScheme
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,  class TDenseSpace>
class ResidualBasedBDFDisplacementScheme
    : public ResidualBasedBDFScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualBasedBDFDisplacementScheme
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedBDFDisplacementScheme );

    /// Base class definition
    typedef Scheme<TSparseSpace,TDenseSpace>                                  BaseType;
    typedef ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace> ImplicitBaseType;
    typedef ResidualBasedBDFScheme<TSparseSpace,TDenseSpace>               BDFBaseType;
    typedef ResidualBasedBDFDisplacementScheme<TSparseSpace, TDenseSpace>    ClassType;

    /// Data type definition
    typedef typename BDFBaseType::TDataType                                  TDataType;
    /// Matrix type definition
    typedef typename BDFBaseType::TSystemMatrixType                  TSystemMatrixType;
    /// Vector type definition
    typedef typename BDFBaseType::TSystemVectorType                  TSystemVectorType;
    /// Local system matrix type definition
    typedef typename BDFBaseType::LocalSystemVectorType          LocalSystemVectorType;
    /// Local system vector type definition
    typedef typename BDFBaseType::LocalSystemMatrixType          LocalSystemMatrixType;

    /// DoF array type definition
    typedef typename BDFBaseType::DofsArrayType                          DofsArrayType;
    /// DoF vector type definition
    typedef typename Element::DofsVectorType                            DofsVectorType;

    /// Nodes containers definition
    typedef ModelPart::NodesContainerType                               NodesArrayType;
    /// Elements containers definition
    typedef ModelPart::ElementsContainerType                         ElementsArrayType;
    /// Conditions containers definition
    typedef ModelPart::ConditionsContainerType                     ConditionsArrayType;

    typedef double              ComponentType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. The BDF method (parameters)
     * @param ThisParameters Parameters with the integration order
     */
    explicit ResidualBasedBDFDisplacementScheme(Parameters ThisParameters)
        : BDFBaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Constructor. The BDF method
     * @param Order The integration order
     * @todo The ideal would be to use directly the dof or the variable itself to identify the type of variable and is derivatives
     */
    explicit ResidualBasedBDFDisplacementScheme(const std::size_t Order = 2)
        :BDFBaseType(Order)
    {
    }

    /** Copy Constructor.
     */
    explicit ResidualBasedBDFDisplacementScheme(ResidualBasedBDFDisplacementScheme& rOther)
        :BDFBaseType(rOther)
    {
    }

    /**
     * Clone
     */
    typename BaseType::Pointer Clone() override
    {
        return Kratos::make_shared<ResidualBasedBDFDisplacementScheme>(*this);
    }

    /** Destructor.
     */
    ~ResidualBasedBDFDisplacementScheme
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

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        BDFBaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Updating time derivatives (nodally for efficiency)
        const auto it_node_begin = rModelPart.Nodes().begin();

        // Getting dimension
        KRATOS_WARNING_IF("ResidualBasedBDFDisplacementScheme", !r_current_process_info.Has(DOMAIN_SIZE)) << "DOMAIN_SIZE not defined. Please define DOMAIN_SIZE. 3D case will be assumed" << std::endl;
        const std::size_t dimension = r_current_process_info.Has(DOMAIN_SIZE) ? r_current_process_info.GetValue(DOMAIN_SIZE) : 3;

        // Getting position
        const int velpos = it_node_begin->HasDofFor(VELOCITY_X) ? static_cast<int>(it_node_begin->GetDofPosition(VELOCITY_X)) : -1;
        const int accelpos = it_node_begin->HasDofFor(ACCELERATION_X) ? static_cast<int>(it_node_begin->GetDofPosition(ACCELERATION_X)) : -1;

        std::array<bool, 3> fixed = {false, false, false};
        const std::array<const Variable<ComponentType>*, 3> disp_components = {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};
        const std::array<const Variable<ComponentType>*, 3> vel_components = {&VELOCITY_X, &VELOCITY_Y, &VELOCITY_Z};
        const std::array<const Variable<ComponentType>*, 3> accel_components = {&ACCELERATION_X, &ACCELERATION_Y, &ACCELERATION_Z};

        block_for_each(rModelPart.Nodes(), fixed, [&](Node<3>& rNode, std::array<bool,3>& rFixedTLS){
            for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                rFixedTLS[i_dim] = false;
            }

            if (accelpos > -1) {
                for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                    if (rNode.GetDof(*accel_components[i_dim], accelpos + i_dim).IsFixed()) {
                        rNode.Fix(*disp_components[i_dim]);
                        rFixedTLS[i_dim] = true;
                    }
                }
            }
            if (velpos > -1) {
                for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                    if (rNode.GetDof(*vel_components[i_dim], velpos + i_dim).IsFixed() && !rFixedTLS[i_dim]) {
                        rNode.Fix(*disp_components[i_dim]);
                    }
                }
            }
        });

        KRATOS_CATCH("ResidualBasedBDFDisplacementScheme.InitializeSolutionStep");
    }

    /**
     * @brief Performing the prediction of the solution
     * @details It predicts the solution for the current step x = xold + vold * Dt
     * @param rModelPart The model of the problem to solve
     * @param rDofSet Set of all primary variables
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

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        const auto it_node_begin = rModelPart.Nodes().begin();

        // Getting position
        KRATOS_ERROR_IF_NOT(it_node_begin->HasDofFor(DISPLACEMENT_X)) << "ResidualBasedBDFDisplacementScheme:: DISPLACEMENT is not added" << std::endl;
        const int disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
        const int velpos = it_node_begin->HasDofFor(VELOCITY_X) ? static_cast<int>(it_node_begin->GetDofPosition(VELOCITY_X)) : -1;
        const int accelpos = it_node_begin->HasDofFor(ACCELERATION_X) ? static_cast<int>(it_node_begin->GetDofPosition(ACCELERATION_X)) : -1;

        // Getting dimension
        KRATOS_WARNING_IF("ResidualBasedBDFDisplacementScheme", !r_current_process_info.Has(DOMAIN_SIZE)) << "DOMAIN_SIZE not defined. Please define DOMAIN_SIZE. 3D case will be assumed" << std::endl;
        const std::size_t dimension = r_current_process_info.Has(DOMAIN_SIZE) ? r_current_process_info.GetValue(DOMAIN_SIZE) : 3;

        // Auxiliar variables
        std::array<bool, 3> predicted = {false, false, false};
        const std::array<const Variable<ComponentType>*, 3> disp_components = {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};
        const std::array<const Variable<ComponentType>*, 3> vel_components = {&VELOCITY_X, &VELOCITY_Y, &VELOCITY_Z};
        const std::array<const Variable<ComponentType>*, 3> accel_components = {&ACCELERATION_X, &ACCELERATION_Y, &ACCELERATION_Z};

        IndexPartition<int>(num_nodes).for_each(predicted, [&](int Index, std::array<bool, 3>& rPredictedTLS){
            auto it_node = it_node_begin + Index;

            for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                rPredictedTLS[i_dim] = false;
            }

            const array_1d<double, 3>& dot2un1 = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double, 3>& dotun1 = it_node->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3>& un1 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double, 3>& dot2un0 = it_node->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3>& dotun0 = it_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& un0 = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            if (accelpos > -1) {
                for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                    if (it_node->GetDof(*accel_components[i_dim], accelpos + i_dim).IsFixed()) {
                        dotun0[i_dim] = dot2un0[i_dim];
                        for (std::size_t i_order = 1; i_order < this->mOrder + 1; ++i_order)
                            dotun0[i_dim] -= this->mBDF[i_order] * it_node->FastGetSolutionStepValue(*vel_components[i_dim], i_order);
                        dotun0[i_dim] /= this->mBDF[i_dim];

                        un0[i_dim] = dotun0[i_dim];
                        for (std::size_t i_order = 1; i_order < this->mOrder + 1; ++i_order)
                            un0[i_dim] -= this->mBDF[i_order] * it_node->FastGetSolutionStepValue(*disp_components[i_dim], i_order);
                        un0[i_dim] /= this->mBDF[i_dim];
                        rPredictedTLS[i_dim] = true;
                    }
                }
            }
            if (velpos > -1) {
                for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                    if (it_node->GetDof(*vel_components[i_dim], velpos + i_dim).IsFixed() && !rPredictedTLS[i_dim]) {
                        un0[i_dim] = dotun0[i_dim];
                        for (std::size_t i_order = 1; i_order < this->mOrder + 1; ++i_order)
                            un0[i_dim] -= this->mBDF[i_order] * it_node->FastGetSolutionStepValue(*disp_components[i_dim], i_order);
                        un0[i_dim] /= this->mBDF[i_dim];
                        rPredictedTLS[i_dim] = true;
                    }
                }
            }
            for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
                if (!it_node->GetDof(*disp_components[i_dim], disppos + i_dim).IsFixed() && !rPredictedTLS[i_dim]) {
                    un0[i_dim] = un1[i_dim] + delta_time * dotun1[i_dim] + 0.5 * std::pow(delta_time, 2) * dot2un1[i_dim];
                }
            }

            // Updating time derivatives
            UpdateFirstDerivative(it_node);
            UpdateSecondDerivative(it_node);
        });

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

        // Check that variables are correctly allocated
        for(auto& rnode : rModelPart.Nodes()) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
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
            "name"               : "bdf_displacement_scheme",
            "integration_order"  : 2
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
        return "bdf_displacement_scheme";
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
        return "ResidualBasedBDFDisplacementScheme";
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
        array_1d<double, 3>& dotun0 = itNode->FastGetSolutionStepValue(VELOCITY);
        noalias(dotun0) = this->mBDF[0] * itNode->FastGetSolutionStepValue(DISPLACEMENT);
        for (std::size_t i_order = 1; i_order < this->mOrder + 1; ++i_order)
            noalias(dotun0) += this->mBDF[i_order] * itNode->FastGetSolutionStepValue(DISPLACEMENT, i_order);
    }

    /**
     * @brief Updating second time derivative (acceleration)
     * @param itNode the node interator
     */
    inline void UpdateSecondDerivative(NodesArrayType::iterator itNode) override
    {
        array_1d<double, 3>& dot2un0 = itNode->FastGetSolutionStepValue(ACCELERATION);
        noalias(dot2un0) = this->mBDF[0] * itNode->FastGetSolutionStepValue(VELOCITY);
        for (std::size_t i_order = 1; i_order < this->mOrder + 1; ++i_order)
            noalias(dot2un0) += this->mBDF[i_order] * itNode->FastGetSolutionStepValue(VELOCITY, i_order);
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
}; /* Class ResidualBasedBDFDisplacementScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BDF_DISPLACEMENT_SCHEME defined */
