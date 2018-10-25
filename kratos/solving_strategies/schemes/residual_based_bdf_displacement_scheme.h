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

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor. The BDF method (parameters)
     * @param ThisParameters Parameters with the integration order
     */
    explicit ResidualBasedBDFDisplacementScheme(Parameters ThisParameters)
        : ResidualBasedBDFDisplacementScheme(ThisParameters.Has("integration_order") ? static_cast<std::size_t>(ThisParameters["integration_order"].GetInt()) : 2)
    {
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

        ProcessInfo& current_process_info = rModelPart.GetProcessInfo();
        const double delta_time = current_process_info[DELTA_TIME];

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );

        #pragma omp parallel for
        for(int i = 0;  i< num_nodes; ++i) {
            auto it_node = rModelPart.Nodes().begin() + i;

            //ATTENTION::: the prediction is performed only on free nodes

            const array_1d<double, 3>& dot2un1 = it_node->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double, 3>& dotun1 = it_node->FastGetSolutionStepValue(VELOCITY,     1);
            const array_1d<double, 3>& un1 = it_node->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double, 3>& dot2un0 = it_node->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3>& dotun0 = it_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& un0 = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            if (it_node->HasDofFor(ACCELERATION_X)) {
                if (it_node -> IsFixed(ACCELERATION_X)) {
                    dotun0[0] = (dot2un0[0] - BDFBaseType::mBDF[1] * dotun1[0])/BDFBaseType::mBDF[0];
                    un0[0] = (dotun0[0] - BDFBaseType::mBDF[1] * un1[0])/BDFBaseType::mBDF[0];
            } } else if (it_node->HasDofFor(VELOCITY_X)) {
                if (it_node -> IsFixed(VELOCITY_X)) {
                    un0[0] = (dotun1[0] - BDFBaseType::mBDF[1] * un1[0])/BDFBaseType::mBDF[0];
            } } else if (it_node -> IsFixed(DISPLACEMENT_X) == false) {
                un0[0] = un1[0] + delta_time * dotun1[0] + 0.5 * std::pow(delta_time, 2) * dot2un1[0];
            }

            if (it_node->HasDofFor(ACCELERATION_Y)) {
                if (it_node -> IsFixed(ACCELERATION_Y)) {
                    dotun0[1] = (dot2un0[1] - BDFBaseType::mBDF[1] * dotun1[1])/BDFBaseType::mBDF[0];
                    un0[1] = (dotun0[1] - BDFBaseType::mBDF[1] * un1[1])/BDFBaseType::mBDF[0];
            } } else if (it_node->HasDofFor(VELOCITY_Y)) {
                if (it_node -> IsFixed(VELOCITY_Y)) {
                    un0[1] = (dotun1[1] - BDFBaseType::mBDF[1] * un1[1])/BDFBaseType::mBDF[0];
            } } else if (it_node -> IsFixed(DISPLACEMENT_Y) == false) {
                un0[1] = un1[1] + delta_time * dotun1[1] + 0.5 * std::pow(delta_time, 2) * dot2un1[1];
            }

            // For 3D cases
            if (it_node -> HasDofFor(DISPLACEMENT_Z)) {
                if (it_node->HasDofFor(ACCELERATION_Z)) {
                    if (it_node -> IsFixed(ACCELERATION_Z)) {
                        dotun0[2] = (dot2un0[2] - BDFBaseType::mBDF[1] * dotun1[2])/BDFBaseType::mBDF[0];
                        un0[2] = (dotun0[2] - BDFBaseType::mBDF[1] * un1[2])/BDFBaseType::mBDF[0];
                } } else if (it_node->HasDofFor(VELOCITY_Y)) {
                    if (it_node -> IsFixed(VELOCITY_Y)) {
                        un0[2] = (dotun1[2] - BDFBaseType::mBDF[1] * un1[2])/BDFBaseType::mBDF[0];
                } } else if (it_node -> IsFixed(DISPLACEMENT_Z) == false) {
                    un0[2] = un1[2] + delta_time * dotun1[2] + 0.5 * std::pow(delta_time, 2) * dot2un1[2];
                }
            }

            for (std::size_t i_order = 2; i_order < BDFBaseType::mOrder + 1; ++i_order) {
                const array_1d<double, 3>& dotun = it_node->FastGetSolutionStepValue(VELOCITY,     i_order);
                const array_1d<double, 3>& un = it_node->FastGetSolutionStepValue(DISPLACEMENT, i_order);

                if (it_node->HasDofFor(ACCELERATION_X)) {
                    if (it_node -> IsFixed(ACCELERATION_X)) {
                        dotun0[0] -= (BDFBaseType::mBDF[i_order] * dotun[0])/BDFBaseType::mBDF[0];
                        un0[0] -= (BDFBaseType::mBDF[i_order] * un[0])/BDFBaseType::mBDF[0];
                } } else if (it_node->HasDofFor(VELOCITY_X)) {
                    if (it_node -> IsFixed(VELOCITY_X)) {
                        un0[0] -= (BDFBaseType::mBDF[i_order] * un[0])/BDFBaseType::mBDF[0];
                } }

                if (it_node->HasDofFor(ACCELERATION_Y)) {
                    if (it_node -> IsFixed(ACCELERATION_Y)) {
                        dotun0[1] -= (BDFBaseType::mBDF[i_order] * dotun[1])/BDFBaseType::mBDF[0];
                        un0[1] -= (BDFBaseType::mBDF[i_order] * un[1])/BDFBaseType::mBDF[0];
                } } else if (it_node->HasDofFor(VELOCITY_Y)) {
                    if (it_node -> IsFixed(VELOCITY_X)) {
                        un0[1] -= (BDFBaseType::mBDF[i_order] * un[1])/BDFBaseType::mBDF[0];
                } }

                // For 3D cases
                if (it_node -> HasDofFor(DISPLACEMENT_Z)) {
                    if (it_node->HasDofFor(ACCELERATION_Z)) {
                        if (it_node -> IsFixed(ACCELERATION_Z)) {
                            dotun0[1] -= (BDFBaseType::mBDF[i_order] * dotun[2])/BDFBaseType::mBDF[0];
                            un0[1] -= (BDFBaseType::mBDF[i_order] * un[2])/BDFBaseType::mBDF[0];
                    } } else if (it_node->HasDofFor(VELOCITY_Y)) {
                        if (it_node -> IsFixed(VELOCITY_X)) {
                            un0[1] -= (BDFBaseType::mBDF[i_order] * un[2])/BDFBaseType::mBDF[0];
                    } }
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
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)

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
        noalias(dotun0) = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(DISPLACEMENT);
        for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
            noalias(dotun0) += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(DISPLACEMENT, i_order);
    }

    /**
     * @brief Updating second time derivative (acceleration)
     * @param itNode the node interator
     */
    inline void UpdateSecondDerivative(NodesArrayType::iterator itNode) override
    {
        array_1d<double, 3>& dot2un0 = itNode->FastGetSolutionStepValue(ACCELERATION);
        noalias(dot2un0) = BDFBaseType::mBDF[0] * itNode->FastGetSolutionStepValue(VELOCITY);
        for (std::size_t i_order = 1; i_order < BDFBaseType::mOrder + 1; ++i_order)
            noalias(dot2un0) += BDFBaseType::mBDF[i_order] * itNode->FastGetSolutionStepValue(VELOCITY, i_order);
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
