// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PENALTY_FRICTIONAL_MORTAR_CRITERIA_H)
#define  KRATOS_PENALTY_FRICTIONAL_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/table_stream_utility.h"
#include "custom_strategies/custom_convergencecriterias/base_mortar_criteria.h"
#include "utilities/color_utilities.h"
#include "custom_utilities/active_set_utilities.h"

namespace Kratos
{
///@addtogroup ContactStructuralMechanicsApplication
///@{

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
 * @class PenaltyFrictionalMortarConvergenceCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom convergence criteria for the mortar condition for frictional case
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class PenaltyFrictionalMortarConvergenceCriteria
    : public  BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PenaltyFrictionalMortarConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION( PenaltyFrictionalMortarConvergenceCriteria );

    /// The base convergence criteria class definition
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > ConvergenceCriteriaBaseType;

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( PRINTING_OUTPUT );
    KRATOS_DEFINE_LOCAL_FLAG( TABLE_IS_INITIALIZED );

    /// The base class definition (and it subclasses)
    typedef BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >          BaseType;
    typedef typename BaseType::TDataType                                       TDataType;
    typedef typename BaseType::DofsArrayType                               DofsArrayType;
    typedef typename BaseType::TSystemMatrixType                       TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType                       TSystemVectorType;

    /// The sparse space used
    typedef TSparseSpace                                                 SparseSpaceType;

    /// The components containers
    typedef ModelPart::NodesContainerType                                 NodesArrayType;
    typedef ModelPart::ConditionsContainerType                       ConditionsArrayType;

    /// The r_table stream definition TODO: Replace by logger
    typedef TableStreamUtility::Pointer                          TablePrinterPointerType;

    /// The index type definition
    typedef std::size_t                                                        IndexType;

    /// The epsilon tolerance definition
    static constexpr double Tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    explicit PenaltyFrictionalMortarConvergenceCriteria(
        const bool PureSlip = false,
        const bool PrintingOutput = false,
        const bool ComputeDynamicFactor = true,
        const bool IODebug = false
        ) : BaseType(ComputeDynamicFactor, IODebug, PureSlip)
    {
        // Set local flags
        BaseType::mOptions.Set(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT, PrintingOutput);
        BaseType::mOptions.Set(PenaltyFrictionalMortarConvergenceCriteria::TABLE_IS_INITIALIZED, false);
    }

    ///Copy constructor
    PenaltyFrictionalMortarConvergenceCriteria( PenaltyFrictionalMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /// Destructor
    ~PenaltyFrictionalMortarConvergenceCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Criterias that need to be called before getting the solution
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        BaseType::PreCriteria(rModelPart, rDofSet, rA, rDx, rb);

        return true;
    }

    /**
     * @brief Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // We call the base class
        BaseType::PostCriteria(rModelPart, rDofSet, rA, rDx, rb);

        // Compute the active set
        const array_1d<std::size_t, 2> is_converged = ActiveSetUtilities::ComputePenaltyFrictionalActiveSet(rModelPart, BaseType::mOptions.Is(BaseType::PURE_SLIP), this->GetEchoLevel());

        // We save to the process info if the active set has converged
        const bool active_set_converged = (is_converged[0] + is_converged[1]) == 0 ? true : false;

        // We get the process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        r_process_info[ACTIVE_SET_CONVERGED] = active_set_converged;

        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
            if (r_process_info.Has(TABLE_UTILITY)) {
                TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                auto& r_table = p_table->GetTable();
                if (is_converged[0] == 0) {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        r_table << BOLDFONT(FGRN("       Achieved"));
                    else
                        r_table << "Achieved";
                } else {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        r_table << BOLDFONT(FRED("   Not achieved"));
                    else
                        r_table << "Not achieved";
                }
                if (is_converged[1] == 0) {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        r_table << BOLDFONT(FGRN("       Achieved"));
                    else
                        r_table << "Achieved";
                } else {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        r_table << BOLDFONT(FRED("   Not achieved"));
                    else
                        r_table << "Not achieved";
                }
            } else {
                if (is_converged[0] == 0) {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        KRATOS_INFO("PenaltyFrictionalMortarConvergenceCriteria") << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    else
                        KRATOS_INFO("PenaltyFrictionalMortarConvergenceCriteria") << "\tActive set convergence is achieved" << std::endl;
                } else {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        KRATOS_INFO("PenaltyFrictionalMortarConvergenceCriteria") << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    else
                        KRATOS_INFO("PenaltyFrictionalMortarConvergenceCriteria") << "\tActive set convergence is not achieved" << std::endl;
                }

                if (is_converged[1] == 0) {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        KRATOS_INFO("PenaltyFrictionalMortarConvergenceCriteria") << BOLDFONT("\tSlip/stick set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    else
                        KRATOS_INFO("PenaltyFrictionalMortarConvergenceCriteria") << "\tSlip/stick set convergence is achieved" << std::endl;
                } else {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        KRATOS_INFO("PenaltyFrictionalMortarConvergenceCriteria") << BOLDFONT("\tSlip/stick set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    else
                        KRATOS_INFO("PenaltyFrictionalMortarConvergenceCriteria") << "\tSlip/stick set  convergence is not achieved" << std::endl;
                }
            }
        }

        return active_set_converged;
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        ConvergenceCriteriaBaseType::mConvergenceCriteriaIsInitialized = true;

        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info.Has(TABLE_UTILITY) && BaseType::mOptions.IsNot(PenaltyFrictionalMortarConvergenceCriteria::TABLE_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& r_table = p_table->GetTable();
            r_table.AddColumn("ACTIVE SET CONV", 15);
            r_table.AddColumn("SLIP/STICK CONV", 15);
            BaseType::mOptions.Set(PenaltyFrictionalMortarConvergenceCriteria::TABLE_IS_INITIALIZED, true);
        }
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
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
     * @brief This method resets the weighted gap in the nodes of the problem
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     */
    void ResetWeightedGap(ModelPart& rModelPart) override
    {
        // Auxiliar zero array
        const array_1d<double, 3> zero_array = ZeroVector(3);

        // We reset the weighted values
        NodesArrayType& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        VariableUtils().SetVariable(WEIGHTED_GAP, 0.0, r_nodes_array);
        VariableUtils().SetVariable(WEIGHTED_SLIP, zero_array, r_nodes_array);
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
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}

}; // Class PenaltyFrictionalMortarConvergenceCriteria

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags PenaltyFrictionalMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::PRINTING_OUTPUT(Kratos::Flags::Create(3));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags PenaltyFrictionalMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::NOT_PRINTING_OUTPUT(Kratos::Flags::Create(3, false));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags PenaltyFrictionalMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(4));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags PenaltyFrictionalMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::NOT_TABLE_IS_INITIALIZED(Kratos::Flags::Create(4, false));

}  // namespace Kratos

#endif /* KRATOS_PENALTY_FRICTIONAL_MORTAR_CRITERIA_H  defined */

