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

#if !defined(KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H)
#define  KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H

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
 * @class ALMFrictionalMortarConvergenceCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom convergence criteria for the mortar condition for frictional case
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class ALMFrictionalMortarConvergenceCriteria
    : public BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ALMFrictionalMortarConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION( ALMFrictionalMortarConvergenceCriteria );

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( PRINTING_OUTPUT );
    KRATOS_DEFINE_LOCAL_FLAG( TABLE_IS_INITIALIZED );

    /// The base convergence criteria class definition
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > ConvergenceCriteriaBaseType;

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
    explicit ALMFrictionalMortarConvergenceCriteria(
        const bool PureSlip = false,
        const bool PrintingOutput = false,
        const bool ComputeDynamicFactor = false,
        const bool IODebug = false
        ) : BaseType(ComputeDynamicFactor, IODebug, PureSlip)
    {
        // Set local flags
        BaseType::mOptions.Set(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT, PrintingOutput);
        BaseType::mOptions.Set(ALMFrictionalMortarConvergenceCriteria::TABLE_IS_INITIALIZED, false);
    }

    ///Copy constructor
    ALMFrictionalMortarConvergenceCriteria( ALMFrictionalMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /// Destructor
    ~ALMFrictionalMortarConvergenceCriteria() override = default;

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

        // Getting process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // Compute the active set
        if (!r_process_info[ACTIVE_SET_COMPUTED]) {
            const array_1d<std::size_t, 2> is_converged = ActiveSetUtilities::ComputeALMFrictionalActiveSet(rModelPart, BaseType::mOptions.Is(BaseType::PURE_SLIP), this->GetEchoLevel());

            // We save to the process info if the active set has converged
            r_process_info[ACTIVE_SET_CONVERGED] = is_converged[0] == 0 ? true : false;
            r_process_info[SLIP_SET_CONVERGED] = is_converged[1] == 0 ? true : false;
            r_process_info[ACTIVE_SET_COMPUTED] = true;
        }

        // Getting converged bools
        const bool active_set_converged = r_process_info[ACTIVE_SET_CONVERGED];
        const bool slip_set_converged = r_process_info[SLIP_SET_CONVERGED];

        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
            if (r_process_info.Has(TABLE_UTILITY)) {
                TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                auto& r_table = p_table->GetTable();
                if (active_set_converged) {
                    if (BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        r_table << BOLDFONT(FGRN("       Achieved"));
                    else
                        r_table << "Achieved";
                } else {
                    if (BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        r_table << BOLDFONT(FRED("   Not achieved"));
                    else
                        r_table << "Not achieved";
                }
                if (BaseType::mOptions.IsNot(BaseType::PURE_SLIP)) {
                    if (slip_set_converged) {
                        if (BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                            r_table << BOLDFONT(FGRN("       Achieved"));
                        else
                            r_table << "Achieved";
                    } else {
                        if (BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                            r_table << BOLDFONT(FRED("   Not achieved"));
                        else
                            r_table << "Not achieved";
                    }
                }
            } else {
                if (active_set_converged) {
                    if (BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        KRATOS_INFO("ALMFrictionalMortarConvergenceCriteria") << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                    else
                        KRATOS_INFO("ALMFrictionalMortarConvergenceCriteria") << "\tActive set convergence is achieved" << std::endl;
                } else {
                    if (BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                        KRATOS_INFO("ALMFrictionalMortarConvergenceCriteria") << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    else
                        KRATOS_INFO("ALMFrictionalMortarConvergenceCriteria") << "\tActive set convergence is not achieved" << std::endl;
                }

                if (BaseType::mOptions.IsNot(BaseType::PURE_SLIP)) {
                    if (slip_set_converged) {
                        if (BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("ALMFrictionalMortarConvergenceCriteria") << BOLDFONT("\tSlip/stick set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        else
                            KRATOS_INFO("ALMFrictionalMortarConvergenceCriteria") << "\tSlip/stick set convergence is achieved" << std::endl;
                    } else {
                        if (BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("ALMFrictionalMortarConvergenceCriteria") << BOLDFONT("\tSlip/stick set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                        else
                            KRATOS_INFO("ALMFrictionalMortarConvergenceCriteria") << "\tSlip/stick set  convergence is not achieved" << std::endl;
                    }
                }
            }
        }

        return (active_set_converged && slip_set_converged);
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        // Calling base criteria
        BaseType::Initialize(rModelPart);

        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info.Has(TABLE_UTILITY) && BaseType::mOptions.IsNot(ALMFrictionalMortarConvergenceCriteria::TABLE_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& r_table = p_table->GetTable();
            r_table.AddColumn("ACTIVE SET CONV", 15);
            if (BaseType::mOptions.IsNot(BaseType::PURE_SLIP)) {
                r_table.AddColumn("SLIP/STICK CONV", 15);
            }
            BaseType::mOptions.Set(ALMFrictionalMortarConvergenceCriteria::TABLE_IS_INITIALIZED, true);
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
        NodesArrayType& nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        VariableUtils().SetVariable(WEIGHTED_GAP, 0.0, nodes_array);
        VariableUtils().SetVariable(WEIGHTED_SLIP, zero_array, nodes_array);
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

}; // Class ALMFrictionalMortarConvergenceCriteria

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags ALMFrictionalMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::PRINTING_OUTPUT(Kratos::Flags::Create(3));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags ALMFrictionalMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(4));

}  // namespace Kratos

#endif /* KRATOS_ALM_FRICTIONAL_MORTAR_CRITERIA_H  defined */
