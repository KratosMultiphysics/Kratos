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

#if !defined(KRATOS_ALM_FRICTIONLESS_MORTAR_CRITERIA_H)
#define  KRATOS_ALM_FRICTIONLESS_MORTAR_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/table_stream_utility.h"
#include "custom_strategies/custom_convergencecriterias/base_mortar_criteria.h"
#include "utilities/color_utilities.h"

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
 * @class ALMFrictionlessMortarConvergenceCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom convergence criteria for the mortar condition for frictionless case with components
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class ALMFrictionlessMortarConvergenceCriteria
    : public  BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ALMFrictionlessMortarConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION( ALMFrictionlessMortarConvergenceCriteria );

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

    /// The table stream definition TODO: Replace by logger
    typedef TableStreamUtility::Pointer                          TablePrinterPointerType;

    /// The index type definition
    typedef std::size_t                                                        IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    explicit ALMFrictionlessMortarConvergenceCriteria(
        const bool PrintingOutput = false,
        const bool GiDIODebug = false
        ) : BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >(GiDIODebug),
        mPrintingOutput(PrintingOutput),
        mTableIsInitialized(false)
    {
    }

    ///Copy constructor
    ALMFrictionlessMortarConvergenceCriteria( ALMFrictionlessMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
      ,mPrintingOutput(rOther.mPrintingOutput)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
    {
    }

    /// Destructor
    ~ALMFrictionlessMortarConvergenceCriteria() override = default;

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

        // Defining the convergence
        IndexType is_converged = 0;

        // We get the process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        // We check the active/inactive set during the first non-linear iteration or for the general semi-smooth case
        if (rModelPart.Is(INTERACTION) || r_process_info[NL_ITERATION_NUMBER] == 1) {
            const double common_epsilon = r_process_info[INITIAL_PENALTY];
            const double scale_factor = r_process_info[SCALE_FACTOR];

            NodesArrayType& nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();

            #pragma omp parallel for reduction(+:is_converged)
            for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
                auto it_node = nodes_array.begin() + i;

                const double epsilon = it_node->Has(INITIAL_PENALTY) ? it_node->GetValue(INITIAL_PENALTY) : common_epsilon;

                const double augmented_normal_pressure = scale_factor * it_node->FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) + epsilon * it_node->FastGetSolutionStepValue(WEIGHTED_GAP);

                it_node->SetValue(AUGMENTED_NORMAL_CONTACT_PRESSURE, augmented_normal_pressure); // NOTE: This value is purely for debugging interest (to see the "effective" pressure)

                if (augmented_normal_pressure < 0.0) { // NOTE: This could be conflictive (< or <=)
                    if (it_node->Is(ACTIVE) == false ) {
                        it_node->Set(ACTIVE, true);
                        is_converged += 1;
                    }
                } else {
                    if (it_node->Is(ACTIVE) == true ) {
                        it_node->Set(ACTIVE, false);
                        is_converged += 1;
                    }
                }
            }
        }

        // We save to the process info if the active set has converged
        const bool active_set_converged = (is_converged == 0 ? true : false);
        r_process_info[ACTIVE_SET_CONVERGED] = active_set_converged;

        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
            if (r_process_info.Has(TABLE_UTILITY)) {
                TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                auto& table = p_table->GetTable();
                if (active_set_converged) {
                    if (mPrintingOutput == false)
                        table << BOLDFONT(FGRN("       Achieved"));
                    else
                        table << "Achieved";
                } else {
                    if (mPrintingOutput == false)
                        table << BOLDFONT(FRED("   Not achieved"));
                    else
                        table << "Not achieved";
                }
            } else {
                if (active_set_converged) {
                    if (mPrintingOutput == false)
                        KRATOS_INFO("ALMFrictionlessMortarConvergenceCriteria")  << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                } else {
                    if (mPrintingOutput == false)
                        KRATOS_INFO("ALMFrictionlessMortarConvergenceCriteria")  << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    else
                        KRATOS_INFO("ALMFrictionlessMortarConvergenceCriteria")  << "\tActive set convergence is not achieved" << std::endl;
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
        if (r_process_info.Has(TABLE_UTILITY) && mTableIsInitialized == false) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& table = p_table->GetTable();
            table.AddColumn("ACTIVE SET CONV", 15);
            mTableIsInitialized = true;
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

    TablePrinterPointerType p_table; /// Pointer to the fancy table
    bool mPrintingOutput;            /// If the colors and bold are printed
    bool mTableIsInitialized;        /// If the table is already initialized

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

}; // Class ALMFrictionlessMortarConvergenceCriteria

///@name Explicit Specializations
///@{

}  // namespace Kratos

#endif /* KRATOS_ALM_FRICTIONLESS_MORTAR_CRITERIA_H  defined */

