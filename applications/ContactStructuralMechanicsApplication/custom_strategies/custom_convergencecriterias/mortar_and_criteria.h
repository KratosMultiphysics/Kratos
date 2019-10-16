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

#if !defined(KRATOS_MORTAR_AND_CRITERIA_H)
#define  KRATOS_MORTAR_AND_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/table_stream_utility.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "utilities/color_utilities.h"
#include "utilities/condition_number_utility.h"

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
 * @class MortarAndConvergenceCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom AND convergence criteria for the mortar condition
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,
         class TDenseSpace
         >
class MortarAndConvergenceCriteria
    : public And_Criteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MortarAndConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION( MortarAndConvergenceCriteria );

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( PRINTING_OUTPUT );
    KRATOS_DEFINE_LOCAL_FLAG( TABLE_IS_INITIALIZED );
    KRATOS_DEFINE_LOCAL_FLAG( CONDITION_NUMBER_IS_INITIALIZED );

    /// The base convergence criteria class definition
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > ConvergenceCriteriaBaseType;

    /// The base class definition (and it subclasses)
    typedef And_Criteria< TSparseSpace, TDenseSpace >                           BaseType;
    typedef typename BaseType::TDataType                                       TDataType;
    typedef typename BaseType::DofsArrayType                               DofsArrayType;
    typedef typename BaseType::TSystemMatrixType                       TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType                       TSystemVectorType;

    /// The sparse space used (and it subclasses)
    typedef TSparseSpace                                                 SparseSpaceType;
    typedef typename TSparseSpace::MatrixType                           SparseMatrixType;
    typedef typename TSparseSpace::VectorType                           SparseVectorType;
    typedef typename TDenseSpace::MatrixType                             DenseMatrixType;
    typedef typename TDenseSpace::VectorType                             DenseVectorType;

    /// The table stream definition TODO: Replace by logger
    typedef TableStreamUtility::Pointer                          TablePrinterPointerType;

    /// The index type definition
    typedef std::size_t                                                        IndexType;

    /// The condition number utility pointer definition
    typedef ConditionNumberUtility::Pointer            ConditionNumberUtilityPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit MortarAndConvergenceCriteria(
        typename ConvergenceCriteriaBaseType::Pointer pFirstCriterion,
        typename ConvergenceCriteriaBaseType::Pointer pSecondCriterion,
        const bool PrintingOutput = false,
        ConditionNumberUtilityPointerType pConditionNumberUtility = nullptr
        )
        :BaseType(pFirstCriterion, pSecondCriterion),
        mpConditionNumberUtility(pConditionNumberUtility)
    {
        // Set local flags
        mOptions.Set(MortarAndConvergenceCriteria::PRINTING_OUTPUT, PrintingOutput);
        mOptions.Set(MortarAndConvergenceCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(MortarAndConvergenceCriteria::CONDITION_NUMBER_IS_INITIALIZED, false);
    }

    /**
     * Copy constructor.
     */
    MortarAndConvergenceCriteria(MortarAndConvergenceCriteria const& rOther)
        :BaseType(rOther)
        ,mOptions(rOther.mOptions)
        ,mpConditionNumberUtility(rOther.mpConditionNumberUtility)
     {
     }

    /** Destructor.
    */
    ~MortarAndConvergenceCriteria () override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Criteria that need to be called after getting the solution
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
        // The process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
            if (r_process_info.Has(TABLE_UTILITY)) {
                TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                p_table->AddToRow<IndexType>(rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER]);
            }
        }

        bool criterion_result = BaseType::PostCriteria(rModelPart, rDofSet, rA, rDx, rb);

        if (mpConditionNumberUtility != nullptr) {
            TSystemMatrixType copy_A(rA); // NOTE: Can not be const, TODO: Change the solvers to const
            const double condition_number = mpConditionNumberUtility->GetConditionNumber(copy_A);

            if (r_process_info.Has(TABLE_UTILITY)) {
                std::cout.precision(4);
                TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                auto& r_table = p_table->GetTable();
                r_table  << condition_number;
            } else {
                if (mOptions.IsNot(MortarAndConvergenceCriteria::PRINTING_OUTPUT))
                    KRATOS_INFO("MortarAndConvergenceCriteria") << "\n" << BOLDFONT("CONDITION NUMBER:") << "\t " << std::scientific << condition_number << std::endl;
                else
                    KRATOS_INFO("MortarAndConvergenceCriteria") << "\n" << "CONDITION NUMBER:" << "\t" << std::scientific << condition_number << std::endl;
            }
        }

        if (criterion_result == true && rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            if (r_process_info.Has(TABLE_UTILITY)) {
                TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                p_table->PrintFooter();
            }

        return criterion_result;
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        // The process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        if (r_process_info.Has(TABLE_UTILITY) && mOptions.IsNot(MortarAndConvergenceCriteria::TABLE_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            (p_table->GetTable()).SetBold(mOptions.IsNot(MortarAndConvergenceCriteria::PRINTING_OUTPUT));
            (p_table->GetTable()).AddColumn("ITER", 4);
        }

        mOptions.Set(MortarAndConvergenceCriteria::TABLE_IS_INITIALIZED, true);
        BaseType::Initialize(rModelPart);

        if (r_process_info.Has(TABLE_UTILITY) && mpConditionNumberUtility != nullptr
            && mOptions.IsNot(MortarAndConvergenceCriteria::CONDITION_NUMBER_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            (p_table->GetTable()).AddColumn("COND.NUM.", 10);
            mOptions.Set(MortarAndConvergenceCriteria::CONDITION_NUMBER_IS_INITIALIZED, true);
        }
    }

    /**
     * @brief This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // The process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
            std::cout.precision(4);
            if (mOptions.IsNot(MortarAndConvergenceCriteria::PRINTING_OUTPUT))
                std::cout << "\n\n" << BOLDFONT("CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[STEP] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            else
                std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[STEP] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;

            if (r_process_info.Has(TABLE_UTILITY)) {
                TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                p_table->PrintHeader();
            }
        }

        BaseType::InitializeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);
    }

    /**
     * @brief This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        BaseType::FinalizeSolutionStep(rModelPart,rDofSet, rA, rDx, rb);
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

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

    Flags mOptions; /// Local flags

    ConditionNumberUtilityPointerType mpConditionNumberUtility; /// The utility to compute the condition number

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

};  // Kratos MortarAndConvergenceCriteria

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags MortarAndConvergenceCriteria<TSparseSpace, TDenseSpace>::PRINTING_OUTPUT(Kratos::Flags::Create(0));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags MortarAndConvergenceCriteria<TSparseSpace, TDenseSpace>::NOT_PRINTING_OUTPUT(Kratos::Flags::Create(0, false));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags MortarAndConvergenceCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(1));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags MortarAndConvergenceCriteria<TSparseSpace, TDenseSpace>::NOT_TABLE_IS_INITIALIZED(Kratos::Flags::Create(1, false));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags MortarAndConvergenceCriteria<TSparseSpace, TDenseSpace>::CONDITION_NUMBER_IS_INITIALIZED(Kratos::Flags::Create(2));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags MortarAndConvergenceCriteria<TSparseSpace, TDenseSpace>::NOT_CONDITION_NUMBER_IS_INITIALIZED(Kratos::Flags::Create(2, false));

}  /* namespace Kratos.*/

#endif /* KRATOS_MORTAR_AND_CRITERIA_H  defined */

