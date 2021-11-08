// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
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

    /// The base class definition
    typedef And_Criteria< TSparseSpace, TDenseSpace >                           BaseType;

    /// The definition of the current class
    typedef MortarAndConvergenceCriteria< TSparseSpace, TDenseSpace >          ClassType;

    /// The dofs array type
    typedef typename BaseType::DofsArrayType                               DofsArrayType;

    /// The sparse matrix type
    typedef typename BaseType::TSystemMatrixType                       TSystemMatrixType;

    /// The dense vector type
    typedef typename BaseType::TSystemVectorType                       TSystemVectorType;

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
     * @brief Default constructor.
     */
    explicit MortarAndConvergenceCriteria()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit MortarAndConvergenceCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

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

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename ConvergenceCriteriaBaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

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

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                        : "mortar_and_criteria",
            "print_convergence_criterion" : false
        })" );

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "mortar_and_criteria";
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
        return "MortarAndConvergenceCriteria";
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

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // Defining condition number utility
        mpConditionNumberUtility = nullptr; // TODO: Define a factory

        // Set local flags
        mOptions.Set(MortarAndConvergenceCriteria::PRINTING_OUTPUT, ThisParameters["print_convergence_criterion"].GetBool());
        mOptions.Set(MortarAndConvergenceCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(MortarAndConvergenceCriteria::CONDITION_NUMBER_IS_INITIALIZED, false);
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
const Kratos::Flags MortarAndConvergenceCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(1));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags MortarAndConvergenceCriteria<TSparseSpace, TDenseSpace>::CONDITION_NUMBER_IS_INITIALIZED(Kratos::Flags::Create(2));

}  /* namespace Kratos.*/

#endif /* KRATOS_MORTAR_AND_CRITERIA_H  defined */
