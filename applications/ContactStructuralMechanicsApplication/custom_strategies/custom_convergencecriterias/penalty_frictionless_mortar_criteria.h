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

#if !defined(KRATOS_PENALTY_FRICTIONLESS_MORTAR_CRITERIA_H)
#define  KRATOS_PENALTY_FRICTIONLESS_MORTAR_CRITERIA_H

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
 * @class PenaltyFrictionlessMortarConvergenceCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom convergence criteria for the mortar condition for frictionless case with components
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace, class TDenseSpace>
class PenaltyFrictionlessMortarConvergenceCriteria
    : public BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PenaltyFrictionlessMortarConvergenceCriteria
    KRATOS_CLASS_POINTER_DEFINITION( PenaltyFrictionlessMortarConvergenceCriteria );

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( PRINTING_OUTPUT );
    KRATOS_DEFINE_LOCAL_FLAG( TABLE_IS_INITIALIZED );

    /// The base convergence criteria class definition
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace >        ConvergenceCriteriaBaseType;

    /// The base class definition
    typedef BaseMortarConvergenceCriteria< TSparseSpace, TDenseSpace >                 BaseType;

    /// The definition of the current class
    typedef PenaltyFrictionlessMortarConvergenceCriteria< TSparseSpace, TDenseSpace > ClassType;

    /// The dofs array type
    typedef typename BaseType::DofsArrayType                                      DofsArrayType;

    /// The sparse matrix type
    typedef typename BaseType::TSystemMatrixType                              TSystemMatrixType;

    /// The dense vector type
    typedef typename BaseType::TSystemVectorType                              TSystemVectorType;

    /// The table stream definition TODO: Replace by logger
    typedef TableStreamUtility::Pointer                                 TablePrinterPointerType;

    /// The index type definition
    typedef std::size_t                                                               IndexType;

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    explicit PenaltyFrictionlessMortarConvergenceCriteria(
        const bool PrintingOutput = false,
        const bool ComputeDynamicFactor = true,
        const bool GiDIODebug = false
        ) : BaseType(ComputeDynamicFactor, GiDIODebug)
    {
        // Set local flags
        BaseType::mOptions.Set(PenaltyFrictionlessMortarConvergenceCriteria::PRINTING_OUTPUT, PrintingOutput);
        BaseType::mOptions.Set(PenaltyFrictionlessMortarConvergenceCriteria::TABLE_IS_INITIALIZED, false);
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit PenaltyFrictionlessMortarConvergenceCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    ///Copy constructor
    PenaltyFrictionlessMortarConvergenceCriteria( PenaltyFrictionlessMortarConvergenceCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /// Destructor
    ~PenaltyFrictionlessMortarConvergenceCriteria() override = default;

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
        const IndexType is_converged = ActiveSetUtilities::ComputePenaltyFrictionlessActiveSet(rModelPart);

        // We save to the process info if the active set has converged
        const bool active_set_converged = (is_converged == 0 ? true : false);

        // We get the process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        r_process_info[ACTIVE_SET_CONVERGED] = active_set_converged;

        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
            if (r_process_info.Has(TABLE_UTILITY)) {
                TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                auto& r_table = p_table->GetTable();
                if (active_set_converged) {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionlessMortarConvergenceCriteria::PRINTING_OUTPUT))
                        r_table << BOLDFONT(FGRN("       Achieved"));
                    else
                        r_table << "Achieved";
                } else {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionlessMortarConvergenceCriteria::PRINTING_OUTPUT))
                        r_table << BOLDFONT(FRED("   Not achieved"));
                    else
                        r_table << "Not achieved";
                }
            } else {
                if (active_set_converged) {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionlessMortarConvergenceCriteria::PRINTING_OUTPUT))
                        KRATOS_INFO("PenaltyFrictionlessMortarConvergenceCriteria")  << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                } else {
                    if (BaseType::mOptions.IsNot(PenaltyFrictionlessMortarConvergenceCriteria::PRINTING_OUTPUT))
                        KRATOS_INFO("PenaltyFrictionlessMortarConvergenceCriteria")  << BOLDFONT("\tActive set") << " convergence is " << BOLDFONT(FRED("not achieved")) << std::endl;
                    else
                        KRATOS_INFO("PenaltyFrictionlessMortarConvergenceCriteria")  << "\tActive set convergence is not achieved" << std::endl;
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
        if (r_process_info.Has(TABLE_UTILITY) && BaseType::mOptions.IsNot(PenaltyFrictionlessMortarConvergenceCriteria::TABLE_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& r_table = p_table->GetTable();
            r_table.AddColumn("ACTIVE SET CONV", 15);
            BaseType::mOptions.Set(PenaltyFrictionlessMortarConvergenceCriteria::TABLE_IS_INITIALIZED, true);
        }
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                        : "penalty_frictionless_mortar_criteria",
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
        return "penalty_frictionless_mortar_criteria";
    }

    ///@}
    ///@name Acces
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
        return "PenaltyFrictionlessMortarConvergenceCriteria";
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
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // Set local flags
        BaseType::mOptions.Set(PenaltyFrictionlessMortarConvergenceCriteria::PRINTING_OUTPUT, ThisParameters["print_convergence_criterion"].GetBool());
        BaseType::mOptions.Set(PenaltyFrictionlessMortarConvergenceCriteria::TABLE_IS_INITIALIZED, false);
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

}; // Class PenaltyFrictionlessMortarConvergenceCriteria

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags PenaltyFrictionlessMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::PRINTING_OUTPUT(Kratos::Flags::Create(3));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags PenaltyFrictionlessMortarConvergenceCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(4));

}  // namespace Kratos

#endif /* KRATOS_PENALTY_FRICTIONLESS_MORTAR_CRITERIA_H  defined */
