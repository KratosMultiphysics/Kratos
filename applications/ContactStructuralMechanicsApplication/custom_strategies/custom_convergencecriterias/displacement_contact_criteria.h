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

#if !defined(KRATOS_DISPLACEMENT_CONTACT_CRITERIA_H)
#define KRATOS_DISPLACEMENT_CONTACT_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/table_stream_utility.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
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

///@name Kratos Classes
///@{

/**
 * @class DisplacementContactCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Convergence criteria for contact problems
 * @details This class implements a convergence control based on nodal displacement (for penalty contact)
 * @author Vicente Mataix Ferrandiz
 */
template<   class TSparseSpace,
            class TDenseSpace >
class DisplacementContactCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of DisplacementContactCriteria
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementContactCriteria );

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( PRINTING_OUTPUT );
    KRATOS_DEFINE_LOCAL_FLAG( TABLE_IS_INITIALIZED );
    KRATOS_DEFINE_LOCAL_FLAG( ROTATION_DOF_IS_CONSIDERED );

    /// The base class definition (and it subclasses)
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;
    typedef typename BaseType::TDataType                    TDataType;
    typedef typename BaseType::DofsArrayType            DofsArrayType;
    typedef typename BaseType::TSystemMatrixType    TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType    TSystemVectorType;

    /// The sparse space used
    typedef TSparseSpace                              SparseSpaceType;

    /// The table stream definition TODO: Replace by logger
    typedef TableStreamUtility::Pointer       TablePrinterPointerType;

    /// The index type definition
    typedef std::size_t                                     IndexType;

    /// The key type definition
    typedef std::size_t                                       KeyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param DispRatioTolerance Relative tolerance for displacement error
     * @param DispAbsTolerance Absolute tolerance for displacement error
     * @param RotRatioTolerance Relative tolerance for rotation error
     * @param RotAbsTolerance Absolute tolerance for rotation error
     * @param pTable The pointer to the output table
     * @param PrintingOutput If the output is going to be printed in a txt file
     */
    explicit DisplacementContactCriteria(
        const TDataType DispRatioTolerance,
        const TDataType DispAbsTolerance,
        const TDataType RotRatioTolerance,
        const TDataType RotAbsTolerance,
        const bool PrintingOutput = false
        )
        : BaseType()
    {
        // Set local flags
        mOptions.Set(DisplacementContactCriteria::PRINTING_OUTPUT, PrintingOutput);
        mOptions.Set(DisplacementContactCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED, false);

        // The displacement solution
        mDispRatioTolerance = DispRatioTolerance;
        mDispAbsTolerance = DispAbsTolerance;

        // The rotation solution
        mRotRatioTolerance = RotRatioTolerance;
        mRotAbsTolerance = RotAbsTolerance;
    }

    /**
     * @brief Default constructor (parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit DisplacementContactCriteria( Parameters ThisParameters = Parameters(R"({})"))
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    // Copy constructor.
    DisplacementContactCriteria( DisplacementContactCriteria const& rOther )
      :BaseType(rOther)
      ,mOptions(rOther.mOptions)
      ,mDispRatioTolerance(rOther.mDispRatioTolerance)
      ,mDispAbsTolerance(rOther.mDispAbsTolerance)
      ,mRotRatioTolerance(rOther.mRotRatioTolerance)
      ,mRotAbsTolerance(rOther.mRotAbsTolerance)
    {
    }

    /// Destructor.
    ~DisplacementContactCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{

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
        if (SparseSpaceType::Size(rDx) != 0) { //if we are solving for something
            // Initialize
            TDataType disp_solution_norm = 0.0, disp_increase_norm = 0.0;
            IndexType disp_dof_num(0);
            TDataType rot_solution_norm = 0.0, rot_increase_norm = 0.0;
            IndexType rot_dof_num(0);

            // First iterator
            const auto it_dof_begin = rDofSet.begin();

            // Auxiliar values
            std::size_t dof_id = 0;
            TDataType dof_value = 0.0, dof_incr = 0.0;

            // Auxiliar displacement DoF check
            const std::function<bool(const VariableData&)> check_without_rot =
            [](const VariableData& rCurrVar) -> bool {return true;};
            const std::function<bool(const VariableData&)> check_with_rot =
            [](const VariableData& rCurrVar) -> bool {return ((rCurrVar == DISPLACEMENT_X) || (rCurrVar == DISPLACEMENT_Y) || (rCurrVar == DISPLACEMENT_Z));};
            const auto* p_check_disp = (mOptions.Is(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED)) ? &check_with_rot : &check_without_rot;

            // Loop over Dofs
            #pragma omp parallel for reduction(+:disp_solution_norm,disp_increase_norm,disp_dof_num,rot_solution_norm,rot_increase_norm,rot_dof_num,dof_id,dof_value,dof_incr)
            for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
                auto it_dof = it_dof_begin + i;

                if (it_dof->IsFree()) {
                    dof_id = it_dof->EquationId();
                    dof_value = it_dof->GetSolutionStepValue(0);
                    dof_incr = rDx[dof_id];

                    const auto& r_curr_var = it_dof->GetVariable();
                    if ((*p_check_disp)(r_curr_var)) {
                        disp_solution_norm += std::pow(dof_value, 2);
                        disp_increase_norm += std::pow(dof_incr, 2);
                        ++disp_dof_num;
                    } else {
                        KRATOS_DEBUG_ERROR_IF_NOT((r_curr_var == ROTATION_X) || (r_curr_var == ROTATION_Y) || (r_curr_var == ROTATION_Z)) << "Variable must be a ROTATION and it is: " << r_curr_var.Name() << std::endl;
                        rot_solution_norm += std::pow(dof_value, 2);
                        rot_increase_norm += std::pow(dof_incr, 2);
                        ++rot_dof_num;
                    }
                }
            }

            if(disp_increase_norm == 0.0) disp_increase_norm = 1.0;
            if(disp_solution_norm == 0.0) disp_solution_norm = 1.0;
            if(rot_increase_norm == 0.0) rot_increase_norm = 1.0;
            if(rot_solution_norm == 0.0) rot_solution_norm = 1.0;

            const TDataType disp_ratio = std::sqrt(disp_increase_norm/disp_solution_norm);

            const TDataType disp_abs = std::sqrt(disp_increase_norm)/ static_cast<TDataType>(disp_dof_num);

            const TDataType rot_ratio = std::sqrt(rot_increase_norm/rot_solution_norm);

            const TDataType rot_abs = std::sqrt(rot_increase_norm)/ static_cast<TDataType>(rot_dof_num);

            // The process info of the model part
            ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

            // We print the results  // TODO: Replace for the new log
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                if (r_process_info.Has(TABLE_UTILITY)) {
                    std::cout.precision(4);
                    TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                    auto& Table = p_table->GetTable();
                    if (mOptions.Is(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                        Table << disp_ratio << mDispRatioTolerance << disp_abs << mDispAbsTolerance << rot_ratio << mRotRatioTolerance << rot_abs << mRotAbsTolerance;
                    } else {
                        Table << disp_ratio << mDispRatioTolerance << disp_abs << mDispAbsTolerance;
                    }
                } else {
                    std::cout.precision(4);
                    if (mOptions.IsNot(DisplacementContactCriteria::PRINTING_OUTPUT)) {
                        KRATOS_INFO("DisplacementContactCriteria") << BOLDFONT("DoF ONVERGENCE CHECK") << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl;
                        KRATOS_INFO("DisplacementContactCriteria") << BOLDFONT("\tDISPLACEMENT: RATIO = ") << disp_ratio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << disp_abs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        if (mOptions.Is(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                            KRATOS_INFO("DisplacementContactCriteria") << BOLDFONT("\tROTATION: RATIO = ") << rot_ratio << BOLDFONT(" EXP.RATIO = ") << mRotRatioTolerance << BOLDFONT(" ABS = ") << rot_abs << BOLDFONT(" EXP.ABS = ") << mRotAbsTolerance << std::endl;
                        }
                    } else {
                        KRATOS_INFO("DisplacementContactCriteria") << "DoF ONVERGENCE CHECK" << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl;
                        KRATOS_INFO("DisplacementContactCriteria") << "\tDISPLACEMENT: RATIO = " << disp_ratio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << disp_abs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        if (mOptions.Is(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                            KRATOS_INFO("DisplacementContactCriteria") << "\tROTATION: RATIO = " << rot_ratio << " EXP.RATIO = " << mRotRatioTolerance << " ABS = " << rot_abs << " EXP.ABS = " << mRotAbsTolerance << std::endl;
                        }
                    }
                }
            }

            // We check if converged
            const bool disp_converged = (disp_ratio <= mDispRatioTolerance || disp_abs <= mDispAbsTolerance);
            const bool rot_converged = mOptions.Is(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED) ? (rot_ratio <= mRotRatioTolerance || rot_abs <= mRotAbsTolerance) : true;

            if (disp_converged && rot_converged) {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& table = p_table->GetTable();
                        if (mOptions.IsNot(DisplacementContactCriteria::PRINTING_OUTPUT))
                            table << BOLDFONT(FGRN("       Achieved"));
                        else
                            table << "Achieved";
                    } else {
                        if (mOptions.IsNot(DisplacementContactCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("DisplacementContactCriteria") << BOLDFONT("\tDoF") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementContactCriteria") << "\tDoF convergence is achieved" << std::endl;
                    }
                }
                return true;
            } else {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& table = p_table->GetTable();
                        if (mOptions.IsNot(DisplacementContactCriteria::PRINTING_OUTPUT))
                            table << BOLDFONT(FRED("   Not achieved"));
                        else
                            table << "Not achieved";
                    } else {
                        if (mOptions.IsNot(DisplacementContactCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("DisplacementContactCriteria") << BOLDFONT("\tDoF") << " convergence is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementContactCriteria") << "\tDoF convergence is not achieved" << std::endl;
                    }
                }
                return false;
            }
        }
        else // In this case all the displacements are imposed!
            return true;
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the contact problem. (unused)
     */
    void Initialize( ModelPart& rModelPart ) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;

        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info.Has(TABLE_UTILITY) && mOptions.IsNot(DisplacementContactCriteria::TABLE_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& r_table = p_table->GetTable();
            r_table.AddColumn("DP RATIO", 10);
            r_table.AddColumn("EXP. RAT", 10);
            r_table.AddColumn("ABS", 10);
            r_table.AddColumn("EXP. ABS", 10);
            if (mOptions.IsNot(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                r_table.AddColumn("RT RATIO", 10);
                r_table.AddColumn("EXP. RAT", 10);
                r_table.AddColumn("ABS", 10);
                r_table.AddColumn("EXP. ABS", 10);
            }
            r_table.AddColumn("CONVERGENCE", 15);
            mOptions.Set(DisplacementContactCriteria::TABLE_IS_INITIALIZED, true);
        }

        // Check rotation dof
        mOptions.Set(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED, ContactUtilities::CheckModelPartHasRotationDoF(rModelPart));
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                            : "displacement_contact_criteria",
            "ensure_contact"                  : false,
            "print_convergence_criterion"     : false,
            "displacement_relative_tolerance" : 1.0e-4,
            "displacement_absolute_tolerance" : 1.0e-9,
            "rotation_relative_tolerance"     : 1.0e-4,
            "rotation_absolute_tolerance"     : 1.0e-9
        })");

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
        return "displacement_contact_criteria";
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
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // The displacement solution
        mDispRatioTolerance = ThisParameters["displacement_relative_tolerance"].GetDouble();
        mDispAbsTolerance = ThisParameters["displacement_absolute_tolerance"].GetDouble();

        // The rotation solution
        mRotRatioTolerance = ThisParameters["rotation_relative_tolerance"].GetDouble();
        mRotAbsTolerance = ThisParameters["rotation_absolute_tolerance"].GetDouble();

        // Set local flags
        mOptions.Set(DisplacementContactCriteria::PRINTING_OUTPUT, ThisParameters["print_convergence_criterion"].GetBool());
        mOptions.Set(DisplacementContactCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(DisplacementContactCriteria::ROTATION_DOF_IS_CONSIDERED, false);
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

    TDataType mDispRatioTolerance; /// The ratio threshold for the norm of the displacement
    TDataType mDispAbsTolerance;   /// The absolute value threshold for the norm of the displacement

    TDataType mRotRatioTolerance; /// The ratio threshold for the norm of the rotation
    TDataType mRotAbsTolerance;   /// The absolute value threshold for the norm of the rotation

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

    ///@name Unaccessible methods
    ///@{
    ///@}
};  // Kratos DisplacementContactCriteria

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementContactCriteria<TSparseSpace, TDenseSpace>::PRINTING_OUTPUT(Kratos::Flags::Create(1));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementContactCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(2));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementContactCriteria<TSparseSpace, TDenseSpace>::ROTATION_DOF_IS_CONSIDERED(Kratos::Flags::Create(3));
}

#endif	/* KRATOS_DISPLACEMENT_CONTACT_CRITERIA_H */
