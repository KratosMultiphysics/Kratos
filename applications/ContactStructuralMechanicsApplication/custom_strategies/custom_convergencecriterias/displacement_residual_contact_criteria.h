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

#if !defined(KRATOS_DISPLACEMENT_RESIDUAL_CONTACT_CRITERIA_H)
#define KRATOS_DISPLACEMENT_RESIDUAL_CONTACT_CRITERIA_H

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
 * @class DisplacementResidualContactCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Convergence criteria for contact problems
 * This class implements a convergence control based on nodal displacement (for penalty contact)
 * @author Vicente Mataix Ferrandiz
 */
template<   class TSparseSpace,
            class TDenseSpace >
class DisplacementResidualContactCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of DisplacementResidualContactCriteria
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementResidualContactCriteria );

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( PRINTING_OUTPUT );
    KRATOS_DEFINE_LOCAL_FLAG( TABLE_IS_INITIALIZED );
    KRATOS_DEFINE_LOCAL_FLAG( ROTATION_DOF_IS_CONSIDERED );
    KRATOS_DEFINE_LOCAL_FLAG( INITIAL_RESIDUAL_IS_SET );

    /// The base class definition (and it subclasses)
    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;
    typedef typename BaseType::TDataType                    TDataType;
    typedef typename BaseType::DofsArrayType            DofsArrayType;
    typedef typename BaseType::TSystemMatrixType    TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType    TSystemVectorType;

    /// The sparse space used
    typedef TSparseSpace                              SparseSpaceType;

    /// The r_table stream definition TODO: Replace by logger
    typedef TableStreamUtility::Pointer       TablePrinterPointerType;

    /// The index type definition
    typedef std::size_t                                     IndexType;

    /// The key type definition
    typedef std::size_t                                       KeyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param DispRatioTolerance Relative tolerance for displacement residual error
     * @param DispAbsTolerance Absolute tolerance for displacement residual error
     * @param RotRatioTolerance Relative tolerance for rotation residual error
     * @param RotAbsTolerance Absolute tolerance for rotation residual error
     * @param pTable The pointer to the output r_table
     * @param PrintingOutput If the output is going to be printed in a txt file
     */
    explicit DisplacementResidualContactCriteria(
        const TDataType DispRatioTolerance,
        const TDataType DispAbsTolerance,
        const TDataType RotRatioTolerance,
        const TDataType RotAbsTolerance,
        const bool PrintingOutput = false
        )
        : BaseType()
    {
        // Set local flags
        mOptions.Set(DisplacementResidualContactCriteria::PRINTING_OUTPUT, PrintingOutput);
        mOptions.Set(DisplacementResidualContactCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED, false);
        mOptions.Set(DisplacementResidualContactCriteria::INITIAL_RESIDUAL_IS_SET, false);

        // The displacement residual
        mDispRatioTolerance = DispRatioTolerance;
        mDispAbsTolerance = DispAbsTolerance;

        // The rotation residual
        mRotRatioTolerance = RotRatioTolerance;
        mRotAbsTolerance = RotAbsTolerance;
    }

    /**
     * @brief Default constructor (parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit DisplacementResidualContactCriteria( Parameters ThisParameters = Parameters(R"({})"))
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    //* Copy constructor.
    DisplacementResidualContactCriteria( DisplacementResidualContactCriteria const& rOther )
      :BaseType(rOther)
      ,mOptions(rOther.mOptions)
      ,mDispRatioTolerance(rOther.mDispRatioTolerance)
      ,mDispAbsTolerance(rOther.mDispAbsTolerance)
      ,mDispInitialResidualNorm(rOther.mDispInitialResidualNorm)
      ,mDispCurrentResidualNorm(rOther.mDispCurrentResidualNorm)
      ,mRotRatioTolerance(rOther.mRotRatioTolerance)
      ,mRotAbsTolerance(rOther.mRotAbsTolerance)
      ,mRotInitialResidualNorm(rOther.mRotInitialResidualNorm)
      ,mRotCurrentResidualNorm(rOther.mRotCurrentResidualNorm)
    {
    }

    /// Destructor.
    ~DisplacementResidualContactCriteria() override = default;

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
        if (SparseSpaceType::Size(rb) != 0) { //if we are solving for something
            // Initialize
            TDataType disp_residual_solution_norm = 0.0;
            IndexType disp_dof_num(0);
            TDataType rot_residual_solution_norm = 0.0;
            IndexType rot_dof_num(0);

            // First iterator
            const auto it_dof_begin = rDofSet.begin();

            // Auxiliar values
            std::size_t dof_id = 0;
            TDataType residual_dof_value = 0.0;

            // Auxiliar displacement DoF check
            const std::function<bool(const VariableData&)> check_without_rot =
            [](const VariableData& rCurrVar) -> bool {return true;};
            const std::function<bool(const VariableData&)> check_with_rot =
            [](const VariableData& rCurrVar) -> bool {return ((rCurrVar == DISPLACEMENT_X) || (rCurrVar == DISPLACEMENT_Y) || (rCurrVar == DISPLACEMENT_Z));};
            const auto* p_check_disp = (mOptions.Is(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED)) ? &check_with_rot : &check_without_rot;

            // Loop over Dofs
            #pragma omp parallel for reduction(+:disp_residual_solution_norm,disp_dof_num,rot_residual_solution_norm,rot_dof_num,dof_id,residual_dof_value)
            for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
                auto it_dof = it_dof_begin + i;

                if (it_dof->IsFree()) {
                    dof_id = it_dof->EquationId();
                    residual_dof_value = rb[dof_id];

                    const auto& r_curr_var = it_dof->GetVariable();
                    if ((*p_check_disp)(r_curr_var)) {
                        disp_residual_solution_norm += std::pow(residual_dof_value, 2);
                        ++disp_dof_num;
                    } else { // We will assume is rotation dof
                        KRATOS_DEBUG_ERROR_IF_NOT((r_curr_var == ROTATION_X) || (r_curr_var == ROTATION_Y) || (r_curr_var == ROTATION_Z)) << "Variable must be a ROTATION and it is: " << r_curr_var.Name() << std::endl;
                        rot_residual_solution_norm += std::pow(residual_dof_value, 2);
                        ++rot_dof_num;
                    }
                }
            }

            mDispCurrentResidualNorm = disp_residual_solution_norm;
            mRotCurrentResidualNorm = rot_residual_solution_norm;

            TDataType residual_disp_ratio = 1.0;
            TDataType residual_rot_ratio = 1.0;

            // We initialize the solution
            if (mOptions.IsNot(DisplacementResidualContactCriteria::INITIAL_RESIDUAL_IS_SET)) {
                mDispInitialResidualNorm = (disp_residual_solution_norm == 0.0) ? 1.0 : disp_residual_solution_norm;
                residual_disp_ratio = 1.0;
                if (mOptions.Is(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                    mRotInitialResidualNorm = (rot_residual_solution_norm == 0.0) ? 1.0 : rot_residual_solution_norm;
                    residual_rot_ratio = 1.0;
                }
                mOptions.Set(DisplacementResidualContactCriteria::INITIAL_RESIDUAL_IS_SET, true);
            }

            // We calculate the ratio of the displacements
            residual_disp_ratio = mDispCurrentResidualNorm/mDispInitialResidualNorm;
            residual_rot_ratio = mRotCurrentResidualNorm/mRotInitialResidualNorm;

            // We calculate the absolute norms
            const TDataType residual_disp_abs = mDispCurrentResidualNorm/disp_dof_num;
            const TDataType residual_rot_abs = mRotCurrentResidualNorm/rot_dof_num;

            // The process info of the model part
            ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

            // We print the results // TODO: Replace for the new log
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                if (r_process_info.Has(TABLE_UTILITY)) {
                    std::cout.precision(4);
                    TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                    auto& r_table = p_table->GetTable();
                    if (mOptions.Is(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                        r_table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance << residual_rot_ratio << mRotRatioTolerance << residual_rot_abs << mRotAbsTolerance;
                    } else {
                        r_table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance;
                    }
                } else {
                    std::cout.precision(4);
                    if (mOptions.IsNot(DisplacementResidualContactCriteria::PRINTING_OUTPUT)) {
                        KRATOS_INFO("DisplacementResidualContactCriteria") << BOLDFONT("RESIDUAL CONVERGENCE CHECK") << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        KRATOS_INFO("DisplacementResidualContactCriteria") << BOLDFONT("\tDISPLACEMENT: RATIO = ") << residual_disp_ratio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << residual_disp_abs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        if (mOptions.Is(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                            KRATOS_INFO("DisplacementResidualContactCriteria") << BOLDFONT("\tDISPLACEMENT: RATIO = ") << residual_rot_ratio << BOLDFONT(" EXP.RATIO = ") << mRotRatioTolerance << BOLDFONT(" ABS = ") << residual_rot_abs << BOLDFONT(" EXP.ABS = ") << mRotAbsTolerance << std::endl;
                        }
                    } else {
                        KRATOS_INFO("DisplacementResidualContactCriteria") << "RESIDUAL CONVERGENCE CHECK" << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        KRATOS_INFO("DisplacementResidualContactCriteria") << "\tDISPLACEMENT: RATIO = " << residual_disp_ratio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << residual_disp_abs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        if (mOptions.Is(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                            KRATOS_INFO("DisplacementResidualContactCriteria") << "\tDISPLACEMENT: RATIO = " << residual_rot_ratio << " EXP.RATIO = " << mRotRatioTolerance << " ABS = " << residual_rot_abs << " EXP.ABS = " << mRotAbsTolerance << std::endl;
                        }
                    }
                }
            }

            r_process_info[CONVERGENCE_RATIO] = residual_disp_ratio;
            r_process_info[RESIDUAL_NORM] = residual_disp_abs;

            // We check if converged
            const bool disp_converged = (residual_disp_ratio <= mDispRatioTolerance || residual_disp_abs <= mDispAbsTolerance);
            const bool rot_converged = (mOptions.Is(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED)) ? (residual_rot_ratio <= mRotRatioTolerance || residual_rot_abs <= mRotAbsTolerance) : true;

            if (disp_converged && rot_converged) {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& r_table = p_table->GetTable();
                        if (mOptions.IsNot(DisplacementResidualContactCriteria::PRINTING_OUTPUT))
                            r_table << BOLDFONT(FGRN("       Achieved"));
                        else
                            r_table << "Achieved";
                    } else {
                        if (mOptions.IsNot(DisplacementResidualContactCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("DisplacementResidualContactCriteria") << BOLDFONT("\tResidual") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementResidualContactCriteria") << "\tResidual convergence is achieved" << std::endl;
                    }
                }
                return true;
            } else {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& r_table = p_table->GetTable();
                        if (mOptions.IsNot(DisplacementResidualContactCriteria::PRINTING_OUTPUT))
                            r_table << BOLDFONT(FRED("   Not achieved"));
                        else
                            r_table << "Not achieved";
                    } else {
                        if (mOptions.IsNot(DisplacementResidualContactCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("DisplacementResidualContactCriteria") << BOLDFONT("\tResidual") << " convergence is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementResidualContactCriteria") << "\tResidual convergence is not achieved" << std::endl;
                    }
                }
                return false;
            }
        } else // In this case all the displacements are imposed!
            return true;
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the contact problem. (unused)
     */
    void Initialize( ModelPart& rModelPart) override
    {
        // Initialize
        BaseType::mConvergenceCriteriaIsInitialized = true;

        // Check rotation dof
        mOptions.Set(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED, ContactUtilities::CheckModelPartHasRotationDoF(rModelPart));

        // Initialize header
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info.Has(TABLE_UTILITY) && mOptions.IsNot(DisplacementResidualContactCriteria::TABLE_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& r_table = p_table->GetTable();
            r_table.AddColumn("DP RATIO", 10);
            r_table.AddColumn("EXP. RAT", 10);
            r_table.AddColumn("ABS", 10);
            r_table.AddColumn("EXP. ABS", 10);
            if (mOptions.Is(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                r_table.AddColumn("RT RATIO", 10);
                r_table.AddColumn("EXP. RAT", 10);
                r_table.AddColumn("ABS", 10);
                r_table.AddColumn("EXP. ABS", 10);
            }
            r_table.AddColumn("CONVERGENCE", 15);
            mOptions.Set(DisplacementResidualContactCriteria::TABLE_IS_INITIALIZED, true);
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
        mOptions.Set(DisplacementResidualContactCriteria::INITIAL_RESIDUAL_IS_SET, false);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                 : "displacement_residual_contact_criteria",
            "ensure_contact"                       : false,
            "print_convergence_criterion"          : false,
            "residual_relative_tolerance"          : 1.0e-4,
            "residual_absolute_tolerance"          : 1.0e-9,
            "rotation_residual_relative_tolerance" : 1.0e-4,
            "rotation_residual_absolute_tolerance" : 1.0e-9
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
        return "displacement_residual_contact_criteria";
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

        // The displacement residual
        mDispRatioTolerance = ThisParameters["residual_relative_tolerance"].GetDouble();
        mDispAbsTolerance = ThisParameters["residual_absolute_tolerance"].GetDouble();

        // The rotation residual
        mRotRatioTolerance = ThisParameters["rotation_residual_relative_tolerance"].GetDouble();
        mRotAbsTolerance = ThisParameters["rotation_residual_absolute_tolerance"].GetDouble();

        // Set local flags
        mOptions.Set(DisplacementResidualContactCriteria::PRINTING_OUTPUT, ThisParameters["print_convergence_criterion"].GetBool());
        mOptions.Set(DisplacementResidualContactCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(DisplacementResidualContactCriteria::ROTATION_DOF_IS_CONSIDERED, false);
        mOptions.Set(DisplacementResidualContactCriteria::INITIAL_RESIDUAL_IS_SET, false);
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

    TDataType mDispRatioTolerance;      /// The ratio threshold for the norm of the displacement residual
    TDataType mDispAbsTolerance;        /// The absolute value threshold for the norm of the displacement residual
    TDataType mDispInitialResidualNorm; /// The reference norm of the displacement residual
    TDataType mDispCurrentResidualNorm; /// The current norm of the displacement residual

    TDataType mRotRatioTolerance;      /// The ratio threshold for the norm of the rotation residual
    TDataType mRotAbsTolerance;        /// The absolute value threshold for the norm of the rotation residual
    TDataType mRotInitialResidualNorm; /// The reference norm of the rotation residual
    TDataType mRotCurrentResidualNorm; /// The current norm of the rotation residual

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
};  // Kratos DisplacementResidualContactCriteria

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementResidualContactCriteria<TSparseSpace, TDenseSpace>::PRINTING_OUTPUT(Kratos::Flags::Create(1));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementResidualContactCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(2));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementResidualContactCriteria<TSparseSpace, TDenseSpace>::ROTATION_DOF_IS_CONSIDERED(Kratos::Flags::Create(3));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementResidualContactCriteria<TSparseSpace, TDenseSpace>::INITIAL_RESIDUAL_IS_SET(Kratos::Flags::Create(4));
}

#endif /* KRATOS_DISPLACEMENT_RESIDUAL_CONTACT_CRITERIA_H */
