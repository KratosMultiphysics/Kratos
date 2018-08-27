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

#if !defined(KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_FRICTIONAL_CONTACT_CRITERIA_H)
#define KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_FRICTIONAL_CONTACT_CRITERIA_H

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
 * @class DisplacementLagrangeMultiplierResidualFrictionalContactCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Convergence criteria for contact problems (only for frictional cases)
 * This class implements a convergence control based on nodal displacement and
 * lagrange multiplier values. The error is evaluated separately for each of them, and
 * relative and absolute tolerances for both must be specified.
 * @author Vicente Mataix Ferrandiz
 */
template<   class TSparseSpace,
            class TDenseSpace >
class DisplacementLagrangeMultiplierResidualFrictionalContactCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of DisplacementLagrangeMultiplierResidualFrictionalContactCriteria
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementLagrangeMultiplierResidualFrictionalContactCriteria );

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
     * @brief Default constructor
     * @param DispRatioTolerance Relative tolerance for displacement residual error
     * @param DispAbsTolerance Absolute tolerance for displacement residual error
     * @param LMRatioTolerance Relative tolerance for lagrange multiplier residual  error
     * @param LMAbsTolerance Absolute tolerance for lagrange multiplier residual error
     * @param EnsureContact To check if the contact is lost
     * @param pTable The pointer to the output table
     * @param PrintingOutput If the output is going to be printed in a txt file
     */
    explicit DisplacementLagrangeMultiplierResidualFrictionalContactCriteria(
        const TDataType DispRatioTolerance,
        const TDataType DispAbsTolerance,
        const TDataType LMNormalRatioTolerance,
        const TDataType LMNormalAbsTolerance,
        const TDataType LMTangentRatioTolerance,
        const TDataType LMTangentAbsTolerance,
        const bool EnsureContact = false,
        const bool PrintingOutput = false
        ) : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mEnsureContact(EnsureContact),
          mPrintingOutput(PrintingOutput),
          mTableIsInitialized(false)
    {
        // The displacement residual
        mDispRatioTolerance = DispRatioTolerance;
        mDispAbsTolerance = DispAbsTolerance;

        // The normal contact residual
        mLMNormalRatioTolerance = LMNormalRatioTolerance;
        mLMNormalAbsTolerance = LMNormalAbsTolerance;

        // The tangent contact residual
        mLMTangentRatioTolerance = LMTangentRatioTolerance;
        mLMTangentAbsTolerance = LMTangentAbsTolerance;

        // We "initialize" the flag-> NOTE: Replace for a ral flag?¿
        mInitialResidualIsSet = false;
    }

    /**
     * @brief Default constructor (parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit DisplacementLagrangeMultiplierResidualFrictionalContactCriteria( Parameters ThisParameters = Parameters(R"({})"))
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >(),
          mTableIsInitialized(false)
    {
        // The default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "ensure_contact"                                 : false,
            "print_convergence_criterion"                    : false,
            "residual_relative_tolerance"                    : 1.0e-4,
            "residual_absolute_tolerance"                    : 1.0e-9,
            "contact_residual_relative_tolerance"            : 1.0e-4,
            "contact_residual_absolute_tolerance"            : 1.0e-9,
            "frictional_contact_residual_relative_tolerance" : 1.0e-4,
            "frictional_contact_residual_absolute_tolerance" : 1.0e-9
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // The displacement residual
        mDispRatioTolerance = ThisParameters["residual_relative_tolerance"].GetDouble();
        mDispAbsTolerance = ThisParameters["residual_absolute_tolerance"].GetDouble();

        // The normal contact residual
        mLMNormalRatioTolerance =  ThisParameters["contact_displacement_absolute_tolerance"].GetDouble();
        mLMNormalAbsTolerance =  ThisParameters["contact_residual_absolute_tolerance"].GetDouble();

        // The tangent contact residual
        mLMTangentRatioTolerance =  ThisParameters["frictional_contact_residual_relative_tolerance"].GetDouble();
        mLMTangentAbsTolerance =  ThisParameters["frictional_contact_residual_absolute_tolerance"].GetDouble();

        // Additional flags -> NOTE: Replace for a ral flag?¿
        mEnsureContact = ThisParameters["ensure_contact"].GetBool();
        mPrintingOutput = ThisParameters["print_convergence_criterion"].GetBool();

        // We "initialize" the flag-> NOTE: Replace for a ral flag?¿
        mInitialResidualIsSet = false;
    }

    //* Copy constructor.
    DisplacementLagrangeMultiplierResidualFrictionalContactCriteria( DisplacementLagrangeMultiplierResidualFrictionalContactCriteria const& rOther )
      :BaseType(rOther)
      ,mInitialResidualIsSet(rOther.mInitialResidualIsSet)
      ,mDispRatioTolerance(rOther.mDispRatioTolerance)
      ,mDispAbsTolerance(rOther.mDispAbsTolerance)
      ,mDispInitialResidualNorm(rOther.mDispInitialResidualNorm)
      ,mDispCurrentResidualNorm(rOther.mDispCurrentResidualNorm)
      ,mLMNormalRatioTolerance(rOther.mLMNormalRatioTolerance)
      ,mLMNormalAbsTolerance(rOther.mLMNormalAbsTolerance)
      ,mLMNormalInitialResidualNorm(rOther.mLMNormalInitialResidualNorm)
      ,mLMNormalCurrentResidualNorm(rOther.mLMNormalCurrentResidualNorm)
      ,mLMTangentRatioTolerance(rOther.mLMTangentRatioTolerance)
      ,mLMTangentAbsTolerance(rOther.mLMTangentAbsTolerance)
      ,mLMTangentInitialResidualNorm(rOther.mLMTangentInitialResidualNorm)
      ,mLMTangentCurrentResidualNorm(rOther.mLMTangentCurrentResidualNorm)
      ,mPrintingOutput(rOther.mPrintingOutput)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
    {
    }

    /// Destructor.
    ~DisplacementLagrangeMultiplierResidualFrictionalContactCriteria() override = default;

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
            TDataType disp_residual_solution_norm = 0.0, normal_lm_residual_solution_norm = 0.0, tangent_lm_residual_solution_norm = 0.0;
            IndexType disp_dof_num(0),lm_dof_num(0);

            // The nodes array
            auto& nodes_array = rModelPart.Nodes();

            // Loop over Dofs
            #pragma omp parallel for reduction(+:disp_residual_solution_norm,normal_lm_residual_solution_norm, tangent_lm_residual_solution_norm,disp_dof_num,lm_dof_num)
            for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
                auto it_dof = rDofSet.begin() + i;

                std::size_t dof_id;
                TDataType residual_dof_value;

                if (it_dof->IsFree()) {
                    // The component of the residual
                    dof_id = it_dof->EquationId();
                    residual_dof_value = rb[dof_id];

                    const auto curr_var = it_dof->GetVariable();
                    if (curr_var == VECTOR_LAGRANGE_MULTIPLIER_X) {
                        // The normal of the node (TODO: how to solve this without accesing all the time to the database?)
                        const auto& it_node = nodes_array.find(it_dof->Id());
                        const array_1d<double, 3>& normal = it_node->FastGetSolutionStepValue(NORMAL);

                        const TDataType normal_comp_residual = residual_dof_value * normal[0];
                        normal_lm_residual_solution_norm += std::pow(normal_comp_residual, 2);
                        tangent_lm_residual_solution_norm += std::pow(residual_dof_value - normal_comp_residual, 2);
                        lm_dof_num++;
                    } else if (curr_var == VECTOR_LAGRANGE_MULTIPLIER_Y) {
                        // The normal of the node (TODO: how to solve this without accesing all the time to the database?)
                        const auto& it_node = nodes_array.find(it_dof->Id());
                        const array_1d<double, 3>& normal = it_node->FastGetSolutionStepValue(NORMAL);

                        const TDataType normal_comp_residual = residual_dof_value * normal[1];
                        normal_lm_residual_solution_norm += std::pow(normal_comp_residual, 2);
                        tangent_lm_residual_solution_norm += std::pow(residual_dof_value - normal_comp_residual, 2);
                        lm_dof_num++;
                    } else if (curr_var == VECTOR_LAGRANGE_MULTIPLIER_Z) {
                        // The normal of the node (TODO: how to solve this without accesing all the time to the database?)
                        const auto& it_node = nodes_array.find(it_dof->Id());
                        const array_1d<double, 3>& normal = it_node->FastGetSolutionStepValue(NORMAL);

                        const TDataType normal_comp_residual = residual_dof_value * normal[2];
                        normal_lm_residual_solution_norm += std::pow(normal_comp_residual, 2);
                        tangent_lm_residual_solution_norm += std::pow(residual_dof_value - normal_comp_residual, 2);
                        lm_dof_num++;
                    } else {
                        disp_residual_solution_norm += residual_dof_value * residual_dof_value;
                        disp_dof_num++;
                    }
                }
            }

            mDispCurrentResidualNorm = disp_residual_solution_norm;
            mLMNormalCurrentResidualNorm = normal_lm_residual_solution_norm;
            mLMTangentCurrentResidualNorm = tangent_lm_residual_solution_norm;

            TDataType residual_disp_ratio = 1.0;
            TDataType residual_normal_lm_ratio = 1.0;
            TDataType residual_tangent_lm_ratio = 1.0;

            // We initialize the solution
            if (mInitialResidualIsSet == false) {
                mDispInitialResidualNorm = (disp_residual_solution_norm == 0.0) ? 1.0 : disp_residual_solution_norm;
                mLMNormalInitialResidualNorm = (normal_lm_residual_solution_norm == 0.0) ? 1.0 : normal_lm_residual_solution_norm;
                mLMTangentInitialResidualNorm = (tangent_lm_residual_solution_norm == 0.0) ? 1.0 : tangent_lm_residual_solution_norm;
                residual_disp_ratio = 1.0;
                residual_normal_lm_ratio = 1.0;
                residual_tangent_lm_ratio = 1.0;
                mInitialResidualIsSet = true;
            }

            // We calculate the ratio of the displacements
            residual_disp_ratio = mDispCurrentResidualNorm/mDispInitialResidualNorm;

            // We calculate the ratio of the LM
            residual_normal_lm_ratio = mLMNormalCurrentResidualNorm/mLMNormalInitialResidualNorm;
            residual_tangent_lm_ratio = mLMTangentCurrentResidualNorm/mLMTangentInitialResidualNorm;

            KRATOS_ERROR_IF(mEnsureContact && residual_normal_lm_ratio == 0.0) << "ERROR::CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?" << std::endl;

            // We calculate the absolute norms
            const TDataType residual_disp_abs = mDispCurrentResidualNorm/disp_dof_num;
            const TDataType residual_normal_lm_abs = mLMNormalCurrentResidualNorm/lm_dof_num;
            const TDataType residual_tangent_lm_abs = mLMTangentCurrentResidualNorm/lm_dof_num;

            // The process info of the model part
            ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

            // We print the results // TODO: Replace for the new log
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                if (r_process_info.Has(TABLE_UTILITY)) {
                    std::cout.precision(4);
                    TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                    auto& Table = p_table->GetTable();
                    Table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance << residual_normal_lm_ratio << mLMNormalRatioTolerance << residual_normal_lm_abs << mLMNormalAbsTolerance << residual_tangent_lm_ratio << mLMTangentRatioTolerance << residual_tangent_lm_abs << mLMTangentAbsTolerance;
                } else {
                    std::cout.precision(4);
                    if (!mPrintingOutput) {
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("RESIDUAL CONVERGENCE CHECK") << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tDISPLACEMENT: RATIO = ") << residual_disp_ratio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << residual_disp_abs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tNORMAL LAGRANGE MUL: RATIO = ") << residual_normal_lm_ratio << BOLDFONT(" EXP.RATIO = ") << mLMNormalRatioTolerance << BOLDFONT(" ABS = ") << residual_normal_lm_abs << BOLDFONT(" EXP.ABS = ") << mLMNormalAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tTANGENT LAGRANGE MUL: RATIO = ") << residual_tangent_lm_ratio << BOLDFONT(" EXP.RATIO = ") << mLMTangentRatioTolerance << BOLDFONT(" ABS = ") << residual_tangent_lm_abs << BOLDFONT(" EXP.ABS = ") << mLMTangentAbsTolerance << std::endl;
                    } else {
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "RESIDUAL CONVERGENCE CHECK" << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tDISPLACEMENT: RATIO = " << residual_disp_ratio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << residual_disp_abs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tNORMAL LAGRANGE MUL: RATIO = " << residual_normal_lm_ratio << " EXP.RATIO = " << mLMNormalRatioTolerance << " ABS = " << residual_normal_lm_abs << " EXP.ABS = " << mLMNormalAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tTANGENT LAGRANGE MUL: RATIO = " << residual_tangent_lm_ratio << " EXP.RATIO = " << mLMTangentRatioTolerance << " ABS = " << residual_tangent_lm_abs << " EXP.ABS = " << mLMTangentAbsTolerance << std::endl;
                    }
                }
            }

            // NOTE: Here we don't include the tangent counter part
            r_process_info[CONVERGENCE_RATIO] = (residual_disp_ratio > residual_normal_lm_ratio) ? residual_disp_ratio : residual_normal_lm_ratio;
            r_process_info[RESIDUAL_NORM] = (residual_normal_lm_abs > mLMNormalAbsTolerance) ? residual_normal_lm_abs : mLMNormalAbsTolerance;

            // We check if converged
            const bool disp_converged = (residual_disp_ratio <= mDispRatioTolerance || residual_disp_abs <= mDispAbsTolerance);
            const bool lm_converged = (!mEnsureContact && residual_normal_lm_ratio == 0.0) ? true : (residual_normal_lm_ratio <= mLMNormalRatioTolerance || residual_normal_lm_abs <= mLMNormalAbsTolerance) && (residual_tangent_lm_ratio <= mLMTangentRatioTolerance || residual_tangent_lm_abs <= mLMTangentAbsTolerance);

            if (disp_converged && lm_converged ) {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& Table = p_table->GetTable();
                        if (!mPrintingOutput)
                            Table << BOLDFONT(FGRN("       Achieved"));
                        else
                            Table << "Achieved";
                    } else {
                        if (!mPrintingOutput)
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tResidual") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tResidual convergence is achieved" << std::endl;
                    }
                }
                return true;
            } else {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& table = p_table->GetTable();
                        if (!mPrintingOutput)
                            table << BOLDFONT(FRED("   Not achieved"));
                        else
                            table << "Not achieved";
                    } else {
                        if (!mPrintingOutput)
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tResidual") << " convergence is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tResidual convergence is not achieved" << std::endl;
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
        BaseType::mConvergenceCriteriaIsInitialized = true;

        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info.Has(TABLE_UTILITY) && mTableIsInitialized == false) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& table = p_table->GetTable();
            table.AddColumn("DP RATIO", 10);
            table.AddColumn("EXP. RAT", 10);
            table.AddColumn("ABS", 10);
            table.AddColumn("EXP. ABS", 10);
            table.AddColumn("N.LM RATIO", 10);
            table.AddColumn("EXP. RAT", 10);
            table.AddColumn("ABS", 10);
            table.AddColumn("EXP. ABS", 10);
            table.AddColumn("T.LM RATIO", 10);
            table.AddColumn("EXP. RAT", 10);
            table.AddColumn("ABS", 10);
            table.AddColumn("EXP. ABS", 10);
            table.AddColumn("CONVERGENCE", 15);
            mTableIsInitialized = true;
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
        mInitialResidualIsSet = false;
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

    bool mInitialResidualIsSet; /// This "flag" is set in order to set that the initial residual is already computed

    bool mEnsureContact; /// This "flag" is used to check that the norm of the LM is always greater than 0 (no contact)

    bool mPrintingOutput;      /// If the colors and bold are printed
    bool mTableIsInitialized;  /// If the table is already initialized

    TDataType mDispRatioTolerance;      /// The ratio threshold for the norm of the displacement residual
    TDataType mDispAbsTolerance;        /// The absolute value threshold for the norm of the displacement residual
    TDataType mDispInitialResidualNorm; /// The reference norm of the displacement residual
    TDataType mDispCurrentResidualNorm; /// The current norm of the displacement residual

    TDataType mLMNormalRatioTolerance;      /// The ratio threshold for the norm of the normal LM residual
    TDataType mLMNormalAbsTolerance;        /// The absolute value threshold for the norm of the normal LM  residual
    TDataType mLMNormalInitialResidualNorm; /// The reference norm of the normal LM residual
    TDataType mLMNormalCurrentResidualNorm; /// The current norm of the normal LM residual

    TDataType mLMTangentRatioTolerance;      /// The ratio threshold for the norm of the tangent LM residual
    TDataType mLMTangentAbsTolerance;        /// The absolute value threshold for the norm of the tangent LM  residual
    TDataType mLMTangentInitialResidualNorm; /// The reference norm of the tangent LM residual
    TDataType mLMTangentCurrentResidualNorm; /// The current norm of the tangent LM residual

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
};

///@} // Kratos classes

///@} // Kratos namespace
}

#endif /* KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_RESIDUAL_FRICTIONAL_CONTACT_CRITERIA_H */

