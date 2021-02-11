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

#if !defined(KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_FRICTIONAL_CONTACT_CRITERIA_H)
#define KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_FRICTIONAL_CONTACT_CRITERIA_H

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/table_stream_utility.h"
#include "utilities/color_utilities.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_utilities/active_set_utilities.h"
#include "utilities/constraint_utilities.h"
#include "custom_utilities/contact_utilities.h"

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
 * @class DisplacementLagrangeMultiplierFrictionalContactCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Convergence criteria for contact problems
 * @details This class implements a convergence control based on nodal displacement and
 * lagrange multiplier values. The error is evaluated separately for each of them, and
 * relative and absolute tolerances for both must be specified.
 * @author Vicente Mataix Ferrandiz
 */
template<   class TSparseSpace,
            class TDenseSpace >
class DisplacementLagrangeMultiplierFrictionalContactCriteria
    : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of DisplacementLagrangeMultiplierFrictionalContactCriteria
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementLagrangeMultiplierFrictionalContactCriteria );

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( ENSURE_CONTACT );
    KRATOS_DEFINE_LOCAL_FLAG( PRINTING_OUTPUT );
    KRATOS_DEFINE_LOCAL_FLAG( TABLE_IS_INITIALIZED );
    KRATOS_DEFINE_LOCAL_FLAG( ROTATION_DOF_IS_CONSIDERED );
    KRATOS_DEFINE_LOCAL_FLAG( PURE_SLIP );

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

    /// The epsilon tolerance definition
    static constexpr double Tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param DispRatioTolerance Relative tolerance for displacement error
     * @param DispAbsTolerance Absolute tolerance for displacement error
     * @param RotRatioTolerance Relative tolerance for rotation error
     * @param RotAbsTolerance Absolute tolerance for rotation error
     * @param LMRatioTolerance Relative tolerance for lagrange multiplier error
     * @param LMAbsTolerance Absolute tolerance for lagrange multiplier error
     * @param NormalTangentRatio Ratio between the normal and tangent that will accepted as converged
     * @param EnsureContact To check if the contact is lost
     * @param pTable The pointer to the output r_table
     * @param PrintingOutput If the output is going to be printed in a txt file
     */
    explicit DisplacementLagrangeMultiplierFrictionalContactCriteria(
        const TDataType DispRatioTolerance,
        const TDataType DispAbsTolerance,
        const TDataType RotRatioTolerance,
        const TDataType RotAbsTolerance,
        const TDataType LMNormalRatioTolerance,
        const TDataType LMNormalAbsTolerance,
        const TDataType LMTangentStickRatioTolerance,
        const TDataType LMTangentStickAbsTolerance,
        const TDataType LMTangentSlipRatioTolerance,
        const TDataType LMTangentSlipAbsTolerance,
        const TDataType NormalTangentRatio,
        const bool EnsureContact = false,
        const bool PureSlip = false,
        const bool PrintingOutput = false
        )
        : BaseType()
    {
        // Set local flags
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::ENSURE_CONTACT, EnsureContact);
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::PRINTING_OUTPUT, PrintingOutput);
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::PURE_SLIP, PureSlip);
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED, false);

        // The displacement solution
        mDispRatioTolerance = DispRatioTolerance;
        mDispAbsTolerance = DispAbsTolerance;

        // The rotation solution
        mRotRatioTolerance = RotRatioTolerance;
        mRotAbsTolerance = RotAbsTolerance;

        // The normal contact solution
        mLMNormalRatioTolerance = LMNormalRatioTolerance;
        mLMNormalAbsTolerance = LMNormalAbsTolerance;

        // The tangent contact solution
        mLMTangentStickRatioTolerance = LMTangentStickRatioTolerance;
        mLMTangentStickAbsTolerance = LMTangentStickAbsTolerance;
        mLMTangentStickRatioTolerance = LMTangentSlipRatioTolerance;
        mLMTangentStickAbsTolerance = LMTangentSlipAbsTolerance;

        // We get the  ratio between the normal and tangent that will accepted as converged
        mNormalTangentRatio = NormalTangentRatio;
    }

    /**
     * @brief Default constructor (parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit DisplacementLagrangeMultiplierFrictionalContactCriteria( Parameters ThisParameters = Parameters(R"({})"))
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    //* Copy constructor.
    DisplacementLagrangeMultiplierFrictionalContactCriteria( DisplacementLagrangeMultiplierFrictionalContactCriteria const& rOther )
      :BaseType(rOther)
      ,mOptions(rOther.mOptions)
      ,mDispRatioTolerance(rOther.mDispRatioTolerance)
      ,mDispAbsTolerance(rOther.mDispAbsTolerance)
      ,mRotRatioTolerance(rOther.mDispRatioTolerance)
      ,mRotAbsTolerance(rOther.mDispAbsTolerance)
      ,mLMNormalRatioTolerance(rOther.mLMNormalRatioTolerance)
      ,mLMNormalAbsTolerance(rOther.mLMNormalAbsTolerance)
      ,mLMTangentStickRatioTolerance(rOther.mLMTangentStickRatioTolerance)
      ,mLMTangentStickAbsTolerance(rOther.mLMTangentStickAbsTolerance)
      ,mLMTangentSlipRatioTolerance(rOther.mLMTangentSlipRatioTolerance)
      ,mLMTangentSlipAbsTolerance(rOther.mLMTangentSlipAbsTolerance)
      ,mNormalTangentRatio(rOther.mNormalTangentRatio)
    {
    }

    /// Destructor.
    ~DisplacementLagrangeMultiplierFrictionalContactCriteria() override = default;

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

            // Getting process info
            ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

            // Initialize
            TDataType disp_solution_norm = 0.0, rot_solution_norm = 0.0, normal_lm_solution_norm = 0.0, tangent_lm_stick_solution_norm = 0.0, tangent_lm_slip_solution_norm = 0.0, disp_increase_norm = 0.0, rot_increase_norm = 0.0, normal_lm_increase_norm = 0.0, tangent_lm_stick_increase_norm = 0.0, tangent_lm_slip_increase_norm = 0.0;
            IndexType disp_dof_num(0), rot_dof_num(0), lm_dof_num(0), lm_stick_dof_num(0), lm_slip_dof_num(0);

            // First iterator
            const auto it_dof_begin = rDofSet.begin();

            // The nodes array
            auto& r_nodes_array = rModelPart.Nodes();

            // Auxiliar values
            std::size_t dof_id = 0;
            TDataType dof_value = 0.0, dof_incr = 0.0;

            // The number of active dofs
            const std::size_t number_active_dofs = rb.size();

            // Auxiliar displacement DoF check
            const std::function<bool(const VariableData&)> check_without_rot =
            [](const VariableData& rCurrVar) -> bool {return true;};
            const std::function<bool(const VariableData&)> check_with_rot =
            [](const VariableData& rCurrVar) -> bool {return ((rCurrVar == DISPLACEMENT_X) || (rCurrVar == DISPLACEMENT_Y) || (rCurrVar == DISPLACEMENT_Z));};
            const auto* p_check_disp = (mOptions.Is(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) ? &check_with_rot : &check_without_rot;

            // Loop over Dofs
            #pragma omp parallel for firstprivate(dof_id, dof_value, dof_incr) reduction(+:disp_solution_norm, rot_solution_norm, normal_lm_solution_norm, tangent_lm_slip_solution_norm, tangent_lm_stick_solution_norm, disp_increase_norm, rot_increase_norm, normal_lm_increase_norm, tangent_lm_slip_increase_norm, tangent_lm_stick_increase_norm, disp_dof_num, rot_dof_num, lm_dof_num, lm_stick_dof_num, lm_slip_dof_num)
            for (int i = 0; i < static_cast<int>(rDofSet.size()); i++) {
                auto it_dof = it_dof_begin + i;

                dof_id = it_dof->EquationId();

                // Check dof id is solved
                if (dof_id < number_active_dofs) {
                    if (mActiveDofs[dof_id] == 1) {
                        dof_value = it_dof->GetSolutionStepValue(0);
                        dof_incr = rDx[dof_id];

                        const auto& r_curr_var = it_dof->GetVariable();
                        if (r_curr_var == VECTOR_LAGRANGE_MULTIPLIER_X || r_curr_var == VECTOR_LAGRANGE_MULTIPLIER_Y || r_curr_var == VECTOR_LAGRANGE_MULTIPLIER_Z) {
                            // The normal of the node (TODO: how to solve this without accesing all the time to the database?)
                            const auto it_node = r_nodes_array.find(it_dof->Id());

                            const double mu = it_node->GetValue(FRICTION_COEFFICIENT);
                            if (mu < std::numeric_limits<double>::epsilon()) {
                                normal_lm_solution_norm += std::pow(dof_value, 2);
                                normal_lm_increase_norm += std::pow(dof_incr, 2);
                            } else {
                                const double normal = it_node->FastGetSolutionStepValue(NORMAL)[r_curr_var.GetComponentIndex()];
                                const TDataType normal_dof_value = dof_value * normal;
                                const TDataType normal_dof_incr = dof_incr * normal;

                                normal_lm_solution_norm += std::pow(normal_dof_value, 2);
                                normal_lm_increase_norm += std::pow(normal_dof_incr, 2);

                                if (it_node->Is(SLIP) || mOptions.Is(DisplacementLagrangeMultiplierFrictionalContactCriteria::PURE_SLIP)) {
                                    tangent_lm_slip_solution_norm += std::pow(dof_value - normal_dof_value, 2);
                                    tangent_lm_slip_increase_norm += std::pow(dof_incr - normal_dof_incr, 2);
                                    ++lm_slip_dof_num;
                                } else {
                                    tangent_lm_stick_solution_norm += std::pow(dof_value - normal_dof_value, 2);
                                    tangent_lm_stick_increase_norm += std::pow(dof_incr - normal_dof_incr, 2);
                                    ++lm_stick_dof_num;
                                }
                            }
                            ++lm_dof_num;
                        } else if ((*p_check_disp)(r_curr_var)) {
                            disp_solution_norm += std::pow(dof_value, 2);
                            disp_increase_norm += std::pow(dof_incr, 2);
                            ++disp_dof_num;
                        } else { // We will assume is rotation dof
                            KRATOS_DEBUG_ERROR_IF_NOT((r_curr_var == ROTATION_X) || (r_curr_var == ROTATION_Y) || (r_curr_var == ROTATION_Z)) << "Variable must be a ROTATION and it is: " << r_curr_var.Name() << std::endl;
                            rot_solution_norm += std::pow(dof_value, 2);
                            rot_increase_norm += std::pow(dof_incr, 2);
                            ++rot_dof_num;
                        }
                    }
                }
            }

            if(disp_increase_norm < Tolerance) disp_increase_norm = 1.0;
            if(rot_increase_norm < Tolerance) rot_increase_norm = 1.0;
            if(normal_lm_increase_norm < Tolerance) normal_lm_increase_norm = 1.0;
            if(tangent_lm_stick_increase_norm < Tolerance) tangent_lm_stick_increase_norm = 1.0;
            if(tangent_lm_slip_increase_norm < Tolerance) tangent_lm_slip_increase_norm = 1.0;
            if(disp_solution_norm < Tolerance) disp_solution_norm = 1.0;

            KRATOS_ERROR_IF(mOptions.Is(DisplacementLagrangeMultiplierFrictionalContactCriteria::ENSURE_CONTACT) && normal_lm_solution_norm < Tolerance) << "WARNING::CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?" << std::endl;

            const TDataType disp_ratio = std::sqrt(disp_increase_norm/disp_solution_norm);
            const TDataType rot_ratio = std::sqrt(rot_increase_norm/rot_solution_norm);
            const TDataType normal_lm_ratio = normal_lm_solution_norm > Tolerance ? std::sqrt(normal_lm_increase_norm/normal_lm_solution_norm) : 0.0;
            const TDataType tangent_lm_stick_ratio = tangent_lm_stick_solution_norm > Tolerance ? std::sqrt(tangent_lm_stick_increase_norm/tangent_lm_stick_solution_norm) : 0.0;
            const TDataType tangent_lm_slip_ratio = tangent_lm_slip_solution_norm > Tolerance ? std::sqrt(tangent_lm_slip_increase_norm/tangent_lm_slip_solution_norm) : 0.0;

            const TDataType disp_abs = std::sqrt(disp_increase_norm)/ static_cast<TDataType>(disp_dof_num);
            const TDataType rot_abs = std::sqrt(rot_increase_norm)/ static_cast<TDataType>(rot_dof_num);
            const TDataType normal_lm_abs = std::sqrt(normal_lm_increase_norm)/ static_cast<TDataType>(lm_dof_num);
            const TDataType tangent_lm_stick_abs = lm_stick_dof_num > 0 ?  std::sqrt(tangent_lm_stick_increase_norm)/ static_cast<TDataType>(lm_stick_dof_num) : 0.0;
            const TDataType tangent_lm_slip_abs = lm_slip_dof_num > 0 ? std::sqrt(tangent_lm_slip_increase_norm)/ static_cast<TDataType>(lm_slip_dof_num) : 0.0;

            const TDataType normal_tangent_stick_ratio = tangent_lm_stick_abs/normal_lm_abs;
            const TDataType normal_tangent_slip_ratio = tangent_lm_slip_abs/normal_lm_abs;

            // We print the results  // TODO: Replace for the new log
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                if (r_process_info.Has(TABLE_UTILITY)) {
                    std::cout.precision(4);
                    TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                    auto& Table = p_table->GetTable();
                    if (mOptions.Is(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                        Table << disp_ratio << mDispRatioTolerance << disp_abs << mDispAbsTolerance << rot_ratio << mRotRatioTolerance << rot_abs << mRotAbsTolerance << normal_lm_ratio << mLMNormalRatioTolerance << normal_lm_abs << mLMNormalAbsTolerance << tangent_lm_stick_ratio << mLMTangentStickRatioTolerance << tangent_lm_stick_abs << mLMTangentStickAbsTolerance << tangent_lm_slip_ratio << mLMTangentSlipRatioTolerance << tangent_lm_slip_abs << mLMTangentSlipAbsTolerance;
                    } else {
                        Table << disp_ratio << mDispRatioTolerance << disp_abs << mDispAbsTolerance << normal_lm_ratio << mLMNormalRatioTolerance << normal_lm_abs << mLMNormalAbsTolerance << tangent_lm_stick_ratio << mLMTangentStickRatioTolerance << tangent_lm_stick_abs << mLMTangentStickAbsTolerance << tangent_lm_slip_ratio << mLMTangentSlipRatioTolerance << tangent_lm_slip_abs << mLMTangentSlipAbsTolerance;
                    }
                } else {
                    std::cout.precision(4);
                    if (mOptions.IsNot(DisplacementLagrangeMultiplierFrictionalContactCriteria::PRINTING_OUTPUT)) {
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << BOLDFONT("DoF ONVERGENCE CHECK") << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << BOLDFONT("\tDISPLACEMENT: RATIO = ") << disp_ratio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << disp_abs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        if (mOptions.Is(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                            KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << BOLDFONT("\tROTATION: RATIO = ") << rot_ratio << BOLDFONT(" EXP.RATIO = ") << mRotRatioTolerance << BOLDFONT(" ABS = ") << rot_abs << BOLDFONT(" EXP.ABS = ") << mRotAbsTolerance << std::endl;
                        }
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << BOLDFONT(" NORMAL LAGRANGE MUL:\tRATIO = ") << normal_lm_ratio << BOLDFONT(" EXP.RATIO = ") << mLMNormalRatioTolerance << BOLDFONT(" ABS = ") << normal_lm_abs << BOLDFONT(" EXP.ABS = ") << mLMNormalAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << BOLDFONT(" STICK LAGRANGE MUL:\tRATIO = ") << tangent_lm_stick_ratio << BOLDFONT(" EXP.RATIO = ") << mLMTangentStickRatioTolerance << BOLDFONT(" ABS = ") << tangent_lm_stick_abs << BOLDFONT(" EXP.ABS = ") << mLMTangentStickAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << BOLDFONT(" SLIP LAGRANGE MUL:\tRATIO = ") << tangent_lm_slip_ratio << BOLDFONT(" EXP.RATIO = ") << mLMTangentSlipRatioTolerance << BOLDFONT(" ABS = ") << tangent_lm_slip_abs << BOLDFONT(" EXP.ABS = ") << mLMTangentSlipAbsTolerance << std::endl;
                    } else {
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << "DoF ONVERGENCE CHECK" << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << "\tDISPLACEMENT: RATIO = " << disp_ratio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << disp_abs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        if (mOptions.Is(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                            KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << "\tROTATION: RATIO = " << rot_ratio << " EXP.RATIO = " << mRotRatioTolerance << " ABS = " << rot_abs << " EXP.ABS = " << mRotAbsTolerance << std::endl;
                        }
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << " NORMAL LAGRANGE MUL:\tRATIO = " << normal_lm_ratio << " EXP.RATIO = " << mLMNormalRatioTolerance << " ABS = " << normal_lm_abs << " EXP.ABS = " << mLMNormalAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << " STICK LAGRANGE MUL:\tRATIO = " << tangent_lm_stick_ratio << " EXP.RATIO = " << mLMTangentStickRatioTolerance << " ABS = " << tangent_lm_stick_abs << " EXP.ABS = " << mLMTangentStickAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << " SLIP LAGRANGE MUL:\tRATIO = " << tangent_lm_slip_ratio << " EXP.RATIO = " << mLMTangentSlipRatioTolerance << " ABS = " << tangent_lm_slip_abs << " EXP.ABS = " << mLMTangentSlipAbsTolerance << std::endl;
                    }
                }
            }

            // We check if converged
            const bool disp_converged = (disp_ratio <= mDispRatioTolerance || disp_abs <= mDispAbsTolerance);
            const bool rot_converged = (mOptions.Is(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) ? (rot_ratio <= mRotRatioTolerance || rot_abs <= mRotAbsTolerance) : true;
            const bool lm_converged = (mOptions.IsNot(DisplacementLagrangeMultiplierFrictionalContactCriteria::ENSURE_CONTACT) && normal_lm_solution_norm < Tolerance) ? true : (normal_lm_ratio <= mLMNormalRatioTolerance || normal_lm_abs <= mLMNormalAbsTolerance) && (tangent_lm_stick_ratio <= mLMTangentStickRatioTolerance || tangent_lm_stick_abs <= mLMTangentStickAbsTolerance || normal_tangent_stick_ratio <= mNormalTangentRatio) && (tangent_lm_slip_ratio <= mLMTangentSlipRatioTolerance || tangent_lm_slip_abs <= mLMTangentSlipAbsTolerance || normal_tangent_slip_ratio <= mNormalTangentRatio);

            if (disp_converged && rot_converged && lm_converged) {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& r_table = p_table->GetTable();
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierFrictionalContactCriteria::PRINTING_OUTPUT))
                            r_table << BOLDFONT(FGRN("       Achieved"));
                        else
                            r_table << "Achieved";
                    } else {
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierFrictionalContactCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << BOLDFONT("\tDoF") << " convergence is " << BOLDFONT(FGRN("achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << "\tDoF convergence is achieved" << std::endl;
                    }
                }
                return true;
            } else {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& r_table = p_table->GetTable();
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierFrictionalContactCriteria::PRINTING_OUTPUT))
                            r_table << BOLDFONT(FRED("   Not achieved"));
                        else
                            r_table << "Not achieved";
                    } else {
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierFrictionalContactCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << BOLDFONT("\tDoF") << " convergence is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementLagrangeMultiplierFrictionalContactCriteria") << "\tDoF convergence is not achieved" << std::endl;
                    }
                }
                return false;
            }
        }
        else // In this case all the displacements are imposed!
            return true;
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart Reference to the ModelPart containing the contact problem. (unused)
     */
    void Initialize( ModelPart& rModelPart ) override
    {
        // Initialize
        BaseType::mConvergenceCriteriaIsInitialized = true;

        // Check rotation dof
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED, ContactUtilities::CheckModelPartHasRotationDoF(rModelPart));

        // Initialize header
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info.Has(TABLE_UTILITY) && mOptions.IsNot(DisplacementLagrangeMultiplierFrictionalContactCriteria::TABLE_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& r_table = p_table->GetTable();
            r_table.AddColumn("DP RATIO", 10);
            r_table.AddColumn("EXP. RAT", 10);
            r_table.AddColumn("ABS", 10);
            r_table.AddColumn("EXP. ABS", 10);
            if (mOptions.Is(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                r_table.AddColumn("RT RATIO", 10);
                r_table.AddColumn("EXP. RAT", 10);
                r_table.AddColumn("ABS", 10);
                r_table.AddColumn("EXP. ABS", 10);
            }
            r_table.AddColumn("N.LM RATIO", 10);
            r_table.AddColumn("EXP. RAT", 10);
            r_table.AddColumn("ABS", 10);
            r_table.AddColumn("EXP. ABS", 10);
            if (mOptions.IsNot(DisplacementLagrangeMultiplierFrictionalContactCriteria::PURE_SLIP)) {
                r_table.AddColumn("STI. RATIO", 10);
                r_table.AddColumn("EXP. RAT", 10);
                r_table.AddColumn("ABS", 10);
                r_table.AddColumn("EXP. ABS", 10);
            }
            r_table.AddColumn("SLIP RATIO", 10);
            r_table.AddColumn("EXP. RAT", 10);
            r_table.AddColumn("ABS", 10);
            r_table.AddColumn("EXP. ABS", 10);
            r_table.AddColumn("CONVERGENCE", 15);
            mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::TABLE_IS_INITIALIZED, true);
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
        // Filling mActiveDofs when MPC exist
        ConstraintUtilities::ComputeActiveDofs(rModelPart, mActiveDofs, rDofSet);
    }

    /**
     * @brief This function finalizes the non-linear iteration
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // Calling base criteria
        BaseType::FinalizeNonLinearIteration(rModelPart, rDofSet, rA, rDx, rb);

        // The current process info
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(ACTIVE_SET_COMPUTED, false);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                                     : "displacement_lagrangemultiplier_frictional_contact_criteria",
            "ensure_contact"                                           : false,
            "pure_slip"                                                : false,
            "print_convergence_criterion"                              : false,
            "displacement_relative_tolerance"                          : 1.0e-4,
            "displacement_absolute_tolerance"                          : 1.0e-9,
            "rotation_relative_tolerance"                              : 1.0e-4,
            "rotation_absolute_tolerance"                              : 1.0e-9,
            "contact_displacement_relative_tolerance"                  : 1.0e-4,
            "contact_displacement_absolute_tolerance"                  : 1.0e-9,
            "frictional_stick_contact_displacement_relative_tolerance" : 1.0e-4,
            "frictional_stick_contact_displacement_absolute_tolerance" : 1.0e-9,
            "frictional_slip_contact_displacement_relative_tolerance"  : 1.0e-4,
            "frictional_slip_contact_displacement_absolute_tolerance"  : 1.0e-9,
            "ratio_normal_tangent_threshold"                           : 1.0e-4
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
        return "displacement_lagrangemultiplier_frictional_contact_criteria";
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

        // The normal contact solution
        mLMNormalRatioTolerance =  ThisParameters["contact_displacement_relative_tolerance"].GetDouble();
        mLMNormalAbsTolerance =  ThisParameters["contact_displacement_absolute_tolerance"].GetDouble();

        // The tangent contact solution
        mLMTangentStickRatioTolerance =  ThisParameters["frictional_stick_contact_displacement_relative_tolerance"].GetDouble();
        mLMTangentStickAbsTolerance =  ThisParameters["frictional_stick_contact_displacement_absolute_tolerance"].GetDouble();
        mLMTangentSlipRatioTolerance =  ThisParameters["frictional_slip_contact_displacement_relative_tolerance"].GetDouble();
        mLMTangentSlipAbsTolerance =  ThisParameters["frictional_slip_contact_displacement_absolute_tolerance"].GetDouble();

        // We get the  ratio between the normal and tangent that will accepted as converged
        mNormalTangentRatio = ThisParameters["ratio_normal_tangent_threshold"].GetDouble();

        // Set local flags
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::ENSURE_CONTACT, ThisParameters["ensure_contact"].GetBool());
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::PRINTING_OUTPUT, ThisParameters["print_convergence_criterion"].GetBool());
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED, false);
        mOptions.Set(DisplacementLagrangeMultiplierFrictionalContactCriteria::PURE_SLIP, ThisParameters["pure_slip"].GetBool());
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

    TDataType mDispRatioTolerance;      /// The ratio threshold for the norm of the displacement
    TDataType mDispAbsTolerance;        /// The absolute value threshold for the norm of the displacement

    TDataType mRotRatioTolerance;      /// The ratio threshold for the norm of the rotation
    TDataType mRotAbsTolerance;        /// The absolute value threshold for the norm of the rotation

    TDataType mLMNormalRatioTolerance;  /// The ratio threshold for the norm of the LM (normal)
    TDataType mLMNormalAbsTolerance;    /// The absolute value threshold for the norm of the LM (normal)

    TDataType mLMTangentStickRatioTolerance; /// The ratio threshold for the norm of the LM (tangent-stick)
    TDataType mLMTangentStickAbsTolerance;   /// The absolute value threshold for the norm of the LM (tangent-stick)
    TDataType mLMTangentSlipRatioTolerance;  /// The ratio threshold for the norm of the LM (tangent-slip)
    TDataType mLMTangentSlipAbsTolerance;    /// The absolute value threshold for the norm of the LM (tangent-slip)

    TDataType mNormalTangentRatio;      /// The ratio to accept a non converged tangent component in case

    std::vector<int> mActiveDofs;       /// This vector contains the dofs that are active

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
};  // Kratos DisplacementLagrangeMultiplierFrictionalContactCriteria

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierFrictionalContactCriteria<TSparseSpace, TDenseSpace>::ENSURE_CONTACT(Kratos::Flags::Create(0));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierFrictionalContactCriteria<TSparseSpace, TDenseSpace>::PRINTING_OUTPUT(Kratos::Flags::Create(1));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierFrictionalContactCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(2));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierFrictionalContactCriteria<TSparseSpace, TDenseSpace>::ROTATION_DOF_IS_CONSIDERED(Kratos::Flags::Create(3));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierFrictionalContactCriteria<TSparseSpace, TDenseSpace>::PURE_SLIP(Kratos::Flags::Create(4));
}

#endif	/* KRATOS_DISPLACEMENT_LAGRANGE_MULTIPLIER_FRICTIONAL_CONTACT_CRITERIA_H */
