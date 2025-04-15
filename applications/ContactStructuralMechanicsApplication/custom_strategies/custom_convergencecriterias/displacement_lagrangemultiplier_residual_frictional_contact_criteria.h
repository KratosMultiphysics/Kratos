// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "utilities/table_stream_utility.h"
#include "custom_strategies/custom_convergencecriterias/base_mortar_criteria.h"
#include "utilities/color_utilities.h"
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

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( ENSURE_CONTACT );
    KRATOS_DEFINE_LOCAL_FLAG( PRINTING_OUTPUT );
    KRATOS_DEFINE_LOCAL_FLAG( TABLE_IS_INITIALIZED );
    KRATOS_DEFINE_LOCAL_FLAG( ROTATION_DOF_IS_CONSIDERED );
    KRATOS_DEFINE_LOCAL_FLAG( PURE_SLIP );
    KRATOS_DEFINE_LOCAL_FLAG( INITIAL_RESIDUAL_IS_SET );
    KRATOS_DEFINE_LOCAL_FLAG( INITIAL_NORMAL_RESIDUAL_IS_SET );
    KRATOS_DEFINE_LOCAL_FLAG( INITIAL_STICK_RESIDUAL_IS_SET );
    KRATOS_DEFINE_LOCAL_FLAG( INITIAL_SLIP_RESIDUAL_IS_SET );

    /// The base class definition
    using BaseType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;

    /// The definition of the current class
    using ClassType = DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>;

    /// The dofs array type
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// The sparse matrix type
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    /// The dense vector type
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    /// The sparse space used
    using SparseSpaceType = TSparseSpace;

    /// The table stream definition TODO: Replace by logger
    using TablePrinterPointerType = TableStreamUtility::Pointer;

    /// The index type definition
    using IndexType = std::size_t;

    /// Zero tolerance definition
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

   /**
     * @brief Default constructor.
     */
    explicit DisplacementLagrangeMultiplierResidualFrictionalContactCriteria()
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit DisplacementLagrangeMultiplierResidualFrictionalContactCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor
     * @param DispRatioTolerance Relative tolerance for displacement residual error
     * @param DispAbsTolerance Absolute tolerance for displacement residual error
     * @param RotRatioTolerance Relative tolerance for rotation residual error
     * @param RotAbsTolerance Absolute tolerance for rotation residual error
     * @param LMRatioTolerance Relative tolerance for lagrange multiplier residual  error
     * @param LMAbsTolerance Absolute tolerance for lagrange multiplier residual error
     * @param EnsureContact To check if the contact is lost
     * @param NormalTangentRatio Ratio between the normal and tangent that will accepted as converged
     * @param pTable The pointer to the output r_table
     * @param PrintingOutput If the output is going to be printed in a txt file
     */
    explicit DisplacementLagrangeMultiplierResidualFrictionalContactCriteria(
        const double DispRatioTolerance,
        const double DispAbsTolerance,
        const double RotRatioTolerance,
        const double RotAbsTolerance,
        const double LMNormalRatioTolerance,
        const double LMNormalAbsTolerance,
        const double LMTangentStickRatioTolerance,
        const double LMTangentStickAbsTolerance,
        const double LMTangentSlipRatioTolerance,
        const double LMTangentSlipAbsTolerance,
        const double NormalTangentRatio,
        const bool EnsureContact = false,
        const bool PureSlip = false,
        const bool PrintingOutput = false
        ) : BaseType()
    {
        // Set local flags
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ENSURE_CONTACT, EnsureContact);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PRINTING_OUTPUT, PrintingOutput);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PURE_SLIP, PureSlip);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_NORMAL_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_STICK_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_SLIP_RESIDUAL_IS_SET, false);

        // The displacement residual
        mDispRatioTolerance = DispRatioTolerance;
        mDispAbsTolerance = DispAbsTolerance;

        // The rotation residual
        mRotRatioTolerance = RotRatioTolerance;
        mRotAbsTolerance = RotAbsTolerance;

        // The normal contact residual
        mLMNormalRatioTolerance = LMNormalRatioTolerance;
        mLMNormalAbsTolerance = LMNormalAbsTolerance;

        // The tangent contact residual
        mLMTangentStickRatioTolerance = LMTangentStickRatioTolerance;
        mLMTangentStickAbsTolerance = LMTangentStickAbsTolerance;
        mLMTangentSlipRatioTolerance = LMTangentSlipRatioTolerance;
        mLMTangentSlipAbsTolerance = LMTangentSlipAbsTolerance;

        // We get the  ratio between the normal and tangent that will accepted as converged
        mNormalTangentRatio = NormalTangentRatio;
    }

    // Copy constructor.
    DisplacementLagrangeMultiplierResidualFrictionalContactCriteria( DisplacementLagrangeMultiplierResidualFrictionalContactCriteria const& rOther )
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
      ,mLMNormalRatioTolerance(rOther.mLMNormalRatioTolerance)
      ,mLMNormalAbsTolerance(rOther.mLMNormalAbsTolerance)
      ,mLMNormalInitialResidualNorm(rOther.mLMNormalInitialResidualNorm)
      ,mLMNormalCurrentResidualNorm(rOther.mLMNormalCurrentResidualNorm)
      ,mLMTangentStickRatioTolerance(rOther.mLMTangentStickRatioTolerance)
      ,mLMTangentStickAbsTolerance(rOther.mLMTangentStickAbsTolerance)
      ,mLMTangentSlipRatioTolerance(rOther.mLMTangentSlipRatioTolerance)
      ,mLMTangentSlipAbsTolerance(rOther.mLMTangentSlipAbsTolerance)
      ,mLMTangentStickInitialResidualNorm(rOther.mLMTangentStickInitialResidualNorm)
      ,mLMTangentStickCurrentResidualNorm(rOther.mLMTangentStickCurrentResidualNorm)
      ,mStickCounter(rOther.mStickCounter)
      ,mSlipCounter(rOther.mSlipCounter)
      ,mNormalTangentRatio(rOther.mNormalTangentRatio)
    {
    }

    /// Destructor.
    ~DisplacementLagrangeMultiplierResidualFrictionalContactCriteria() override = default;

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
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
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
        if (SparseSpaceType::Size(rb) != 0) { //if we are solving for something

            // Getting process info
            ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

            // Initialize
            double disp_residual_solution_norm = 0.0, rot_residual_solution_norm = 0.0,normal_lm_residual_solution_norm = 0.0, tangent_lm_stick_residual_solution_norm = 0.0, tangent_lm_slip_residual_solution_norm = 0.0;
            IndexType disp_dof_num(0), rot_dof_num(0), lm_dof_num(0), lm_stick_dof_num(0), lm_slip_dof_num(0);

            // The nodes array
            auto& r_nodes_array = rModelPart.Nodes();

            // Auxiliary values
            struct AuxValues {
                std::size_t dof_id = 0;
                double residual_dof_value = 0.0;
            };
            const bool pure_slip = mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PURE_SLIP);

            // The number of active dofs
            const std::size_t number_active_dofs = rb.size();

            // Auxiliary displacement DoF check
            const std::function<bool(const VariableData&)> check_without_rot =
            [](const VariableData& rCurrVar) -> bool {return true;};
            const std::function<bool(const VariableData&)> check_with_rot =
            [](const VariableData& rCurrVar) -> bool {return ((rCurrVar == DISPLACEMENT_X) || (rCurrVar == DISPLACEMENT_Y) || (rCurrVar == DISPLACEMENT_Z));};
            const auto* p_check_disp = (mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) ? &check_with_rot : &check_without_rot;

            // Loop over Dofs
            using TenReduction = CombinedReduction<SumReduction<double>, SumReduction<double>, SumReduction<double>, SumReduction<double>, SumReduction<double>, SumReduction<IndexType>, SumReduction<IndexType>, SumReduction<IndexType>, SumReduction<IndexType>, SumReduction<IndexType>>;
            std::tie(disp_residual_solution_norm, rot_residual_solution_norm, normal_lm_residual_solution_norm, tangent_lm_slip_residual_solution_norm, tangent_lm_stick_residual_solution_norm, disp_dof_num, rot_dof_num, lm_dof_num, lm_slip_dof_num, lm_stick_dof_num) = block_for_each<TenReduction>(rDofSet, AuxValues(), [this,&number_active_dofs,p_check_disp,&pure_slip,&r_nodes_array,&rb](Dof<double>& rDof, AuxValues& aux_values) {
                aux_values.dof_id = rDof.EquationId();

                // Check dof id is solved
                if (aux_values.dof_id < number_active_dofs) {
                    if (mActiveDofs[aux_values.dof_id] == 1) {
                        // The component of the residual
                        aux_values.residual_dof_value = rb[aux_values.dof_id];

                        const auto& r_curr_var = rDof.GetVariable();
                        if (r_curr_var == VECTOR_LAGRANGE_MULTIPLIER_X || r_curr_var == VECTOR_LAGRANGE_MULTIPLIER_Y || r_curr_var == VECTOR_LAGRANGE_MULTIPLIER_Z) {
                            // The normal of the node (TODO: how to solve this without accesing all the time to the database?)
                            const auto it_node = r_nodes_array.find(rDof.Id());
                            const double mu = it_node->GetValue(FRICTION_COEFFICIENT);

                            if (mu < ZeroTolerance) {
                                return std::make_tuple(0.0,0.0,std::pow(aux_values.residual_dof_value, 2),0.0,0.0,0,0,1,0,0);
                            } else {
                                const double normal = it_node->FastGetSolutionStepValue(NORMAL)[r_curr_var.GetComponentIndex()];
                                const double normal_comp_residual = aux_values.residual_dof_value * normal;
                                if (it_node->Is(SLIP) || pure_slip) {
                                    return std::make_tuple(0.0,0.0,std::pow(normal_comp_residual, 2),std::pow(aux_values.residual_dof_value - normal_comp_residual, 2),0.0,0,0,1,1,0);
                                } else {
                                    return std::make_tuple(0.0,0.0,std::pow(normal_comp_residual, 2),0.0,std::pow(aux_values.residual_dof_value - normal_comp_residual, 2),0,0,1,0,1);
                                }
                            }
                            return std::make_tuple(0.0,0.0,0.0,0.0,0.0,0,0,0,0,0);
                        } else if ((*p_check_disp)(r_curr_var)) {
                            return std::make_tuple(std::pow(aux_values.residual_dof_value, 2),0.0,0.0,0.0,0.0,1,0,0,0,0);
                        } else { // We will assume is rotation dof
                            KRATOS_DEBUG_ERROR_IF_NOT((r_curr_var == ROTATION_X) || (r_curr_var == ROTATION_Y) || (r_curr_var == ROTATION_Z)) << "Variable must be a ROTATION and it is: " << r_curr_var.Name() << std::endl;
                            return std::make_tuple(0.0,std::pow(aux_values.residual_dof_value, 2),0.0,0.0,0.0,0,1,0,0,0);
                        }
                    }
                }
                return std::make_tuple(0.0,0.0,0.0,0.0,0.0,0,0,0,0,0);
            });

            // Auxiliary dofs counters
            if (mStickCounter > 0) {
                if (lm_stick_dof_num == 0) {
                    mStickCounter = 0;
                    mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_STICK_RESIDUAL_IS_SET, false);
                }
            } else {
                if (lm_stick_dof_num > 0) {
                    mStickCounter = lm_stick_dof_num;
                    mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_STICK_RESIDUAL_IS_SET, false);
                }
            }
            if (mSlipCounter > 0) {
                if (lm_slip_dof_num == 0) {
                    mSlipCounter = 0;
                    mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_SLIP_RESIDUAL_IS_SET, false);
                }
            } else {
                if (lm_slip_dof_num > 0) {
                    mSlipCounter = lm_slip_dof_num;
                    mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_SLIP_RESIDUAL_IS_SET, false);
                }
            }

            mDispCurrentResidualNorm = disp_residual_solution_norm;
            mRotCurrentResidualNorm = rot_residual_solution_norm;
            mLMNormalCurrentResidualNorm = normal_lm_residual_solution_norm;
            mLMTangentStickCurrentResidualNorm = tangent_lm_stick_residual_solution_norm;
            mLMTangentSlipCurrentResidualNorm = tangent_lm_slip_residual_solution_norm;

            double residual_disp_ratio = 1.0;
            double residual_rot_ratio = 1.0;
            double residual_normal_lm_ratio = 1.0;
            double residual_tangent_lm_stick_ratio = 1.0;
            double residual_tangent_lm_slip_ratio = 1.0;

            // We initialize the solution
            if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_RESIDUAL_IS_SET)) {
                mDispInitialResidualNorm = (disp_residual_solution_norm < ZeroTolerance) ? 1.0 : disp_residual_solution_norm;
                residual_disp_ratio = 1.0;
                if (mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                    mRotInitialResidualNorm = (rot_residual_solution_norm < ZeroTolerance) ? 1.0 : rot_residual_solution_norm;
                    residual_rot_ratio = 1.0;
                }
                mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_RESIDUAL_IS_SET, true);
            }

            // We calculate the ratio of the displacements
            residual_disp_ratio = mDispCurrentResidualNorm/mDispInitialResidualNorm;
            residual_rot_ratio = mRotCurrentResidualNorm/mRotInitialResidualNorm;

            // We initialize the solution
            if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_NORMAL_RESIDUAL_IS_SET)) {
                mLMNormalInitialResidualNorm = (normal_lm_residual_solution_norm < ZeroTolerance) ? 1.0 : normal_lm_residual_solution_norm;
                residual_normal_lm_ratio = 1.0;
                mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_NORMAL_RESIDUAL_IS_SET, true);
            }

            // We calculate the ratio of the normal LM
            residual_normal_lm_ratio = mLMNormalCurrentResidualNorm/mLMNormalInitialResidualNorm;

            // We initialize the solution
            if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_STICK_RESIDUAL_IS_SET) && lm_stick_dof_num > 0) {
                mLMTangentStickInitialResidualNorm = (tangent_lm_stick_residual_solution_norm < ZeroTolerance) ? 1.0 : tangent_lm_stick_residual_solution_norm;
                residual_tangent_lm_stick_ratio = 1.0;
                mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_STICK_RESIDUAL_IS_SET, true);
            }
            if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_SLIP_RESIDUAL_IS_SET) && lm_slip_dof_num > 0) {
                mLMTangentSlipInitialResidualNorm = (tangent_lm_slip_residual_solution_norm < ZeroTolerance) ? 1.0 : tangent_lm_slip_residual_solution_norm;
                residual_tangent_lm_slip_ratio = 1.0;
                mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_SLIP_RESIDUAL_IS_SET, true);
            }

            // We calculate the ratio of the tangent LM
            if (lm_stick_dof_num > 0) {
                residual_tangent_lm_stick_ratio = mLMTangentStickCurrentResidualNorm/mLMTangentStickInitialResidualNorm;
            } else {
                residual_tangent_lm_stick_ratio = 0.0;
            }
            if (lm_slip_dof_num > 0) {
                residual_tangent_lm_slip_ratio = mLMTangentSlipCurrentResidualNorm/mLMTangentSlipInitialResidualNorm;
            } else {
                residual_tangent_lm_slip_ratio = 0.0;
            }

            KRATOS_ERROR_IF(mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ENSURE_CONTACT) && residual_normal_lm_ratio < ZeroTolerance) << "ERROR::CONTACT LOST::ARE YOU SURE YOU ARE SUPPOSED TO HAVE CONTACT?" << std::endl;

            // We calculate the absolute norms
            const double residual_disp_abs = mDispCurrentResidualNorm/static_cast<double>(disp_dof_num);
            const double residual_rot_abs = mRotCurrentResidualNorm/static_cast<double>(rot_dof_num);
            const double residual_normal_lm_abs = mLMNormalCurrentResidualNorm/static_cast<double>(lm_dof_num);
            const double residual_tangent_lm_stick_abs = lm_stick_dof_num > 0 ? mLMTangentStickCurrentResidualNorm/static_cast<double>(lm_dof_num) : 0.0;
//             const double residual_tangent_lm_stick_abs = lm_stick_dof_num > 0 ? mLMTangentStickCurrentResidualNorm/static_cast<double>(lm_stick_dof_num) : 0.0;
            const double residual_tangent_lm_slip_abs = lm_slip_dof_num > 0 ? mLMTangentSlipCurrentResidualNorm/static_cast<double>(lm_dof_num) : 0.0;
//             const double residual_tangent_lm_slip_abs = lm_slip_dof_num > 0 ? mLMTangentSlipCurrentResidualNorm/static_cast<double>(lm_slip_dof_num) : 0.0;
            const double normal_tangent_stick_ratio = residual_tangent_lm_stick_abs/residual_normal_lm_abs;
            const double normal_tangent_slip_ratio = residual_tangent_lm_slip_abs/residual_normal_lm_abs;

            // We print the results // TODO: Replace for the new log
            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                if (r_process_info.Has(TABLE_UTILITY)) {
                    std::cout.precision(4);
                    TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                    auto& r_table = p_table->GetTable();
                    if (mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PURE_SLIP)) {
                            r_table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance << residual_rot_ratio << mRotRatioTolerance << residual_rot_abs << mRotAbsTolerance << residual_normal_lm_ratio << mLMNormalRatioTolerance << residual_normal_lm_abs << mLMNormalAbsTolerance << residual_tangent_lm_stick_ratio << mLMTangentStickRatioTolerance << residual_tangent_lm_stick_abs << mLMTangentStickAbsTolerance << residual_tangent_lm_slip_ratio << mLMTangentSlipRatioTolerance << residual_tangent_lm_slip_abs << mLMTangentSlipAbsTolerance;
                        } else {
                            r_table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance << residual_rot_ratio << mRotRatioTolerance << residual_rot_abs << mRotAbsTolerance << residual_normal_lm_ratio << mLMNormalRatioTolerance << residual_normal_lm_abs << mLMNormalAbsTolerance << residual_tangent_lm_slip_ratio << mLMTangentSlipRatioTolerance << residual_tangent_lm_slip_abs << mLMTangentSlipAbsTolerance;
                        }
                    } else {
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PURE_SLIP)) {
                            r_table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance << residual_normal_lm_ratio << mLMNormalRatioTolerance << residual_normal_lm_abs << mLMNormalAbsTolerance << residual_tangent_lm_stick_ratio << mLMTangentStickRatioTolerance << residual_tangent_lm_stick_abs << mLMTangentStickAbsTolerance << residual_tangent_lm_slip_ratio << mLMTangentSlipRatioTolerance << residual_tangent_lm_slip_abs << mLMTangentSlipAbsTolerance;
                        } else {
                            r_table << residual_disp_ratio << mDispRatioTolerance << residual_disp_abs << mDispAbsTolerance << residual_normal_lm_ratio << mLMNormalRatioTolerance << residual_normal_lm_abs << mLMNormalAbsTolerance << residual_tangent_lm_slip_ratio << mLMTangentSlipRatioTolerance << residual_tangent_lm_slip_abs << mLMTangentSlipAbsTolerance;
                        }
                    }
                } else {
                    std::cout.precision(4);
                    if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PRINTING_OUTPUT)) {
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("RESIDUAL CONVERGENCE CHECK") << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tDISPLACEMENT: RATIO = ") << residual_disp_ratio << BOLDFONT(" EXP.RATIO = ") << mDispRatioTolerance << BOLDFONT(" ABS = ") << residual_disp_abs << BOLDFONT(" EXP.ABS = ") << mDispAbsTolerance << std::endl;
                        if (mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tROTATION: RATIO = ") << residual_rot_ratio << BOLDFONT(" EXP.RATIO = ") << mRotRatioTolerance << BOLDFONT(" ABS = ") << residual_rot_abs << BOLDFONT(" EXP.ABS = ") << mRotAbsTolerance << std::endl;
                        }
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tNORMAL LAGRANGE MUL: RATIO = ") << residual_normal_lm_ratio << BOLDFONT(" EXP.RATIO = ") << mLMNormalRatioTolerance << BOLDFONT(" ABS = ") << residual_normal_lm_abs << BOLDFONT(" EXP.ABS = ") << mLMNormalAbsTolerance << std::endl;
                        KRATOS_INFO_IF("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria", mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PURE_SLIP)) << BOLDFONT("\tSTICK LAGRANGE MUL: RATIO = ") << residual_tangent_lm_stick_ratio << BOLDFONT(" EXP.RATIO = ") << mLMTangentStickRatioTolerance << BOLDFONT(" ABS = ") << residual_tangent_lm_stick_abs << BOLDFONT(" EXP.ABS = ") << mLMTangentStickAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tSLIP LAGRANGE MUL: RATIO = ") << residual_tangent_lm_slip_ratio << BOLDFONT(" EXP.RATIO = ") << mLMTangentSlipRatioTolerance << BOLDFONT(" ABS = ") << residual_tangent_lm_slip_abs << BOLDFONT(" EXP.ABS = ") << mLMTangentSlipAbsTolerance << std::endl;
                    } else {
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "RESIDUAL CONVERGENCE CHECK" << "\tSTEP: " << r_process_info[STEP] << "\tNL ITERATION: " << r_process_info[NL_ITERATION_NUMBER] << std::endl << std::scientific;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tDISPLACEMENT: RATIO = " << residual_disp_ratio << " EXP.RATIO = " << mDispRatioTolerance << " ABS = " << residual_disp_abs << " EXP.ABS = " << mDispAbsTolerance << std::endl;
                        if (mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tROTATION: RATIO = " << residual_rot_ratio << " EXP.RATIO = " << mRotRatioTolerance << " ABS = " << residual_rot_abs << " EXP.ABS = " << mRotAbsTolerance << std::endl;
                        }
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tNORMAL LAGRANGE MUL: RATIO = " << residual_normal_lm_ratio << " EXP.RATIO = " << mLMNormalRatioTolerance << " ABS = " << residual_normal_lm_abs << " EXP.ABS = " << mLMNormalAbsTolerance << std::endl;
                        KRATOS_INFO_IF("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria", mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PURE_SLIP)) << "\tSTICK LAGRANGE MUL: RATIO = " << residual_tangent_lm_stick_ratio << " EXP.RATIO = " << mLMTangentStickRatioTolerance << " ABS = " << residual_tangent_lm_stick_abs << " EXP.ABS = " << mLMTangentStickAbsTolerance << std::endl;
                        KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tSLIP LAGRANGE MUL: RATIO = " << residual_tangent_lm_slip_ratio << " EXP.RATIO = " << mLMTangentSlipRatioTolerance << " ABS = " << residual_tangent_lm_slip_abs << " EXP.ABS = " << mLMTangentSlipAbsTolerance << std::endl;
                    }
                }
            }

            // NOTE: Here we don't include the tangent counter part
            r_process_info[CONVERGENCE_RATIO] = (residual_disp_ratio > residual_normal_lm_ratio) ? residual_disp_ratio : residual_normal_lm_ratio;
            r_process_info[RESIDUAL_NORM] = (residual_normal_lm_abs > mLMNormalAbsTolerance) ? residual_normal_lm_abs : mLMNormalAbsTolerance;

            // We check if converged
            const bool disp_converged = (residual_disp_ratio <= mDispRatioTolerance || residual_disp_abs <= mDispAbsTolerance);
            const bool rot_converged = (mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) ? (residual_rot_ratio <= mRotRatioTolerance || residual_rot_abs <= mRotAbsTolerance) : true;
            const bool lm_converged = (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ENSURE_CONTACT) && residual_normal_lm_ratio == 0.0) ? true : (residual_normal_lm_ratio <= mLMNormalRatioTolerance || residual_normal_lm_abs <= mLMNormalAbsTolerance) && (residual_tangent_lm_stick_ratio <= mLMTangentStickRatioTolerance || residual_tangent_lm_stick_abs <= mLMTangentStickAbsTolerance || normal_tangent_stick_ratio <= mNormalTangentRatio) && (residual_tangent_lm_slip_ratio <= mLMTangentSlipRatioTolerance || residual_tangent_lm_slip_abs <= mLMTangentSlipAbsTolerance || normal_tangent_slip_ratio <= mNormalTangentRatio);

            if (disp_converged && rot_converged && lm_converged ) {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0) {
                    if (r_process_info.Has(TABLE_UTILITY)) {
                        TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
                        auto& r_table = p_table->GetTable();
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PRINTING_OUTPUT))
                            r_table << BOLDFONT(FGRN("       Achieved"));
                        else
                            r_table << "Achieved";
                    } else {
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PRINTING_OUTPUT))
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
                        auto& r_table = p_table->GetTable();
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PRINTING_OUTPUT))
                            r_table << BOLDFONT(FRED("   Not achieved"));
                        else
                            r_table << "Not achieved";
                    } else {
                        if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PRINTING_OUTPUT))
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << BOLDFONT("\tResidual") << " convergence is " << BOLDFONT(FRED(" not achieved")) << std::endl;
                        else
                            KRATOS_INFO("DisplacementLagrangeMultiplierResidualFrictionalContactCriteria") << "\tResidual convergence is not achieved" << std::endl;
                    }
                }
                return false;
            }
        } else { // In this case all the displacements are imposed!
            return true;
        }
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
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED, ContactUtilities::CheckModelPartHasRotationDoF(rModelPart));

        // Initialize header
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        if (r_process_info.Has(TABLE_UTILITY) && mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::TABLE_IS_INITIALIZED)) {
            TablePrinterPointerType p_table = r_process_info[TABLE_UTILITY];
            auto& r_table = p_table->GetTable();
            r_table.AddColumn("DP RATIO", 10);
            r_table.AddColumn("EXP. RAT", 10);
            r_table.AddColumn("ABS", 10);
            r_table.AddColumn("EXP. ABS", 10);
            if (mOptions.Is(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED)) {
                r_table.AddColumn("RT RATIO", 10);
                r_table.AddColumn("EXP. RAT", 10);
                r_table.AddColumn("ABS", 10);
                r_table.AddColumn("EXP. ABS", 10);
            }
            r_table.AddColumn("N.LM RATIO", 10);
            r_table.AddColumn("EXP. RAT", 10);
            r_table.AddColumn("ABS", 10);
            r_table.AddColumn("EXP. ABS", 10);
            if (mOptions.IsNot(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PURE_SLIP)) {
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
            mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::TABLE_IS_INITIALIZED, true);
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
        // Initialize flags
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_NORMAL_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_STICK_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_SLIP_RESIDUAL_IS_SET, false);

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
            "name"                                                 : "displacement_lagrangemultiplier_ressidual_frictional_contact_criteria",
            "ensure_contact"                                       : false,
            "pure_slip"                                            : false,
            "print_convergence_criterion"                          : false,
            "residual_relative_tolerance"                          : 1.0e-4,
            "residual_absolute_tolerance"                          : 1.0e-9,
            "rotation_residual_relative_tolerance"                 : 1.0e-4,
            "rotation_residual_absolute_tolerance"                 : 1.0e-9,
            "contact_residual_relative_tolerance"                  : 1.0e-4,
            "contact_residual_absolute_tolerance"                  : 1.0e-9,
            "frictional_stick_contact_residual_relative_tolerance" : 1.0e-4,
            "frictional_stick_contact_residual_absolute_tolerance" : 1.0e-9,
            "frictional_slip_contact_residual_relative_tolerance"  : 1.0e-4,
            "frictional_slip_contact_residual_absolute_tolerance"  : 1.0e-9
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
        return "displacement_lagrangemultiplier_ressidual_frictional_contact_criteria";
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
        return "DisplacementLagrangeMultiplierResidualFrictionalContactCriteria";
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

        // The normal contact residual
        mLMNormalRatioTolerance =  ThisParameters["contact_displacement_absolute_tolerance"].GetDouble();
        mLMNormalAbsTolerance =  ThisParameters["contact_residual_absolute_tolerance"].GetDouble();

        // The tangent contact residual
        mLMTangentStickRatioTolerance =  ThisParameters["frictional_stick_contact_residual_relative_tolerance"].GetDouble();
        mLMTangentStickAbsTolerance =  ThisParameters["frictional_stick_contact_residual_absolute_tolerance"].GetDouble();
        mLMTangentSlipRatioTolerance =  ThisParameters["frictional_slip_contact_residual_relative_tolerance"].GetDouble();
        mLMTangentSlipAbsTolerance =  ThisParameters["frictional_slip_contact_residual_absolute_tolerance"].GetDouble();

        // We get the  ratio between the normal and tangent that will accepted as converged
        mNormalTangentRatio = ThisParameters["ratio_normal_tangent_threshold"].GetDouble();

        // Set local flags
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ENSURE_CONTACT, ThisParameters["ensure_contact"].GetBool());
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PRINTING_OUTPUT, ThisParameters["print_convergence_criterion"].GetBool());
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::TABLE_IS_INITIALIZED, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::ROTATION_DOF_IS_CONSIDERED, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::PURE_SLIP, ThisParameters["pure_slip"].GetBool());
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_NORMAL_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_STICK_RESIDUAL_IS_SET, false);
        mOptions.Set(DisplacementLagrangeMultiplierResidualFrictionalContactCriteria::INITIAL_SLIP_RESIDUAL_IS_SET, false);
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Flags mOptions; /// Local flags

    double mDispRatioTolerance;                /// The ratio threshold for the norm of the displacement residual
    double mDispAbsTolerance;                  /// The absolute value threshold for the norm of the displacement residual
    double mDispInitialResidualNorm;           /// The reference norm of the displacement residual
    double mDispCurrentResidualNorm;           /// The current norm of the displacement residual

    double mRotRatioTolerance;                /// The ratio threshold for the norm of the rotation residual
    double mRotAbsTolerance;                  /// The absolute value threshold for the norm of the rotation residual
    double mRotInitialResidualNorm;           /// The reference norm of the rotation residual
    double mRotCurrentResidualNorm;           /// The current norm of the rotation residual

    double mLMNormalRatioTolerance;            /// The ratio threshold for the norm of the normal LM residual
    double mLMNormalAbsTolerance;              /// The absolute value threshold for the norm of the normal LM  residual
    double mLMNormalInitialResidualNorm;       /// The reference norm of the normal LM residual
    double mLMNormalCurrentResidualNorm;       /// The current norm of the normal LM residual

    double mLMTangentStickRatioTolerance;      /// The ratio threshold for the norm of the tangent LM residual (stick)
    double mLMTangentStickAbsTolerance;        /// The absolute value threshold for the norm of the tangent LM  residual (stick)
    double mLMTangentSlipRatioTolerance;       /// The ratio threshold for the norm of the tangent LM residual (slip)
    double mLMTangentSlipAbsTolerance;         /// The absolute value threshold for the norm of the tangent LM  residual (slip)
    double mLMTangentStickInitialResidualNorm; /// The reference norm of the tangent LM residual (stick)
    double mLMTangentStickCurrentResidualNorm; /// The current norm of the tangent LM residual (stick)
    double mLMTangentSlipInitialResidualNorm;  /// The reference norm of the tangent LM residual (slip)
    double mLMTangentSlipCurrentResidualNorm;  /// The current norm of the tangent LM residual (slip)

    std::size_t mStickCounter = 0;                /// This is an auxiliary counter for stick dofs
    std::size_t mSlipCounter = 0;                 /// This is an auxiliary counter for slip dofs

    double mNormalTangentRatio;                /// The ratio to accept a non converged tangent component in case

    std::vector<int> mActiveDofs;                 /// This vector contains the dofs that are active

    ///@}
};  // Kratos DisplacementLagrangeMultiplierResidualFrictionalContactCriteria

///@name Local flags creation
///@{

/// Local Flags
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::ENSURE_CONTACT(Kratos::Flags::Create(0));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::PRINTING_OUTPUT(Kratos::Flags::Create(1));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::TABLE_IS_INITIALIZED(Kratos::Flags::Create(2));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::ROTATION_DOF_IS_CONSIDERED(Kratos::Flags::Create(3));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::PURE_SLIP(Kratos::Flags::Create(4));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::INITIAL_RESIDUAL_IS_SET(Kratos::Flags::Create(5));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::INITIAL_NORMAL_RESIDUAL_IS_SET(Kratos::Flags::Create(6));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::INITIAL_STICK_RESIDUAL_IS_SET(Kratos::Flags::Create(7));
template<class TSparseSpace, class TDenseSpace>
const Kratos::Flags DisplacementLagrangeMultiplierResidualFrictionalContactCriteria<TSparseSpace, TDenseSpace>::INITIAL_SLIP_RESIDUAL_IS_SET(Kratos::Flags::Create(8));
}