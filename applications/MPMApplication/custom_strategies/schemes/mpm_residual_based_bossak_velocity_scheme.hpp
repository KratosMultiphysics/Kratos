//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


#if !defined(KRATOS_MPM_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME )
#define      KRATOS_MPM_RESIDUAL_BASED_BOSSAK_VELOCITY_SCHEME

//#pragma once //??

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residual_based_implicit_time_scheme.h"
//#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "custom_utilities/mpm_boundary_rotation_utility.h"
#include "includes/checks.h"

namespace Kratos
{

/**
 * @class MPMResidualBasedBossakScheme
 * @ingroup KratosMPM
 * @brief Bossak integration scheme (for linear and nonlinear dynamic problems) for displacements adjusted for Material Point Method
 * @details This is an implicit scheme based of the Bossak algorithm for displacements suitable for quasi-static and dynamic problems.
 * Furthermore, this scheme has been adjusted for mixed formulation where pressure is also solved as one of the DoFs.
 * The parameter Alpha of Bossak introduces damping, the value of Bossak is from 0 to -0.5 (negative)
 * Implementation according to: "An alpha modification of Newmark's method; W.L. Wood, M. Bossak, O.C. Zienkiewicz;
 * Numerical Methods in Engineering; 1980"
 * MPM implementation according to: "Implicit time integration for the material point method"; J. E. Guilkey, J. A. Weiss
 * The parameter Alpha of Bossak introduces damping, the value of Bossak is from 0 to -0.5 (negative)
 */
template<class TSparseSpace,  class TDenseSpace >
class MPMResidualBasedBossakVelocityScheme
    : public ResidualBasedImplicitTimeScheme<TSparseSpace,TDenseSpace> //°?
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( MPMResidualBasedBossakVelocityScheme );

   /// Base type for the scheme
    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    /// Implicit base type for the scheme
    using ImplicitBaseType = ResidualBasedImplicitTimeScheme<TSparseSpace, TDenseSpace>;

    /// Class type for the scheme
    using ClassType = MPMResidualBasedBossakVelocityScheme<TSparseSpace, TDenseSpace>;

    /// Data type used within the ImplicitBaseType
    using TDataType = typename ImplicitBaseType::TDataType;

    /// Array type for degrees of freedom within ImplicitBaseType
    using DofsArrayType = typename ImplicitBaseType::DofsArrayType;

    /// Vector type for degrees of freedom within an Element
    using DofsVectorType = typename Element::DofsVectorType;

    /// Type for the system matrix within ImplicitBaseType
    using TSystemMatrixType = typename ImplicitBaseType::TSystemMatrixType;

    /// Type for the system vector within ImplicitBaseType
    using TSystemVectorType = typename ImplicitBaseType::TSystemVectorType;

    /// Type for local system vectors within ImplicitBaseType
    using LocalSystemVectorType = typename ImplicitBaseType::LocalSystemVectorType;

    /// Type for local system matrices within ImplicitBaseType
    using LocalSystemMatrixType = typename ImplicitBaseType::LocalSystemMatrixType;

    /// Iterator for nodes in a ModelPart
    using NodeIterator = ModelPart::NodeIterator;

    /// Container type for nodes in a ModelPart
    using NodesArrayType = ModelPart::NodesContainerType;

    /// Container type for elements in a ModelPart
    using ElementsArrayType = ModelPart::ElementsContainerType;

    /// Container type for conditions in a ModelPart
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    /// Pointer type for the BaseType
    using BaseTypePointer = typename BaseType::Pointer;

    /// Component type as 'double'
    using ComponentType = double;

    using ImplicitBaseType::mMatrix;


    ///@}
    ///@name Life Cycle
    ///@{

    /** @brief Construct from a @ref Parameters object.
     *  @param ThisParameters The parameters containing the configuration.
     */
    //explicit MPMResidualBasedBossakVelocityScheme(Parameters ThisParameters)
    //    : ImplicitBaseType()
    //{
    //    // Validate and assign defaults
    //    ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
    //    this->AssignSettings(ThisParameters);
//
    //    // For pure Newmark Scheme
    //    mNewmark.gamma = 0.5;
//
    //    AuxiliarInitializeBossak();
    //}
    
//    explicit MPMResidualBasedBossakVelocityScheme(Parameters ThisParameters, ModelPart& rGridModelPart, unsigned int DomainSize,
//        unsigned int BlockSize, double Alpha = 0.0,
//        double NewmarkBeta = 0.25, bool IsDynamic = true)
//            : ImplicitBaseType(),
//            mGridModelPart(rGridModelPart), mRotationTool(DomainSize, BlockSize, IS_STRUCTURE)
//    {
//        // Validate and assign defaults
//        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
//        this->AssignSettings(ThisParameters);
//
//        // For pure Newmark Scheme
//        mNewmark.gamma = 0.5;
//
//        AuxiliarInitializeBossak();
//
//        // To distinguish quasi-static and dynamic
//        mIsDynamic = IsDynamic;
//
//        // For Rotation Utility
//        mDomainSize = DomainSize;
//        mBlockSize  = BlockSize;
//    }

//    explicit MPMResidualBasedBossakVelocityScheme(ModelPart& rGridModelPart, unsigned int DomainSize,
//        unsigned int BlockSize, Parameters ThisParameters)
//            : ImplicitBaseType()
//            ,mGridModelPart(rGridModelPart), mRotationTool(DomainSize, BlockSize, IS_STRUCTURE)
//    {
// 
//        // For Rotation Utility
//        mDomainSize = DomainSize;
//        mBlockSize  = BlockSize;
// 
//        // Validate and assign defaults
//        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
//        this->AssignSettings(ThisParameters);
// 
//        // For pure Newmark Scheme
//        mNewmark.gamma = 0.5;
// 
//        AuxiliarInitializeBossak();
//    }

    /** @brief Constructor from a Bossak parameter.
     *  @param Alpha is the Bossak parameter. Default value is 0, which is the Newmark method.
     *  @note The Newmark beta parameter is set to 0.25, for mean constant acceleration.
     */

//    explicit MPMResidualBasedBossakVelocityScheme(ModelPart& rGridModelPart, unsigned int DomainSize,
//        unsigned int BlockSize, const double Alpha = 0.0, double NewmarkBeta = 0.25, bool IsDynamic = true)
//            : ImplicitBaseType(),
//            mGridModelPart(rGridModelPart), mRotationTool(DomainSize, BlockSize, IS_STRUCTURE)
//    {
//        // Validate and assign defaults
//        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
//        this->AssignSettings(ThisParameters);
//
//        // For pure Newmark Scheme
//        mNewmark.gamma = 0.5;
//
//        AuxiliarInitializeBossak();
//
//        // To distinguish quasi-static and dynamic
//        mIsDynamic = IsDynamic;
//
//        // For Rotation Utility
//        mDomainSize = DomainSize;
//        mBlockSize  = BlockSize;
//    }

//    explicit MPMResidualBasedBossakVelocityScheme(ModelPart& rGridModelPart, unsigned int DomainSize,
//        unsigned int BlockSize, const double Alpha = 0.0, const double NewmarkBeta = 0.25, bool IsDynamic = true)
//            : MPMResidualBasedBossakVelocityScheme(rGridModelPart, DomainSize, BlockSize, Alpha, NewmarkBeta)
//    {
//        // To distinguish quasi-static and dynamic
//        mIsDynamic = IsDynamic;
//    }

    /** @brief Constructor.
     *  @param Alpha The Bossak parameter.
     *  @param NewarkBeta The Newmark parameter.
     *  @note A Bossak parameter (@a Alpha) of 0 results in the Newmark method.
     *  @note A Newmark parameter (@a NewmarkBeta) of 0.25 results in mean constant acceleration.
     */
//    explicit MPMResidualBasedBossakVelocityScheme(const double Alpha, const double NewmarkBeta,
//        ModelPart& rGridModelPart, unsigned int DomainSize,
//        unsigned int BlockSize, bool IsDynamic = true)
//            :ImplicitBaseType(),
//            mGridModelPart(rGridModelPart), mRotationTool(DomainSize, BlockSize, IS_STRUCTURE)
//    {
//        // For pure Newmark Scheme
//        mBossak.alpha = Alpha;
//        mNewmark.beta = NewmarkBeta;
//        mNewmark.gamma = 0.5;
//
//        AuxiliarInitializeBossak();
//
//        // To distinguish quasi-static and dynamic
//        mIsDynamic = IsDynamic;
//
//        // For Rotation Utility
//        mDomainSize = DomainSize;
//        mBlockSize  = BlockSize;
//    }
//    explicit MPMResidualBasedBossakVelocityScheme(ModelPart& rGridModelPart, unsigned int DomainSize,
//        unsigned int BlockSize, const double Alpha, const double NewmarkBeta)
//            :ImplicitBaseType()
//            ,mGridModelPart(rGridModelPart), mRotationTool(DomainSize, BlockSize, IS_STRUCTURE)
//    {
//
//        // For Rotation Utility
//        mDomainSize = DomainSize;
//        mBlockSize  = BlockSize;
//
//        // For pure Newmark Scheme
//        mBossak.alpha = Alpha;
//        mNewmark.beta = NewmarkBeta;
//        mNewmark.gamma = 0.5;
//
//        AuxiliarInitializeBossak();
//    }

    explicit MPMResidualBasedBossakVelocityScheme(ModelPart& rGridModelPart, unsigned int DomainSize,
        unsigned int BlockSize, bool IsDynamic, const double Alpha = 0.0, const double NewmarkBeta = 0.25)
            : ImplicitBaseType(),
            mGridModelPart(rGridModelPart), mRotationTool(DomainSize, BlockSize, IS_STRUCTURE)
    {
        // To distinguish quasi-static and dynamic
        mIsDynamic = IsDynamic;

        // For Rotation Utility
        mDomainSize = DomainSize;
        mBlockSize  = BlockSize;

        // For pure Newmark Scheme
        //mBossak.alpha = Alpha; //perchè Alpha non è zero?
        mBossak.alpha = 0.0;
        mNewmark.beta = NewmarkBeta;
        mNewmark.gamma = 0.5;

        AuxiliarInitializeBossak();
    }
    
    /**
     * @brief Copy Constructor.
     */
    // ° ?
    explicit MPMResidualBasedBossakVelocityScheme(MPMResidualBasedBossakVelocityScheme& rOther)
        :ImplicitBaseType(rOther)
        ,mBossak(rOther.mBossak)
        ,mNewmark(rOther.mNewmark)
        ,mVector(rOther.mVector)
        ,mGridModelPart(rOther.mGridModelPart)
        ,mRotationTool(rOther.mDomainSize,rOther.mBlockSize,IS_STRUCTURE)
    {
    }

    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new MPMResidualBasedBossakVelocityScheme(*this) );
    }

    ~MPMResidualBasedBossakVelocityScheme
    () override {}

    ///@}
    ///@name Operations
    ///@{

    /** @brief Construct a dynamically allocated new scheme from a @ref Parameters object.
     *  @param ThisParameters The configuration parameters.
     */
    //typename BaseType::Pointer Create(Parameters ThisParameters) const override //°a cosa serve? Come sistemarlo?
    //{
    //    //return Kratos::make_shared<ClassType>(ThisParameters);
    //    return Kratos::make_shared<ImplicitBaseType>(ThisParameters);
    //}

    /// @brief Recalculate the Newmark coefficients, taking the alpha parameters into account.
    void CalculateBossakCoefficients()
    {
        mBossak.beta  = (1.0 - mBossak.alpha) * (1.0 - mBossak.alpha) * mNewmark.beta;
        mBossak.gamma = mNewmark.gamma  - mBossak.alpha;
    }

    /** @brief Update state variables within a newton iteration at the end of the time step.
     *  @details @f[ u_{n+1}^{k+1} = u_{n+1}^k+ \Delta u @f]
     *  @param rModelPart @ref ModelPart to update.
     *  @param rDofSet Set of all primary variables.
     *  @param rA Left hand side matrix.
     *  @param rDx Primary variable updates.
     *  @param rb Right hand side vector.
     */
    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        // Rotate the current VELOCITY to the modified coordinate system since rDx is currently at the modified coordinate system
        mRotationTool.RotateVelocities(rModelPart);

        // Update of VELOCITY (by DOF)
        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        // Rotate the VELOCITY back to the original coordinate system to calculate the velocity and acceleration
        mRotationTool.RecoverVelocities(rModelPart);

        // Updating time derivatives (nodally for efficiency)
        const int num_nodes = static_cast<int>( rModelPart.Nodes().size() );
        const auto it_node_begin = rModelPart.Nodes().begin();

        #pragma omp parallel for
        for(int i = 0;  i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            // In MPM the delta_displacement is the current displacement as the previous_displacement is always reset //° ??
            array_1d<double, 3 > & r_current_displacement = it_node->FastGetSolutionStepValue(DISPLACEMENT);

            array_1d<double, 3>& r_current_velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& r_previous_velocity = it_node->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3>& r_current_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION);
            const array_1d<double, 3>& r_previous_acceleration = it_node->FastGetSolutionStepValue(ACCELERATION, 1);

            if (mIsDynamic){
                UpdateDisplacement(r_current_displacement, r_current_velocity, r_previous_velocity, r_previous_acceleration);
                UpdateAcceleration(r_current_acceleration, r_current_velocity, r_previous_velocity, r_previous_acceleration);
            }
            // add update displacement for static case
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Apply the predictor.
     * @details @f[ x_{k+1} = x_{k} + v_{k} \cdot \Delta t @f]
     * @param rModelPart @ref ModelPart to update.
     * @param rDofSet Set of all primary variables.
     * @param rA Left hand side matrix.
     * @param rDx Primary variable updates.
     * @param rb Right hand side vector.
     */
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY;

		#pragma omp parallel for
		for(int iter = 0; iter < static_cast<int>(rModelPart.Nodes().size()); ++iter)
		{
			auto i = rModelPart.NodesBegin() + iter;
            const array_1d<double, 3 > & r_previous_displacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
			const array_1d<double, 3 > & r_previous_velocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);
            const array_1d<double, 3 > & r_previous_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

            array_1d<double, 3 > & r_current_velocity  = (i)->FastGetSolutionStepValue(VELOCITY);

            // Velocity prediction for implicit MPM
            if (!(i->pGetDof(VELOCITY_X)->IsFixed()))
                r_current_velocity[0] = 0.0;
            else
                r_current_velocity[0]  = r_previous_velocity[0];

            if (!(i->pGetDof(VELOCITY_Y)->IsFixed()))
                r_current_velocity[1] = 0.0;
            else
                r_current_velocity[1]  = r_previous_velocity[1];

            if (i->HasDofFor(VELOCITY_Z))
            {
                if (!(i->pGetDof(VELOCITY_Z)->IsFixed()))
                    r_current_velocity[2] = 0.0;
                else
                    r_current_velocity[2]  = r_previous_velocity[2];
            }

            // Pressure prediction for implicit MPM
            if (i->HasDofFor(PRESSURE))
            {
                double& r_current_pressure        = (i)->FastGetSolutionStepValue(PRESSURE);
                const double& r_previous_pressure = (i)->FastGetSolutionStepValue(PRESSURE, 1);

                if (!(i->pGetDof(PRESSURE))->IsFixed())
                    r_current_pressure = r_previous_pressure;
            }

            // Updating time derivatives
            array_1d<double, 3 > & current_displacement    = (i)->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & current_acceleration     = (i)->FastGetSolutionStepValue(ACCELERATION);

            if (mIsDynamic){
                UpdateDisplacement(current_displacement, r_current_velocity, r_previous_velocity, r_previous_acceleration);
                UpdateAcceleration(current_acceleration, r_current_velocity, r_previous_velocity, r_previous_acceleration);
            } //else {
            //    update displcement
            //}

		}

        KRATOS_CATCH( "" );
    }

    /** @brief Prepare state variables for a new solution step.
     *  @note This function should be called before solving a solution step,
     *        or restarting a failed one.
     *  @param rModelPart @ref ModelPart to update.
     *  @param A Left hand side matrix.
     *  @param Dx Primary variable updates.
     *  @param b Right hand side vector.
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY

        // Loop over the grid nodes performed to clear all nodal information
		#pragma omp parallel for
		for(int iter = 0; iter < static_cast<int>(mGridModelPart.Nodes().size()); ++iter)
		{
			auto i = mGridModelPart.NodesBegin() + iter;

            // Variables to be cleaned
            double & r_nodal_mass     = (i)->FastGetSolutionStepValue(NODAL_MASS);
            array_1d<double, 3 > & r_nodal_momentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
            array_1d<double, 3 > & r_nodal_inertia  = (i)->FastGetSolutionStepValue(NODAL_INERTIA);

            array_1d<double, 3 > & r_nodal_displacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > & r_nodal_velocity     = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & r_nodal_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);

            double & r_nodal_old_pressure = (i)->FastGetSolutionStepValue(PRESSURE,1);
            double & r_nodal_pressure = (i)->FastGetSolutionStepValue(PRESSURE);

            // Clear
            r_nodal_mass = 0.0;
            r_nodal_momentum.clear();
            r_nodal_inertia.clear();

            r_nodal_displacement.clear();
            r_nodal_velocity.clear();
            r_nodal_acceleration.clear();
            r_nodal_old_pressure = 0.0;
            r_nodal_pressure = 0.0;

            // Other additional variables
            if ((i)->SolutionStepsDataHas(NODAL_AREA)){
                double & r_nodal_area = (i)->FastGetSolutionStepValue(NODAL_AREA);
                r_nodal_area          = 0.0;
            }
            if(i->SolutionStepsDataHas(NODAL_MPRESSURE)) {
                double & r_nodal_mpressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                r_nodal_mpressure = 0.0;
            }
		}

        // Extrapolate from Material Point Elements and Conditions
        ImplicitBaseType::InitializeSolutionStep(rModelPart,rA,rDx,rb);

        // Assign nodal variables after extrapolation
        #pragma omp parallel for
        for(int iter = 0; iter < static_cast<int>(mGridModelPart.Nodes().size()); ++iter)
        {
            auto i = mGridModelPart.NodesBegin() + iter;
            const double & r_nodal_mass = (i)->FastGetSolutionStepValue(NODAL_MASS);

            if (r_nodal_mass > std::numeric_limits<double>::epsilon())
            {
                const array_1d<double, 3 > & r_nodal_momentum   = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                const array_1d<double, 3 > & r_nodal_inertia    = (i)->FastGetSolutionStepValue(NODAL_INERTIA);

                array_1d<double, 3 > & r_nodal_velocity_1     = (i)->FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double, 3 > & r_nodal_acceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
                double & r_nodal_pressure = (i)->FastGetSolutionStepValue(PRESSURE,1);

                double delta_nodal_pressure = 0.0;

                // For mixed formulation
                if (i->HasDofFor(PRESSURE) && i->SolutionStepsDataHas(NODAL_MPRESSURE))
                {
                    double & nodal_mpressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                    delta_nodal_pressure = nodal_mpressure/r_nodal_mass;
                }

                const array_1d<double, 3 > delta_nodal_velocity = r_nodal_momentum/r_nodal_mass;
                const array_1d<double, 3 > delta_nodal_acceleration = r_nodal_inertia/r_nodal_mass;

                r_nodal_velocity_1 += delta_nodal_velocity;
                r_nodal_acceleration += delta_nodal_acceleration;

                r_nodal_pressure += delta_nodal_pressure;
            }
        }

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const double delta_time = r_current_process_info[DELTA_TIME];

        // Initializing Bossak constants
        mBossak.c0 = ( 1.0 / ( delta_time * mBossak.gamma) );
        mBossak.c1 = ( (1.0 - mBossak.gamma) / mBossak.gamma );
        mBossak.c2 = ( delta_time * mBossak.beta / mBossak.gamma );
        mBossak.c3 = ( delta_time * (1.0 - mBossak.beta / mBossak.gamma) );
        mBossak.c4 = ( delta_time * delta_time * ( 0.5 - mBossak.beta / mBossak.gamma) );
        mBossak.c5 = ( 2 / delta_time );

        KRATOS_CATCH( "" )
    }

    /** @brief Check whether the scheme and the provided @ref ModelPart are configured correctly.
     *  @details Checks may be expensive as the function is designed to catch user errors.
     *  @param rModelPart @ref ModelPart to check.
     *  @return 0 if ok, nonzero otherwise.
     */
    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY;

        const int err = ImplicitBaseType::Check(rModelPart);
        if(err != 0) return err;

        // Check that variables are correctly allocated
        for (const auto& rnode : rModelPart.Nodes()) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z, rnode)
        }

        // Check for minimum value of the buffer index
        // Verify buffer size
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2)
            << "Insufficient buffer size. Buffer size should be greater than 2. Current size is: "
            << rModelPart.GetBufferSize() << std::endl;

        // Check for admissible value of the AlphaBossak
         KRATOS_ERROR_IF(mBossak.alpha > 0.0 || mBossak.alpha < -0.5) << "Value not admissible for "
            << "AlphaBossak. Admissible values are between 0.0 and -0.5\nCurrent value is: "
            << mBossak.alpha << std::endl;

        static const double epsilon = 1e-12;
        KRATOS_ERROR_IF_NOT(std::abs(mNewmark.beta - 0.0)   < epsilon ||
                            std::abs(mNewmark.beta - 0.167) < epsilon ||
                            std::abs(mNewmark.beta - 0.25)  < epsilon)
            << "Value not admissible for NewmarkBeta. Admissible values are:\n"
            << "0.0 for central-differencing\n"
            << "0.25 for mean-constant-acceleration\n"
            << "0.167 for linear-acceleration\n"
            << "Current value is: " << mNewmark.beta << std::endl;

        return 0;
        KRATOS_CATCH( "" );
    }

    /// @brief Release dynamic memory allocated by this instance.
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    /// @brief This function returns the default @ref Parameters to help avoiding conflicts between different constructors.
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"          : "bossak_scheme",
            "damp_factor_m" : -0.3,
            "newmark_beta"  : 0.25
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = ImplicitBaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /// @brief Return the name of the class as used in the settings (@a snake_case).
    static std::string Name()
    {
        return "bossak_scheme";
    }

    ///@}
    ///@name Input and output
    ///@{

    /// @brief Return information as a string.
    std::string Info() const override
    {
        return "MPMResidualBasedBossakVelocityScheme";
    }

    /// @brief Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// @brief Print the instance's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}

// da qui in poi: mpm residual displacement
    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rCurrentElement The element to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& rConstElemRef = rCurrentElement;
        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        rConstElemRef.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            rCurrentElement.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            //rCurrentElement.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);
            //AddDynamicsToLHS(LHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
            //AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
            AddDynamicsToLHS(LHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
            AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);

        }

        // If there is a slip condition, apply it on a rotated system of coordinates
        mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement.GetGeometry());
        mRotationTool.ElementApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement.GetGeometry());

        KRATOS_CATCH( "" )
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentElement The element to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {

        KRATOS_TRY

        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_elem_ref = rCurrentElement;

        // Basic operations for the element considered
        rCurrentElement.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
        r_const_elem_ref.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            rCurrentElement.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            //rCurrentElement.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);
            //AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
            AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
        }

        // If there is a slip condition, apply it on a rotated system of coordinates
        mRotationTool.RotateRHS(RHS_Contribution,rCurrentElement.GetGeometry());
        mRotationTool.ElementApplySlipCondition(RHS_Contribution,rCurrentElement.GetGeometry());

        KRATOS_CATCH( "" )
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCurrentCondition The condition to compute
     * @param rLHSContribution The LHS matrix contribution
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {

        KRATOS_TRY

        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCurrentCondition;

        rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        r_const_cond_ref.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            rCurrentCondition.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            //rCurrentCondition.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);
            //AddDynamicsToLHS(LHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
            //AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
            AddDynamicsToLHS(LHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
            AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);
        }

        // Rotate contributions (to match coordinates for slip conditions)
        mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition.GetGeometry());
        mRotationTool.ConditionApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition.GetGeometry());

        KRATOS_CATCH( "" )
    }

    /**
     * @brief Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition The condition to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const IndexType this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCurrentCondition;
        rCurrentCondition.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
        r_const_cond_ref.EquationIdVector(EquationId,rCurrentProcessInfo);

        if(mIsDynamic)
        {
            rCurrentCondition.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);
            //rCurrentCondition.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);
            //AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);
            AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mMatrix.M[this_thread], rCurrentProcessInfo);        
        }

        // Rotate contributions (to match coordinates for slip conditions)
        mRotationTool.RotateRHS(RHS_Contribution,rCurrentCondition.GetGeometry());
        mRotationTool.ConditionApplySlipCondition(RHS_Contribution,rCurrentCondition.GetGeometry());

        KRATOS_CATCH( "" )
    }



protected:
    /// @todo Move to @ref ImplicitBaseType
    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    /// @brief Bossak Alpha parameters.
    struct BossakAlphaMethod
    {
        /// @brief Bossak Alpha.
        double alpha;

        /// @brief Bossak Beta.
        double beta;

        /// @brief Bossak Gamma.
        double gamma;

        /// @brief System constants
        double c0, c1, c2, c3, c4, c5;
    };

    /// @brief Newmark parameters used for integration.
    struct NewmarkMethod
    {
        /// @brief Newmark Beta.
        double beta;

        /// @brief Newmark Gamma.
        double gamma;
    };

    /// @brief Velocities and accelerations used for integration.
    struct GeneralVectors
    {
        /// @brief Velocity.
        std::vector< Vector > v;

        /// @brief Acceleration.
        std::vector< Vector > a;

        /// @brief Previous acceleration.
        std::vector< Vector > ap;
    };

    /// @brief Bossak Alpha parameters.
    BossakAlphaMethod mBossak;

    /// @brief Newmark Beta parameters.
    NewmarkMethod mNewmark;

    /// @brief Aggregate struct for velocities and accelerations.
    GeneralVectors mVector;

    ///@name Protected Operations
    ///@{

    /** @brief Update displacement.
     *  @param rCurrentDisplacement Current velocity.
     *  @param rDeltaDisplacement Displacement increment.
     *  @param rPreviousVelocity Previous velocity.
     *  @param rPreviousAcceleration Previous acceleration.
     */
    inline void UpdateDisplacement(
        array_1d<double, 3>& rCurrentDisplacement,
        array_1d<double, 3>& rCurrentVelocity,
        const array_1d<double, 3>& rPreviousVelocity,
        const array_1d<double, 3>& rPreviousAcceleration
        )
    {
        noalias(rCurrentDisplacement) = mBossak.c2 * rCurrentVelocity + mBossak.c3 * rPreviousVelocity + mBossak.c4 * rPreviousAcceleration;
    }

    /** @brief Update the second time derivative.
     *  @param rCurrentAcceleration Current velocity.
     *  @param rDeltaDisplacement Displacement increment.
     *  @param rPreviousVelocity Previous velocity.
     *  @param rPreviousAcceleration Previous acceleration.
     */
    inline void UpdateAcceleration(
        array_1d<double, 3>& rCurrentAcceleration,
        array_1d<double, 3>& rCurrentVelocity,
        const array_1d<double, 3>& rPreviousVelocity,
        const array_1d<double, 3>& rPreviousAcceleration
        )
    {
        noalias(rCurrentAcceleration) = mBossak.c0 * (rCurrentVelocity - rPreviousVelocity) - mBossak.c1 * rPreviousAcceleration;
    }

    /** @brief Add dynamic left hand side contribution from @ref Element s.
     *  @details @f[ M \cdot c_0 + D \cdot c_1 + K @f]
     *  @param LHS_Contribution Dynamic contribution for the left hand side.
     *  @param D Damping matrix.
     *  @param M Mass matrix.
     *  @param rCurrentProcessInfo Current @ref ProcessInfo.
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        //LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        // Adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
            //noalias(LHS_Contribution) += M * (1.0 - mBossak.alpha) * mBossak.c0; //°?
            noalias(LHS_Contribution) += M * mBossak.c0;
        // Adding  damping contribution
        //if (D.size1() != 0) // if D matrix declared
        //    noalias(LHS_Contribution) += D * mBossak.c1;
    }

    /** @brief Add dynamic right hand side contribution from an @ref Element.
     *  @details @f[ b - (1 - \alpha) \cdot M \cdot a_{n+1} - \alpha \cdot M \cdot a_n - D \cdot v_n @f]
     *  @param rElement @ref Element to compute the contribution from.
     *  @param RHS_Contribution Dynamic contribution for the right hand side.
     *  @param D Damping matrix.
     *  @param M Mass matrix.
     *  @param rCurrentProcessInfo Current @ref ProcessInfo.
     */
    void AddDynamicsToRHS( //°incompleto
        Element& rElement,
        LocalSystemVectorType& RHS_Contribution,
        //LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        const auto& r_const_elem_ref = rElement;
        // Adding inertia contribution
        if (M.size1() != 0) {

            //r_const_elem_ref.GetSecondDerivativesVector(mVector.a[this_thread], 0);
            //mVector.a[this_thread] *= (1.00 - mBossak.alpha);
//
            //r_const_elem_ref.GetSecondDerivativesVector(mVector.ap[this_thread], 1);
            //noalias(mVector.a[this_thread]) += mBossak.alpha * mVector.ap[this_thread];
//
            //noalias(RHS_Contribution) -= prod(M, mVector.a[this_thread]);
        }

        // Adding damping contribution
        //if (D.size1() != 0) {
        //    r_const_elem_ref.GetFirstDerivativesVector(mVector.v[this_thread], 0);
        //    noalias(RHS_Contribution) -= prod(D, mVector.v[this_thread]);
        //}
    }

    /** @brief Add dynamic right hand side contribution of a @ref Condition.
     *  @details @f[ b - (1-alpha)*M*a_n+1 - alpha*M*a_n - D*v_n @f]
     *  @param rCondition @ref Condition to compute the contribution from.
     *  @param RHS_Contribution Dynamic contribution for the right hand side.
     *  @param D Damping matrix.
     *  @param M Mass matrix.
     *  @param rCurrentProcessInfo Current @ref ProcessInfo.
     */
    void AddDynamicsToRHS( //°incompleto
        Condition& rCondition,
        LocalSystemVectorType& RHS_Contribution,
        //LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        const std::size_t this_thread = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCondition;

        // Adding inertia contribution
        //if (M.size1() != 0) {
        //    r_const_cond_ref.GetSecondDerivativesVector(mVector.a[this_thread], 0);
        //    mVector.a[this_thread] *= (1.00 - mBossak.alpha);
//
        //    r_const_cond_ref.GetSecondDerivativesVector(mVector.ap[this_thread], 1);
        //    noalias(mVector.a[this_thread]) += mBossak.alpha * mVector.ap[this_thread];
//
        //    noalias(RHS_Contribution) -= prod(M, mVector.a[this_thread]);
        //}

        // Adding damping contribution
        // Damping contribution
        //if (D.size1() != 0) {
        //    r_const_cond_ref.GetFirstDerivativesVector(mVector.v[this_thread], 0);
//
        //    noalias(RHS_Contribution) -= prod(D, mVector.v[this_thread]);
        //}
    }

    /// @brief Assign member variables from @ref Parameters.
    void AssignSettings(const Parameters ThisParameters) override
    {
        ImplicitBaseType::AssignSettings(ThisParameters);
        mBossak.alpha = ThisParameters["damp_factor_m"].GetDouble();
        mNewmark.beta = ThisParameters["newmark_beta"].GetDouble();
    }

    ///@}

    // da qui in poi: mpm residual displacement
    // MPM Background Grid
    ModelPart& mGridModelPart;

    // To distinguish quasi-static and dynamic
    bool mIsDynamic;

    // For Rotation Utility
    unsigned int mDomainSize;
    unsigned int mBlockSize;
    MPMBoundaryRotationUtility<LocalSystemMatrixType,LocalSystemVectorType> mRotationTool;

    void ClearReaction() const
    {
        #pragma omp parallel for
        for (int iter = 0; iter < static_cast<int>(mGridModelPart.Nodes().size()); ++iter) {
            auto i = mGridModelPart.NodesBegin() + iter;
            (i)->FastGetSolutionStepValue(REACTION).clear();
        }
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    /// @brief Utility function for initializing some members.
    void AuxiliarInitializeBossak()
    {
        // Initialize Bossak coefficients
        CalculateBossakCoefficients();

        // Allocate auxiliary memory
        const std::size_t num_threads = ParallelUtilities::GetNumThreads();

        mVector.v.resize(num_threads);
        mVector.a.resize(num_threads);
        mVector.ap.resize(num_threads);

        KRATOS_DETAIL("MECHANICAL SCHEME: The Bossak Time Integration Scheme ") << "[alpha_m= " << mBossak.alpha << " beta= " << mNewmark.beta << " gamma= " << mNewmark.gamma << "]" <<std::endl;
    }

    ///@}
}; // Class ResidualBasedBossakDisplacementScheme

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME defined */

