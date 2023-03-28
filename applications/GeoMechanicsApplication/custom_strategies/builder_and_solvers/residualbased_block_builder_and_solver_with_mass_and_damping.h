//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborators:   Vicente Mataix
//
//

#pragma once

/* System includes */
#include <unordered_set>

/* External includes */
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "utilities/timer.h"
#include "utilities/variable_utils.h"
#include "includes/kratos_flags.h"
#include "includes/lock_object.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/builtin_timer.h"
#include "utilities/atomic_utilities.h"
#include "spaces/ublas_space.h"

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
 * @class ResidualBasedBlockBuilderAndSolverWithMassAndDamping
 * @ingroup KratosCore
 * @brief Current class provides an implementation for standard builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similar to the calculation of the total residual
 * @tparam TSparseSpace The sparse system considered
 * @tparam TDenseSpace The dense system considered
 * @tparam TLinearSolver The linear solver considered
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockBuilderAndSolverWithMassAndDamping
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the flags
    KRATOS_DEFINE_LOCAL_FLAG( SILENT_WARNINGS );

    /// Definition of the pointer
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithMassAndDamping);

    /// Definition of the base class
    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// The definition of the current class
    typedef ResidualBasedBlockBuilderAndSolverWithMassAndDamping<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /// Additional definitions
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedBlockBuilderAndSolverWithMassAndDamping() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ResidualBasedBlockBuilderAndSolverWithMassAndDamping(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolverWithMassAndDamping(typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverWithMassAndDamping() override
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param b The RHS vector
     */
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemMatrixType mass_contribution = LocalSystemMatrixType(0, 0);
        LocalSystemMatrixType damping_contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);


        //
        // auto p_mass_matrix = TSparseSpace::CreateEmptyMatrixPointer();
        // auto p_damping_matrix = TSparseSpace::CreateEmptyMatrixPointer();
        //
        // mMassMatrix = *p_mass_matrix;
        // mDampingMatrix = *p_damping_matrix;

        mMassMatrix.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);

        ConstructMatrixStructure(pScheme, mMassMatrix, rModelPart);

        mDampingMatrix.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
        ConstructMatrixStructure(pScheme, mDampingMatrix, rModelPart);

        TSparseSpace::SetToZero(mMassMatrix);
        TSparseSpace::SetToZero(mDampingMatrix);


        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // Assemble all elements
        const auto timer = BuiltinTimer();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution,mass_contribution,damping_contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++) {
                auto it_elem = el_begin + k;

                if (it_elem->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_elem, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    
                    // assemble mass and damping matrix
                    it_elem->CalculateMassMatrix(mass_contribution, CurrentProcessInfo);
                    it_elem->CalculateDampingMatrix(damping_contribution, CurrentProcessInfo);

                    if (mass_contribution.size1() != 0)
                    {
                        AssembleLHS(mMassMatrix, mass_contribution, EquationId);
                    }
                    if (damping_contribution.size1() != 0)
                    {
                        AssembleLHS(mDampingMatrix, damping_contribution, EquationId);
                    }
                    

                    // Assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
                }

            }

            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++) {
                auto it_cond = cond_begin + k;

                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_cond, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    // assemble mass and damping matrix
                    it_cond->CalculateMassMatrix(mass_contribution, CurrentProcessInfo);
                    it_cond->CalculateDampingMatrix(damping_contribution, CurrentProcessInfo);

                    if (mass_contribution.size1() != 0)
                    {
                        AssembleLHS(mMassMatrix, mass_contribution, EquationId);
                    }
                    if (damping_contribution.size1() != 0)
                    {
                        AssembleLHS(mDampingMatrix, damping_contribution, EquationId);
                    }

                    // Assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >= 1) << "Build time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }


    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve
     * just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

        if(rModelPart.MasterSlaveConstraints().size() != 0) {
            const auto timer_constraints = BuiltinTimer();
            Timer::Start("ApplyConstraints");
            ApplyConstraints(pScheme, rModelPart, A, b);
            Timer::Stop("ApplyConstraints");
            KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >=1) << "Constraints build time: " << timer_constraints.ElapsedSeconds() << std::endl;
        }

        ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);
        ApplyDirichletConditions(pScheme, rModelPart, mMassMatrix, Dx, b);
        ApplyDirichletConditions(pScheme, rModelPart, mDampingMatrix, Dx, b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >=1) << "System solve time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        KRATOS_CATCH("")
    }


    /**
 * @brief Corresponds to the previews, but the System's matrix is considered already built and only the RHS is built again
 * @param pScheme The integration scheme considered
 * @param rModelPart The model part of the problem to solve
 * @param rA The LHS matrix
 * @param rDx The Unknowns vector
 * @param rb The RHS vector
 */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    ) override
    {
        KRATOS_TRY

            BuildRHS(pScheme, rModelPart, rb);

        if (rModelPart.MasterSlaveConstraints().size() != 0) {
            Timer::Start("ApplyRHSConstraints");
            ApplyRHSConstraints(pScheme, rModelPart, rb);
            Timer::Stop("ApplyRHSConstraints");
        }

        ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >= 1) << "System solve time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
    * @brief Function to perform the build of the RHS.
    * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
    * @param pScheme The integration scheme considered
    * @param rModelPart The model part of the problem to solve
    */
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

            Timer::Start("BuildRHS");

        BuildRHSNoDirichlet(pScheme, rModelPart, b);

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        block_for_each(BaseType::mDofSet, [&](Dof<double>& rDof) {
            const std::size_t i = rDof.EquationId();

            if (rDof.IsFixed())
                b[i] = 0.0;
            });

        Timer::Stop("BuildRHS");

        KRATOS_CATCH("")
    }


    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                 : "block_builder_and_solver_with_mass_and_damping",
            "block_builder"                        : true,
            "diagonal_values_for_dirichlet_dofs"   : "use_max_diagonal",
            "silent_warnings"                      : false
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
        return "block_builder_and_solver_with_mass_and_damping";
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
        return "ResidualBasedBlockBuilderAndSolverWithMassAndDamping";
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

    // TSystemMatrixType mT;                             /// This is matrix containing the global relation for the constraints
    // TSystemVectorType mConstantVector;                /// This is vector containing the rigid movement of the constraint
    // std::vector<IndexType> mSlaveIds;                 /// The equation ids of the slaves
    // std::vector<IndexType> mMasterIds;                /// The equation ids of the master
    // std::unordered_set<IndexType> mInactiveSlaveDofs; /// The set containing the inactive slave dofs
    TSystemMatrixType mMassMatrix;
    TSystemMatrixType mDampingMatrix;
    // double mScaleFactor = 1.0;                        /// The manually set scale factor

    // SCALING_DIAGONAL mScalingDiagonal = SCALING_DIAGONAL::NO_SCALING; /// We identify the scaling considered for the dirichlet dofs
    // Flags mOptions;                                                   /// Some flags used internally

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


    void GetFirstAndSecondDeriviativeVector(TSystemVectorType& rFirstDerivativeVector, TSystemVectorType& rSecondDerivativeVector, ModelPart& rModelPart)
    {
        NodesArrayType& rNodes = rModelPart.Nodes();
        const int n_nodes = static_cast<int>(rNodes.size());

        if (rFirstDerivativeVector.size() != BaseType::mEquationSystemSize) {
            rFirstDerivativeVector.resize(BaseType::mEquationSystemSize, false);
        }
        TSparseSpace::SetToZero(rFirstDerivativeVector);

        if (rSecondDerivativeVector.size() != BaseType::mEquationSystemSize) {
            rSecondDerivativeVector.resize(BaseType::mEquationSystemSize, false);
        }
        TSparseSpace::SetToZero(rSecondDerivativeVector);

        std::size_t id_x;
        std::size_t id_y;
        std::size_t id_z;

		#pragma omp parallel firstprivate(n_nodes, id_x, id_y, id_z)
        {
			#pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_nodes; i++) {
                typename NodesArrayType::iterator it = rNodes.begin() + i;
                // If the node is active
                if (it->IsActive()) {

                    id_x = it->GetDof(DISPLACEMENT_X).EquationId();
                    id_y = it->GetDof(DISPLACEMENT_Y).EquationId();

                    rFirstDerivativeVector[id_x] = it->FastGetSolutionStepValue(VELOCITY_X);
                    rFirstDerivativeVector[id_y] = it->FastGetSolutionStepValue(VELOCITY_Y);

                    rSecondDerivativeVector[id_x] = it->FastGetSolutionStepValue(ACCELERATION_X);
                    rSecondDerivativeVector[id_y] = it->FastGetSolutionStepValue(ACCELERATION_Y);

                    if (it->HasDofFor(DISPLACEMENT_Z))
                    {
                        id_z = it->GetDof(DISPLACEMENT_Z).EquationId();

                        rFirstDerivativeVector[id_z] = it->FastGetSolutionStepValue(VELOCITY_Z);
                        rSecondDerivativeVector[id_z] = it->FastGetSolutionStepValue(ACCELERATION_Z);
                    }
                }
            }
        }
    }


    void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = rModelPart.Conditions();

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationIds;

        // assemble all elements

        const int nelements = static_cast<int>(pElements.size());
        #pragma omp parallel firstprivate(nelements, RHS_Contribution, EquationIds)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i=0; i<nelements; i++) {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                // If the element is active
                if(it->IsActive()) {

                    //calculate elemental Right Hand Side Contribution
                    it->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
                    it->EquationIdVector(EquationIds, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationIds);
                }
            }

            RHS_Contribution.resize(0, false);

            // assemble all conditions
            const int nconditions = static_cast<int>(ConditionsArray.size());
            #pragma omp for schedule(guided, 512)
            for (int i = 0; i<nconditions; i++) {
                auto it = ConditionsArray.begin() + i;
                // If the condition is active
                if(it->IsActive()) {


                    it->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
                    it->EquationIdVector(EquationIds, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationIds);

                }
            }
        }

        TSystemVectorType firstDerivativeVector;
        TSystemVectorType secondDerivativeVector;

        GetFirstAndSecondDeriviativeVector(firstDerivativeVector, secondDerivativeVector, rModelPart);

        TSystemVectorType mass_contribution = prod(mMassMatrix, secondDerivativeVector);
        TSystemVectorType damping_contribution = prod(mDampingMatrix, firstDerivativeVector);

        b -= damping_contribution + mass_contribution;

        KRATOS_CATCH("")

    }


private:

}; /* Class ResidualBasedBlockBuilderAndSolver */

///@}

///@name Type Definitions
///@{

// Here one should use the KRATOS_CREATE_LOCAL_FLAG, but it does not play nice with template parameters
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
const Kratos::Flags ResidualBasedBlockBuilderAndSolverWithMassAndDamping<TSparseSpace, TDenseSpace, TLinearSolver>::SILENT_WARNINGS(Kratos::Flags::Create(0));

///@}

} /* namespace Kratos.*/
