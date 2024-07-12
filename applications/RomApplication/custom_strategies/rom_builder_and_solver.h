//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Riccardo Rossi
//                  Raul Bravo
//

#if !defined(KRATOS_ROM_BUILDER_AND_SOLVER)
#define KRATOS_ROM_BUILDER_AND_SOLVER

/* System includes */

/* External includes */
#include "concurrentqueue/concurrentqueue.h"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "utilities/builtin_timer.h"
#include "utilities/reduction_utilities.h"

/* Application includes */
#include "rom_application_variables.h"
#include "custom_utilities/rom_auxiliary_utilities.h"

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

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ROMBuilderAndSolver : public BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    //TODO: UPDATE THIS
    /**
     * This struct is used in the component wise calculation only
     * is defined here and is used to declare a member variable in the component wise builder and solver
     * private pointers can only be accessed by means of set and get functions
     * this allows to set and not copy the Element_Variables and Condition_Variables
     * which will be asked and set by another strategy object
     */

    ///@name Type Definitions
    ///@{

    // Class pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ROMBuilderAndSolver);

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// The definition of the current class
    typedef ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    /// Definition of the classes from the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /// Additional definitions
    typedef typename ModelPart::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;
    typedef LocalSystemMatrixType RomSystemMatrixType;
    typedef LocalSystemVectorType RomSystemVectorType;

    /// DoF types definition
    typedef Node NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;
    typedef moodycamel::ConcurrentQueue<DofType::Pointer> DofQueue;

    ///@}
    ///@name Life cycle
    ///@{

    explicit ROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters): BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = this->ValidateAndAssignParameters(this_parameters_copy, this->GetDefaultParameters());
        this->AssignSettings(this_parameters_copy);
    }

    explicit ROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
    }

    ~ROMBuilderAndSolver() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(pNewLinearSystemSolver,ThisParameters);
    }

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 1)) << "Setting up the dofs" << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of threads" << ParallelUtilities::GetNumThreads() << "\n" << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        // Get model part data
        if (mHromWeightsInitialized == false) {
            InitializeHROMWeights(rModelPart);
        }

        auto dof_queue = ExtractDofSet(pScheme, rModelPart);

        // Fill a sorted auxiliary array of with the DOFs set
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;
        auto dof_array = SortAndRemoveDuplicateDofs(dof_queue);

        // Update base builder and solver DOFs array and set corresponding flag
        BaseType::GetDofSet().swap(dof_array);
        BaseType::SetDofSetIsInitializedFlag(true);

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(BaseType::GetDofSet().size() == 0) << "No degrees of freedom!" << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::GetDofSet().size() << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished setting up the dofs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        if (BaseType::GetCalculateReactionsFlag())
        {
            for (const auto& r_dof: BaseType::GetDofSet())
            {
                KRATOS_ERROR_IF_NOT(r_dof.HasReaction())
                    << "Reaction variable not set for the following :\n"
                    << "Node : " << r_dof.Id() << '\n'
                    << "Dof  : " << r_dof      << '\n'
                    << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif
        KRATOS_CATCH("");
    }

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        TSparseSpace::SetToZero(rb);

        //refresh RHS to have the correct reactions (with contributions on Dirichlet BCs)
        BuildRHSNoDirichlet(rModelPart, rb);

        //NOTE: dofs are assumed to be numbered consecutively in the BuilderAndSolver
        block_for_each(BaseType::mDofSet, [&](Dof<double>& rDof){
            const std::size_t i = rDof.EquationId();

            rDof.GetSolutionStepReactionValue() = -rb[i];
        });
    }

    void SetUpSystem(ModelPart &rModelPart) override
    {
        auto& r_dof_set = BaseType::GetDofSet();
        BaseType::mEquationSystemSize = r_dof_set.size();
        IndexPartition<IndexType>(r_dof_set.size()).for_each([&](IndexType Index)
        {
            auto dof_iterator = r_dof_set.begin() + Index;
            dof_iterator->SetEquationId(Index);
        });
    }

    // Vector ProjectToReducedBasis(
	// 	const TSystemVectorType& rX,
	// 	ModelPart::NodesContainerType& rNodes
	// )
    // {
    //     Vector rom_unknowns = ZeroVector(GetNumberOfROMModes());
    //     for(const auto& node : rNodes)
    //     {
    //         unsigned int node_aux_id = node.GetValue(AUX_ID);
    //         const auto& nodal_rom_basis = node.GetValue(ROM_BASIS);
	// 			for (int i = 0; i < GetNumberOfROMModes(); ++i) {
	// 				for (int j = 0; j < mNodalDofs; ++j) {
	// 					rom_unknowns[i] += nodal_rom_basis(j, i)*rX(node_aux_id*mNodalDofs + j);
	// 				}
	// 			}
    //     }
    //     return rom_unknowns;
	// }

    SizeType GetNumberOfROMModes() const noexcept
    {
        return mNumberOfRomModes;
    } 

    void ProjectToFineBasis(
        const TSystemVectorType& rRomUnkowns,
        const ModelPart& rModelPart,
        TSystemVectorType& rDx) const
    {
        const auto& r_dof_set = BaseType::GetDofSet();
        block_for_each(r_dof_set, [&](const DofType& r_dof)
        {
            const auto& r_node = rModelPart.GetNode(r_dof.Id());
            const Matrix& r_rom_nodal_basis = r_node.GetValue(ROM_BASIS);
            const Matrix::size_type row_id = mMapPhi.at(r_dof.GetVariable().Key());
            rDx[r_dof.EquationId()] = inner_prod(row(r_rom_nodal_basis, row_id), rRomUnkowns);
        });
    }

    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        // Call the base B&S InitializeSolutionStep
        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Reset the ROM solution increment in the root modelpart database
        auto& r_root_mp = rModelPart.GetRootModelPart();
        r_root_mp.GetValue(ROM_CURRENT_SOLUTION_TOTAL) = ZeroVector(GetNumberOfROMModes());
    }

    
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        RomSystemMatrixType Arom = ZeroMatrix(GetNumberOfROMModes(), GetNumberOfROMModes());
        RomSystemVectorType brom = ZeroVector(GetNumberOfROMModes());

        BuildROM(pScheme, rModelPart, Arom, brom);
        SolveROM(rModelPart, Arom, brom, Dx);

        KRATOS_CATCH("")
    }

    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType &pA,
        TSystemVectorPointerType &pDx,
        TSystemVectorPointerType &pb,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY

        // If not initialized, initalize the system arrays to an empty vector/matrix
        if (!pA) {
            TSystemMatrixPointerType p_new_A = Kratos::make_shared<TSystemMatrixType>(0, 0);
            pA.swap(p_new_A);
        }
        if (!pDx) {
            TSystemVectorPointerType p_new_Dx = Kratos::make_shared<TSystemVectorType>(0);
            pDx.swap(p_new_Dx);
        }
        if (!pb) {
            TSystemVectorPointerType p_new_b = Kratos::make_shared<TSystemVectorType>(0);
            pb.swap(p_new_b);
        }

        TSystemVectorType& r_Dx = *pDx;
        if (r_Dx.size() != BaseType::GetEquationSystemSize()) {
            r_Dx.resize(BaseType::GetEquationSystemSize(), false);
        }

        TSystemVectorType& r_b = *pb;
        if (r_b.size() != BaseType::GetEquationSystemSize()) {
            r_b.resize(BaseType::GetEquationSystemSize(), false);
        }

        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "rom_bns_settings" : {}
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    static std::string Name()
    {
        return "rom_builder_and_solver";
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
    virtual std::string Info() const override
    {
        return "ROMBuilderAndSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@}
    ///@name Protected static member variables
    ///@{


    ///@}
    ///@name Protected member variables
    ///@{

    SizeType mNodalDofs;
    std::unordered_map<Kratos::VariableData::KeyType, Matrix::size_type> mMapPhi;

    ElementsArrayType mSelectedElements;
    ConditionsArrayType mSelectedConditions;

    bool mHromSimulation = false;
    bool mHromWeightsInitialized = false;

    ///@}
    ///@name Protected operators
    ///@{
    
    void BuildRHSNoDirichlet(
        ModelPart& rModelPart,
        TSystemVectorType& rb)
    {
        KRATOS_TRY

        // Get ProcessInfo from main model part
        const auto& r_current_process_info = rModelPart.GetProcessInfo();

        auto& r_elements = mHromSimulation ? mSelectedElements : rModelPart.Elements();
        if(!r_elements.empty())
        {
            block_for_each(r_elements, Kratos::Vector(), [&](Element& r_element, Kratos::Vector& r_rhs_elem)
            {
                DofsVectorType dofs;

                r_element.CalculateRightHandSide(r_rhs_elem, r_current_process_info);
                r_element.GetDofList(dofs, r_current_process_info);
                for (IndexType i = 0; i < dofs.size(); ++i){
                    double& r_bi = rb[dofs[i]->EquationId()];
                    AtomicAdd(r_bi, r_rhs_elem[i]); // Building RHS.
                }
            });
        }

        auto& r_conditions = mHromSimulation ? mSelectedConditions : rModelPart.Conditions();
        if(!r_conditions.empty())
        {
            block_for_each(r_conditions, Kratos::Vector(), [&](Condition& r_condition, Kratos::Vector& r_rhs_cond)
            {
                DofsVectorType dofs = {};
                r_condition.CalculateRightHandSide(r_rhs_cond, r_current_process_info);
                r_condition.GetDofList(dofs, r_current_process_info);
                for (IndexType i = 0; i < dofs.size(); ++i){
                    double& r_bi = rb[dofs[i]->EquationId()];
                    AtomicAdd(r_bi, r_rhs_cond[i]); // Building RHS.
                }
            });
        }   

        KRATOS_CATCH("")

    }

    ///@}
    ///@name Protected operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // Set member variables
        mNodalDofs = ThisParameters["nodal_unknowns"].size();
        mNumberOfRomModes = ThisParameters["number_of_rom_dofs"].GetInt();

        // Set up a map with key the variable key and value the correct row in ROM basis
        IndexType k = 0;
        for (const auto& r_var_name : ThisParameters["nodal_unknowns"].GetStringArray()) {
            if(KratosComponents<Variable<double>>::Has(r_var_name)) {
                const auto& var = KratosComponents<Variable<double>>::Get(r_var_name);
                mMapPhi[var.Key()] = k++;
            } else {
                KRATOS_ERROR << "Variable \""<< r_var_name << "\" not valid" << std::endl;
            }
        }
    }


    void InitializeHROMWeights(ModelPart& rModelPart)
    {
        KRATOS_TRY

        using ElementQueue = moodycamel::ConcurrentQueue<Element::Pointer>;
        using ConditionQueue = moodycamel::ConcurrentQueue<Condition::Pointer>;

        // Inspecting elements
        ElementQueue element_queue;
        block_for_each(rModelPart.Elements().GetContainer(),
            [&](Element::Pointer p_element)
        {
            if (p_element->Has(HROM_WEIGHT)) {
                element_queue.enqueue(std::move(p_element));
            } else {
                p_element->SetValue(HROM_WEIGHT, 1.0);
            }
        });

        // Inspecting conditions
        ConditionQueue condition_queue;
        block_for_each(rModelPart.Conditions().GetContainer(),
            [&](Condition::Pointer p_condition)
        {
            if (p_condition->Has(HROM_WEIGHT)) {
                condition_queue.enqueue(std::move(p_condition));
            } else {
                p_condition->SetValue(HROM_WEIGHT, 1.0);
            }
        });

        // Dequeueing elements
        std::size_t err_id;
        mSelectedElements.reserve(element_queue.size_approx());
        Element::Pointer p_element;
        while ( (err_id = element_queue.try_dequeue(p_element)) != 0) {
            mSelectedElements.push_back(std::move(p_element));
        }

        // Dequeueing conditions
        mSelectedConditions.reserve(condition_queue.size_approx());
        Condition::Pointer p_condition;
        while ( (err_id = condition_queue.try_dequeue(p_condition)) != 0) {
            mSelectedConditions.push_back(std::move(p_condition));
        }

        // Wrap-up
        mHromSimulation = !(mSelectedElements.empty() && mSelectedConditions.empty());
        mHromWeightsInitialized = true;

        KRATOS_CATCH("")
    }
    
    static DofQueue ExtractDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart)
    {
        KRATOS_TRY

        DofQueue dof_queue;

        // Emulates ConcurrentQueue::enqueue_bulk adding move semantics to avoid atomic ops
        const auto enqueue_bulk_move = [](DofQueue& r_queue, DofsVectorType& r_dof_list) {
            for(auto& p_dof: r_dof_list) {
                r_queue.enqueue(std::move(p_dof));
            }
            r_dof_list.clear();
        };

        // Inspecting elements
        DofsVectorType tls_dof_list; // Preallocation
        block_for_each(rModelPart.Elements(), tls_dof_list,
            [&](const Element& r_element, DofsVectorType& r_dof_list)
        {
            pScheme->GetDofList(r_element, r_dof_list, rModelPart.GetProcessInfo());
            enqueue_bulk_move(dof_queue, r_dof_list);
        });

        // Inspecting conditions
        block_for_each(rModelPart.Conditions(), tls_dof_list,
            [&](const Condition& r_condition, DofsVectorType& r_dof_list)
        {
            pScheme->GetDofList(r_condition, r_dof_list, rModelPart.GetProcessInfo());
            enqueue_bulk_move(dof_queue, r_dof_list);
        });

        // Inspecting master-slave constraints
        std::pair<DofsVectorType, DofsVectorType> tls_ms_dof_lists; // Preallocation
        block_for_each(rModelPart.MasterSlaveConstraints(), tls_ms_dof_lists,
            [&](const MasterSlaveConstraint& r_constraint, std::pair<DofsVectorType, DofsVectorType>& r_dof_lists)
        {
            r_constraint.GetDofList(r_dof_lists.first, r_dof_lists.second, rModelPart.GetProcessInfo());

            enqueue_bulk_move(dof_queue, r_dof_lists.first);
            enqueue_bulk_move(dof_queue, r_dof_lists.second);
        });

        return dof_queue;

        KRATOS_CATCH("")
    }

    static DofsArrayType SortAndRemoveDuplicateDofs(DofQueue& rDofQueue)
    {
        KRATOS_TRY

        DofsArrayType dof_array;
        dof_array.reserve(rDofQueue.size_approx());
        DofType::Pointer p_dof;
        std::size_t err_id;
        while ( (err_id = rDofQueue.try_dequeue(p_dof)) != 0) {
            dof_array.push_back(std::move(p_dof));
        }

        dof_array.Unique(); // Sorts internally

        return dof_array;

        KRATOS_CATCH("")
    }

    /**
    * Thread Local Storage containing dynamically allocated structures to avoid reallocating each iteration.
    */
    struct AssemblyTLS
    {
        AssemblyTLS(SizeType NRomModes)
            : romA(ZeroMatrix(NRomModes, NRomModes)),
              romB(ZeroVector(NRomModes))
        { }
        AssemblyTLS() = delete;

        Matrix phiE = {};                // Elemental Phi
        LocalSystemMatrixType lhs = {};  // Elemental LHS
        LocalSystemVectorType rhs = {};  // Elemental RHS
        EquationIdVectorType eq_id = {}; // Elemental equation ID vector
        DofsVectorType dofs = {};        // Elemental dof vector
        RomSystemMatrixType romA;        // reduced LHS
        RomSystemVectorType romB;        // reduced RHS
        RomSystemMatrixType aux = {};    // Auxiliary: romA = phi.t * (LHS * phi) := phi.t * aux
    };

    
    /**
     * Class to sum-reduce matrices and vectors.
     */
    template<typename T>
    struct NonTrivialSumReduction
    {
        typedef T value_type;
        typedef T return_type;

        T mValue;
        bool mInitialized = false;

        void Init(const value_type& first_value)
        {
            mValue = first_value;
            mInitialized = true;
        }

        /// access to reduced value
        return_type GetValue() const
        {
            return mValue;
        }

        void LocalReduce(const value_type& value)
        {
            if(!mInitialized) {
                Init(value);
            } else {
                noalias(mValue) += value;
            }
        }

        void ThreadSafeReduce(const NonTrivialSumReduction& rOther)
        {
            if(!rOther.mInitialized) return;

            const std::lock_guard<LockObject> scope_lock(ParallelUtilities::GetGlobalLock());
            LocalReduce(rOther.mValue);
        }
    };

    /**
     * Resizes a Matrix if it's not the right size
     */
    template<typename TMatrix>
    static void ResizeIfNeeded(TMatrix& mat, const SizeType rows, const SizeType cols)
    {
        if(mat.size1() != rows || mat.size2() != cols) {
            mat.resize(rows, cols, false);
        }
    }

    /**
     * Builds the reduced system of equations on rank 0 
     */
    virtual void BuildROM(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        RomSystemMatrixType &rA,
        RomSystemVectorType &rb)
    {
        KRATOS_TRY

        // Define a dense matrix to hold the reduced problem
        rA = ZeroMatrix(GetNumberOfROMModes(), GetNumberOfROMModes());
        rb = ZeroVector(GetNumberOfROMModes());

        // Build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Get ProcessInfo from main model part
        const auto& r_current_process_info = rModelPart.GetProcessInfo();


        // Assemble all entities
        const auto assembling_timer = BuiltinTimer();

        using SystemSumReducer = CombinedReduction<NonTrivialSumReduction<RomSystemMatrixType>, NonTrivialSumReduction<RomSystemVectorType>>;
        AssemblyTLS assembly_tls_container(GetNumberOfROMModes());

        auto& elements = mHromSimulation ? mSelectedElements : rModelPart.Elements();
        if(!elements.empty())
        {
            std::tie(rA, rb) =
            block_for_each<SystemSumReducer>(elements, assembly_tls_container, 
                [&](Element& r_element, AssemblyTLS& r_thread_prealloc)
            {
                return CalculateLocalContribution(r_element, r_thread_prealloc, *pScheme, r_current_process_info);
            });
        }

        auto& conditions = mHromSimulation ? mSelectedConditions : rModelPart.Conditions();
        if(!conditions.empty())
        {
            RomSystemMatrixType Aconditions;
            RomSystemVectorType bconditions;

            std::tie(Aconditions, bconditions) =
            block_for_each<SystemSumReducer>(conditions, assembly_tls_container, 
                [&](Condition& r_condition, AssemblyTLS& r_thread_prealloc)
            {
                return CalculateLocalContribution(r_condition, r_thread_prealloc, *pScheme, r_current_process_info);
            });

            rA += Aconditions;
            rb += bconditions;
        }

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Build time: " << assembling_timer.ElapsedSeconds() << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Solves reduced system of equations and broadcasts it
     */
    virtual void SolveROM(
        ModelPart &rModelPart,
        RomSystemMatrixType &rA,
        RomSystemVectorType &rb,
        TSystemVectorType &rDx)
    {
        KRATOS_TRY

        RomSystemVectorType dxrom(GetNumberOfROMModes());
        
        const auto solving_timer = BuiltinTimer();
        MathUtils<double>::Solve(rA, dxrom, rb);
        // KRATOS_WATCH(dxrom)
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Solve reduced system time: " << solving_timer.ElapsedSeconds() << std::endl;

        // Save the ROM solution increment in the root modelpart database
        auto& r_root_mp = rModelPart.GetRootModelPart();
        noalias(r_root_mp.GetValue(ROM_CURRENT_SOLUTION_TOTAL)) += dxrom;

        // project reduced solution back to full order model
        const auto backward_projection_timer = BuiltinTimer();
        ProjectToFineBasis(dxrom, rModelPart, rDx);
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Project to fine basis time: " << backward_projection_timer.ElapsedSeconds() << std::endl;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Protected access
    ///@{


    ///@}
    ///@name Protected inquiry
    ///@{


    ///@}
    ///@name Protected life cycle
    ///@{

private:

    SizeType mNumberOfRomModes;

    ///@}
    ///@name Private operations 
    ///@{

    /**
     * Computes the local contribution of an element or condition
     */
    template<typename TEntity>
    std::tuple<LocalSystemMatrixType, LocalSystemVectorType> CalculateLocalContribution(
        TEntity& rEntity,
        AssemblyTLS& rPreAlloc,
        TSchemeType& rScheme,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE))
        {
            rPreAlloc.romA = ZeroMatrix(GetNumberOfROMModes(), GetNumberOfROMModes());
            rPreAlloc.romB = ZeroVector(GetNumberOfROMModes());
            return std::tie(rPreAlloc.romA, rPreAlloc.romB);
        }

        // Calculate elemental contribution
        rScheme.CalculateSystemContributions(rEntity, rPreAlloc.lhs, rPreAlloc.rhs, rPreAlloc.eq_id, rCurrentProcessInfo);
        rEntity.GetDofList(rPreAlloc.dofs, rCurrentProcessInfo);

        const SizeType ndofs = rPreAlloc.dofs.size();
        ResizeIfNeeded(rPreAlloc.phiE, ndofs, GetNumberOfROMModes());
        ResizeIfNeeded(rPreAlloc.aux, ndofs, GetNumberOfROMModes());

        const auto &r_geom = rEntity.GetGeometry();
        RomAuxiliaryUtilities::GetPhiElemental(rPreAlloc.phiE, rPreAlloc.dofs, r_geom, mMapPhi);

        const double h_rom_weight = mHromSimulation ? rEntity.GetValue(HROM_WEIGHT) : 1.0;

        noalias(rPreAlloc.aux) = prod(rPreAlloc.lhs, rPreAlloc.phiE);
        noalias(rPreAlloc.romA) = prod(trans(rPreAlloc.phiE), rPreAlloc.aux) * h_rom_weight;
        noalias(rPreAlloc.romB) = prod(trans(rPreAlloc.phiE), rPreAlloc.rhs) * h_rom_weight;

        return std::tie(rPreAlloc.romA, rPreAlloc.romB);
    }


    ///@}
}; /* Class ROMBuilderAndSolver */

///@}
///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_ROM_BUILDER_AND_SOLVER  defined */