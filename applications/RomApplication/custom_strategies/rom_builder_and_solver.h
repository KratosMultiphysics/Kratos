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

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;
    typedef moodycamel::ConcurrentQueue<DofType::Pointer> DofQueue;

    ///@}
    ///@name Life cycle
    ///@{

    explicit ROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters)
        : BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = this->ValidateAndAssignParameters(this_parameters_copy, this->GetDefaultParameters());
        this->AssignSettings(this_parameters_copy);
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

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the dofs" << std::endl;
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
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        if (BaseType::GetCalculateReactionsFlag())
        {
            for (const auto& r_dof: BaseType::mDofSet)
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

    void SetUpSystem(ModelPart &rModelPart) override
    {
        auto& r_dof_set = BaseType::GetDofSet();
        BaseType::mEquationSystemSize = r_dof_set.size();
        int ndofs = static_cast<int>(r_dof_set.size());

        #pragma omp parallel for firstprivate(ndofs)
        for (int i = 0; i < ndofs; i++){
            auto dof_iterator = r_dof_set.begin() + i;
            dof_iterator->SetEquationId(i);
        }
    }

    // Vector ProjectToReducedBasis(
	// 	const TSystemVectorType& rX,
	// 	ModelPart::NodesContainerType& rNodes
	// )
    // {
    //     Vector rom_unknowns = ZeroVector(mNumberOfRomModes);
    //     for(const auto& node : rNodes)
    //     {
    //         unsigned int node_aux_id = node.GetValue(AUX_ID);
    //         const auto& nodal_rom_basis = node.GetValue(ROM_BASIS);
	// 			for (int i = 0; i < mNumberOfRomModes; ++i) {
	// 				for (int j = 0; j < mNodalDofs; ++j) {
	// 					rom_unknowns[i] += nodal_rom_basis(j, i)*rX(node_aux_id*mNodalDofs + j);
	// 				}
	// 			}
    //     }
    //     return rom_unknowns;
	// }

    void ProjectToFineBasis(
        const TSystemVectorType& rRomUnkowns,
        const ModelPart& rModelPart,
        TSystemVectorType& rDx) const
    {
        const auto& r_dof_set = BaseType::mDofSet;
        block_for_each(r_dof_set, [&](const DofType& r_dof)
        {
            const NodeType& node = rModelPart.GetNode(r_dof.Id());
            const Matrix& rom_nodal_basis = node.GetValue(ROM_BASIS);
            const Matrix::size_type row_id = mMapPhi.at(r_dof.GetVariable().Key());
            rDx[r_dof.EquationId()] = inner_prod(row(rom_nodal_basis, row_id), rRomUnkowns);
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
        r_root_mp.GetValue(ROM_SOLUTION_INCREMENT) = ZeroVector(mNumberOfRomModes);
    }

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        // Define a dense matrix to hold the reduced problem
        Matrix Arom = ZeroMatrix(mNumberOfRomModes, mNumberOfRomModes);
        Vector brom = ZeroVector(mNumberOfRomModes);
        // TSystemVectorType x(Dx.size());

        // const auto forward_projection_timer = BuiltinTimer();
        // Vector xrom = ZeroVector(mNumberOfRomModes);
        //this->ProjectToReducedBasis(x, rModelPart.Nodes(),xrom);
        // KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)) << "Project to reduced basis time: " << forward_projection_timer.ElapsedSeconds() << std::endl;

        // Build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Get ProcessInfo from main model part
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Assemble all entities
        const auto assembling_timer = BuiltinTimer();

        // Preallocating
        Element::EquationIdVectorType EquationId;
        Matrix PhiElemental;
        Matrix tempA = ZeroMatrix(mNumberOfRomModes,mNumberOfRomModes);
        Vector tempb = ZeroVector(mNumberOfRomModes);
        Matrix aux;

        auto& elements = mHromSimulation ? mSelectedElements : rModelPart.Elements();
        for (auto& r_element: elements)
        {
            // Detect if the element is active or not. If the user did not make any choice the element is active by default
            bool element_is_active = true;
            if (r_element.IsDefined(ACTIVE)) {
                element_is_active = r_element.Is(ACTIVE);
            }

            // Calculate elemental contribution
            if (element_is_active){
                pScheme->CalculateSystemContributions(r_element, LHS_Contribution, RHS_Contribution, EquationId, r_current_process_info);
                Element::DofsVectorType dofs;
                r_element.GetDofList(dofs, r_current_process_info);
                const auto &geom = r_element.GetGeometry();
                if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mNumberOfRomModes) {
                    PhiElemental.resize(dofs.size(), mNumberOfRomModes,false);
                }
                if(aux.size1() != dofs.size() || aux.size2() != mNumberOfRomModes) {
                    aux.resize(dofs.size(), mNumberOfRomModes,false);
                }
                RomAuxiliaryUtilities::GetPhiElemental(PhiElemental, dofs, geom, mMapPhi);
                noalias(aux) = prod(LHS_Contribution, PhiElemental);
                double h_rom_weight = r_element.GetValue(HROM_WEIGHT);
                noalias(tempA) += prod(trans(PhiElemental), aux) * h_rom_weight;
                noalias(tempb) += prod(trans(PhiElemental), RHS_Contribution) * h_rom_weight;
            }
        }

        auto& conditions = mHromSimulation ? mSelectedConditions : rModelPart.Conditions();
        for (auto& r_condition: conditions){

            // Detect if the element is active or not. If the user did not make any choice the condition is active by default
            bool condition_is_active = true;
            if (r_condition.IsDefined(ACTIVE)) {
                condition_is_active = r_condition.Is(ACTIVE);
            }

            // Calculate condition contribution
            if (condition_is_active) {
                Condition::DofsVectorType dofs;
                r_condition.GetDofList(dofs, r_current_process_info);
                pScheme->CalculateSystemContributions(r_condition, LHS_Contribution, RHS_Contribution, EquationId, r_current_process_info);
                const auto &geom = r_condition.GetGeometry();
                if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mNumberOfRomModes) {
                    PhiElemental.resize(dofs.size(), mNumberOfRomModes,false);
                }
                if(aux.size1() != dofs.size() || aux.size2() != mNumberOfRomModes) {
                    aux.resize(dofs.size(), mNumberOfRomModes,false);
                }
                RomAuxiliaryUtilities::GetPhiElemental(PhiElemental, dofs, geom, mMapPhi);
                noalias(aux) = prod(LHS_Contribution, PhiElemental);
                double h_rom_weight = r_condition.GetValue(HROM_WEIGHT);
                noalias(tempA) += prod(trans(PhiElemental), aux) * h_rom_weight;
                noalias(tempb) += prod(trans(PhiElemental), RHS_Contribution) * h_rom_weight;
            }
        }

        noalias(Arom) += tempA;
        noalias(brom) += tempb;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << assembling_timer.ElapsedSeconds() << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        //solve for the rom unkowns dunk = Arom^-1 * brom
        Vector dxrom(mNumberOfRomModes);
        const auto solving_timer = BuiltinTimer();
        MathUtils<double>::Solve(Arom, dxrom, brom);
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)) << "Solve reduced system time: " << solving_timer.ElapsedSeconds() << std::endl;

        // Save the ROM solution increment in the root modelpart database
        // This can be used later on to recover the solution in a visualization submodelpart
        auto& r_root_mp = rModelPart.GetRootModelPart();
        noalias(r_root_mp.GetValue(ROM_SOLUTION_INCREMENT)) += dxrom;

        // //update database
        // noalias(xrom) += dxrom;

        // project reduced solution back to full order model
        const auto backward_projection_timer = BuiltinTimer();
        ProjectToFineBasis(dxrom, rModelPart, Dx);
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)) << "Project to fine basis time: " << backward_projection_timer.ElapsedSeconds() << std::endl;
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
            "number_of_rom_dofs" : 10
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
    SizeType mNumberOfRomModes;
    std::unordered_map<Kratos::VariableData::KeyType, Matrix::size_type> mMapPhi;

    ElementsArrayType mSelectedElements;
    ConditionsArrayType mSelectedConditions;

    bool mHromSimulation = false;
    bool mHromWeightsInitialized = false;

    ///@}
    ///@name Protected operators
    ///@{


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
        DofsVectorType tmp_dof_list; // Preallocation
        block_for_each(rModelPart.Elements().GetContainer(), tmp_dof_list,
            [&](Element::Pointer p_element, DofsVectorType& dof_list)
        {
            if (p_element->Has(HROM_WEIGHT)) {
                element_queue.enqueue(std::move(p_element));
            } else {
                p_element->SetValue(HROM_WEIGHT, 1.0);
            }
        });

        // Inspecting conditions
        ConditionQueue condition_queue;
        block_for_each(rModelPart.Conditions().GetContainer(), tmp_dof_list,
            [&](Condition::Pointer p_condition, DofsVectorType& dof_list)
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
        const auto enqueue_bulk_move = [](DofQueue& queue, DofsVectorType& dof_list) {
            for(auto& p_dof: dof_list) {
                queue.enqueue(std::move(p_dof));
            }
            dof_list.clear();
        };

        // Inspecting elements
        DofsVectorType tmp_dof_list; // Preallocation
        block_for_each(rModelPart.Elements().GetContainer(), tmp_dof_list,
            [&](Element::Pointer p_element, DofsVectorType& dof_list)
        {
            pScheme->GetDofList(*p_element, dof_list, rModelPart.GetProcessInfo());
            enqueue_bulk_move(dof_queue, dof_list);
        });

        // Inspecting conditions
        block_for_each(rModelPart.Conditions().GetContainer(), tmp_dof_list,
            [&](Condition::Pointer p_condition, DofsVectorType& dof_list)
        {
            pScheme->GetDofList(*p_condition, dof_list, rModelPart.GetProcessInfo());
            enqueue_bulk_move(dof_queue, dof_list);
        });

        // Inspecting master-slave constraints
        std::pair<DofsVectorType, DofsVectorType> tmp_ms_dof_lists; // Preallocation
        block_for_each(rModelPart.MasterSlaveConstraints(), tmp_ms_dof_lists,
            [&](MasterSlaveConstraint& r_constraint, std::pair<DofsVectorType, DofsVectorType>& dof_lists)
        {
            r_constraint.GetDofList(dof_lists.first, dof_lists.second, rModelPart.GetProcessInfo());

            enqueue_bulk_move(dof_queue, dof_lists.first);
            enqueue_bulk_move(dof_queue, dof_lists.second);
        });

        return dof_queue;

        KRATOS_CATCH("")
    }

    static DofsArrayType SortAndRemoveDuplicateDofs(DofQueue& dof_queue)
    {
        KRATOS_TRY

        DofsArrayType dof_array;
        dof_array.reserve(dof_queue.size_approx());
        DofType::Pointer p_dof;
        std::size_t err_id;
        while ( (err_id = dof_queue.try_dequeue(p_dof)) != 0) {
            dof_array.push_back(std::move(p_dof));
        }

        dof_array.Unique(); // Sorts internally

        return dof_array;

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


    ///@}
}; /* Class ROMBuilderAndSolver */

///@}
///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_ROM_BUILDER_AND_SOLVER  defined */
