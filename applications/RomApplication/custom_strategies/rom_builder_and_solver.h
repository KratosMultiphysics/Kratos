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
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;
    typedef typename std::unordered_set<DofPointerType, DofPointerHasher> DofSetType;

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
        const auto& r_elements_array = rModelPart.Elements();
        const auto& r_conditions_array = rModelPart.Conditions();
        const auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
        const int number_of_elements = static_cast<int>(r_elements_array.size());
        const int number_of_conditions = static_cast<int>(r_conditions_array.size());
        const int number_of_constraints = static_cast<int>(r_constraints_array.size());
        const auto& r_current_process_info = rModelPart.GetProcessInfo();

        DofsVectorType dof_list;
        DofsVectorType second_dof_list; // NOTE: The second dof list is only used on constraints to include master/slave relations

        DofSetType dof_global_set;
        dof_global_set.reserve(number_of_elements * 20);

        if (mHromWeightsInitialized == false){
            int number_of_hrom_entities = 0;
            #pragma omp parallel firstprivate(dof_list, second_dof_list) reduction(+:number_of_hrom_entities)
            {
                // We create the temporal set and we reserve some space on them
                DofSetType dofs_tmp_set;
                dofs_tmp_set.reserve(20000);

                // Loop the array of elements
                ElementsArrayType selected_elements_private;
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_elements; ++i) {
                    auto it_elem = r_elements_array.begin() + i;

                    // Detect whether the element has an hyperreduced weight (H-ROM simulation) or not (ROM simulation)
                    if ((it_elem)->Has(HROM_WEIGHT)){
                        selected_elements_private.push_back(*it_elem.base());
                        number_of_hrom_entities++;
                    } else {
                        it_elem->SetValue(HROM_WEIGHT, 1.0);
                    }

                    // Gets list of DOF involved on every element
                    pScheme->GetDofList(*it_elem, dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                }

                // Loop the array of conditions
                ConditionsArrayType selected_conditions_private;
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_conditions; ++i) {
                    auto it_cond = r_conditions_array.begin() + i;

                    // Detect whether the condition has an hyperreduced weight (H-ROM simulation) or not (ROM simulation)
                    // Note that those conditions used for displaying results are to be ignored as they will not be assembled
                    if (it_cond->Has(HROM_WEIGHT)){
                        selected_conditions_private.push_back(*it_cond.base());
                        number_of_hrom_entities++;
                    } else {
                        it_cond->SetValue(HROM_WEIGHT, 1.0);
                    }

                    // Gets list of DOF involved on every condition
                    pScheme->GetDofList(*it_cond, dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                }

                // Loop the array of constraints
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_constraints; ++i) {
                    auto it_const = r_constraints_array.begin() + i;

                    // Gets list of Dof involved on every constraint
                    it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                    dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
                }

                #pragma omp critical
                {
                    // Collect the elements and conditions belonging to the H-ROM mesh
                    // These are those that feature a weight and are to be assembled
                    for (auto &r_cond : selected_conditions_private){
                        mSelectedConditions.push_back(&r_cond);
                    }
                    for (auto &r_elem : selected_elements_private){
                        mSelectedElements.push_back(&r_elem);
                    }

                    // We merge all the sets in one thread
                    dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
                }
            }

            // Update H-ROM flags
            if (number_of_hrom_entities) {
                mHromSimulation = true;
            }
            mHromWeightsInitialized = true;

        } else {
            #pragma omp parallel firstprivate(dof_list, second_dof_list)
            {
                // We create the temporal set and we reserve some space on them
                DofSetType dofs_tmp_set;
                dofs_tmp_set.reserve(20000);

                // Loop the array of elements
                ElementsArrayType selected_elements_private;
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_elements; ++i) {
                    auto it_elem = r_elements_array.begin() + i;

                    // Gets list of DOF involved on every element
                    pScheme->GetDofList(*it_elem, dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                }

                // Loop the array of conditions
                ConditionsArrayType selected_conditions_private;
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_conditions; ++i) {
                    auto it_cond = r_conditions_array.begin() + i;

                    // Gets list of DOF involved on every condition
                    pScheme->GetDofList(*it_cond, dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                }

                // Loop the array of constraints
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_constraints; ++i) {
                    auto it_const = r_constraints_array.begin() + i;

                    // Gets list of Dof involved on every constraint
                    it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                    dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
                }

                #pragma omp critical
                {
                    // We merge all the sets in one thread
                    dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
                }
            }
        }

        // Fill a sorted auxiliary array of with the DOFs set
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;
        DofsArrayType Doftemp;
        Doftemp.reserve(dof_global_set.size());
        for (auto it = dof_global_set.begin(); it != dof_global_set.end(); it++) {
            Doftemp.push_back(*it);
        }
        Doftemp.Sort();

        // Update base builder and solver DOFs array and set corresponding flag
        BaseType::GetDofSet() = Doftemp;
        BaseType::SetDofSetIsInitializedFlag(true);

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(BaseType::GetDofSet().size() == 0) << "No degrees of freedom!" << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        if (BaseType::GetCalculateReactionsFlag()) {
            for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator) {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction())
                    << "Reaction variable not set for the following : " << std::endl
                    << "Node : " << dof_iterator->Id() << std::endl
                    << "Dof : " << (*dof_iterator) << std::endl
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
        ModelPart& rModelPart,
        TSystemVectorType& rDx)
    {
        auto& r_dof_set = BaseType::GetDofSet();
        const int dofs_number = r_dof_set.size();
        const auto dofs_begin = r_dof_set.begin();

        #pragma omp parallel firstprivate(dofs_begin, dofs_number)
        {
            const Matrix *pcurrent_rom_nodal_basis = nullptr;
            unsigned int old_dof_id;
            #pragma omp for nowait
            for (int k = 0; k < dofs_number; k++) {
                auto dof = dofs_begin + k;
                if (pcurrent_rom_nodal_basis == nullptr) {
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS));
                    old_dof_id = dof->Id();
                } else if (dof->Id() != old_dof_id ) {
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS));
                    old_dof_id = dof->Id();
                }
                rDx[dof->EquationId()] = inner_prod(row(*pcurrent_rom_nodal_basis, mMapPhi[dof->GetVariable().Key()]), rRomUnkowns);
            }
        }
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

        // Getting the elements from the model
        // Only selected conditions and elements are used for the calculation in an H-ROM simulation
        const auto el_begin = mHromSimulation ? mSelectedElements.begin() : rModelPart.ElementsBegin();
        const int nelements = mHromSimulation ? mSelectedElements.size() : rModelPart.NumberOfElements();
        const auto cond_begin = mHromSimulation ? mSelectedConditions.begin() : rModelPart.ConditionsBegin();
        const int nconditions = mHromSimulation ? mSelectedConditions.size() : rModelPart.NumberOfConditions();

        // Get ProcessInfo from main model part
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        // Assemble all entities
        const auto assembling_timer = BuiltinTimer();
        #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
        {
            Matrix PhiElemental;
            Matrix tempA = ZeroMatrix(mNumberOfRomModes,mNumberOfRomModes);
            Vector tempb = ZeroVector(mNumberOfRomModes);
            Matrix aux;

            #pragma omp for nowait
            for (int k = 0; k < static_cast<int>(nelements); k++) {
                auto it_el = el_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool element_is_active = true;
                if ((it_el)->IsDefined(ACTIVE)) {
                    element_is_active = (it_el)->Is(ACTIVE);
                }

                // Calculate elemental contribution
                if (element_is_active){
                    pScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, r_current_process_info);
                    Element::DofsVectorType dofs;
                    it_el->GetDofList(dofs, r_current_process_info);
                    const auto &geom = it_el->GetGeometry();
                    if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mNumberOfRomModes) {
                        PhiElemental.resize(dofs.size(), mNumberOfRomModes,false);
                    }
                    if(aux.size1() != dofs.size() || aux.size2() != mNumberOfRomModes) {
                        aux.resize(dofs.size(), mNumberOfRomModes,false);
                    }
                    RomAuxiliaryUtilities::GetPhiElemental(PhiElemental, dofs, geom, mMapPhi);
                    noalias(aux) = prod(LHS_Contribution, PhiElemental);
                    double h_rom_weight = it_el->GetValue(HROM_WEIGHT);
                    noalias(tempA) += prod(trans(PhiElemental), aux) * h_rom_weight;
                    noalias(tempb) += prod(trans(PhiElemental), RHS_Contribution) * h_rom_weight;
                }
            }

            #pragma omp for nowait
            for (int k = 0; k < static_cast<int>(nconditions); k++){
                auto it = cond_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the condition is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE)) {
                    condition_is_active = (it)->Is(ACTIVE);
                }

                // Calculate condition contribution
                if (condition_is_active) {
                    Condition::DofsVectorType dofs;
                    it->GetDofList(dofs, r_current_process_info);
                    pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, r_current_process_info);
                    const auto &geom = it->GetGeometry();
                    if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mNumberOfRomModes) {
                        PhiElemental.resize(dofs.size(), mNumberOfRomModes,false);
                    }
                    if(aux.size1() != dofs.size() || aux.size2() != mNumberOfRomModes) {
                        aux.resize(dofs.size(), mNumberOfRomModes,false);
                    }
                    RomAuxiliaryUtilities::GetPhiElemental(PhiElemental, dofs, geom, mMapPhi);
                    noalias(aux) = prod(LHS_Contribution, PhiElemental);
                    double h_rom_weight = it->GetValue(HROM_WEIGHT);
                    noalias(tempA) += prod(trans(PhiElemental), aux) * h_rom_weight;
                    noalias(tempb) += prod(trans(PhiElemental), RHS_Contribution) * h_rom_weight;
                }
            }

            #pragma omp critical
            {
                noalias(Arom) += tempA;
                noalias(brom) += tempb;
            }

        }

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
