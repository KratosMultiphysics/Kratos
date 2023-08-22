//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Sebastian Ares de Parga Regalado
//

#pragma once

/* System includes */

/* External includes */
#include "concurrentqueue/concurrentqueue.h"
#include <Eigen/Core>
#include <Eigen/Dense>

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/global_rom_builder_and_solver.h"
#include "utilities/builtin_timer.h"
#include "utilities/reduction_utilities.h"
#include "utilities/dense_householder_qr_decomposition.h"
#include "processes/find_nodal_neighbours_process.h"

/* Application includes */
#include "rom_application_variables.h"
#include "custom_utilities/rom_auxiliary_utilities.h"
#include "custom_utilities/rom_residuals_utility.h"

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
 * @class GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver
 * @ingroup RomApplication
 * @brief Current class provides an implementation for LeastSquaresPetrovGalerkinROM builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual).
 * The LHS is constituted by multiplying the Jacobian or its approximation with the ROM RIGHT BASIS, 
 * yielding a rectangular system (FOM size) that is then solved using the QR decomposition.
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet (as for the FOM).
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information (as for the FOM).
 * @tparam TSparseSpace The sparse system considered
 * @tparam TDenseSpace The dense system considered
 * @tparam TLinearSolver The linear solver considered
 * @author SebastianAres
 */

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver : public ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    // Class pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver);

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// The definition of the current class
    typedef GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    /// Definition of the classes from the base class
    using BaseBuilderAndSolverType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    typedef GlobalROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
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

    /// Definition of the classes from the ROM base class
    // typedef ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> ROMBuilderAndSolverType;

    /// Additional definitions
    typedef typename ModelPart::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    // Non-distributed, dense:
    typedef LocalSystemMatrixType RomSystemMatrixType;
    typedef LocalSystemVectorType RomSystemVectorType;

    //Distributed, dense
    typedef RomSystemMatrixType LSPGSystemMatrixType; 
    typedef RomSystemVectorType LSPGSystemVectorType;
    //      ^ Change this to a distributed dense type

    /// DoF types definition
    typedef Node NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;
    typedef moodycamel::ConcurrentQueue<DofType::Pointer> DofQueue;

    ///@}
    ///@name Life cycle
    ///@{

    explicit GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters) : BaseType(pNewLinearSystemSolver) 
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = this->ValidateAndAssignParameters(this_parameters_copy, this->GetDefaultParameters());
        this->AssignSettings(this_parameters_copy);
    } 

    ~GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 1)) << "Setting up the dofs" << std::endl;
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of threads" << ParallelUtilities::GetNumThreads() << "\n" << std::endl;
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        // Get model part data
        if (this->mHromWeightsInitialized == false) {
            this->InitializeHROMWeights(rModelPart);
        }

        // Compute the complementary mesh for HROM
        if (this->mHromSimulation){
            FindNeighbouringElementsAndConditions(rModelPart);
            // ComputeComplementaryElementsAndConditions(rModelPart);
        }

        auto dof_queue = this->ExtractDofSet(pScheme, rModelPart);
        
        // Fill a sorted auxiliary array of with the DOFs set
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;
        auto dof_array = this->SortAndRemoveDuplicateDofs(dof_queue);

        // Update base builder and solver DOFs array and set corresponding flag
        BaseType::GetDofSet().swap(dof_array);
        BaseType::SetDofSetIsInitializedFlag(true);

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(BaseType::GetDofSet().size() == 0) << "No degrees of freedom!" << std::endl;
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::GetDofSet().size() << std::endl;
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished setting up the dofs" << std::endl;

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

    void FindNeighbouringElementsAndConditions(ModelPart& rModelPart) 
    {
        mNeighbouringAndSelectedElements.clear();
        mNeighbouringAndSelectedConditions.clear();

        // Create sets for storing Ids
        std::set<int> selected_element_ids;
        std::set<int> selected_condition_ids;

        FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType> find_nodal_elements_neighbours_process(rModelPart);
        find_nodal_elements_neighbours_process.Execute();
        FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType> find_nodal_conditions_neighbours_process(rModelPart);
        find_nodal_conditions_neighbours_process.Execute();

        for (auto it_elem = this->mSelectedElements.ptr_begin(); it_elem != this->mSelectedElements.ptr_end(); ++it_elem) {
            mNeighbouringAndSelectedElements.push_back(*it_elem);
            selected_element_ids.insert((*it_elem)->Id());

            const auto& r_geom = (*it_elem)->GetGeometry();
            const SizeType n_nodes = r_geom.PointsNumber();
            for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
                NodeType::Pointer p_node = r_geom(i_node);
                auto& neighbour_elements = p_node->GetValue(NEIGHBOUR_ELEMENTS);
                for (auto& neighbor_elem : neighbour_elements.GetContainer()) {
                    if(selected_element_ids.find(neighbor_elem->Id()) == selected_element_ids.end()) {
                        mNeighbouringAndSelectedElements.push_back(&*neighbor_elem);
                        selected_element_ids.insert(neighbor_elem->Id());
                    }
                }

                // Check and add neighbouring conditions
                auto& neighbour_conditions = p_node->GetValue(NEIGHBOUR_CONDITIONS);
                for (auto& neighbor_cond : neighbour_conditions.GetContainer()) {
                    if(selected_condition_ids.find(neighbor_cond->Id()) == selected_condition_ids.end()) {
                        mNeighbouringAndSelectedConditions.push_back(&*neighbor_cond);
                        selected_condition_ids.insert(neighbor_cond->Id());
                    }
                }
            }
        }

        for (auto it_cond = this->mSelectedConditions.ptr_begin(); it_cond != this->mSelectedConditions.ptr_end(); ++it_cond) {
            mNeighbouringAndSelectedConditions.push_back(*it_cond);
            selected_condition_ids.insert((*it_cond)->Id());

            const auto& r_geom = (*it_cond)->GetGeometry();
            const SizeType n_nodes = r_geom.PointsNumber();
            for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
                NodeType::Pointer p_node = r_geom(i_node);
                auto& neighbour_conditions = p_node->GetValue(NEIGHBOUR_CONDITIONS);
                for (auto& neighbor_cond : neighbour_conditions.GetContainer()) {
                    if(selected_condition_ids.find(neighbor_cond->Id()) == selected_condition_ids.end()) {
                        mNeighbouringAndSelectedConditions.push_back(&*neighbor_cond);
                        selected_condition_ids.insert(neighbor_cond->Id());
                    }
                }

                // Check and add neighbouring elements
                auto& neighbour_elements = p_node->GetValue(NEIGHBOUR_ELEMENTS);
                for (auto& neighbor_elem : neighbour_elements.GetContainer()) {
                    if(selected_element_ids.find(neighbor_elem->Id()) == selected_element_ids.end()) {
                        mNeighbouringAndSelectedElements.push_back(&*neighbor_elem);
                        selected_element_ids.insert(neighbor_elem->Id());
                    }
                }
            }
        }

        std::unordered_map<int, Element::Pointer> unique_elements_map;
        for (const auto& elem : mNeighbouringAndSelectedElements) {
            unique_elements_map[elem.Id()] = intrusive_ptr<Element>(const_cast<Element*>(&elem));
        }
        mNeighbouringAndSelectedElements.clear();
        for (const auto& pair : unique_elements_map) {
            mNeighbouringAndSelectedElements.push_back(pair.second);
        }

        std::unordered_map<int, Condition::Pointer> unique_conditions_map;
        for (const auto& cond : mNeighbouringAndSelectedConditions) {
            unique_conditions_map[cond.Id()] = intrusive_ptr<Condition>(const_cast<Condition*>(&cond));
        }
        mNeighbouringAndSelectedConditions.clear();
        for (const auto& pair : unique_conditions_map) {
            mNeighbouringAndSelectedConditions.push_back(pair.second);
        }


        for (auto elem : this->mSelectedElements) {
        std::cout << elem.Id() << " ";
        }
        std::cout << std::endl;

        for (auto elem : this->mNeighbouringAndSelectedElements) {
            std::cout << elem.Id() << " ";
        }
        std::cout << std::endl;

        for (auto cond : this->mSelectedConditions) {
            std::cout << cond.Id() << " ";
        }
        std::cout << std::endl;

        for (auto cond : this->mNeighbouringAndSelectedConditions) {
            std::cout << cond.Id() << " ";
        }
        std::cout << std::endl;
    }


    void CheckForDuplicates() 
    {
        std::set<int> element_ids;
        std::set<int> condition_ids;

        for(auto& elem : mNeighbouringAndSelectedElements) {
            element_ids.insert(elem.Id());
        }

        for(auto& cond : mNeighbouringAndSelectedConditions) {
            condition_ids.insert(cond.Id());
        }

        if(element_ids.size() != mNeighbouringAndSelectedElements.size()) {
            std::cout << "There are duplicate elements in mNeighbouringAndSelectedElements." << std::endl;
        }

        if(condition_ids.size() != mNeighbouringAndSelectedConditions.size()) {
            std::cout << "There are duplicate conditions in mNeighbouringAndSelectedConditions." << std::endl;
        }
    }



    // TODO: Parallel Utilities
    void ComputeComplementaryElementsAndConditions(ModelPart& rModelPart) 
    {
        // Creating a standard vector of IDs from the selected elements and conditions
        std::vector<unsigned int> selectedElementIds;
        for (const auto& element: this->mSelectedElements) {
            selectedElementIds.push_back(element.Id());
        }

        std::vector<unsigned int> selectedConditionIds;
        for (const auto& condition: this->mSelectedConditions) {
            selectedConditionIds.push_back(condition.Id());
        }

        // Compute complementary elements
        for (auto& p_element : rModelPart.Elements().GetContainer())
        {
            auto found_element = std::find(selectedElementIds.begin(), selectedElementIds.end(), p_element->Id());

            // If the element is not found in selectedElementIds, add it to mComplementaryElements.
            if (found_element == selectedElementIds.end()) {
                mComplementaryElements.push_back(p_element);;
            }
        }

        // Compute complementary conditions
        for (auto& p_condition : rModelPart.Conditions().GetContainer())
        {
            auto found_condition = std::find(selectedConditionIds.begin(), selectedConditionIds.end(), p_condition->Id());

            // If the condition is not found in selectedConditionIds, add it to mComplementaryConditions.
            if (found_condition == selectedConditionIds.end()) {
                mComplementaryConditions.push_back(p_condition);
            }
        }
    }

    /**
     * Builds and projects the reduced system of equations 
     */
    virtual void BuildAndProjectROM(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb,
        TSystemVectorType &rDx)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        const auto assembling_timer = BuiltinTimer();

        if (rA.size1() != BaseType::mEquationSystemSize || rA.size2() != BaseType::mEquationSystemSize) {
            rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            BaseType::ConstructMatrixStructure(pScheme, rA, rModelPart);
        }

        this->Build(pScheme, rModelPart, rA, rb);

        BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        ProjectROM(rModelPart, rA, rb);

        double time = assembling_timer.ElapsedSeconds();
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Build and project time: " << time << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Projects the reduced system of equations 
     */
    virtual void ProjectROM(
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb)
    {
        KRATOS_TRY

        if (mRightRomBasisInitialized==false){
            mPhiGlobal = ZeroMatrix(BaseBuilderAndSolverType::GetEquationSystemSize(), GetNumberOfROMModes());
            mRightRomBasisInitialized = true;
        }

        BuildRightROMBasis(rModelPart);

        auto a_wrapper = UblasWrapper<double>(rA);
        const auto& eigen_rA = a_wrapper.matrix();
        Eigen::Map<EigenDynamicVector> eigen_rb(rb.data().begin(), rb.size());
        Eigen::Map<EigenDynamicMatrix> eigen_mPhiGlobal(mPhiGlobal.data().begin(), mPhiGlobal.size1(), mPhiGlobal.size2());
        
        EigenDynamicMatrix eigen_rA_times_mPhiGlobal = eigen_rA * eigen_mPhiGlobal; //TODO: Make it in parallel.

        // Compute the matrix multiplication
        mEigenRomA = eigen_rA_times_mPhiGlobal.transpose() * eigen_rA_times_mPhiGlobal; //TODO: Make it in parallel.
        mEigenRomB = eigen_rA_times_mPhiGlobal.transpose() * eigen_rb; //TODO: Make it in parallel.
        
        KRATOS_CATCH("")
    }

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        BuildAndProjectROM(pScheme, rModelPart, A, b, Dx);
        
        this->SolveROM(rModelPart, A, b, Dx);

        KRATOS_CATCH("")
    }

    // void BuildAndSolve(
    //     typename TSchemeType::Pointer pScheme,
    //     ModelPart &rModelPart,
    //     TSystemMatrixType &A,
    //     TSystemVectorType &Dx,
    //     TSystemVectorType &b) override
    // {
    //     KRATOS_TRY
    //     LSPGSystemMatrixType Arom = ZeroMatrix(BaseType::GetEquationSystemSize(), this->GetNumberOfROMModes());
    //     LSPGSystemVectorType brom = ZeroVector(BaseType::GetEquationSystemSize());
    //     BuildROM(pScheme, rModelPart, Arom, brom);

    //     // Obtain the assembled residuals vector (To build a basis for Petrov-Galerkin)
    //     if (mTrainPetrovGalerkinFlag){
    //         TSystemVectorType r_residual;
    //         GetAssembledResiduals(pScheme, rModelPart, r_residual);
    //     }

    //     SolveROM(rModelPart, Arom, brom, Dx);


    //     KRATOS_CATCH("")
    // }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "lspg_rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "train_petrov_galerkin" : false
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    static std::string Name() 
    {
        return "lspg_rom_builder_and_solver";
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
        return "GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's 
    virtual void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    double time_rom_system_contruction = 0;
    double time_rom_system_solving = 0;
    double time_rom_projection_back = 0;

    void ResetTimeMeasures(){
        time_rom_system_contruction = 0;
        time_rom_system_solving = 0;
        time_rom_projection_back = 0;
    }

    double Get_contruction_time(){
        return  time_rom_system_contruction;
    }

    double Get_solving_time(){
        return time_rom_system_solving;
    }

    double Get_projection_time(){
        return time_rom_projection_back;
    }


    ///@}
protected:

    ///@}
    ///@name Protected static member variables
    ///@{


    ///@}
    ///@name Protected member variables
    ///@{

    SizeType mNodalDofs;
    ElementsArrayType mComplementaryElements;
    ConditionsArrayType mComplementaryConditions;
    ElementsArrayType mNeighbouringAndSelectedElements;
    ConditionsArrayType mNeighbouringAndSelectedConditions;
    bool mTrainPetrovGalerkinFlag = false;

    ///@}
    ///@name Protected operators
    ///@{


    ///@}
    ///@name Protected operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // // Set member variables
        mTrainPetrovGalerkinFlag = ThisParameters["train_petrov_galerkin"].GetBool();
    }

    /**
    * Thread Local Storage containing dynamically allocated structures to avoid reallocating each iteration.
    */
    struct AssemblyTLS
    {
        Matrix phiE = {};                // Elemental Phi
        LocalSystemMatrixType lhs = {};  // Elemental LHS
        EquationIdVectorType eq_id = {}; // Elemental equation ID vector
        DofsVectorType dofs = {};        // Elemental dof vector
        RomSystemMatrixType romA;        // reduced LHS
        RomSystemVectorType romB;        // reduced RHS
    };

    /**
     * Resizes a Matrix if it's not the right size
     */
    template<typename TMatrix>
    static void ResizeIfNeeded(TMatrix& rMat, const SizeType Rows, const SizeType Cols)

    {
        if(rMat.size1() != Rows || rMat.size2() != Cols) {
            rMat.resize(Rows, Cols, false);
        }
    };


    /**
     * Builds the reduced system of equations on rank 0 
     */
    void BuildROM(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        LSPGSystemMatrixType &rA,
        LSPGSystemVectorType &rb) override
    {
        KRATOS_TRY
        // Define a dense matrix to hold the reduced problem
        rA = ZeroMatrix(BaseType::GetEquationSystemSize(), this->GetNumberOfROMModes());
        rb = ZeroVector(BaseType::GetEquationSystemSize());

        // Build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Get ProcessInfo from main model part
        const auto& r_current_process_info = rModelPart.GetProcessInfo();


        // Assemble all entities
        const auto assembling_timer = BuiltinTimer();

        AssemblyTLS assembly_tls_container;

        // const auto& r_elements = rModelPart.Elements();
        const auto& r_elements = this->mHromSimulation ? this->mSelectedElements : rModelPart.Elements();

        if(!r_elements.empty())
        {
            block_for_each(r_elements, assembly_tls_container, 
                [&](Element& r_element, AssemblyTLS& r_thread_prealloc)
            {
                CalculateLocalContributionLSPG(r_element, rA, rb, r_thread_prealloc, *pScheme, r_current_process_info, true);
            });
        }


        // const auto& r_conditions = rModelPart.Conditions();
        const auto& r_conditions = this->mHromSimulation ? this->mSelectedConditions : rModelPart.Conditions();

        if(!r_conditions.empty())
        {
            block_for_each(r_conditions, assembly_tls_container, 
                [&](Condition& r_condition, AssemblyTLS& r_thread_prealloc)
            {
                CalculateLocalContributionLSPG(r_condition, rA, rb, r_thread_prealloc, *pScheme, r_current_process_info, true);
            });
        }

        //////////
        //BUILD COMPLETE OPERATOR FOR HROM 
        Matrix rA_left;
        Vector rb_left;
        if (this->mHromSimulation){
            // rA_left = rA;
            // rb_left = rb;  
            AssemblyTLS assembly_tls_container_left;
            rA_left = ZeroMatrix(BaseType::GetEquationSystemSize(), this->GetNumberOfROMModes());
            rb_left = ZeroVector(BaseType::GetEquationSystemSize());
            // const auto& r_elements_left = rModelPart.Elements();
            auto& r_elements_left = mNeighbouringAndSelectedElements;
            // auto& r_elements_left = this->mHromSimulation ? this->mSelectedElements : rModelPart.Elements();
            // auto& r_elements_left = mComplementaryElements;

            if(!r_elements_left.empty())
            {
                block_for_each(r_elements_left, assembly_tls_container_left, 
                    [&](Element& r_element, AssemblyTLS& r_thread_prealloc_left)
                {
                    CalculateLocalContributionLSPG(r_element, rA_left, rb_left, r_thread_prealloc_left, *pScheme, r_current_process_info, false);
                });
            }

            // const auto& r_conditions_left = rModelPart.Conditions();
            auto& r_conditions_left = mNeighbouringAndSelectedConditions;
            // auto& r_conditions_left = mComplementaryConditions;
            // auto& r_conditions_left = this->mHromSimulation ? this->mSelectedConditions : rModelPart.Conditions();
            if(!r_conditions_left.empty())
            {
                block_for_each(r_conditions_left, assembly_tls_container_left, 
                    [&](Condition& r_condition, AssemblyTLS& r_thread_prealloc_left)
                {
                    CalculateLocalContributionLSPG(r_condition, rA_left, rb_left, r_thread_prealloc_left, *pScheme, r_current_process_info, false);
                });
            }

            // Initialize the mask vector with zeros
            Vector hrom_dof_mask_vector = ZeroVector(BaseType::GetEquationSystemSize());

            // Build the mask vector for selected elements and conditions
            BuildHromDofMaskVector(hrom_dof_mask_vector, r_current_process_info);

            // Zero out rows in the matrix that correspond to zero in the mask vector
            ApplyMaskToMatrixRows(rA_left, hrom_dof_mask_vector);
        }

        const auto projection_timer = BuiltinTimer();

        using EigenDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

        // Create the Eigen matrix using the buffer
        Eigen::Map<EigenDynamicMatrix> eigen_matrix(rA.data().begin(), rA.size1(), rA.size2());
        Eigen::Map<EigenDynamicVector> eigen_vector(rb.data().begin(), rb.size());

        if (this->mHromSimulation){
            Eigen::Map<EigenDynamicMatrix> eigen_matrix_left(rA_left.data().begin(), rA_left.size1(), rA_left.size2());
            // KRATOS_WATCH(eigen_matrix)
            // KRATOS_WATCH(eigen_matrix_left)
            mA_eigen = eigen_matrix_left.transpose() * eigen_matrix;
            mb_eigen = eigen_matrix_left.transpose() * eigen_vector;
        }
        else{
            // Compute the matrix multiplication
            mA_eigen = eigen_matrix.transpose() * eigen_matrix;
            mb_eigen = eigen_matrix.transpose() * eigen_vector;
        }
        
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Projection time: " << projection_timer.ElapsedSeconds() << std::endl;
        time_rom_system_contruction += assembling_timer.ElapsedSeconds();
        double time = assembling_timer.ElapsedSeconds();
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Build time: " << time << std::endl;
        
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    void ApplyMaskToMatrixRows(Matrix& rA_left, 
        const Vector& hrom_dof_mask_vector)
    {
        if(rA_left.size1() != hrom_dof_mask_vector.size()) 
        {
            throw std::invalid_argument("Matrix and vector sizes do not match");
        }
        
        for(std::size_t i = 0; i < hrom_dof_mask_vector.size(); ++i)
        {
            if(hrom_dof_mask_vector[i] == 0)
            {
                row(rA_left, i) = zero_vector<double>(rA_left.size2());
            }
        }
    }


    void BuildHromDofMaskVector(
        Vector& rHromDofMaskVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const auto& r_hrom_elements = this->mSelectedElements;
        const auto& r_hrom_conditions = this->mSelectedConditions;

        // Ensuring the vector has the correct type for atomic operations.
        std::vector<std::atomic<int>> atomicHromDofMaskVector(rHromDofMaskVector.size());

        block_for_each(r_hrom_elements, [&](Element& r_element)
        {
            Element::DofsVectorType hrom_dofs;
            r_element.GetDofList(hrom_dofs, rCurrentProcessInfo);
            for(std::size_t i = 0; i < hrom_dofs.size(); ++i)
            {
                const Dof<double>& r_dof = *hrom_dofs[i];
                std::atomic_fetch_or(&atomicHromDofMaskVector[r_dof.EquationId()], 1);
            }
        });

        block_for_each(r_hrom_conditions, [&](Condition& r_condition)
        {
            Condition::DofsVectorType hrom_dofs;
            r_condition.GetDofList(hrom_dofs, rCurrentProcessInfo);
            for(std::size_t i = 0; i < hrom_dofs.size(); ++i)
            {
                const Dof<double>& r_dof = *hrom_dofs[i];
                std::atomic_fetch_or(&atomicHromDofMaskVector[r_dof.EquationId()], 1);
            }
        });

        // Copy back the atomic vector to the original one after all threads have finished their work.
        for(std::size_t i = 0; i < atomicHromDofMaskVector.size(); ++i)
        {
            rHromDofMaskVector[i] = atomicHromDofMaskVector[i].load();
        }
    }


    /**
     * Solves reduced system of equations and broadcasts it
     */
    void SolveROM(
        ModelPart &rModelPart,
        LSPGSystemMatrixType &rA,
        LSPGSystemVectorType &rb,
        TSystemVectorType &rDx) override
    {
        KRATOS_TRY

        LSPGSystemVectorType dxrom(this->GetNumberOfROMModes());
        
        const auto solving_timer = BuiltinTimer();

        using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
        Eigen::Map<EigenDynamicVector> dxrom_eigen(dxrom.data().begin(), dxrom.size());
        dxrom_eigen = mA_eigen.colPivHouseholderQr().solve(mb_eigen);
        // KRATOS_WATCH(mA_eigen)
        // KRATOS_WATCH(dxrom_eigen)
        time_rom_system_solving+= solving_timer.ElapsedSeconds();
        double time = solving_timer.ElapsedSeconds();
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Solve reduced system time: " << time << std::endl;
        

        // Save the ROM solution increment in the root modelpart database
        auto& r_root_mp = rModelPart.GetRootModelPart();
        noalias(r_root_mp.GetValue(ROM_SOLUTION_INCREMENT)) += dxrom;

        // project reduced solution back to full order model
        const auto backward_projection_timer = BuiltinTimer();
        this->ProjectToFineBasis(dxrom, rModelPart, rDx);
        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Project to fine basis time: " << backward_projection_timer.ElapsedSeconds() << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Builds the reduced system of equations on rank 0 
     */
    void GetAssembledResiduals(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        LSPGSystemVectorType &rb) 
    {
        KRATOS_TRY
        const auto residual_writing_timer = BuiltinTimer();
        // Define a dense matrix to hold the reduced problem
        rb = ZeroVector(BaseType::GetEquationSystemSize());

        // Build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Get ProcessInfo from main model part
        const auto& r_current_process_info = rModelPart.GetProcessInfo();

        AssemblyTLS assembly_tls_container;

        const auto& r_elements = rModelPart.Elements();

        if(!r_elements.empty())
        {
            block_for_each(r_elements, assembly_tls_container, 
                [&](Element& r_element, AssemblyTLS& r_thread_prealloc)
            {
                CalculateAssembledResiduals(r_element, rb, r_thread_prealloc, *pScheme, r_current_process_info);
            });
        }


        const auto& r_conditions = rModelPart.Conditions();

        if(!r_conditions.empty())
        {
            block_for_each(r_conditions, assembly_tls_container, 
                [&](Condition& r_condition, AssemblyTLS& r_thread_prealloc)
            {
                CalculateAssembledResiduals(r_condition, rb, r_thread_prealloc, *pScheme, r_current_process_info);
            });
        }

        std::stringstream matrix_market_vector_name;
        matrix_market_vector_name << "R_" << rModelPart.GetProcessInfo()[TIME] << "_" << rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER] <<  ".res.mm";
        SparseSpaceType::WriteMatrixMarketVector((matrix_market_vector_name.str()).c_str(), rb);

        KRATOS_INFO_IF("GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Write residuals to train Petrov Galerkin time: " << residual_writing_timer.ElapsedSeconds() << std::endl;
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
    bool mRightRomBasisInitialized = false;
    Eigen::MatrixXd mA_eigen;
    Eigen::VectorXd mb_eigen;
    Matrix mPhiGlobal;

    ///@}
    ///@name Private operations 
    ///@{

    /**
     * Computes the local contribution of an element or condition for LSPG
     */
    template<typename TEntity>
    void CalculateLocalContributionLSPG(
        TEntity& rEntity,
        LSPGSystemMatrixType& rAglobal,
        LSPGSystemVectorType& rBglobal,
        AssemblyTLS& rPreAlloc,
        TSchemeType& rScheme,
        const ProcessInfo& rCurrentProcessInfo,
        bool rHromWeightFlag)
    {
        if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE)) return;

        // Calculate elemental contribution
        rScheme.CalculateSystemContributions(rEntity, rPreAlloc.lhs, rPreAlloc.romB, rPreAlloc.eq_id, rCurrentProcessInfo);
        if (rHromWeightFlag){
            const double hrom_weight = this->mHromSimulation ? rEntity.GetValue(HROM_WEIGHT) : 1.0;
            rPreAlloc.lhs *= hrom_weight;
            rPreAlloc.romB *= hrom_weight;
        }
        
        rEntity.GetDofList(rPreAlloc.dofs, rCurrentProcessInfo);

        const SizeType ndofs = rPreAlloc.dofs.size();
        ResizeIfNeeded(rPreAlloc.phiE, ndofs, this->GetNumberOfROMModes());
        ResizeIfNeeded(rPreAlloc.romA, ndofs, this->GetNumberOfROMModes());

        const auto &r_geom = rEntity.GetGeometry();
        RomAuxiliaryUtilities::GetPhiElemental(rPreAlloc.phiE, rPreAlloc.dofs, r_geom, this->mMapPhi);

        // if (rHromWeightFlag){
        //     const double hrom_weight = this->mHromSimulation ? rEntity.GetValue(HROM_WEIGHT) : 1.0;
        //     noalias(rPreAlloc.romA) = prod(rPreAlloc.lhs, rPreAlloc.phiE) * hrom_weight;
        // }
        // else{
        //     noalias(rPreAlloc.romA) = prod(rPreAlloc.lhs, rPreAlloc.phiE);
        // }

        noalias(rPreAlloc.romA) = prod(rPreAlloc.lhs, rPreAlloc.phiE);

        // Assembly
        for(SizeType row=0; row < ndofs; ++row)
        {
            const SizeType global_row = rPreAlloc.eq_id[row];

            double& r_bi = rBglobal(global_row);
            AtomicAdd(r_bi, rPreAlloc.romB[row]);

            if(rPreAlloc.dofs[row]->IsFree())
            {
                for(SizeType col=0; col < this->GetNumberOfROMModes(); ++col)
                {
                    // const SizeType global_col = rPreAlloc.eq_id[col];
                    const SizeType global_col = col;
                    double& r_aij = rAglobal(global_row, global_col);
                    AtomicAdd(r_aij, rPreAlloc.romA(row, col));
                }
            }
        }
    }

    /**
     * Computes the global assembled residual 
     */
    template<typename TEntity>
    void CalculateAssembledResiduals(
        TEntity& rEntity,
        LSPGSystemVectorType& rBglobal,
        AssemblyTLS& rPreAlloc,
        TSchemeType& rScheme,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE)) return;

        // Calculate elemental contribution
        rScheme.CalculateRHSContribution(rEntity, rPreAlloc.romB, rPreAlloc.eq_id, rCurrentProcessInfo);
        rEntity.GetDofList(rPreAlloc.dofs, rCurrentProcessInfo);

        const SizeType ndofs = rPreAlloc.dofs.size();

        // Assembly
        for(SizeType row=0; row < ndofs; ++row)
        {
            const SizeType global_row = rPreAlloc.eq_id[row];

            if(rPreAlloc.dofs[row]->IsFree())
            {
                double& r_bi = rBglobal(global_row);
                AtomicAdd(r_bi, rPreAlloc.romB[row]);
            }
        }
    }




    ///@}
}; /* Class GlobalLeastSquaresPetrovGalerkinROMBuilderAndSolver */

///@}
///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/