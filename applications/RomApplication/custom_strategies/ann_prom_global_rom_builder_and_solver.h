//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Nicolás Sibuet, Sebastian Ares de Parga
//                  
//

#pragma once

// System includes

// External includes
#include "concurrentqueue/concurrentqueue.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "utilities/builtin_timer.h"
#include "utilities/reduction_utilities.h"
#include "custom_utilities/ublas_wrapper.h"

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

/**
 * @class AnnPromGlobalROMBuilderAndSolver
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class AnnPromGlobalROMBuilderAndSolver : public ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
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
    KRATOS_CLASS_POINTER_DEFINITION(AnnPromGlobalROMBuilderAndSolver);

    // The size_t types
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    /// The definition of the current class
    using ClassType = AnnPromGlobalROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Definition of the classes from the base class
    using BaseBuilderAndSolverType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    using BaseType = ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TSchemeType = typename BaseBuilderAndSolverType::TSchemeType;
    using DofsArrayType = typename BaseBuilderAndSolverType::DofsArrayType;
    using TSystemMatrixType = typename BaseBuilderAndSolverType::TSystemMatrixType;
    using TSystemVectorType = typename BaseBuilderAndSolverType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseBuilderAndSolverType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseBuilderAndSolverType::LocalSystemMatrixType;
    using TSystemMatrixPointerType = typename BaseBuilderAndSolverType::TSystemMatrixPointerType;
    using TSystemVectorPointerType = typename BaseBuilderAndSolverType::TSystemVectorPointerType;
    using ElementsArrayType = typename BaseBuilderAndSolverType::ElementsArrayType;
    using ConditionsArrayType = typename BaseBuilderAndSolverType::ConditionsArrayType;

    /// Additional definitions
    using MasterSlaveConstraintContainerType = typename ModelPart::MasterSlaveConstraintContainerType;
    using EquationIdVectorType = typename Element::EquationIdVectorType;
    using DofsVectorType = typename Element::DofsVectorType;
    using CompressedMatrixType = boost::numeric::ublas::compressed_matrix<double>;
    using RomSystemMatrixType = LocalSystemMatrixType;
    using RomSystemVectorType = LocalSystemVectorType;
    using EigenDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    using EigenSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

    /// DoF types definition
    using DofType = typename Node::DofType;
    using DofPointerType = typename DofType::Pointer;
    using DofQueue = moodycamel::ConcurrentQueue<DofType::Pointer>;

    ///@}
    ///@name Life cycle
    ///@{

    explicit AnnPromGlobalROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters): BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = this->ValidateAndAssignParameters(this_parameters_copy, this->GetDefaultParameters());
        this->AssignSettings(this_parameters_copy);
    }

    explicit AnnPromGlobalROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
    }

    ~AnnPromGlobalROMBuilderAndSolver() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseBuilderAndSolverType::Pointer Create(
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

        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 1)) << "Setting up the dofs" << std::endl;
        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of threads" << ParallelUtilities::GetNumThreads() << "\n" << std::endl;
        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        // Get model part data
        if (mHromWeightsInitialized == false) {
            InitializeHROMWeights(rModelPart);
        }

        auto dof_queue = ExtractDofSet(pScheme, rModelPart);

        // Fill a sorted auxiliary array of with the DOFs set
        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;
        auto dof_array = SortAndRemoveDuplicateDofs(dof_queue);

        // Update base builder and solver DOFs array and set corresponding flag
        BaseBuilderAndSolverType::GetDofSet().swap(dof_array);
        BaseBuilderAndSolverType::SetDofSetIsInitializedFlag(true);

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(BaseBuilderAndSolverType::GetDofSet().size() == 0) << "No degrees of freedom!" << std::endl;
        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseBuilderAndSolverType::GetDofSet().size() << std::endl;
        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Finished setting up the dofs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        if (BaseBuilderAndSolverType::GetCalculateReactionsFlag())
        {
            for (const auto& r_dof: BaseBuilderAndSolverType::GetDofSet())
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
        auto& r_dof_set = BaseBuilderAndSolverType::GetDofSet();
        BaseBuilderAndSolverType::mEquationSystemSize = r_dof_set.size();
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

    void SetNumberOfROMModes(SizeType numberOfRomModes)
    {
        mNumberOfRomModes = numberOfRomModes;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "SetNumberOfROMModes: " << mNumberOfRomModes << std::endl;
    }

    void SetNumberOfNNLayers(SizeType numberOfNNLayers)
    {
        mNNLayers = vector<Matrix> (numberOfNNLayers);
        mGradientLayers = vector<Matrix> (numberOfNNLayers);
        mLayersOut = vector<Vector> (numberOfNNLayers+1);
        mIntermediateGradients = vector<Matrix> (numberOfNNLayers);
    }

    void SetNNLayer(
        SizeType layerIndex,
        Matrix nnLayer
    )
    {
        mNNLayers[layerIndex]= nnLayer;
        mGradientLayers[layerIndex] = Matrix (nnLayer.size1(), nnLayer.size2());
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "mNNLayers[layerIndex] size: " << mNNLayers[layerIndex].size1() << mNNLayers[layerIndex].size2() << std::endl;
    }

    void InitializeLastItDecoderOutAndGradientsAndAuxiliaries()
    {
        mIntermediateGradients[0] = Matrix (mNNLayers[0].size1(), mNNLayers[0].size2());
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "mIntermediateGradients[0] size: " << mIntermediateGradients[0].size1() << " " << mIntermediateGradients[0].size2() << std::endl;


        for (int i = 0; i < mGradientLayers.size()-1; i++)
        {
            mLayersOut[i] = Vector (mNNLayers[i].size1());
            mIntermediateGradients[i+1] = Matrix (mNNLayers[0].size1(), mNNLayers[i+1].size2());

            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "mLayersOut[i] size: " << mLayersOut[i].size() << std::endl;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "mIntermediateGradients[i+1] size: " << mIntermediateGradients[i+1].size1() << " " << mIntermediateGradients[i+1].size2() << std::endl;
        }

        mLayersOut[mLayersOut.size()-2] = Vector (mNNLayers[mNNLayers.size()-1].size1());
        mLayersOut[mLayersOut.size()-1] = Vector (mNNLayers[mNNLayers.size()-1].size2());

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "mLayersOut[end-1] size: " << mLayersOut[mLayersOut.size()-2].size() << std::endl;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "mLayersOut[end] size: " << mLayersOut[mLayersOut.size()-1].size() << std::endl;
    

        // Vector x;
        // GetXAndDecoderGradient(xrom, mPhiGlobal, x);

        // mLastItDecoderOut = x;

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "LastItDecoderOut: " << mLastItDecoderOut << std::endl;
    }

    void SetPhiMatrices(
        Matrix svdPhiInf,
        Matrix svdPhiSup
    )
    {
        mSVDPhiMatrices = vector<Matrix>(2);
        mSVDPhiMatrices[0] = svdPhiInf;
        mSVDPhiMatrices[1] = svdPhiSup;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "SetPhiMatrices: " << mSVDPhiMatrices << std::endl;
    }

    void SetTestMatrix(Matrix testMatrix)
    {
        mTestMatrix = testMatrix;
    }

    void SetTestVector(Vector testVector)
    {
        mTestVector = testVector;
    }

    bool GetMonotonicityPreservingFlag() const noexcept
    {
        return mMonotonicityPreservingFlag;
    }

    void BuildRightROMBasis(
        const ModelPart& rModelPart,
        Matrix& rPhiGlobal) 
    {

        const auto& r_dof_set = BaseBuilderAndSolverType::GetDofSet();
        block_for_each(r_dof_set, [&](const DofType& r_dof)
        {
            const auto& r_node = rModelPart.GetNode(r_dof.Id());
            const Matrix& r_rom_nodal_basis = r_node.GetValue(ROM_BASIS);
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildRightROMBasis: r_rom_nodal_basis:" << r_rom_nodal_basis << std::endl;
            const Matrix::size_type row_id = mMapPhi.at(r_dof.GetVariable().Key());
            if (r_dof.IsFixed())
            {
                noalias(row(rPhiGlobal, r_dof.EquationId())) = ZeroVector(r_rom_nodal_basis.size2());
            }
            else{
                noalias(row(rPhiGlobal, r_dof.EquationId())) = row(r_rom_nodal_basis, row_id);
            }
        });
    }

    void MonotonicityPreserving(
        TSystemMatrixType& rA,
        TSystemVectorType& rB
    ) 
    {
        const auto& r_dof_set = BaseType::GetDofSet();
        Vector dofs_values = ZeroVector(r_dof_set.size());

        block_for_each(r_dof_set, [&](Dof<double>& rDof){
            const std::size_t id = rDof.EquationId();
            dofs_values[id] = rDof.GetSolutionStepValue();
        });
        double *values_vector = rA.value_data().begin();
        std::size_t *index1_vector = rA.index1_data().begin();
        std::size_t *index2_vector = rA.index2_data().begin();

        IndexPartition<std::size_t>(rA.size1()).for_each(
            [&](std::size_t i)
            {
                for (std::size_t k = index1_vector[i]; k < index1_vector[i + 1]; k++) {
                    const double value = values_vector[k];
                    if (value > 0.0) {
                        const auto j = index2_vector[k];
                        if (j > i) {
                            // TODO: Partition in blocks to gain efficiency by avoiding thread locks.
                            // Values conflicting with other threads
                            auto& r_aij = rA(i,j).ref();
                            AtomicAdd(r_aij, -value);
                            auto& r_aji = rA(j,i).ref();
                            AtomicAdd(r_aji, -value);
                            auto& r_aii = rA(i,i).ref();
                            AtomicAdd(r_aii, value);
                            auto& r_ajj = rA(j,j).ref();
                            AtomicAdd(r_ajj, value);
                            auto& r_bi = rB[i];
                            AtomicAdd(r_bi, value*dofs_values[j] - value*dofs_values[i]);
                            auto& r_bj = rB[j];
                            AtomicAdd(r_bj, value*dofs_values[i] - value*dofs_values[j]);
                        }
                    }
                }
            }
        );
    }

    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        // Call the base B&S InitializeSolutionStep
        BaseBuilderAndSolverType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        auto& r_root_mp = rModelPart.GetRootModelPart();
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "InitializeSolutionStep: rom_solution_increment: " << r_root_mp.GetValue(ROM_SOLUTION_INCREMENT) << std::endl;


        // Reset the ROM solution increment in the root modelpart database
        // auto& r_root_mp = rModelPart.GetRootModelPart();
        // r_root_mp.GetValue(ROM_SOLUTION_INCREMENT) = ZeroVector(GetNumberOfROMModes());
    }

    
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndSolve: " << "StartedBuildAndSolveStep" << std::endl;

        BuildAndProjectROM(pScheme, rModelPart, A, b, Dx);

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndSolve: " << "Between BuildAndProject and Solve" << std::endl;
        
        SolveROM(rModelPart, mEigenRomA, mEigenRomB, Dx);

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
        if (r_Dx.size() != BaseBuilderAndSolverType::GetEquationSystemSize()) {
            r_Dx.resize(BaseBuilderAndSolverType::GetEquationSystemSize(), false);
        }

        TSystemVectorType& r_b = *pb;
        if (r_b.size() != BaseBuilderAndSolverType::GetEquationSystemSize()) {
            r_b.resize(BaseBuilderAndSolverType::GetEquationSystemSize(), false);
        }

        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "ann_prom_global_rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "rom_bns_settings" : {
                "monotonicity_preserving": false
            }
        })");
        default_parameters.AddMissingParameters(BaseBuilderAndSolverType::GetDefaultParameters());

        return default_parameters;
    }

    static std::string Name()
    {
        return "ann_prom_global_rom_builder_and_solver";
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
        return "AnnPromGlobalROMBuilderAndSolver";
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
    bool mRightRomBasisInitialized = false;

    ///@}
    ///@name Protected operators
    ///@{

    ///@}
    ///@name Protected operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseBuilderAndSolverType::AssignSettings(ThisParameters);

        // Set member variables
        mNodalDofs = ThisParameters["nodal_unknowns"].size();
        mNumberOfRomModes = ThisParameters["number_of_rom_dofs"].GetInt();
        mMonotonicityPreservingFlag = ThisParameters["rom_bns_settings"]["monotonicity_preserving"].GetBool();

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

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rA (init):" << rA.size1() << rA.size2() << std::endl;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rb (init):" << rb.size() << std::endl;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rDx (init):" << rDx.size() << std::endl;

        if (rA.size1() != BaseType::mEquationSystemSize || rA.size2() != BaseType::mEquationSystemSize) {
            rA.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            BaseType::ConstructMatrixStructure(pScheme, rA, rModelPart);
        }

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rA (structured):" << rA.size1() << rA.size2() << std::endl;

        Build(pScheme, rModelPart, rA, rb);

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rA (final):" << rA.size1() << rA.size2() << std::endl;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rb (final):" << rb.size() << std::endl;

        if (mMonotonicityPreservingFlag) {
            BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);
            MonotonicityPreserving(rA, rb);
        }

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rA (monoton):" << rA.size1() << rA.size2() << std::endl;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rb (monoton):" << rb.size() << std::endl;

        auto& r_root_mp = rModelPart.GetRootModelPart();
        RomSystemVectorType& rRomUnknowns = r_root_mp.GetValue(ROM_SOLUTION_INCREMENT);


        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: rRomUnknowns (init):" << rRomUnknowns << std::endl;
        
        // Vector testDiff = mTestVector - rDx;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "BuildAndProjectROM: error in rDx (init): " << testDiff << std::endl;

        ProjectROM(rModelPart, rA, rb);

        double time = assembling_timer.ElapsedSeconds();
        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Build and project time: " << time << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the LHS and RHS multiplied by its corresponding hrom weight. 
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

        // // Getting the elements from the model
        const int nelements = mHromSimulation ? mSelectedElements.size() : rModelPart.Elements().size();

        // // Getting the array of the conditions
        const int nconditions = mHromSimulation ? mSelectedConditions.size() : rModelPart.Conditions().size();

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = mHromSimulation ? mSelectedElements.begin() : rModelPart.Elements().begin();
        ModelPart::ConditionsContainerType::iterator cond_begin = mHromSimulation ? mSelectedConditions.begin() : rModelPart.Conditions().begin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // Assemble all elements
        const auto timer = BuiltinTimer();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++) {
                auto it_elem = el_begin + k;

                if (it_elem->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_elem, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //Get HROM weight and multiply it by its contribution
                    const double h_rom_weight = mHromSimulation ? it_elem->GetValue(HROM_WEIGHT) : 1.0;
                    LHS_Contribution *= h_rom_weight;
                    RHS_Contribution *= h_rom_weight;

                    // Assemble the elemental contribution
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
                }

            }

            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++) {
                auto it_cond = cond_begin + k;

                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_cond, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //Get HROM weight and multiply it by its contribution
                    const double h_rom_weight = mHromSimulation ? it_cond->GetValue(HROM_WEIGHT) : 1.0;
                    LHS_Contribution *= h_rom_weight;
                    RHS_Contribution *= h_rom_weight;

                    // Assemble the elemental contribution
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_INFO_IF("GlobalROMResidualBasedBlockBuilderAndSolver", this->GetEchoLevel() >= 1) << "Build time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("GlobalROMResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

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
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "mPhiGlobal: " << mPhiGlobal << std::endl;
            auto& r_root_mp = rModelPart.GetRootModelPart();
            Vector& xrom = r_root_mp.GetValue(ROM_SOLUTION_INCREMENT);
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "xrom: " << xrom << std::endl;
            // Vector x = Vector (BaseBuilderAndSolverType::GetEquationSystemSize());
            mLastItDecoderOut = ZeroVector(BaseBuilderAndSolverType::GetEquationSystemSize());
            GetXAndDecoderGradient(xrom, mPhiGlobal, mLastItDecoderOut);
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "After GetGradients" << std::endl;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "After LastItDecoderOut update" << std::endl;
            mRightRomBasisInitialized = true;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "LastItDecoderOut: " << mLastItDecoderOut << std::endl;
        }

        WriteElementalBasis(mPhiGlobal, rModelPart);

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "Inside ProjectROM, Right basis already initialized" << std::endl;

        BuildRightROMBasis(rModelPart, mPhiGlobal);

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "ProjectROM: mPhiGlobal (updated):" << mPhiGlobal << std::endl;

        auto a_wrapper = UblasWrapper<double>(rA);
        const auto& eigen_rA = a_wrapper.matrix();
        Eigen::Map<EigenDynamicVector> eigen_rb(rb.data().begin(), rb.size());
        Eigen::Map<EigenDynamicMatrix> eigen_mPhiGlobal(mPhiGlobal.data().begin(), mPhiGlobal.size1(), mPhiGlobal.size2());

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "ProjectROM: eigen_mPhiGlobal:" << eigen_mPhiGlobal << std::endl;
        
        EigenDynamicMatrix eigen_rA_times_mPhiGlobal = eigen_rA * eigen_mPhiGlobal; //TODO: Make it in parallel.

        // Compute the matrix multiplication
        mEigenRomA = eigen_mPhiGlobal.transpose() * eigen_rA_times_mPhiGlobal; //TODO: Make it in parallel.
        mEigenRomB = eigen_mPhiGlobal.transpose() * eigen_rb; //TODO: Make it in parallel.
        
        KRATOS_CATCH("")
    }

    /**
     * Solves reduced system of equations and broadcasts it
     */
    virtual void SolveROM(
        ModelPart &rModelPart,
        EigenDynamicMatrix &rEigenRomA,
        EigenDynamicVector &rEigenRomB,
        TSystemVectorType &rDx)
    {
        KRATOS_TRY

        RomSystemVectorType dxrom(GetNumberOfROMModes());
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "SolveROM: dxrom (init): " << dxrom << std::endl;
        
        const auto solving_timer = BuiltinTimer();

        using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
        Eigen::Map<EigenDynamicVector> dxrom_eigen(dxrom.data().begin(), dxrom.size());
        dxrom_eigen = rEigenRomA.colPivHouseholderQr().solve(rEigenRomB);
        
        double time = solving_timer.ElapsedSeconds();
        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Solve reduced system time: " << time << std::endl;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "SolveROM: Solve reduced system time: " << time << std::endl;
        
        // Save the ROM solution increment in the root modelpart database
        auto& r_root_mp = rModelPart.GetRootModelPart();

        RomSystemVectorType& xrom = r_root_mp.GetValue(ROM_SOLUTION_INCREMENT);
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "SolveROM: xrom (init)" << xrom << std::endl;

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "SolveROM: dxrom (updated): " << dxrom << std::endl;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "SolveROM: xrom (updated)" << xrom << std::endl;

        // project reduced solution back to full order model
        const auto backward_projection_timer = BuiltinTimer();
        Vector x (rDx.size());
        GetXAndDecoderGradient(xrom+dxrom, mPhiGlobal,x);
        ProjectToFineBasis(x, rDx);
        KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Project to fine basis time: " << backward_projection_timer.ElapsedSeconds() << std::endl;

        noalias(xrom) += dxrom;

        KRATOS_CATCH("")
    }

    void ProjectToFineBasis(
        const Vector& rx,
        TSystemVectorType& rDx)
    {
        rDx = rx - mLastItDecoderOut;
        mLastItDecoderOut = rx;
    }

    void ProjectToFineBasis_old(
        const TSystemVectorType& rRomUnknowns,
        const RomSystemVectorType& dxrom,
        TSystemVectorType& rDx)
    {
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "ProjectToFineBasis: PostProcess(NN(q_inf))" << NNForwardPass(rRomUnknowns) << std::endl;

        Eigen::Map<EigenDynamicVector> eigen_rDx(rDx.data().begin(), rDx.size());
        eigen_rDx = NNForwardPass(rRomUnknowns+dxrom) - NNForwardPass(rRomUnknowns);
    }

    void GetXAndDecoderGradient(
        const RomSystemVectorType& rRomUnknowns,
        Matrix& rPhiGlobal,
        Vector& rx)
    {   
        mGradientLayers = mNNLayers;
        
        mLayersOut[0] = rRomUnknowns;

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetXAndDecoderGradient: mLayersOut[0]:" << mLayersOut[0] << std::endl;

        int i=0;
        for(const auto& layers_item : mNNLayers)
        {
            mLayersOut[i+1] = boost::numeric::ublas::prod(boost::numeric::ublas::trans(layers_item), mLayersOut[i]);
            if(i<mNNLayers.size()-1){
                IndexPartition<IndexType>(mLayersOut[i+1].size()).for_each([&](IndexType Index)
                {   
                    if (mLayersOut[i+1][Index]<0.0){
                        for (int j = 0; j < mGradientLayers[i].size1(); j++)
                        {
                            mGradientLayers[i](j,Index) = mGradientLayers[i](j,Index)*exp(mLayersOut[i+1][Index]);
                        }
                        mLayersOut[i+1][Index] = exp(mLayersOut[i+1][Index]) - 1.0;
                    }
                });
            }
            i++;
        }

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: mLayersOut[end]:" << mLayersOut[mLayersOut.size()-1] << std::endl;

        mIntermediateGradients[0] = mGradientLayers[0];

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: mIntermediateGradients[0]:" << mIntermediateGradients[0] << std::endl;


        for (int i = 0; i < mGradientLayers.size()-1; i++)
        {
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: i:" << i << std::endl;
            mIntermediateGradients[i+1] = boost::numeric::ublas::prod(mIntermediateGradients[i],mGradientLayers[i+1]);
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: mIntermediateGradients[i+1] size:" << mIntermediateGradients[i+1].size1() << " " << mIntermediateGradients[i+1].size2() << std::endl;
        }

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: mIntermediateGradients[end]:" << mIntermediateGradients[mIntermediateGradients.size()-1] << std::endl;

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: After Intermediate Calculations:" << std::endl;
        PostProcessGradient(mIntermediateGradients[mIntermediateGradients.size()-1], rPhiGlobal);
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: AlmostFinished" << std::endl;
        PostProcessQSup(mLayersOut[mLayersOut.size()-1], rRomUnknowns, rx);
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: Finished" << std::endl;

    }

    void GetDecoderGradient_old(
        const RomSystemVectorType& rRomUnknowns,
        Matrix& rPhiGlobal)
    {   
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: rRomUnknowns:" << rRomUnknowns << std::endl;

        vector<Matrix> grad_layers = mNNLayers;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: mNNLayers size: " << mNNLayers.size() << std::endl;
        Vector rRomUnknowns_temp = rRomUnknowns;
        Eigen::Map<EigenDynamicVector> eigen_rRomUnknowns_temp(rRomUnknowns_temp.data().begin(), rRomUnknowns_temp.size());
        EigenDynamicVector layer_in = eigen_rRomUnknowns_temp;
        // Eigen::Map<EigenDynamicVector> eigen_layer_in(layer_in.data().begin(), layer_in.size());
        
        EigenDynamicVector layer_out;
        int i=0;
        for(const auto& layers_item : mNNLayers)
        {
            // noalias(layer_out) = prod(layer,layer_in);
            Matrix layer = layers_item;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: layer " << i << ": " << layer << std::endl;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: layer_in " << i << ": " << layer_in << std::endl;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: grad_layers " << i << " (init): " << grad_layers[i] << std::endl;
            Eigen::Map<EigenDynamicMatrix> eigen_layer(layer.data().begin(), layer.size1(), layer.size2());
            // layer_out = eigen_layer.transpose() * eigen_layer_in;
            layer_out = eigen_layer.transpose() * layer_in;
            if(i<mNNLayers.size()-1)
            {
                IndexPartition<IndexType>(layer_out.size()).for_each([&](IndexType Index)
                {   
                    if (layer_out[Index]<0.0){
                        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: layer_out element " << Index << ": " << layer_out[Index] << std::endl;
                        for (int j = 0; j < grad_layers[i].size1(); j++)
                        {
                            // We will need to check the dimensionspf the matrices and determine if we have to multiply rows or columns
                            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: grad_layers row " << j << "," << Index << "(init): " << grad_layers[i](j,Index) << std::endl;
                            grad_layers[i](j,Index) = grad_layers[i](j,Index)*exp(layer_out[Index]);
                            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: grad_layers row " << j << "," << Index << "(updated): " << grad_layers[i](j,Index) << std::endl;
                        }
                    }
                
                    // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: grad_layers " << i << " (init): " << grad_layers[i] << std::endl;
                });
                // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: grad_layers " << i << " (final): " << grad_layers[i] << std::endl;
                // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: layer_out " << i << " (init): " << layer_out << std::endl;
                ELUActivation(layer_out);
            }
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: layer_out " << i << " (final): " << layer_out << std::endl;
            layer_in = layer_out;
            i++;
        }

        EigenDynamicMatrix intermediate_gradients;
        Eigen::Map<EigenDynamicMatrix> eigen_gradients_in_temp(grad_layers[0].data().begin(), grad_layers[0].size1(), grad_layers[0].size2());
        EigenDynamicMatrix gradients_in = eigen_gradients_in_temp.transpose();
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: gradients_in[0]: " << gradients_in << std::endl;

        for (int i = 0; i < grad_layers.size()-1; i++)
        {
            // noalias(intermediate_gradients) = prod(grad_layers[i+1],gradients_in);Matrix layer = layers_item;
            Eigen::Map<EigenDynamicMatrix> eigen_grad_layer(grad_layers[i+1].data().begin(), grad_layers[i+1].size1(), grad_layers[i+1].size2());
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: eigen_grad_layer " << i+1 << ": " << eigen_grad_layer << std::endl;
            intermediate_gradients = eigen_grad_layer.transpose() * gradients_in;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: intermediate_gradients " << i+1 << ": " << intermediate_gradients << std::endl;
            gradients_in = intermediate_gradients;
        }

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: intermediate_gradients (qsup): " << intermediate_gradients << std::endl;

        PostProcessGradient(intermediate_gradients, rPhiGlobal);
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: rPhiGlobal: " << rPhiGlobal << std::endl;

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: mTestMatrix " << mTestMatrix << std::endl;

        // Eigen::Map<EigenDynamicMatrix> eigen_test_matrix(mTestMatrix.data().begin(), mTestMatrix.size1(), mTestMatrix.size2());
        // EigenDynamicMatrix testDiff = eigen_test_matrix - rPhiGlobal;
        // Matrix testDiff = mTestMatrix - rPhiGlobal;

        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "GetDecoderGradient: error in gradient (final): " << testDiff << std::endl;


    }

    EigenDynamicVector NNForwardPass(
        const Vector& nn_input)
    {
        Vector nn_input_temp = nn_input;
        Eigen::Map<EigenDynamicVector> eigen_nn_input_temp(nn_input_temp.data().begin(), nn_input_temp.size());
        EigenDynamicVector layer_in = eigen_nn_input_temp;

        EigenDynamicVector layer_out;

        int i=0;
        for(const auto& layers_item : mNNLayers)
        {
            Matrix layer = layers_item;
            Eigen::Map<EigenDynamicMatrix> eigen_layer(layer.data().begin(), layer.size1(), layer.size2());
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "NNForwardPass: layer_in " << i << ": " << layer_in << std::endl;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "NNForwardPass: layer " << i << ": " << layer << std::endl;
            layer_out = eigen_layer.transpose() * layer_in;
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "NNForwardPass: layer_out " << i << ": " << layer_out << std::endl;
            if (i<mNNLayers.size()-1)
            {
                ELUActivation(layer_out);
            }
            layer_in = layer_out;
            i++;
        }
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "NNForwardPass: NN(q_inf)" << layer_out << std::endl;

        return PostProcessQSup(layer_out, nn_input);
    }

    void ELUActivation(
        EigenDynamicVector& layer_out)
    {
        IndexPartition<IndexType>(layer_out.size()).for_each([&](IndexType Index)
        {
            if (layer_out[Index]<0.0){
                layer_out[Index] = exp(layer_out[Index]) - 1.0;
            }
        });
    }

    void PostProcessQSup(
        const Vector& q_sup,
        const Vector& nn_input,
        Vector& rx)
    {
        rx = boost::numeric::ublas::prod(mSVDPhiMatrices[0], nn_input) + boost::numeric::ublas::prod(mSVDPhiMatrices[1], q_sup);
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "PostProcessQSup: rx:"  << rx << std::endl;

    }

    EigenDynamicVector PostProcessQSup_old(
        const EigenDynamicVector& q_sup,
        const Vector& nn_input)
    {
        Vector q_inf = nn_input;
        Eigen::Map<EigenDynamicVector> eigen_q_inf(q_inf.data().begin(), q_inf.size());
        Matrix SVD_Phi_inf = mSVDPhiMatrices[0];
        Matrix SVD_Phi_sup = mSVDPhiMatrices[1];
        Eigen::Map<EigenDynamicMatrix> eigen_SVD_Phi_inf(SVD_Phi_inf.data().begin(), SVD_Phi_inf.size1(), SVD_Phi_inf.size2());
        Eigen::Map<EigenDynamicMatrix> eigen_SVD_Phi_sup(SVD_Phi_sup.data().begin(), SVD_Phi_sup.size1(), SVD_Phi_sup.size2());
        EigenDynamicVector nn_out = eigen_SVD_Phi_inf * eigen_q_inf;
        nn_out += eigen_SVD_Phi_sup * q_sup;
        return nn_out;
    }

    void PostProcessGradient(
        const Matrix& intermediate_gradients,
        Matrix& rPhiGlobal)
    {
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "PostProcessGradient: mSVDPhiMatrices[0] size:" << mSVDPhiMatrices[0].size1() << " " << mSVDPhiMatrices[0].size2() << std::endl;
        // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "PostProcessGradient: mSVDPhiMatrices[1] size:" << mSVDPhiMatrices[1].size1() << " " << mSVDPhiMatrices[1].size2() << std::endl;

        rPhiGlobal=mSVDPhiMatrices[0];
        rPhiGlobal += boost::numeric::ublas::prod(mSVDPhiMatrices[1], boost::numeric::ublas::trans(intermediate_gradients));
    }

    void PostProcessGradient_old(
        const EigenDynamicMatrix& intermediate_gradients,
        Matrix& rPhiGlobal)
    {
        rPhiGlobal=mSVDPhiMatrices[0];
        Matrix SVD_Phi_sup = mSVDPhiMatrices[1];
        Eigen::Map<EigenDynamicMatrix> eigen_SVD_Phi_sup(SVD_Phi_sup.data().begin(), SVD_Phi_sup.size1(), SVD_Phi_sup.size2());
        Eigen::Map<EigenDynamicMatrix> eigen_rPhiGlobal(rPhiGlobal.data().begin(), rPhiGlobal.size1(), rPhiGlobal.size2());
        eigen_rPhiGlobal += eigen_SVD_Phi_sup * intermediate_gradients;
    }

    void WriteElementalBasis(
        const Matrix& rPhiGlobal,
        ModelPart& rModelPart) 
    {
        block_for_each(rModelPart.Nodes().GetContainer(),
            [&](Node::Pointer p_node)
        {
            Matrix& r_rom_nodal_basis = p_node->GetValue(ROM_BASIS);
            // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "WriteElementalBasis: r_rom_nodal_basis " << p_node->Id() << ": " << r_rom_nodal_basis << std::endl;
            
            for(int i=0; i < mNodalDofs; i++)
            {
                
                // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "WriteElementalBasis: rPhiGlobal row node " << p_node->Id() << " dof " << i << ": " << row(rPhiGlobal, (p_node->Id()-1)*mNodalDofs+i) << std::endl;
                // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "WriteElementalBasis: r_rom_nodal_basis row node " << p_node->Id() << " dof " << i << ": " << row(r_rom_nodal_basis, i) << std::endl;
                row(r_rom_nodal_basis,i)=row(rPhiGlobal, (p_node->Id()-1)*mNodalDofs+i);
                // KRATOS_INFO_IF("AnnPromGlobalROMBuilderAndSolver", (this->GetEchoLevel() >= 0)) << "WriteElementalBasis: r_rom_nodal_basis row node " << p_node->Id() << " dof " << i << ": " << row(r_rom_nodal_basis, i) << std::endl;
            }
            p_node->SetValue(ROM_BASIS, r_rom_nodal_basis);
            
        });
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
    ///@name Private member variables 
    ///@{
        
    SizeType mNumberOfRomModes;
    EigenDynamicMatrix mEigenRomA;
    EigenDynamicVector mEigenRomB;
    Matrix mPhiGlobal;
    bool mMonotonicityPreservingFlag;
    vector<Matrix> mNNLayers;
    vector<Matrix> mGradientLayers;
    vector<Matrix> mIntermediateGradients;
    vector<Vector> mLayersOut;
    Vector mLastItDecoderOut;
    vector<Matrix> mSVDPhiMatrices;
    Matrix mTestMatrix;
    Vector mTestVector;

    ///@}
    ///@name Private operations 
    ///@{

    ///@}
}; /* Class AnnPromGlobalROMBuilderAndSolver */

///@}
///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/