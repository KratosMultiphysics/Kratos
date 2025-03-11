//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes
#include "Epetra_MpiComm.h"

// Project includes
#include "containers/model.h"
#include "processes/levelset_convection_process.h"

// Application includes
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"

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

/// Short class definition.
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and convects a level set distance
 * on the top of it
*/
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver >
class TrilinosLevelSetConvectionProcess
    : public LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    KRATOS_DEFINE_LOCAL_FLAG(PERFORM_STEP1);
    KRATOS_DEFINE_LOCAL_FLAG(DO_EXPENSIVE_CHECKS);

    ///@name Type Definitions
    ///@{

    typedef LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename TLinearSolver::Pointer LinearSolverPointerType;
    typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverPointerType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    /// Pointer definition of TrilinosLevelSetConvectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosLevelSetConvectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    TrilinosLevelSetConvectionProcess(
        Epetra_MpiComm& rEpetraCommunicator,
        Model& rModel,
        typename TLinearSolver::Pointer pLinearSolver,
        Parameters ThisParameters)
        : TrilinosLevelSetConvectionProcess(
            rEpetraCommunicator,
            rModel.GetModelPart(ThisParameters["model_part_name"].GetString()),
            pLinearSolver,
            ThisParameters)
    {
    }

    TrilinosLevelSetConvectionProcess(
        Epetra_MpiComm& rEpetraCommunicator,
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        Parameters ThisParameters)
        : BaseType(
            rBaseModelPart,
            ThisParameters)
        , mrEpetraCommunicator(rEpetraCommunicator)
    {
        KRATOS_TRY

        const int row_size_guess = (TDim == 2 ? 15 : 40);
        auto p_builder_and_solver = Kratos::make_shared< TrilinosBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>>(
            mrEpetraCommunicator,
            row_size_guess,
            pLinearSolver);
        InitializeConvectionStrategy(p_builder_and_solver);

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    TrilinosLevelSetConvectionProcess(TrilinosLevelSetConvectionProcess const& rOther) = delete;

    /// Destructor.
    ~TrilinosLevelSetConvectionProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
    std::string Info() const override {
        return "TrilinosLevelSetConvectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "TrilinosLevelSetConvectionProcess";
    }

    // /// Print object's data.
    // void PrintData(std::ostream& rOStream) const override {
    // }

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

    Epetra_MpiComm& mrEpetraCommunicator;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ReGenerateConvectionModelPart(ModelPart& rBaseModelPart) override {

        KRATOS_TRY

        // Check buffer size
        const auto base_buffer_size = rBaseModelPart.GetBufferSize();
        KRATOS_ERROR_IF(base_buffer_size < 2) <<
            "Base model part buffer size is " << base_buffer_size << ". Set it to a minimum value of 2." << std::endl;

        if(rBaseModelPart.GetModel().HasModelPart("DistanceConvectionPart"))
            rBaseModelPart.GetModel().DeleteModelPart("DistanceConvectionPart");

        BaseType::mpDistanceModelPart= &(rBaseModelPart.GetModel().CreateModelPart("DistanceConvectionPart"));

        // Generate
        BaseType::mpDistanceModelPart->Nodes().clear();
        BaseType::mpDistanceModelPart->Conditions().clear();
        BaseType::mpDistanceModelPart->Elements().clear();

        BaseType::mpDistanceModelPart->SetProcessInfo(rBaseModelPart.pGetProcessInfo());
        BaseType::mpDistanceModelPart->SetBufferSize(base_buffer_size);
        BaseType::mpDistanceModelPart->SetProperties(rBaseModelPart.pProperties());
        BaseType::mpDistanceModelPart->Tables() = rBaseModelPart.Tables();

        // Assigning the nodes to the new model part
        BaseType::mpDistanceModelPart->Nodes() = rBaseModelPart.Nodes();

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof<Variable<double>>(*BaseType::mpLevelSetVar, rBaseModelPart);

        // Copy communicator data
        Communicator& r_base_comm = rBaseModelPart.GetCommunicator();
        Communicator::Pointer p_new_comm = r_base_comm.Create();

        p_new_comm->SetNumberOfColors(r_base_comm.GetNumberOfColors());
        p_new_comm->NeighbourIndices() = r_base_comm.NeighbourIndices();
        p_new_comm->LocalMesh().SetNodes(r_base_comm.LocalMesh().pNodes());
        p_new_comm->InterfaceMesh().SetNodes(r_base_comm.InterfaceMesh().pNodes());
        p_new_comm->GhostMesh().SetNodes(r_base_comm.GhostMesh().pNodes());
        for (unsigned int i = 0; i < r_base_comm.GetNumberOfColors(); ++i){
            p_new_comm->pInterfaceMesh(i)->SetNodes(r_base_comm.pInterfaceMesh(i)->pNodes());
            p_new_comm->pLocalMesh(i)->SetNodes(r_base_comm.pLocalMesh(i)->pNodes());
            p_new_comm->pGhostMesh(i)->SetNodes(r_base_comm.pGhostMesh(i)->pNodes());
        }

        BaseType::mpDistanceModelPart->SetCommunicator(p_new_comm);

        // Generating the elements
        (BaseType::mpDistanceModelPart->Elements()).reserve(rBaseModelPart.NumberOfElements());
        KRATOS_ERROR_IF(BaseType::mpConvectionFactoryElement == nullptr) << "Convection factory element has not been set yet." << std::endl;
        for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem){
            // Create the new element from the factory registered one
            auto p_element = BaseType::mpConvectionFactoryElement->Create(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            (BaseType::mpDistanceModelPart->Elements()).push_back(p_element);
            (BaseType::mpDistanceModelPart->GetCommunicator()).LocalMesh().Elements().push_back(p_element);
        }

        // Initialize the nodal and elemental databases
        BaseType::InitializeDistanceModelPartDatabases();

        // Resize the arrays
        const auto n_nodes = BaseType::mpDistanceModelPart->NumberOfNodes();
        (this->mVelocity).resize(n_nodes);
        (this->mVelocityOld).resize(n_nodes);
        (this->mMeshVelocity).resize(n_nodes);
        (this->mMeshVelocityOld).resize(n_nodes);
        (this->mOldDistance).resize(n_nodes);

        if (this->mIsBfecc){
            (this->mError).resize(n_nodes);
            (this->mSigmaPlus).resize(n_nodes);
            (this->mSigmaMinus).resize(n_nodes);
            (this->mLimiter).resize(n_nodes);
        }

        (this->mDistancePartIsInitialized) = true;

        KRATOS_CATCH("")
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void InitializeConvectionStrategy(BuilderSolverPointerType pBuilderAndSolver)
    {
        KRATOS_TRY

        // Get auxiliary member variables from base class
        auto& r_base_model_part = BaseType::mrBaseModelPart;
        const auto& r_level_set_var = *BaseType::mpLevelSetVar;
        const auto& r_convect_var = *BaseType::mpConvectVar;
        const auto& r_mesh_convect_var = *BaseType::mpMeshConvectVar;

        // Check the nodal database of the current partition
        VariableUtils().CheckVariableExists<Variable<double>>(r_level_set_var, r_base_model_part.Nodes());
        VariableUtils().CheckVariableExists<Variable<array_1d<double,3>>>(r_convect_var, r_base_model_part.Nodes());
        VariableUtils().CheckVariableExists<Variable<array_1d<double,3>>>(r_mesh_convect_var, r_base_model_part.Nodes());

        // Check if the modelpart is globally empty
        KRATOS_ERROR_IF(r_base_model_part.GetCommunicator().GlobalNumberOfNodes() == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(r_base_model_part.GetCommunicator().GlobalNumberOfElements() == 0) << "The model has no elements." << std::endl;

        // Check if any partition has incorrect elements
        if constexpr (TDim == 2){
            bool has_incorrect_elems = r_base_model_part.NumberOfElements() ? r_base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle : false;
            KRATOS_ERROR_IF(has_incorrect_elems) << "In 2D the element type is expected to be a triangle" << std::endl;
        } else if constexpr (TDim == 3) {
            bool has_incorrect_elems = r_base_model_part.NumberOfElements() ? r_base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra : false;
            KRATOS_ERROR_IF(has_incorrect_elems) << "In 3D the element type is expected to be a tetrahedra" << std::endl;
        }

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        ReGenerateConvectionModelPart(r_base_model_part);

        // Generate a linear strategy
        const bool calculate_reactions = false;
        const bool reform_dof_at_each_iteration = false;
        const bool calculate_norm_Dx_flag = false;
        auto p_scheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();
        (this->mpSolvingStrategy) = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            *BaseType::mpDistanceModelPart,
            p_scheme,
            pBuilderAndSolver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_Dx_flag);

        //TODO: check flag DO_EXPENSIVE_CHECKS
        this->mpSolvingStrategy->SetEchoLevel(0);
        this->mpSolvingStrategy->Check();
        this->mpSolvingStrategy->Initialize();

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    TrilinosLevelSetConvectionProcess& operator=(TrilinosLevelSetConvectionProcess const& rOther);

    ///@}
}; // Class TrilinosLevelSetConvectionProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}
}  // namespace Kratos.
