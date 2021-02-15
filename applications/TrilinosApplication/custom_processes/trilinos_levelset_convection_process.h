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

#if !defined(KRATOS_TRILINOS_LEVELSET_CONVECTION_PROCESS_INCLUDED )
#define  KRATOS_TRILINOS_LEVELSET_CONVECTION_PROCESS_INCLUDED

// System includes

// External includes
#include "Epetra_MpiComm.h"

// Project includes
#include "containers/model.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "processes/levelset_convection_process.h"
#include "utilities/pointer_communicator.h"
#include "utilities/pointer_map_communicator.h"

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
    typedef typename BaseType::SchemeType::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverPointerType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    /// Pointer definition of TrilinosLevelSetConvectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosLevelSetConvectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     */
    TrilinosLevelSetConvectionProcess(
        Epetra_MpiComm& rEpetraCommunicator,
        Variable<double>& rLevelSetVar,
        Variable< array_1d< double, 3 > >& rConvectVar,
        ModelPart& rBaseModelPart,
        LinearSolverPointerType pLinearSolver,
        const double MaxCFL = 1.0,
        const double CrossWindStabilizationFactor = 0.7,
        const unsigned int MaxSubSteps = 0,
        const bool IsBFECC = false,
        const bool PartialDt = false)
        : LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(
            rLevelSetVar,
            rConvectVar,
            rBaseModelPart,
            MaxCFL,
            MaxSubSteps,
            IsBFECC,
            PartialDt),
        mrEpetraCommunicator(rEpetraCommunicator)
    {
        KRATOS_TRY

        // Check that there is at least one element and node in the model
        int n_nodes = rBaseModelPart.NumberOfNodes();
        int n_elems = rBaseModelPart.NumberOfElements();

        if (n_nodes > 0){
            VariableUtils().CheckVariableExists< Variable< double > >(rLevelSetVar, rBaseModelPart.Nodes());
            VariableUtils().CheckVariableExists< Variable< array_1d < double, 3 > > >(rConvectVar, rBaseModelPart.Nodes());
        }

        n_nodes = rBaseModelPart.GetCommunicator().GetDataCommunicator().SumAll(n_nodes);
        n_elems = rBaseModelPart.GetCommunicator().GetDataCommunicator().SumAll(n_elems);

        // Check if the modelpart is globaly empty
        KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;

        // Check if any partition has incorrect elements
        if(TDim == 2){
            int has_incorrect_elems = rBaseModelPart.NumberOfElements() ? rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle : 0;
            has_incorrect_elems = rBaseModelPart.GetCommunicator().GetDataCommunicator().SumAll(has_incorrect_elems);
            KRATOS_ERROR_IF(has_incorrect_elems) << "In 2D the element type is expected to be a triangle" << std::endl;
        } else if(TDim == 3) {
            int has_incorrect_elems = rBaseModelPart.NumberOfElements() ? rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra : 0;
            has_incorrect_elems = rBaseModelPart.GetCommunicator().GetDataCommunicator().SumAll(has_incorrect_elems);
            KRATOS_ERROR_IF(has_incorrect_elems) << "In 3D the element type is expected to be a tetrahedra" << std::endl;
        }

        // Allocate if needed the variable DYNAMIC_TAU of the process info, and if it does not exist, set it to zero
        if( rBaseModelPart.GetProcessInfo().Has(DYNAMIC_TAU) == false){
            rBaseModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU,0.0);
        }

        // Allocate if needed the variable CONVECTION_DIFFUSION_SETTINGS of the process info, and create it if it does not exist
        if(rBaseModelPart.GetProcessInfo().Has(CONVECTION_DIFFUSION_SETTINGS) == false){
            ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
            rBaseModelPart.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
            p_conv_diff_settings->SetUnknownVariable(rLevelSetVar);
            p_conv_diff_settings->SetConvectionVariable(rConvectVar);
        }

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        (this->mDistancePartIsInitialized) = false;
        ReGenerateConvectionModelPart(rBaseModelPart);

        // Generate a linear strategy
        SchemePointerType p_scheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();

        const int row_size_guess = (TDim == 2 ? 15 : 40);
        BuilderSolverPointerType p_builder_and_solver = Kratos::make_shared< TrilinosBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver > >(
            mrEpetraCommunicator,
            row_size_guess,
            pLinearSolver);

        const bool calculate_reactions = false;
        const bool reform_dof_at_each_iteration = false;
        const bool calculate_norm_Dx_flag = false;

        (this->mpSolvingStrategy) = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            *BaseType::mpDistanceModelPart,
            p_scheme,
            p_builder_and_solver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_Dx_flag);

        (this->mpSolvingStrategy)->SetEchoLevel(0);

        rBaseModelPart.GetProcessInfo().SetValue(CROSS_WIND_STABILIZATION_FACTOR, CrossWindStabilizationFactor);

        //TODO: check flag DO_EXPENSIVE_CHECKS
        (this->mpSolvingStrategy)->Check();

        KRATOS_CATCH("")
    }

    TrilinosLevelSetConvectionProcess(
        Epetra_MpiComm& rEpetraCommunicator,
        Variable<double>& rLevelSetVar,
        ModelPart& rBaseModelPart,
        LinearSolverPointerType pLinearSolver,
        const double MaxCFL = 1.0,
        const double CrossWindStabilizationFactor = 0.7,
        const unsigned int MaxSubSteps = 0)
        :   TrilinosLevelSetConvectionProcess(
                rEpetraCommunicator,
                rLevelSetVar,
                VELOCITY,
                rBaseModelPart,
                pLinearSolver,
                MaxCFL,
                CrossWindStabilizationFactor,
                MaxSubSteps,
                false,
                false) {}

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
        VariableUtils().AddDof< Variable < double> >(this->mrLevelSetVar, rBaseModelPart);

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
        for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem){
            Element::Pointer p_element = Kratos::make_intrusive< LevelSetConvectionElementSimplex < TDim, TDim+1 > >(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            // Assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = it_elem->pGetGeometry();

            (BaseType::mpDistanceModelPart->Elements()).push_back(p_element);
            (BaseType::mpDistanceModelPart->GetCommunicator()).LocalMesh().Elements().push_back(p_element);
        }

        // Resize the arrays
        const auto n_nodes = BaseType::mpDistanceModelPart->NumberOfNodes();
        (this->mVelocity).resize(n_nodes);
        (this->mVelocityOld).resize(n_nodes);
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


    void EvaluateLimiter() override
    {
        // ****************************************************************************************
        // ****************************************************************************************
        // Calculating nodal limiter using \beta_ij = 1 (works fine on symmetric structural meshes)
        // D. Kuzmin et al. / Comput. Methods Appl. Mech. Engrg. 322 (2017) 23–41
        const double epsilon = 1.0e-15;
        const double power_bfecc = 2.0;

        //mProjectedGradientProcess.Execute();

        DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();

        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(this->mpDistanceModelPart->NumberOfNodes()); ++i_node){
            auto it_node = this->mpDistanceModelPart->NodesBegin() + i_node;

            it_node->SetValue(this->mrLevelSetVar, it_node->FastGetSolutionStepValue(this->mrLevelSetVar)); //Store mrLevelSetVar

            const auto X_i = it_node->Coordinates();
            const auto grad_i = it_node->GetValue(DISTANCE_GRADIENT);

            double S_plus = 0.0;
            double S_minus = 0.0;

            GlobalPointersVector< Node<3 > >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);
            std::unique_ptr<GlobalPointerCommunicator< Node<3 > > > pointer_comm(new
                GlobalPointerCommunicator< Node<3 > >(r_default_comm, global_pointer_list));

            auto coordinate_proxy = pointer_comm->Apply(
                [](GlobalPointer<Node<3> >& global_pointer) -> Point::CoordinatesArrayType
                {
                    return global_pointer->Coordinates();
                }
            );

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {
                // if (it_node->Id() == j_node->Id())
                //     continue;

                auto& global_pointer = global_pointer_list(j);
                auto X_j = coordinate_proxy.Get(global_pointer);

                S_plus += std::max(0.0, inner_prod(grad_i, X_i-X_j));
                S_minus += std::min(0.0, inner_prod(grad_i, X_i-X_j));
            }

            this->mSigmaPlus[i_node] = std::min(1.0, (std::abs(S_minus)+epsilon)/(S_plus+epsilon));
            this->mSigmaMinus[i_node] = std::min(1.0, (S_plus+epsilon)/(std::abs(S_minus)+epsilon));

            //pointer_comm.~GlobalPointerCommunicator();
        }

        //Calculating beta_ij in a way that the linearity is preserved on non-symmetrical meshes
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(this->mpDistanceModelPart->NumberOfNodes()); ++i_node){
            auto it_node = this->mpDistanceModelPart->NodesBegin() + i_node;
            const double distance_i = it_node->FastGetSolutionStepValue(this->mrLevelSetVar);
            const auto X_i = it_node->Coordinates();
            const auto grad_i = it_node->GetValue(DISTANCE_GRADIENT);

            double numerator = 0.0;
            double denominator = 0.0;

            GlobalPointersVector< Node<3 > >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);
            //GlobalPointerCommunicator< Node<3 > > pointer_comm(r_default_comm, global_pointer_list);
            std::unique_ptr<GlobalPointerCommunicator< Node<3 > > > pointer_comm(new
                GlobalPointerCommunicator< Node<3 > >(r_default_comm, global_pointer_list));

            auto coordinate_proxy = pointer_comm->Apply(
                [](GlobalPointer<Node<3> >& global_pointer) -> Point::CoordinatesArrayType
                {
                    return global_pointer->Coordinates();
                }
            );

            auto distance_proxy = pointer_comm->Apply(
                [&](GlobalPointer<Node<3> >& global_pointer) -> double
                {
                    return global_pointer->FastGetSolutionStepValue(this->mrLevelSetVar); //TODO: Get mrLevelSetVar
                }
            );

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {
                // if (it_node->Id() == j_node->Id())
                //     continue;

                auto& global_pointer = global_pointer_list(j);
                auto X_j = coordinate_proxy.Get(global_pointer);
                const double distance_j = distance_proxy.Get(global_pointer);

                double beta_ij = 1.0;

                if (inner_prod(grad_i, X_i-X_j) > 0)
                    beta_ij = this->mSigmaPlus[i_node];
                else if (inner_prod(grad_i, X_i-X_j) < 0)
                    beta_ij = this->mSigmaMinus[i_node];

                numerator += beta_ij*(distance_i - distance_j);
                denominator += beta_ij*std::abs(distance_i - distance_j);
            }

            const double fraction = (std::abs(numerator)/*  + epsilon */) / (denominator + epsilon);
            this->mLimiter[i_node] = 1.0 - std::pow(fraction, power_bfecc);

            //pointer_comm.~GlobalPointerCommunicator();
        }
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

    /// Copy constructor.
    //TrilinosLevelSetConvectionProcess(TrilinosLevelSetConvectionProcess const& rOther);

    ///@}
}; // Class TrilinosLevelSetConvectionProcess

// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver > const Kratos::Flags TrilinosLevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::PERFORM_STEP1(Kratos::Flags::Create(0));
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver > const Kratos::Flags TrilinosLevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::DO_EXPENSIVE_CHECKS(Kratos::Flags::Create(1));

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}
}  // namespace Kratos.

#endif // KRATOS_TRILINOS_LEVELSET_CONVECTION_PROCESS_INCLUDED  defined
