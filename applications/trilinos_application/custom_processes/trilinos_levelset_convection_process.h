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
#include "mpi.h"
#include "Epetra_MpiComm.h"

// Project includes
#include "containers/model.h"
#include "includes/communicator.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "processes/levelset_convection_process.h"
// #include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

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
        ModelPart& rBaseModelPart,
        LinearSolverPointerType pLinearSolver,
        const double MaxCFL = 1.0,
        const double CrossWindStabilizationFactor = 0.7,
        const unsigned int MaxSubSteps = 0)
        : LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>(
            rLevelSetVar, 
            rBaseModelPart, 
            MaxCFL, 
            MaxSubSteps),
        mrEpetraCommunicator(rEpetraCommunicator)
    {
        KRATOS_TRY
        
        // Check that there is at least one element and node in the model
        int n_nodes = rBaseModelPart.NumberOfNodes();
        int n_elems = rBaseModelPart.NumberOfElements();

        if (n_nodes > 0){
            VariableUtils().CheckVariableExists< Variable< double > >(rLevelSetVar, rBaseModelPart.Nodes());
            VariableUtils().CheckVariableExists< Variable< array_1d < double, 3 > > >(VELOCITY, rBaseModelPart.Nodes());
        }

        if(TDim == 2){
            KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle) << 
                "In 2D the element type is expected to be a triangle" << std::endl;
        } else if(TDim == 3) {
            KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra) << 
                "In 3D the element type is expected to be a tetrahedra" << std::endl;
        }

        rBaseModelPart.GetCommunicator().SumAll(n_nodes);
        rBaseModelPart.GetCommunicator().SumAll(n_elems);

        KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;
        
        // Allocate if needed the variable DYNAMIC_TAU of the process info, and if it does not exist, set it to zero
        if( rBaseModelPart.GetProcessInfo().Has(DYNAMIC_TAU) == false){
            rBaseModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU,0.0);
        }

        // Allocate if needed the variable CONVECTION_DIFFUSION_SETTINGS of the process info, and create it if it does not exist
        if(rBaseModelPart.GetProcessInfo().Has(CONVECTION_DIFFUSION_SETTINGS) == false){
            ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
            rBaseModelPart.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
            p_conv_diff_settings->SetUnknownVariable(rLevelSetVar);
            p_conv_diff_settings->SetConvectionVariable(VELOCITY);
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

        ModelPart& r_distance_model_part = rBaseModelPart.GetOwnerModel().GetModelPart("DistanceConvectionPart");
        (this->mpSolvingStrategy) = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            r_distance_model_part,
            p_scheme,
            pLinearSolver,
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

        if(rBaseModelPart.GetOwnerModel().HasModelPart("DistanceConvectionPart"))
            rBaseModelPart.GetOwnerModel().DeleteModelPart("DistanceConvectionPart");

        ModelPart& r_distance_model_part = rBaseModelPart.GetOwnerModel().CreateModelPart("DistanceConvectionPart");

        // Generate

        r_distance_model_part.Nodes().clear();
        r_distance_model_part.Conditions().clear();
        r_distance_model_part.Elements().clear();

        r_distance_model_part.SetProcessInfo(rBaseModelPart.pGetProcessInfo());
        r_distance_model_part.SetBufferSize(base_buffer_size);
        r_distance_model_part.SetProperties(rBaseModelPart.pProperties());
        r_distance_model_part.Tables() = rBaseModelPart.Tables();

        // Assigning the nodes to the new model part
        r_distance_model_part.Nodes() = rBaseModelPart.Nodes();

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

        r_distance_model_part.SetCommunicator(p_new_comm);

        // Generating the elements
        (r_distance_model_part.Elements()).reserve(rBaseModelPart.NumberOfElements());
        for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem){
            Element::Pointer p_element = Kratos::make_shared< LevelSetConvectionElementSimplex < TDim, TDim+1 > >(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            // Assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = it_elem->pGetGeometry();
            
            (r_distance_model_part.Elements()).push_back(p_element);
            (r_distance_model_part.GetCommunicator()).LocalMesh().Elements().push_back(p_element);
        }
       
        // Resize the arrays
        const auto n_nodes = r_distance_model_part.NumberOfNodes();
        (this->mVelocity).resize(n_nodes);
        (this->mVelocityOld).resize(n_nodes);
        (this->mOldDistance).resize(n_nodes);

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


