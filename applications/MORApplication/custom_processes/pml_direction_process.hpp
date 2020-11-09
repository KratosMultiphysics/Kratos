//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt

#if !defined(KRATOS_PML_DIRECTION_PROCESS_H_INCLUDED)
#define KRATOS_PML_DIRECTION_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "mor_application_variables.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"



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

/**
 * @class PMLDirectionProcess
 *
 * @ingroup MORApplication
 *
 * @brief This method provides PML damping vector field.
 * @details Based on solving a convection problem from the inner interface of the PML to the outer boundary, the PML damping direction is calculated for every node in the PML.
*/
template< class TSparseSpace, class TDenseSpace, class TLinearSolver >
class PMLDirectionProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;

    // typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    // typedef typename BaseType::TSchemeType TSchemeType;

    // typedef typename BaseType::DofsArrayType DofsArrayType;

    // typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    // typedef typename BaseType::TSystemVectorType TSystemVectorType;

    // typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    // typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    // typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    // typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    // typedef typename ComplexSparseSpaceType::MatrixType TSolutionMatrixType;


    /// Pointer definition of MonolithicMappingProcess
    KRATOS_CLASS_POINTER_DEFINITION(PMLDirectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PMLDirectionProcess(
        ModelPart& rModelPartPML,
        ModelPart& rModelPartInterface,
        ModelPart& rModelPartBoundary,
        typename TLinearSolver::Pointer plinear_solver
        ) : mrModelPartPML(rModelPartPML),
        mrModelPartInterface(rModelPartInterface),
        mrModelPartBoundary(rModelPartBoundary)
    {
        KRATOS_TRY
        
        // Check that there is at least one element and node in the model
        const auto n_nodes = mrModelPartPML.NumberOfNodes();
        const auto n_elems = mrModelPartPML.NumberOfElements();

        KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;

        // Set build level

        mrModelPartPML.GetProcessInfo()[BUILD_LEVEL] = 401; 

        // Generate a linear strategy
        typename SchemeType::Pointer pScheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();
        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        BuilderSolverTypePointer pBuilderSolver = Kratos::make_shared< ResidualBasedBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver > >(plinear_solver);
        mpSolvingStrategy = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            mrModelPartPML,
            pScheme,
            plinear_solver,
            pBuilderSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);
        


        mpSolvingStrategy->SetEchoLevel(0);



        // mpK = TSparseSpace::CreateEmptyMatrixPointer();
        
        // mpRHS = TSparseSpace::CreateEmptyVectorPointer();
        // mpDx = TSparseSpace::CreateEmptyVectorPointer();



        //TODO: check flag DO_EXPENSIVE_CHECKS
        mpSolvingStrategy->Check();

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~PMLDirectionProcess() override = default;


    ///@}
    ///@name Operators
    ///@{

    void operator()(){
        Execute();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& rPMLNodes = mrModelPartPML.Nodes();
        ModelPart::NodesContainerType& rInterfaceNodes = mrModelPartInterface.Nodes();
        ModelPart::NodesContainerType& rBoundaryNodes = mrModelPartBoundary.Nodes();
        
        for ( ModelPart::NodesContainerType::iterator it_node = rInterfaceNodes.begin();
                it_node != rInterfaceNodes.end() ; ++it_node)
        {
            it_node->SetValue(PRESCRIBED_POTENTIAL, -1);
        }
        
        for ( ModelPart::NodesContainerType::iterator it_node = rBoundaryNodes.begin();
                it_node != rBoundaryNodes.end() ; ++it_node)
        {
            it_node->SetValue(PRESCRIBED_POTENTIAL, 1);
        }


        mpSolvingStrategy->Solve();
        // TSystemMatrixType& r_K  = *mpK;
        // TSystemVectorType& r_RHS  = *mpRHS;
        // TSystemVectorType& r_Dx = *mpDx;

        // //setting up the list of the DOFs to be solved
        // BuiltinTimer setup_dofs_time;
        // pBuilderSolver->SetUpDofSet(p_scheme, rModelPartPML);
        // KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
        //     << setup_dofs_time.ElapsedSeconds() << std::endl;


        // //shaping correctly the system
        // BuiltinTimer setup_system_time;
        // pBuilderSolver->SetUpSystem(rModelPartPML);
        // KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
        //     << setup_system_time.ElapsedSeconds() << std::endl;

        // //setting up the Vectors involved to the correct size
        // BuiltinTimer system_matrix_resize_time;
        // pBuilderSolver->ResizeAndInitializeVectors(pScheme, r_K, r_RHS, r_Dx, rModelPartPML);

        // if (mpDx->size() != pBuilderSolver->GetEquationSystemSize())
        //     mpDx->resize(pBuilderSolver->GetEquationSystemSize(), false);

        // KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
        //     << system_matrix_resize_time.ElapsedSeconds() << std::endl;
     

        // //set up system matrices

        // //set up the stiffness matrix and rhs
        // TSparseSpace::SetToZero(r_tmp_RHS);
        // pBuilderSolver->Build(pScheme, rModelPartPML, r_K, r_RHS);
        // //DirichletUtility::ApplyDirichletConditions<TSparseSpace>(r_K, r_tmp_RHS, fixed_dofs, 1.0);












        KRATOS_CATCH("")
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
        return "PMLDirectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PMLDirectionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
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


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

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
    ModelPart& mrModelPartPML;
    ModelPart& mrModelPartInterface;
    ModelPart& mrModelPartBoundary;

    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

    // TSystemVectorPointerType mpRHS; /// The RHS vector
    // TSystemVectorPointerType mpDx; /// The solution vector
    // TSystemMatrixPointerType mpK; /// The stiffness matrix (real part)
    

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
    ///@na
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MonolithicMappingProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MonolithicMappingProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MonolithicMappingProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MonolithicMappingProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

};

}
#endif /* KRATOS_MONOLITHIC_MAPPING_PROCESS_H_INCLUDED defined */
