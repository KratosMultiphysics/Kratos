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
#include "processes/compute_nodal_gradient_process.h"
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
    typedef std::size_t SizeType;

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;

    /// The definitionof the geometry
    typedef Geometry<NodeType> GeometryType;

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

        
        //  for ( ModelPart::NodesContainerType::iterator it_node = rPMLNodes.begin();
        //         it_node != rPMLNodes.end() ; ++it_node)
        // {
        //     std::cout<<" Potential "<< it_node->FastGetSolutionStepValue(PRESSURE) <<std::endl;
        // }


        typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable> GradientType;
        GradientType process = GradientType(mrModelPartPML, PRESSURE, PML_IMAG_DISTANCE, NODAL_AREA, false);
        process.Execute();

        for ( ModelPart::NodesContainerType::iterator it_node = rPMLNodes.begin();
                    it_node != rPMLNodes.end() ; ++it_node)
            {
                std::cout<<" PML direction "<< it_node->GetValue(PML_IMAG_DISTANCE) <<std::endl;
            }








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
