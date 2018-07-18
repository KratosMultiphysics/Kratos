#ifndef KRATOS_TRILINOS_STOKES_INITIALIZATION_PROCESS_H
#define KRATOS_TRILINOS_STOKES_INITIALIZATION_PROCESS_H

// System includes
#include <string>
#include <iostream>


// External includes
#include "Epetra_MpiComm.h"

// Project includes
#include "includes/mpi_communicator.h"

// Application includes
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
//#include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver_periodic.h"
#include "custom_utilities/parallel_fill_communicator.h"

#include "../FluidDynamicsApplication/custom_processes/stokes_initialization_process.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

/// A process to provide initial values for Navier-Stokes problems.
/**
 * A stationary stokes problem (using the same boundary conditions as the full Navier-Stokes problem)
 * is solver in order to produce a divergence-free velocity distribution and the corresponding pressure
 * field, which can be used as an initial condition for the Navier-Stokes problem.
 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class  TrilinosStokesInitializationProcess : public StokesInitializationProcess<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  TrilinosStokesInitializationProcess
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosStokesInitializationProcess);

    typedef StokesInitializationProcess<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    TrilinosStokesInitializationProcess(Epetra_MpiComm& rComm,
                                        ModelPart& rModelPart,
                                        typename TLinearSolver::Pointer pLinearSolver,
                                        unsigned int DomainSize,
                                        const Variable<int>& PeriodicPairIndicesVar):
        BaseType(rModelPart,pLinearSolver,DomainSize,this),
        mrComm(rComm),
        mrPeriodicVar(PeriodicPairIndicesVar)
    {
        KRATOS_TRY;

        ModelPart& rReferenceModelPart = BaseType::mrReferenceModelPart;
        typename TLinearSolver::Pointer& pLinearSolver = BaseType::mpLinearSolver;
        unsigned int DomainSize = BaseType::mDomainSize;
        auto& pStokesModelPart = BaseType::mpStokesModelPart;
        typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer& pSolutionStrategy = BaseType::mpSolutionStrategy;

        // Initialize new model part (same nodes, new elements, no conditions)
        auto tmp  = Kratos::make_unique<ModelPart>("StokesModelPart");
        pStokesModelPart.swap(tmp);
        pStokesModelPart->GetNodalSolutionStepVariablesList() = rReferenceModelPart.GetNodalSolutionStepVariablesList();
        pStokesModelPart->SetBufferSize(1);
        pStokesModelPart->SetNodes( rReferenceModelPart.pNodes() );
        pStokesModelPart->SetProcessInfo( rReferenceModelPart.pGetProcessInfo() );
        pStokesModelPart->SetProperties( rReferenceModelPart.pProperties() );

        // Create a communicator for the new model part and copy the partition information about nodes.
        Communicator& rReferenceComm = rReferenceModelPart.GetCommunicator();
        typename Communicator::Pointer pStokesMPIComm = Kratos::make_shared<MPICommunicator>( &(rReferenceModelPart.GetNodalSolutionStepVariablesList()) );
        pStokesMPIComm->SetNumberOfColors( rReferenceComm.GetNumberOfColors() ) ;
        pStokesMPIComm->NeighbourIndices() = rReferenceComm.NeighbourIndices();
        pStokesMPIComm->LocalMesh().SetNodes( rReferenceComm.LocalMesh().pNodes() );
        pStokesMPIComm->InterfaceMesh().SetNodes( rReferenceComm.InterfaceMesh().pNodes() );
        pStokesMPIComm->GhostMesh().SetNodes( rReferenceComm.GhostMesh().pNodes() );
        for (unsigned int i = 0; i < rReferenceComm.GetNumberOfColors(); i++)
        {
            pStokesMPIComm->pInterfaceMesh(i)->SetNodes( rReferenceComm.pInterfaceMesh(i)->pNodes() );
            pStokesMPIComm->pLocalMesh(i)->SetNodes( rReferenceComm.pLocalMesh(i)->pNodes() );
            pStokesMPIComm->pGhostMesh(i)->SetNodes( rReferenceComm.pGhostMesh(i)->pNodes() );
        }
        pStokesModelPart->SetCommunicator( pStokesMPIComm );

        // Retrieve Stokes element model
        std::string ElementName;
        if (DomainSize == 2)
            ElementName = std::string("StationaryStokes2D");
        else
            ElementName = std::string("StationaryStokes3D");

        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

        // Generate Stokes elements
        for (ModelPart::ElementsContainerType::iterator itElem = rReferenceModelPart.ElementsBegin(); itElem != rReferenceModelPart.ElementsEnd(); itElem++)
        {
            Element::Pointer pElem = rReferenceElement.Create(itElem->Id(), itElem->GetGeometry(), itElem->pGetProperties() );
            pStokesModelPart->Elements().push_back(pElem);
            pStokesModelPart->GetCommunicator().LocalMesh().Elements().push_back(pElem);
        }

        // Solution scheme: Linear static scheme
        auto pScheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > >();

        // Builder and solver
        int guess_row_size;
        if(DomainSize == 2) guess_row_size = 15;
        else guess_row_size = 40;

        auto pBuildAndSolver = Kratos::make_shared<TrilinosBlockBuilderAndSolverPeriodic<TSparseSpace,TDenseSpace,TLinearSolver> >(mrComm,guess_row_size,pLinearSolver,mrPeriodicVar);

        // Strategy
        bool ReactionFlag = false;
        bool ReformDofSetFlag = false;
        bool CalculateNormDxFlag = false;
        bool MoveMeshFlag = false;
        pSolutionStrategy = Kratos::make_shared< ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(*pStokesModelPart,
                                                                                                                         pScheme,
                                                                                                                         pLinearSolver,
                                                                                                                         pBuildAndSolver,
                                                                                                                         ReactionFlag,
                                                                                                                         ReformDofSetFlag,
                                                                                                                         CalculateNormDxFlag,
                                                                                                                         MoveMeshFlag);


        pSolutionStrategy->SetEchoLevel(0);
        pSolutionStrategy->Check();

        BaseType::mIsCleared = false;

        KRATOS_CATCH("");
    }

    /// Destructor.

    virtual ~TrilinosStokesInitializationProcess()
    {
        BaseType::mpSolutionStrategy->Clear();
    }


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

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << " TrilinosStokesInitializationProcess";
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " TrilinosStokesInitializationProcess";
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

    Epetra_MpiComm& mrComm;

    const Variable<int>& mrPeriodicVar;

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

    TrilinosStokesInitializationProcess & operator=( TrilinosStokesInitializationProcess const& rOther)
    {
        return *this;
    }

    /// Copy constructor.

    TrilinosStokesInitializationProcess( TrilinosStokesInitializationProcess const& rOther)
    { }


    ///@}

}; // Class  TrilinosStokesInitializationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
inline std::istream & operator >>(std::istream& rIStream,
                                  TrilinosStokesInitializationProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    return rIStream;
}

/// output stream function

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const  TrilinosStokesInitializationProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TRILINOS_STOKES_INITIALIZATION_PROCESS_H
