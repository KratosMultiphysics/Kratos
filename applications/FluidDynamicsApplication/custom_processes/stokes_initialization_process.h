#ifndef KRATOS_STOKES_INITIALIZATION_PROCESS_H
#define KRATOS_STOKES_INITIALIZATION_PROCESS_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/residualbased_block_builder_and_solver_periodic.h"

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
class  StokesInitializationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of  StokesInitializationProcess
    KRATOS_CLASS_POINTER_DEFINITION( StokesInitializationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    StokesInitializationProcess(ModelPart& rModelPart,
                                typename TLinearSolver::Pointer pLinearSolver,
                                unsigned int DomainSize,
                                const Variable<int>& PeriodicPairIndicesVar):
        Process(),
        mrReferenceModelPart(rModelPart),
        mpLinearSolver(pLinearSolver),
        mDomainSize(DomainSize)
    {
        KRATOS_TRY;

        if(!mrReferenceModelPart.GetModel().HasModelPart("StokesModelPart"))
            mrReferenceModelPart.GetModel().DeleteModelPart("StokesModelPart");

        ModelPart& r_stokes_part = mrReferenceModelPart.GetModel().CreateModelPart("StokesModelPart");

        r_stokes_part.GetNodalSolutionStepVariablesList() = mrReferenceModelPart.GetNodalSolutionStepVariablesList();
        r_stokes_part.SetBufferSize(1);
        r_stokes_part.SetNodes( mrReferenceModelPart.pNodes() );
        r_stokes_part.SetProcessInfo(mrReferenceModelPart.pGetProcessInfo());
        r_stokes_part.SetProperties(mrReferenceModelPart.pProperties());

        // Retrieve Stokes element model
        std::string ElementName;
        if (mDomainSize == 2)
            ElementName = std::string("StationaryStokes2D");
        else
            ElementName = std::string("StationaryStokes3D");

        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

        // Generate Stokes elements
        for (ModelPart::ElementsContainerType::const_iterator itElem = mrReferenceModelPart.ElementsBegin(); itElem != mrReferenceModelPart.ElementsEnd(); itElem++)
        {
            Element::Pointer pElem = rReferenceElement.Create(itElem->Id(), itElem->GetGeometry(), itElem->pGetProperties() );
            r_stokes_part.Elements().push_back(pElem);
        }

        // pointer types for the solution strategy construcion
        typedef typename Scheme< TSparseSpace, TDenseSpace >::Pointer SchemePointerType;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
        typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

        // Solution scheme: Linear static scheme
        SchemePointerType pScheme = SchemePointerType( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > () );

        // Builder and solver
        BuilderSolverTypePointer pBuildAndSolver;
        pBuildAndSolver = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolverPeriodic<TSparseSpace, TDenseSpace, TLinearSolver > (mpLinearSolver,
                                                                                                                                              PeriodicPairIndicesVar));


        // Strategy
        bool ReactionFlag = false;
        bool ReformDofSetFlag = false;
        bool CalculateNormDxFlag = false;
        bool MoveMeshFlag = false;
        mpSolutionStrategy = StrategyPointerType( 
            new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(r_stokes_part,
                                                                                                                            pScheme,
                                                                                                                            mpLinearSolver,
                                                                                                                            pBuildAndSolver,
                                                                                                                            ReactionFlag,
                                                                                                                            ReformDofSetFlag,
                                                                                                                            CalculateNormDxFlag,
                                                                                                                            MoveMeshFlag) );
        mpSolutionStrategy->SetEchoLevel(0);
        mpSolutionStrategy->Check();

        mIsCleared = false;

        KRATOS_CATCH("");
    }

    /// Destructor.

    ~StokesInitializationProcess() override
    {
        mrReferenceModelPart.GetModel().DeleteModelPart("StokesModelPart");
       // mpSolutionStrategy->Clear();
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Solve a stationary Stokes problem
    void Execute() override
    {
        KRATOS_TRY;

        if (!mIsCleared)
        {
            // Solve Stokes problem.
            this->Solve();

            // Destroy the auxiliary ModelPart and solution structures.
            this->Clear();
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Stokes initialization process was called twice","");
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    void SetConditions(ModelPart& rStokesPart, ModelPart::ConditionsContainerType::Pointer pConditions)
    {
        
        ModelPart& r_stokes_part = mrReferenceModelPart.GetModel().GetModelPart("StokesModelPart");    
        rStokesPart.SetConditions(pConditions);
        r_stokes_part.GetCommunicator().LocalMesh().SetConditions(pConditions);
    }

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
        buffer << " StokesInitializationProcess";
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " StokesInitializationProcess";
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

    ModelPart& mrReferenceModelPart;

    typename TLinearSolver::Pointer mpLinearSolver;

    unsigned int mDomainSize;

    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mpSolutionStrategy;

    bool mIsCleared;

    ///@}
    ///@name Protected Life Cycle
    ///@{

    /// Protected constructor to be used by derived classes
    StokesInitializationProcess(ModelPart& rModelPart,
                                typename TLinearSolver::Pointer pLinearSolver,
                                unsigned int DomainSize,
                                const StokesInitializationProcess* pThis):
        Process(),
        mrReferenceModelPart(rModelPart),
        mpLinearSolver(pLinearSolver),
        mDomainSize(DomainSize)
    {}

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{



    // Solve Stokes problem.
    virtual void Solve()
    {
        KRATOS_TRY;

        mpSolutionStrategy->Solve();

        KRATOS_CATCH("");
    }


    // Liberate memory.
    virtual void Clear()
    {
        mpSolutionStrategy->Clear();
        mIsCleared = true;
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

    StokesInitializationProcess & operator=( StokesInitializationProcess const& rOther)
    {
        return *this;
    }

    /// Copy constructor.

    StokesInitializationProcess( StokesInitializationProcess const& rOther)
    { }


    ///@}

}; // Class  StokesInitializationProcess

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
                                  StokesInitializationProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    return rIStream;
}

/// output stream function

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const  StokesInitializationProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.


#endif // KRATOS_STOKES_INITIALIZATION_PROCESS_H
