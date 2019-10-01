//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2011-07-22 17:06:00 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_TRILINOS_SPALART_ALLMARAS_H_INCLUDED )
#define  KRATOS_TRILINOS_SPALART_ALLMARAS_H_INCLUDED

// System includes

// External includes
#include "Epetra_MpiComm.h"

// Project includes
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"

// Application includes
#include "custom_strategies/builder_and_solvers/trilinos_elimination_builder_and_solver.h"
#include "custom_strategies/schemes/trilinos_residualbased_incremental_aitken_static_scheme.h"
#include "mpi/includes/mpi_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"

#include "../FluidDynamicsApplication/custom_processes/spalart_allmaras_turbulence_model.h"

namespace Kratos
{
///@addtogroup TrilinosApplication
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

/// Trilinos implementation of the Spalart-Allmaras turbulence model.
/**  To be used in combination with Fluid or VMS elements. See SpalartAllmarasTurbulenceModel
 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class TrilinosSpalartAllmarasTurbulenceModel : public SpalartAllmarasTurbulenceModel<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TrilinosSpalartAllmarasTurbulenceModel
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosSpalartAllmarasTurbulenceModel);

    typedef SpalartAllmarasTurbulenceModel<TSparseSpace,TDenseSpace,TLinearSolver> BaseSpAlType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TrilinosSpalartAllmarasTurbulenceModel(Epetra_MpiComm& rComm,
                                           ModelPart& rModelPart,
                                           typename TLinearSolver::Pointer pLinearSolver,
                                           unsigned int DomainSize,
                                           double NonLinearTol,
                                           unsigned int MaxIter,
                                           bool ReformDofSet,
                                           unsigned int TimeOrder)
        :
        BaseSpAlType(rModelPart)
    {
        KRATOS_TRY;

        BaseSpAlType::mdomain_size = DomainSize;
        BaseSpAlType::mtol = NonLinearTol;
        BaseSpAlType::mmax_it = MaxIter;
        BaseSpAlType::mtime_order = TimeOrder;
        BaseSpAlType::madapt_for_fractional_step = false;

        //************************************************************************************************
        //check that the variables needed are in the model part
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(DISTANCE)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", DISTANCE);
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(VELOCITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", VELOCITY);
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(MOLECULAR_VISCOSITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", MOLECULAR_VISCOSITY);
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(TURBULENT_VISCOSITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", TURBULENT_VISCOSITY);
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", MESH_VELOCITY);
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(VISCOSITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", VISCOSITY);
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", NODAL_AREA);
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(TEMP_CONV_PROJ)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", TEMP_CONV_PROJ);
        if (!(BaseSpAlType::mr_model_part.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", PARTITION_INDEX);

        if (BaseSpAlType::mr_model_part.GetBufferSize() < 3)
            KRATOS_THROW_ERROR(std::logic_error, "insufficient buffer size for BDF2, currently buffer size is ", BaseSpAlType::mr_model_part.GetBufferSize());

        //************************************************************************************************
        //construct a new auxiliary model part
        BaseSpAlType::mrSpalartModelPart.SetBufferSize(3);
        BaseSpAlType::mrSpalartModelPart.GetNodalSolutionStepVariablesList() = BaseSpAlType::mr_model_part.GetNodalSolutionStepVariablesList();
        BaseSpAlType::mrSpalartModelPart.SetBufferSize(BaseSpAlType::mr_model_part.GetBufferSize());
        BaseSpAlType::mrSpalartModelPart.SetNodes(BaseSpAlType::mr_model_part.pNodes());
        BaseSpAlType::mrSpalartModelPart.SetProcessInfo(BaseSpAlType::mr_model_part.pGetProcessInfo());
        BaseSpAlType::mrSpalartModelPart.SetProperties(BaseSpAlType::mr_model_part.pProperties());

        // Create a communicator for the new model part and copy the partition information about nodes.
        Communicator& rReferenceComm = BaseSpAlType::mr_model_part.GetCommunicator();
        typename Communicator::Pointer pSpalartMPIComm = typename Communicator::Pointer( new MPICommunicator( &(BaseSpAlType::mr_model_part.GetNodalSolutionStepVariablesList()) ) );
        pSpalartMPIComm->SetNumberOfColors( rReferenceComm.GetNumberOfColors() ) ;
        pSpalartMPIComm->NeighbourIndices() = rReferenceComm.NeighbourIndices();
        pSpalartMPIComm->LocalMesh().SetNodes( rReferenceComm.LocalMesh().pNodes() );
        pSpalartMPIComm->InterfaceMesh().SetNodes( rReferenceComm.InterfaceMesh().pNodes() );
        pSpalartMPIComm->GhostMesh().SetNodes( rReferenceComm.GhostMesh().pNodes() );
        for (unsigned int i = 0; i < rReferenceComm.GetNumberOfColors(); i++)
        {
            pSpalartMPIComm->pInterfaceMesh(i)->SetNodes( rReferenceComm.pInterfaceMesh(i)->pNodes() );
            pSpalartMPIComm->pLocalMesh(i)->SetNodes( rReferenceComm.pLocalMesh(i)->pNodes() );
            pSpalartMPIComm->pGhostMesh(i)->SetNodes( rReferenceComm.pGhostMesh(i)->pNodes() );
        }
        BaseSpAlType::mrSpalartModelPart.SetCommunicator( pSpalartMPIComm );

        std::string ElementName;
        if (BaseSpAlType::mdomain_size == 2)
            ElementName = std::string("SpalartAllmaras2D");
        else
            ElementName = std::string("SpalartAllmaras3D");

        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

        //generating the elements
        for (ModelPart::ElementsContainerType::iterator iii = BaseSpAlType::mr_model_part.ElementsBegin(); iii != BaseSpAlType::mr_model_part.ElementsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();
            Element::Pointer p_element = rReferenceElement.Create(iii->Id(), iii->GetGeometry(), properties);
            BaseSpAlType::mrSpalartModelPart.Elements().push_back(p_element);
        }

        std::string ConditionName;
        if (BaseSpAlType::mdomain_size == 2)
            ConditionName = std::string("Condition2D");
        else
            ConditionName = std::string("Condition3D");
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(ConditionName);

        for (ModelPart::ConditionsContainerType::iterator iii = BaseSpAlType::mr_model_part.ConditionsBegin(); iii != BaseSpAlType::mr_model_part.ConditionsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();
            Condition::Pointer p_condition = rReferenceCondition.Create(iii->Id(), iii->GetGeometry(), properties);
            BaseSpAlType::mrSpalartModelPart.Conditions().push_back(p_condition);
        }

        // Create a communicator for the new model part
        ParallelFillCommunicator CommunicatorGeneration(BaseSpAlType::mrSpalartModelPart);
        CommunicatorGeneration.Execute();
        //CommunicatorGeneration.PrintDebugInfo()

        // pointer types for the solution strategy construcion
        typedef typename Scheme< TSparseSpace, TDenseSpace >::Pointer SchemePointerType;
        typedef typename ConvergenceCriteria< TSparseSpace, TDenseSpace >::Pointer ConvergenceCriteriaPointerType;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
        typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

        // Solution scheme: Aitken iterations
        const double DefaultAitkenOmega = 1.0;
        SchemePointerType pScheme = SchemePointerType( new TrilinosResidualBasedIncrementalAitkenStaticScheme< TSparseSpace, TDenseSpace > (DefaultAitkenOmega) );

        // Convergence criteria
        const double NearlyZero = 1.0e-20;
        ConvergenceCriteriaPointerType pConvCriteria = ConvergenceCriteriaPointerType( new ResidualCriteria<TSparseSpace,TDenseSpace>(NonLinearTol,NearlyZero) );

        //definitions for trilinos
        int guess_row_size;
        if(BaseSpAlType::mdomain_size == 2) guess_row_size = 15;
        else guess_row_size = 40;

        // Builder and solver
        BuilderSolverTypePointer pBuildAndSolver = BuilderSolverTypePointer(new TrilinosResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> (rComm,guess_row_size,pLinearSolver));

        // Strategy
        bool CalculateReactions = false;
        bool MoveMesh = false;

        BaseSpAlType::mpSolutionStrategy = StrategyPointerType( new ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(BaseSpAlType::mrSpalartModelPart,pScheme,pLinearSolver,pConvCriteria,pBuildAndSolver,MaxIter,CalculateReactions,ReformDofSet,MoveMesh));
        BaseSpAlType::mpSolutionStrategy->SetEchoLevel(0);
        BaseSpAlType::mpSolutionStrategy->Check();

        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~TrilinosSpalartAllmarasTurbulenceModel()
    {}

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
        buffer << "TrilinosSpalartAllmarasTurbulenceModel";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TrilinosSpalartAllmarasTurbulenceModel";
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

    /// Initialize a Trilinos Solution Strategy for the Spalart Allmaras turbulence model
    virtual void SolutionStrategyConfiguration(ModelPart& rSpalartModelPart,
            typename TLinearSolver::Pointer pNewLinearSolver,
            bool reform_dofset)
    {

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
    TrilinosSpalartAllmarasTurbulenceModel & operator=(TrilinosSpalartAllmarasTurbulenceModel const& rOther)
    {
        BaseSpAlType::operator=(rOther);
        return *this;
    }

    /// Copy constructor.
    TrilinosSpalartAllmarasTurbulenceModel(TrilinosSpalartAllmarasTurbulenceModel const& rOther)
        : BaseSpAlType(rOther)
    {
    }


    ///@}

}; // Class SpalartAllmarasTurbulenceModel

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
                                  TrilinosSpalartAllmarasTurbulenceModel<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const TrilinosSpalartAllmarasTurbulenceModel<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TRILINOS_SPALART_ALLMARAS_H_INCLUDED  defined


