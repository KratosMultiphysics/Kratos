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
#include <string>
#include <iostream>


// External includes
#include "Epetra_MpiComm.h"

// Project includes
#include "includes/define.h"
#include "includes/communicator.h"
#include "includes/mpi_communicator.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver.h"
#include "custom_utilities/parallel_fill_communicator.h"

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

        typedef SpalartAllmarasTurbulenceModel<TSparseSpace,TDenseSpace,TLinearSolver> BaseTurbulenceModelType;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor
        TrilinosSpalartAllmarasTurbulenceModel(Epetra_MpiComm& rComm,
                                               ModelPart& rModelPart,
                                               typename TLinearSolver::Pointer pNewLinearSolver,
                                               unsigned int DomainSize,
                                               double NonLinearTol,
                                               unsigned int MaxIter,
                                               bool ReformDofSet,
                                               unsigned int TimeOrder)
        :
        BaseTurbulenceModelType(rModelPart)
        {
            BaseTurbulenceModelType::mdomain_size = DomainSize;
            BaseTurbulenceModelType::mtol = NonLinearTol;
            BaseTurbulenceModelType::mmax_it = MaxIter;
            BaseTurbulenceModelType::mtime_order = TimeOrder;
            BaseTurbulenceModelType::madapt_for_fractional_step = false;

            //************************************************************************************************
            //check that the variables needed are in the model part
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(DISTANCE)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", DISTANCE);
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(VELOCITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", VELOCITY);
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(MOLECULAR_VISCOSITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", MOLECULAR_VISCOSITY);
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(TURBULENT_VISCOSITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", TURBULENT_VISCOSITY);
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", MESH_VELOCITY);
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(VISCOSITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", VISCOSITY);
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", NODAL_AREA);
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(TEMP_CONV_PROJ)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", TEMP_CONV_PROJ);
            if (!(BaseTurbulenceModelType::mr_model_part.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", PARTITION_INDEX);

            if (BaseTurbulenceModelType::mr_model_part.GetBufferSize() < 3)
                KRATOS_ERROR(std::logic_error, "insufficient buffer size for BDF2, currently buffer size is ", BaseTurbulenceModelType::mr_model_part.GetBufferSize())

            //************************************************************************************************
            //construct a new auxiliary model part
            BaseTurbulenceModelType::mspalart_model_part.SetBufferSize(3);
            BaseTurbulenceModelType::mspalart_model_part.SetBufferSize(BaseTurbulenceModelType::mr_model_part.GetBufferSize());
            BaseTurbulenceModelType::mspalart_model_part.SetNodes(BaseTurbulenceModelType::mr_model_part.pNodes());
            BaseTurbulenceModelType::mspalart_model_part.SetProcessInfo(BaseTurbulenceModelType::mr_model_part.pGetProcessInfo());
            BaseTurbulenceModelType::mspalart_model_part.SetProperties(BaseTurbulenceModelType::mr_model_part.pProperties());

            typename Communicator::Pointer pSpalartMPIComm = typename Communicator::Pointer( new MPICommunicator() );
            BaseTurbulenceModelType::mspalart_model_part.SetCommunicator( pSpalartMPIComm );

            std::string ElementName;
            if (BaseTurbulenceModelType::mdomain_size == 2)
                ElementName = std::string("SpalartAllmaras2D");
            else
                ElementName = std::string("SpalartAllmaras3D");

            const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

            //generating the elements
            for (ModelPart::ElementsContainerType::iterator iii = BaseTurbulenceModelType::mr_model_part.ElementsBegin(); iii != BaseTurbulenceModelType::mr_model_part.ElementsEnd(); iii++)
            {
                Properties::Pointer properties = iii->pGetProperties();
                Element::Pointer p_element = rReferenceElement.Create(iii->Id(), iii->GetGeometry(), properties);
                BaseTurbulenceModelType::mspalart_model_part.Elements().push_back(p_element);
            }

            const Condition& rReferenceCondition = KratosComponents<Condition>::Get("Condition2D");
            for (ModelPart::ConditionsContainerType::iterator iii = BaseTurbulenceModelType::mr_model_part.ConditionsBegin(); iii != BaseTurbulenceModelType::mr_model_part.ConditionsEnd(); iii++)
            {
                Properties::Pointer properties = iii->pGetProperties();
                Condition::Pointer p_condition = rReferenceCondition.Create(iii->Id(), iii->GetGeometry(), properties);
                BaseTurbulenceModelType::mspalart_model_part.Conditions().push_back(p_condition);
            }

            // Create a communicator for the new model part
            ParallelFillCommunicator CommunicatorGeneration(BaseTurbulenceModelType::mspalart_model_part);
            CommunicatorGeneration.Execute();
            //CommunicatorGeneration.PrintDebugInfo()
            
            //initializing solution strategy
            typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
            typename SchemeType::Pointer pScheme = typename SchemeType::Pointer(new TrilinosResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

            bool CalculateReactions = false;
            bool CalculateNormDxFlag = true;

            // Builder and Solver definition
            typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
            typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

            //definitions for trilinos
            int guess_row_size;
            if(BaseTurbulenceModelType::mdomain_size == 2) guess_row_size = 15;
            else guess_row_size = 40;

            BuilderSolverTypePointer pBaS = BuilderSolverTypePointer(
                    new TrilinosResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> (rComm,guess_row_size,pNewLinearSolver));

            BaseTurbulenceModelType::mpstep1 = StrategyPointerType(
                    new ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver> (BaseTurbulenceModelType::mspalart_model_part,pScheme,pNewLinearSolver,pBaS,CalculateReactions,ReformDofSet,CalculateNormDxFlag));
            BaseTurbulenceModelType::mpstep1->SetEchoLevel(0);
            BaseTurbulenceModelType::mpstep1->Check();
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
        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "TrilinosSpalartAllmarasTurbulenceModel";
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "TrilinosSpalartAllmarasTurbulenceModel";
        }

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
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
            BaseTurbulenceModelType::operator=(rOther);
            return *this;
        }

        /// Copy constructor.
        TrilinosSpalartAllmarasTurbulenceModel(TrilinosSpalartAllmarasTurbulenceModel const& rOther)
        : BaseTurbulenceModelType(rOther)
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


