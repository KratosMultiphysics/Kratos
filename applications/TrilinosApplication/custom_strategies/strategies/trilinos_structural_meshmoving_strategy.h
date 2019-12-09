//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

#if !defined(KRATOS_TRILINOS_STRUCTURAL_MESHMOVING_STRATEGY)
#define KRATOS_TRILINOS_STRUCTURAL_MESHMOVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/model_part.h"
#include "containers/model.h"
#include "includes/mesh_moving_variables.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "mpi/utilities/parallel_fill_communicator.h"

/* Trilinos includes */
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"

namespace Kratos
{
/**@name Kratos Globals */
/*@{ */

/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */

/**@name  Enum's */
/*@{ */

/*@} */
/**@name  Functions */
/*@{ */

/*@} */
/**@name Kratos Classes */
/*@{ */

/// Short class definition.
/**   Detail class definition.

 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class TrilinosStructuralMeshMovingStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosStructuralMeshMovingStrategy);

    typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef Variable<array_1d<double, 3>> VariableWithComponentsType;

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    TrilinosStructuralMeshMovingStrategy(Epetra_MpiComm& Communicator,
                                         ModelPart& model_part,
                                         typename TLinearSolver::Pointer pNewLinearSolver,
                                         int TimeOrder = 2,
                                         bool ReformDofSetAtEachStep = false,
                                         bool ComputeReactions = false,
                                         bool CalculateMeshVelocities = true,
                                         int EchoLevel = 0)
        :
        SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part),
        mrReferenceModelPart(model_part)
    {
        KRATOS_TRY

        if(mrReferenceModelPart.GetModel().HasModelPart(mrReferenceModelPart.Name()+"_StructuralMeshMovingPart"))
            KRATOS_ERROR << "StructuralMeshMovingPart already existing when constructing TrilinosLaplacianMeshMovingStrategy";

        // Passed variables
        m_reform_dof_set_at_each_step = ReformDofSetAtEachStep;
        m_compute_reactions = ComputeReactions;
        m_calculate_mesh_velocities = CalculateMeshVelocities;
        m_echo_level = EchoLevel;
        m_time_order = TimeOrder;
        bool calculate_norm_dx_flag = false;

        // Definitions for trilinos
        int guess_row_size;
        guess_row_size = 15;

        // Generating Mesh Part
        GenerateMeshPart();

        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(
            new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>());

        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
        BuilderSolverTypePointer builderSolver = BuilderSolverTypePointer(
            new TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(
                Communicator, guess_row_size, pNewLinearSolver));

        mstrategy = typename BaseType::Pointer(
            new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
                *mpmesh_model_part,
                pscheme,
                pNewLinearSolver,
                builderSolver,
                m_compute_reactions,
                m_reform_dof_set_at_each_step,
                calculate_norm_dx_flag));

        mstrategy->SetEchoLevel(m_echo_level);

        KRATOS_CATCH("")
    }

    /** Destructor.
     */
    virtual ~TrilinosStructuralMeshMovingStrategy()
    {
        mrReferenceModelPart.GetModel().DeleteModelPart(mrReferenceModelPart.Name()+"_StructuralMeshMovingPart");
    }

    /** Destructor.
     */

    void Initialize() override
    {}

    double Solve() override
    {
        KRATOS_TRY;

        // Setting mesh to initial configuration
        for (ModelPart::NodeIterator i = (*mpmesh_model_part).NodesBegin();
             i != (*mpmesh_model_part).NodesEnd();
             ++i)
        {
            (i)->X() = (i)->X0();
            (i)->Y() = (i)->Y0();
            (i)->Z() = (i)->Z0();
        }

        // Solve for mesh movement
        mstrategy->Solve();

        // Clearing the system if needed
        if (m_reform_dof_set_at_each_step == true)
            mstrategy->Clear();

        return 0.0;

        KRATOS_CATCH("")
    }

    /*@} */
    /**@name Operators
     */
    /*@{ */

    /*@} */
    /**@name Operations */
    /*@{ */

    /*@} */
    /**@name Access */
    /*@{ */

    /*@} */
    /**@name Inquiry */
    /*@{ */

    /*@} */
    /**@name Friends */
    /*@{ */

    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /*@} */
    /**@name Protected Inquiry */
    /*@{ */

    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */

    /*@} */
    /**@name Member Variables */
    /*@{ */
    ModelPart& mrReferenceModelPart;
    ModelPart* mpmesh_model_part;

    typename BaseType::Pointer mstrategy;

    int m_echo_level;
    int m_time_order;
    bool m_reform_dof_set_at_each_step;
    bool m_compute_reactions;
    bool m_calculate_mesh_velocities;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */

    void GenerateMeshPart()
    {
        if(!mrReferenceModelPart.GetModel().HasModelPart(mrReferenceModelPart.Name()+"_StructuralMeshMovingPart"))
            mrReferenceModelPart.GetModel().DeleteModelPart(mrReferenceModelPart.Name()+"_StructuralMeshMovingPart");

        mpmesh_model_part  = &mrReferenceModelPart.GetModel().CreateModelPart(mrReferenceModelPart.Name()+"_StructuralMeshMovingPart");

        // Initializing mesh nodes
        mpmesh_model_part->Nodes() = BaseType::GetModelPart().Nodes();

        // Creating mesh elements
        ModelPart::ElementsContainerType& MeshElems = mpmesh_model_part->Elements();
        // this works for multiple geometries because we use the
        // element's Create() member function that accepts a geometry
        // pointer
        const Element& rElem =
            KratosComponents<Element>::Get("StructuralMeshMovingElement3D4N");

        for (auto it = BaseType::GetModelPart().ElementsBegin();
             it != BaseType::GetModelPart().ElementsEnd();
             ++it)
            MeshElems.push_back(rElem.Create(
                it->Id(), it->pGetGeometry(), it->pGetProperties()));

        // Optimize communicaton plan
        ParallelFillCommunicator CommunicatorGeneration(*mpmesh_model_part);
        CommunicatorGeneration.Execute();
    }

    /*@} */
    /**@name Private  Access */
    /*@{ */

    /*@} */
    /**@name Private Inquiry */
    /*@{ */

    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /** Copy constructor.
     */
    TrilinosStructuralMeshMovingStrategy(const TrilinosStructuralMeshMovingStrategy& Other);

    /*@} */

}; /* Class TrilinosStructuralMeshMovingStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */
}
/* namespace Kratos.*/

#endif /* KRATOS_STRUCTURAL_MESHMOVING_STRATEGY  defined */
