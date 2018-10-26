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

/* Trilinos includes */
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_utilities/parallel_fill_communicator.h"

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

        if(mrReferenceModelPart.GetOwnerModel().HasModelPart("StructuralMeshMovingPart"))
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
        mrReferenceModelPart.GetOwnerModel().DeleteModelPart("StructuralMeshMovingPart");
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

        // Update FEM database
        if (m_calculate_mesh_velocities == true)
            CalculateMeshVelocities();

        MoveMesh();

        // Clearing the system if needed
        if (m_reform_dof_set_at_each_step == true)
            mstrategy->Clear();

        return 0.0;

        KRATOS_CATCH("")
    }

    void CalculateMeshVelocities()
    {
        KRATOS_TRY;

        double delta_time = BaseType::GetModelPart().GetProcessInfo()[DELTA_TIME];

        KRATOS_ERROR_IF(delta_time <= 0.0)<< "Invalid DELTA_TIME." << std::endl;

        double coeff = 1 / delta_time;
        if (m_time_order == 1) // mesh velocity calculated as (x(n+1)-x(n))/Dt
        {
            for (ModelPart::NodeIterator i = (*mpmesh_model_part).NodesBegin();
                 i != (*mpmesh_model_part).NodesEnd();
                 ++i)
            {
                array_1d<double, 3>& mesh_v = (i)->FastGetSolutionStepValue(MESH_VELOCITY);
                array_1d<double, 3>& disp =
                    (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
                array_1d<double, 3>& dispold =
                    (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
                noalias(mesh_v) = disp - dispold;
                mesh_v *= coeff;
            }
        }
        else if (m_time_order == 2) // mesh velocity calculated as
        // (3*x(n+1)-4*x(n)+x(n-1))/(2*Dt)
        {
            double c1 = 1.50 * coeff;
            double c2 = -2.0 * coeff;
            double c3 = 0.50 * coeff;

            for (ModelPart::NodeIterator i = (*mpmesh_model_part).NodesBegin();
                 i != (*mpmesh_model_part).NodesEnd();
                 ++i)
            {
                array_1d<double, 3>& mesh_v = (i)->FastGetSolutionStepValue(MESH_VELOCITY);
                noalias(mesh_v) = c1 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT);
                noalias(mesh_v) +=
                    c2 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 1);
                noalias(mesh_v) +=
                    c3 * (i)->FastGetSolutionStepValue(MESH_DISPLACEMENT, 2);
            }
        }
        else{
        KRATOS_ERROR << "Wrong TimeOrder: Acceptable values are: 1 and 2"
                   << std::endl;
        }

        KRATOS_CATCH("")
    }


    void MoveMesh() override
    {
        for (ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin();
             i != BaseType::GetModelPart().NodesEnd();
             ++i)
        {
            (i)->X() = (i)->X0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_X);
            (i)->Y() = (i)->Y0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_Y);
            (i)->Z() = (i)->Z0() + i->GetSolutionStepValue(MESH_DISPLACEMENT_Z);
        }
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
        if(!mrReferenceModelPart.GetOwnerModel().HasModelPart("StructuralMeshMovingPart"))
            mrReferenceModelPart.GetOwnerModel().DeleteModelPart("StructuralMeshMovingPart");

        mpmesh_model_part  = &mrReferenceModelPart.GetOwnerModel().CreateModelPart("StructuralMeshMovingPart");

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
