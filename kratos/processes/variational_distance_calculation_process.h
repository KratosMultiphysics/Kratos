//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//


#if !defined(KRATOS_VARIATIONAL_DISTANCE_CALCULATION_PROCESS_INCLUDED )
#define  KRATOS_VARIATIONAL_DISTANCE_CALCULATION_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"



// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"

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
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and recomputes a signed distance function
mantaining as much as possible the position of the zero of the function prior to the call.

This is achieved by minimizing the function  ( 1 - norm( gradient( distance ) )**2
with the restriction that "distance" is a finite elment function
*/

template< unsigned int TDim,
          class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver >
class VariationalDistanceCalculationProcess
    : public Process
{
public:
    
    KRATOS_DEFINE_LOCAL_FLAG(PERFORM_STEP1);
    KRATOS_DEFINE_LOCAL_FLAG(DO_EXPENSIVE_CHECKS);
    
    ///@name Type Definitions
    ///@{

    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;

///@}
    ///@name Pointer Definitions
    /// Pointer definition of VariationalDistanceCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(VariationalDistanceCalculationProcess);




    ///@}
    ///@name Life Cycle
    ///@{

    /**This process recomputed the distance function mantaining the zero of the existing distance distribution
     * for this reason the DISTANCE should be initialized to values distinct from zero in at least some portions of the domain
     * alternatively, the DISTANCE shall be fixed to zero at least on some nodes, and the process will compute a positive distance
     * respecting that zero
     * @param base_model_parr - is the model part on the top of which the calculation will be performed
     * @param plinear_solver  - linear solver to be used internally
     * @max_iterations        - maximum number of iteration to be employed in the nonlinear optimization process.
     *                        - can also be set to 0 if a (very) rough approximation is enough
     *
     * EXAMPLE OF USAGE FROM PYTHON:
     *
     class distance_linear_solver_settings:
         solver_type = "AMGCL"
         tolerance = 1E-3
         max_iteration = 200
         scaling = False
         krylov_type = "CG"
         smoother_type = "SPAI0"
         verbosity = 0

     import linear_solver_factory
     distance_linear_solver = linear_solver_factory.ConstructSolver(distance_linear_solver_settings)

     max_iterations=1
     distance_calculator = VariationalDistanceCalculationProcess2D(fluid_model_part, distance_linear_solver, max_iterations)
     distance_calculator.Execute()
     */
    VariationalDistanceCalculationProcess(ModelPart& base_model_part,
                                          typename TLinearSolver::Pointer plinear_solver,
                                          unsigned int max_iterations = 10
                                         )
        :mr_base_model_part(base_model_part)
    {
        KRATOS_TRY


        mmax_iterations = max_iterations;

        mdistance_part_is_initialized = false; //this will be set to true upon completing ReGenerateDistanceModelPart

        //check that there is at least one element and node in the model
        if(base_model_part.Nodes().size() == 0) KRATOS_THROW_ERROR(std::logic_error, "the model has no Nodes","");
        if(base_model_part.Elements().size() == 0) KRATOS_THROW_ERROR(std::logic_error, "the model has no Elements","");
        if(base_model_part.NodesBegin()->SolutionStepsDataHas(DISTANCE) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing DISTANCE variable on solution step data","");
        if(TDim == 2)
        {
            if(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle)
                KRATOS_THROW_ERROR(std::logic_error, "In 2D the element type is expected to be a triangle","");
        }
        else if(TDim == 3)
        {
            if(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra)
                KRATOS_THROW_ERROR(std::logic_error, "In 3D the element type is expected to be a tetrahedra","");
        }

        //generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        ReGenerateDistanceModelPart(base_model_part);

        //generate a linear strategy


        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace >() );
        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        BuilderSolverTypePointer pBuilderSolver = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
        mp_solving_strategy = typename SolvingStrategyType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(*mp_distance_model_part,pscheme,plinear_solver,pBuilderSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag) );

        //TODO: check flag DO_EXPENSIVE_CHECKS
        mp_solving_strategy->Check();

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~VariationalDistanceCalculationProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY;

        if(mdistance_part_is_initialized == false)
            ReGenerateDistanceModelPart(mr_base_model_part);

        //TODO: check flag    PERFORM_STEP1
        //step1 - solve a poisson problem with a source term which depends on the sign of the existing distance function
        mp_distance_model_part->pGetProcessInfo()->SetValue(FRACTIONAL_STEP,1);
        
        //unfix the distances
        const int nnodes = static_cast<int>(mp_distance_model_part->Nodes().size());
        #pragma omp parallel for
        for(int iii=0; iii<nnodes; iii++)
        {
            auto it = mp_distance_model_part->NodesBegin() + iii;
            it->Free(DISTANCE);
            
            
            double& d = it->FastGetSolutionStepValue(DISTANCE);
            it->SetValue(DISTANCE, d); //save the distances
              
            if(d == 0){
                it->Fix(DISTANCE);
                d = 1e-15;
            }
            else{
                if(d > 0.0)
                    d = 1.0e9;
                else
                    d = -1.0e9;
            }
            
            
        }
        
        
        const int nelem = static_cast<int>(mp_distance_model_part->Elements().size());
        #pragma omp parallel for
        for(int iii=0; iii<nelem; iii++)
        {
            auto it = mp_distance_model_part->ElementsBegin() + iii;
            array_1d<double,TDim+1> distances;
            auto& geom = it->GetGeometry();
            
            for(unsigned int i=0; i<TDim+1; i++)
                distances[i] = geom[i].GetValue(DISTANCE);
            
            array_1d<double,TDim+1> original_distances = distances;
            
            unsigned int positives = 0, negatives=0;
            for(unsigned int i=0; i<TDim+1; i++)
            {
                if(distances[i] >= 0) positives++;
                else negatives++;
            }
        
            if(positives> 0  && negatives>0) //the element is cut by the interface
            {
                if(TDim==3)
                    GeometryUtils::CalculateTetrahedraDistances(geom, distances);
                else
                    GeometryUtils::CalculateTriangleDistances(geom,distances);
                
                //assign the sign 
                for(unsigned int i=0; i<TDim+1; i++)
                {
                    if(original_distances[i] < 0) distances[i] = -distances[i];
                }
               
                for(unsigned int i=0; i<TDim+1; i++)
                {
                    
                    double& d = geom[i].FastGetSolutionStepValue(DISTANCE);
                    if(std::abs(d) > std::abs(distances[i])) d = distances[i];
                    geom[i].Fix(DISTANCE);
                }
            }
        }
        
        //
        if(mp_distance_model_part->GetCommunicator().TotalProcesses() != 1) //MPI case
        {
            #pragma omp parallel for
            for(int iii=0; iii<nnodes; iii++)
            {
                auto it = mp_distance_model_part->NodesBegin() + iii;
                it->FastGetSolutionStepValue(DISTANCE) = std::abs(it->FastGetSolutionStepValue(DISTANCE));
            }
            
            mp_distance_model_part->GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);
            
            #pragma omp parallel for
            for(int iii=0; iii<nnodes; iii++)
            {
                auto it = mp_distance_model_part->NodesBegin() + iii;
                if(it->GetValue(DISTANCE) < 0.0)
                    it->FastGetSolutionStepValue(DISTANCE) = -it->FastGetSolutionStepValue(DISTANCE);
            }
        }
        
        
        mp_solving_strategy->Solve();


        //step2 - minimize the target residual
        mp_distance_model_part->pGetProcessInfo()->SetValue(FRACTIONAL_STEP,2);
        for(unsigned int it = 0; it<mmax_iterations; it++)
        {
             mp_solving_strategy->Solve();
        }
        
        //unfix the distances
        for(auto it = mp_distance_model_part->NodesBegin(); it!=mp_distance_model_part->NodesEnd(); ++it)
            it->Free(DISTANCE);

        KRATOS_CATCH("")
    }

    virtual void Clear()
    {
        mp_distance_model_part->Nodes().clear();
        mp_distance_model_part->Conditions().clear();
        mp_distance_model_part->Elements().clear();
//        mp_distance_model_part->GetProcessInfo().clear();
        mdistance_part_is_initialized = false;

        mp_solving_strategy->Clear();

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
        return "VariationalDistanceCalculationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "VariationalDistanceCalculationProcess";
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

    /// Minimal constructor for derived classes
    VariationalDistanceCalculationProcess(ModelPart &base_model_part, unsigned int max_iterations):
        mr_base_model_part(base_model_part)
    {
        mdistance_part_is_initialized = false;
        mmax_iterations = max_iterations;
    }


    ///@}
    ///@name Protected member Variables
    ///@{
    bool mdistance_part_is_initialized;
    unsigned int mmax_iterations;
    ModelPart::Pointer mp_distance_model_part;
    ModelPart& mr_base_model_part;

    typename SolvingStrategyType::Pointer mp_solving_strategy;


    ///@}
    ///@name Protected Operators
    ///@{
    virtual void ReGenerateDistanceModelPart(ModelPart& base_model_part)
    {
        KRATOS_TRY

        //generate
        ModelPart::Pointer pAuxModelPart = ModelPart::Pointer( new ModelPart("DistancePart",1) );
        mp_distance_model_part.swap(pAuxModelPart);

        mp_distance_model_part->Nodes().clear();
        mp_distance_model_part->Conditions().clear();
        mp_distance_model_part->Elements().clear();

        mp_distance_model_part->SetProcessInfo(  base_model_part.pGetProcessInfo() );
        mp_distance_model_part->SetBufferSize(base_model_part.GetBufferSize());
        mp_distance_model_part->SetProperties(base_model_part.pProperties());
        mp_distance_model_part->Tables() = base_model_part.Tables();

        //assigning the nodes to the new model part
        mp_distance_model_part->Nodes() = base_model_part.Nodes();

        //ensure that the nodes have distance as a DOF
        for (ModelPart::NodesContainerType::iterator iii = base_model_part.NodesBegin(); iii != base_model_part.NodesEnd(); iii++)
        {
            iii->AddDof(DISTANCE);
        }

        //generating the elements
        mp_distance_model_part->Elements().reserve(base_model_part.Elements().size());
        for (ModelPart::ElementsContainerType::iterator iii = base_model_part.ElementsBegin(); iii != base_model_part.ElementsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();
            Element::Pointer p_element = Element::Pointer(new DistanceCalculationElementSimplex<TDim>(
                                             iii->Id(),
                                             iii->pGetGeometry(),
                                             iii->pGetProperties() ) );

            //assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = iii->pGetGeometry();
            
            mp_distance_model_part->Elements().push_back(p_element);
        }


        //using the conditions to mark the boundary with the flag boundary
        //note that we DO NOT add the conditions to the model part
        for (ModelPart::NodesContainerType::iterator iii = mp_distance_model_part->NodesBegin(); iii != mp_distance_model_part->NodesEnd(); iii++)
        {
            iii->Set(BOUNDARY,false);
        }
        for (ModelPart::ConditionsContainerType::iterator iii = base_model_part.ConditionsBegin(); iii != base_model_part.ConditionsEnd(); iii++)
        {
            Geometry< Node<3> >& geom = iii->GetGeometry();
            for(unsigned int i=0; i<geom.size(); i++) geom[i].Set(BOUNDARY,true);
        }

        //next is for mpi (but mpi would also imply calling an mpi strategy)
        Communicator::Pointer pComm = base_model_part.GetCommunicator().Create();
        mp_distance_model_part->SetCommunicator(pComm);

        mdistance_part_is_initialized = true;

        KRATOS_CATCH("")
    }


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
    VariationalDistanceCalculationProcess& operator=(VariationalDistanceCalculationProcess const& rOther);

    /// Copy constructor.
    //VariationalDistanceCalculationProcess(VariationalDistanceCalculationProcess const& rOther);


    ///@}

}; // Class VariationalDistanceCalculationProcess

//avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::PERFORM_STEP1(Kratos::Flags::Create(0));

template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::DO_EXPENSIVE_CHECKS(Kratos::Flags::Create(1));



///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator >> (std::istream& rIStream,
                                  VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>& rThis);

/// output stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const VariationalDistanceCalculationProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_VARIATIONAL_DISTANCE_CALCULATION_PROCESS_INCLUDED  defined 


