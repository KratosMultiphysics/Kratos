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
//                   Ruben Zorrilla
//                    
//


#if !defined(KRATOS_VARIATIONAL_DISTANCE_CALCULATION_PROCESS_INCLUDED )
#define  KRATOS_VARIATIONAL_DISTANCE_CALCULATION_PROCESS_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "elements/distance_calculation_element_simplex.h"
#include "linear_solvers/linear_solver.h"
#include "processes/process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"

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
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver >
class VariationalDistanceCalculationProcess : public Process
{
public:
    
    KRATOS_DEFINE_LOCAL_FLAG(PERFORM_STEP1);
    KRATOS_DEFINE_LOCAL_FLAG(DO_EXPENSIVE_CHECKS);
    
    ///@name Type Definitions
    ///@{

    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef typename SchemeType::Pointer SchemePointerType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverPointerType;
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

    VariationalDistanceCalculationProcess(
        ModelPart& base_model_part,
        typename TLinearSolver::Pointer plinear_solver,
        unsigned int max_iterations = 10)
        :mr_base_model_part(base_model_part)
    {
        KRATOS_TRY

        mmax_iterations = max_iterations;
        mdistance_part_is_initialized = false; //this will be set to true upon completing ReGenerateDistanceModelPart

        // Check that there is at least one element and node in the model
        KRATOS_ERROR_IF(base_model_part.NumberOfNodes() == 0) << "The model part has no nodes." << std::endl;
        KRATOS_ERROR_IF(base_model_part.NumberOfElements() == 0) << "The model Part has no elements." << std::endl;

        // Check if nodes have DISTANCE variable
        VariableUtils().CheckVariableExists<Variable<double > >(DISTANCE, base_model_part.Nodes());
        VariableUtils().CheckVariableExists<Variable<double > >(FLAG_VARIABLE, base_model_part.Nodes());

        if(TDim == 2){
            KRATOS_ERROR_IF(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle) <<
                "In 2D the element type is expected to be a triangle." << std::endl;
        } else if(TDim == 3) {
            KRATOS_ERROR_IF(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra) <<
                "In 3D the element type is expected to be a tetrahedron" << std::endl;
        }

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        ReGenerateDistanceModelPart(base_model_part);

        // Generate a linear strategy
        SchemePointerType pscheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;
        BuilderSolverPointerType pBuilderSolver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> >(plinear_solver);

        mp_solving_strategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver> >(
            *mp_distance_model_part,
            pscheme,
            plinear_solver,
            pBuilderSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        // TODO: check flag DO_EXPENSIVE_CHECKS
        mp_solving_strategy->Check();

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~VariationalDistanceCalculationProcess() override {};

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

        if(mdistance_part_is_initialized == false){
            ReGenerateDistanceModelPart(mr_base_model_part);
        }

        // TODO: check flag    PERFORM_STEP1
        // Step1 - solve a poisson problem with a source term which depends on the sign of the existing distance function
        mp_distance_model_part->pGetProcessInfo()->SetValue(FRACTIONAL_STEP,1);
        
        // Unfix the distances
        const int nnodes = static_cast<int>(mp_distance_model_part->NumberOfNodes());
        #pragma omp parallel for
        for(int i_node = 0; i_node < nnodes; ++i_node){
            auto it_node = mp_distance_model_part->NodesBegin() + i_node;
            double& d = it_node->FastGetSolutionStepValue(DISTANCE);
            double& fix_flag = it_node->FastGetSolutionStepValue(FLAG_VARIABLE);

            // Free the DISTANCE values
            fix_flag = 1.0; 
            it_node->Free(DISTANCE);

            // Save the distances
            it_node->SetValue(DISTANCE, d); 
            
            if(d == 0){
                d = 1.0e-15;
                fix_flag = -1.0;
                it_node->Fix(DISTANCE);
            } else {
                if(d > 0.0){
                    d = 1.0e15; // Set to a large number, to make sure that that the minimal distance is computed according to CaculateTetrahedraDistances
                } else {
                    d = -1.0e15;
                }
            } 
        }

        const int nelem = static_cast<int>(mp_distance_model_part->NumberOfElements());

        #pragma omp parallel for
        for(int i_elem = 0; i_elem < nelem; ++i_elem){
            auto it_elem = mp_distance_model_part->ElementsBegin() + i_elem;
            array_1d<double,TDim+1> distances;
            auto& geom = it_elem->GetGeometry();
            
            for(unsigned int i=0; i<TDim+1; i++){
                distances[i] = geom[i].GetValue(DISTANCE);
            }
            
            const array_1d<double,TDim+1> original_distances = distances;
        
            // The element is cut by the interface
            if(this->IsSplit(distances)){
                // Compute the unsigned distance using GeometryUtils
                if(TDim==3){
                    GeometryUtils::CalculateTetrahedraDistances(geom, distances);
                } else {
                    GeometryUtils::CalculateTriangleDistances(geom, distances);
                }
                
                // Assign the sign using the original distance values
                for(unsigned int i = 0; i < TDim+1; ++i){
                    if(original_distances[i] < 0){
                        distances[i] = -distances[i];
                    }
                }
               
                for(unsigned int i = 0; i < TDim+1; ++i){
                    double &d = geom[i].FastGetSolutionStepValue(DISTANCE);
                    double &fix_flag = geom[i].FastGetSolutionStepValue(FLAG_VARIABLE);
                    geom[i].SetLock();
                    if(std::abs(d) > std::abs(distances[i])){
                        d = distances[i];
                    }
                    fix_flag = -1.0;
                    geom[i].Fix(DISTANCE);
                    geom[i].UnSetLock();
                }
            }
        }

        // SHALL WE SYNCHRONIZE SOMETHING IN HERE?¿?¿??¿ WE'VE CHANGED THE NODAL DISTANCE VALUES FROM THE ELEMENTS...
        this->SynchronizeFixity();
        this->SynchronizeDistance();

        // Compute the maximum and minimum distance for the fixed nodes
        double max_dist = 0.0;
        double min_dist = 0.0;
        for(int i_node = 0; i_node < nnodes; ++i_node){
            auto it_node = mp_distance_model_part->NodesBegin() + i_node;
            if(it_node->IsFixed(DISTANCE)){
                const double& d = it_node->FastGetSolutionStepValue(DISTANCE);
                if(d > max_dist){
                    max_dist = d;
                }
                if(d < min_dist){
                    min_dist = d;
                }
            }
        }

        // Synchronize the maximum and minimum distance values
        auto &r_communicator = mp_distance_model_part->GetCommunicator();
        r_communicator.MaxAll(max_dist);
        r_communicator.MinAll(min_dist);

        // Assign the max dist to all of the non-fixed positive nodes 
        // and the minimum one to the non-fixed negatives
        #pragma omp parallel for
        for(int i_node = 0; i_node < nnodes; ++i_node){
            auto it_node = mp_distance_model_part->NodesBegin() + i_node;
            if(!it_node->IsFixed(DISTANCE)){
                double& d = it_node->FastGetSolutionStepValue(DISTANCE);
                if(d>0){
                    d = max_dist;
                } else {
                    d = min_dist;
                }
            }
        }
        
        mp_solving_strategy->Solve();

        // Step2 - minimize the target residual
        mp_distance_model_part->pGetProcessInfo()->SetValue(FRACTIONAL_STEP,2);
        for(unsigned int it = 0; it<mmax_iterations; it++){
             mp_solving_strategy->Solve();
        }

        // Unfix the distances
        #pragma omp parallel for
        for(int i_node = 0; i_node < nnodes; ++i_node){
            auto it_node = (mp_distance_model_part->NodesBegin()) + i_node;
            it_node->Free(DISTANCE);
        }

        KRATOS_CATCH("")
    }

    virtual void Clear()
    {
        mp_distance_model_part->Nodes().clear();
        mp_distance_model_part->Conditions().clear();
        mp_distance_model_part->Elements().clear();
        // mp_distance_model_part->GetProcessInfo().clear();
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
    VariationalDistanceCalculationProcess(
        ModelPart &base_model_part, 
        unsigned int max_iterations)
        : mr_base_model_part(base_model_part)
    {
        mdistance_part_is_initialized = false;
        mmax_iterations = max_iterations;
    }

    ///@}
    ///@name Protected member Variables
    ///@{

    bool mdistance_part_is_initialized;
    unsigned int mmax_iterations;
    Kratos::unique_ptr<ModelPart> mp_distance_model_part;
    ModelPart& mr_base_model_part;

    typename SolvingStrategyType::UniquePointer mp_solving_strategy;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void ReGenerateDistanceModelPart(ModelPart& base_model_part)
    {
        KRATOS_TRY

        // Generate
        mp_distance_model_part = Kratos::make_unique<ModelPart>("DistancePart",1);
 
        mp_distance_model_part->Nodes().clear();
        mp_distance_model_part->Conditions().clear();
        mp_distance_model_part->Elements().clear();

        mp_distance_model_part->SetProcessInfo(  base_model_part.pGetProcessInfo() );
        mp_distance_model_part->SetBufferSize(base_model_part.GetBufferSize());
        mp_distance_model_part->SetProperties(base_model_part.pProperties());
        mp_distance_model_part->Tables() = base_model_part.Tables();

        // Assigning the nodes to the new model part
        mp_distance_model_part->Nodes() = base_model_part.Nodes();

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof<Variable<double> >(DISTANCE, base_model_part);

        // Generating the elements
        mp_distance_model_part->Elements().reserve(base_model_part.Elements().size());
        for (auto it_elem = base_model_part.ElementsBegin(); it_elem != base_model_part.ElementsEnd(); ++it_elem){
            Properties::Pointer properties = it_elem->pGetProperties();
            Element::Pointer p_element = Kratos::make_shared<DistanceCalculationElementSimplex<TDim> >(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            // Assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = it_elem->pGetGeometry();
            
            mp_distance_model_part->Elements().push_back(p_element);
        }

        // Using the conditions to mark the boundary with the flag boundary
        // Note that we DO NOT add the conditions to the model part
        VariableUtils().SetFlag<ModelPart::NodesContainerType>(BOUNDARY, false, mp_distance_model_part->Nodes());
        // Note that above we have assigned the same geometry. Thus the flag is 
        // set in the distance model part despite we are iterating the base one
        for (auto it_cond = base_model_part.ConditionsBegin(); it_cond != base_model_part.ConditionsEnd(); ++it_cond){
            Geometry< Node<3> >& geom = it_cond->GetGeometry();
            for(unsigned int i=0; i<geom.size(); i++){
                geom[i].Set(BOUNDARY,true);
            }
        }

        mdistance_part_is_initialized = true;

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

    bool IsSplit(const array_1d<double,TDim+1> &rDistances){
        unsigned int positives = 0, negatives = 0;

        for(unsigned int i = 0; i < TDim+1; ++i){
            if(rDistances[i] >= 0){
                ++positives;
            } else {
                ++negatives;
            }
        }

        if (positives > 0 && negatives > 0){
            return true;
        }

        return false;
    }

    void SynchronizeDistance(){

        auto &r_communicator = mp_distance_model_part->GetCommunicator();

        // Only required in the MPI case
        if(r_communicator.TotalProcesses() != 1){
            int nnodes = static_cast<int>(mp_distance_model_part->NumberOfNodes());

            // Set the distance absolute value
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = mp_distance_model_part->NodesBegin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = std::abs(it_node->FastGetSolutionStepValue(DISTANCE));
            }
            
            // Synchronize the unsigned value to minimum
            r_communicator.SynchronizeCurrentDataToMin(DISTANCE);
            
            // Set the distance sign again by retrieving it from the non-historical database
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = mp_distance_model_part->NodesBegin() + i_node;
                if(it_node->GetValue(DISTANCE) < 0.0){
                    it_node->FastGetSolutionStepValue(DISTANCE) = -it_node->FastGetSolutionStepValue(DISTANCE);
                }
            }
        }
    }

    void SynchronizeFixity(){

        auto &r_communicator = mp_distance_model_part->GetCommunicator();

        // Only required in the MPI case
        if(r_communicator.TotalProcesses() != 1){
            int nnodes = static_cast<int>(mp_distance_model_part->NumberOfNodes());

            // Synchronize the fixity flag variable to minium 
            // (-1.0 means fixed and 1.0 means free)
            r_communicator.SynchronizeCurrentDataToMin(FLAG_VARIABLE);

            // Set the fixity according to the synchronized flag
            #pragma omp parallel for
            for(int i_node = 0; i_node < nnodes; ++i_node){
                auto it_node = mp_distance_model_part->NodesBegin() + i_node;
                const double &r_fix_flag = it_node->FastGetSolutionStepValue(FLAG_VARIABLE);
                if (r_fix_flag == -1.0){
                    it_node->Fix(DISTANCE);
                }
            }
        }
    }

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


