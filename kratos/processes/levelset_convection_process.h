//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//

#if !defined(KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED )
#define  KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/convection_diffusion_settings.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "elements/levelset_convection_element_simplex.h"
#include "geometries/geometry_data.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
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
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and convects a level set distance
 * on the top of it
*/
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver >
class LevelSetConvectionProcess
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
    ///@{

    /// Pointer definition of LevelSetConvectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(LevelSetConvectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     */
    LevelSetConvectionProcess(
        Variable<double>& rLevelSetVar,
        ModelPart& base_model_part,
        typename TLinearSolver::Pointer plinear_solver,
        const double max_cfl = 1.0,
        const double cross_wind_stabilization_factor = 0.7,
        const unsigned int max_substeps = 0)
        : mr_base_model_part(base_model_part), 
        mrLevelSetVar(rLevelSetVar), 
        mmax_allowed_cfl(max_cfl), 
        mMaxSubsteps(max_substeps)
    {
        KRATOS_TRY
        
        // Check that there is at least one element and node in the model
        const auto n_nodes = base_model_part.NumberOfNodes();
        const auto n_elems = base_model_part.NumberOfElements();

        KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;

        VariableUtils().CheckVariableExists< Variable< double > >(rLevelSetVar, base_model_part.Nodes());
        VariableUtils().CheckVariableExists< Variable< array_1d < double, 3 > > >(VELOCITY, base_model_part.Nodes());

        if(TDim == 2){
            KRATOS_ERROR_IF(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle) << 
                "In 2D the element type is expected to be a triangle" << std::endl;
        } else if(TDim == 3) {
            KRATOS_ERROR_IF(base_model_part.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra) << 
                "In 3D the element type is expected to be a tetrahedra" << std::endl;
        }
        
        // Allocate if needed the variable DYNAMIC_TAU of the process info, and if it does not exist, set it to zero
        if( base_model_part.GetProcessInfo().Has(DYNAMIC_TAU) == false){
            base_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU,0.0);
        }

        // Allocate if needed the variable CONVECTION_DIFFUSION_SETTINGS of the process info, and create it if it does not exist
        if( base_model_part.GetProcessInfo().Has(CONVECTION_DIFFUSION_SETTINGS) == false){
            ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
            base_model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
            p_conv_diff_settings->SetUnknownVariable(rLevelSetVar);
            p_conv_diff_settings->SetConvectionVariable(VELOCITY);
        }

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        mdistance_part_is_initialized = false;
        ReGenerateConvectionModelPart(base_model_part);

        // Generate a linear strategy
        typename SchemeType::Pointer pscheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();
        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        BuilderSolverTypePointer pBuilderSolver = Kratos::make_shared< ResidualBasedBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver > >(plinear_solver);
        mp_solving_strategy = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            *mp_distance_model_part,
            pscheme,
            plinear_solver,
            pBuilderSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        mp_solving_strategy->SetEchoLevel(0);
        
        base_model_part.GetProcessInfo().SetValue(CROSS_WIND_STABILIZATION_FACTOR, cross_wind_stabilization_factor);
        
        //TODO: check flag DO_EXPENSIVE_CHECKS
        mp_solving_strategy->Check();

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~LevelSetConvectionProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    void operator()(){
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY;

        if(mdistance_part_is_initialized == false){
            ReGenerateConvectionModelPart(mr_base_model_part);
        }

        // Evaluate steps needed to achieve target max_cfl
        const auto n_substep = EvaluateNumberOfSubsteps();

        // Save the variables to be employed so that they can be restored after the solution
        ProcessInfo& rCurrentProcessInfo = mp_distance_model_part->GetProcessInfo();
        const auto & r_previous_var = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->GetUnknownVariable();
        const double previous_delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);

        // Save current level set value and current and previous step velocity values
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mp_distance_model_part->NumberOfNodes()); ++i_node){
            const auto it_node = mp_distance_model_part->NodesBegin() + i_node;
            mv[i_node] = it_node->FastGetSolutionStepValue(VELOCITY);
            mvold[i_node] = it_node->FastGetSolutionStepValue(VELOCITY,1);
            mold_dist[i_node] = it_node->FastGetSolutionStepValue(mrLevelSetVar,1);
        }

        const double dt = previous_delta_time / static_cast<double>(n_substep);
        rCurrentProcessInfo.SetValue(DELTA_TIME, dt);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(mrLevelSetVar);

        const int rank = mr_base_model_part.GetCommunicator().MyPID();

        for(unsigned int step = 1; step <= n_substep; ++step){

            KRATOS_INFO_IF("LevelSetConvectionProcess", mp_solving_strategy->GetEchoLevel() > 0 && rank == 0) << 
                "Doing step "<< step << " of " << n_substep << std::endl;

            // Compute shape functions of old and new step
            const double Nold = 1.0 - static_cast<double>(step) / static_cast<double>(n_substep);
            const double Nnew = 1.0 - Nold;
            
            const double Nold_before = 1.0 - static_cast<double>(step-1) / static_cast<double>(n_substep);
            const double Nnew_before = 1.0 - Nold_before;            
            
            // Emulate clone time step by copying the new distance onto the old one
            #pragma omp parallel for
            for (int i_node = 0; i_node < static_cast<int>(mp_distance_model_part->NumberOfNodes()); ++i_node){
                auto it_node = mp_distance_model_part->NodesBegin() + i_node;

                const array_1d<double,3>& v = mv[i_node];
                const array_1d<double,3>& v_old = mvold[i_node];

                it_node->FastGetSolutionStepValue(VELOCITY) = Nold * v_old + Nnew * v;
                it_node->FastGetSolutionStepValue(VELOCITY, 1) = Nold_before * v_old + Nnew_before * v;
                it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
            }
            
            mp_solving_strategy->Solve();
        }

        // Reset the processinfo to the original settings 
        rCurrentProcessInfo.SetValue(DELTA_TIME, previous_delta_time);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(r_previous_var);
        
        // Reset the velocities and levelset values to the one saved before the solution process
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mp_distance_model_part->NumberOfNodes()); ++i_node){
            auto it_node = mp_distance_model_part->NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(VELOCITY) = mv[i_node];
            it_node->FastGetSolutionStepValue(VELOCITY,1) = mvold[i_node];
            it_node->FastGetSolutionStepValue(mrLevelSetVar,1) = mold_dist[i_node];
        }
        
        KRATOS_CATCH("")
    }

    virtual void Clear(){
        mp_distance_model_part->Nodes().clear();
        mp_distance_model_part->Conditions().clear();
        mp_distance_model_part->Elements().clear();
        // mp_distance_model_part->GetProcessInfo().clear();
        mdistance_part_is_initialized = false;

        mp_solving_strategy->Clear();
        
        mv.clear();
        mvold.clear();
        mold_dist.clear();
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
    std::string Info() const override {
        return "LevelSetConvectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "LevelSetConvectionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
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

    ModelPart& mr_base_model_part;
    ModelPart::UniquePointer mp_distance_model_part;

    Variable<double>& mrLevelSetVar;

    const double mmax_allowed_cfl;
    
    bool mdistance_part_is_initialized;

	const unsigned int mMaxSubsteps;

    std::vector< double > mold_dist;
    std::vector< array_1d<double,3> > mv, mvold;

    typename SolvingStrategyType::UniquePointer mp_solving_strategy;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ReGenerateConvectionModelPart(ModelPart& base_model_part){
        KRATOS_TRY

        // Generate
        const auto base_buffer_size = base_model_part.GetBufferSize();
        ModelPart::UniquePointer p_aux_model_part = Kratos::make_unique<ModelPart>("DistancePart", base_buffer_size);
        mp_distance_model_part.swap(p_aux_model_part);

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
        VariableUtils().AddDof< Variable < double> >(mrLevelSetVar, base_model_part);

        // Generating the elements
        mp_distance_model_part->Elements().reserve(base_model_part.NumberOfElements());
        for (auto it_elem = base_model_part.ElementsBegin(); it_elem != base_model_part.ElementsEnd(); ++it_elem){
            Element::Pointer p_element = Kratos::make_shared< LevelSetConvectionElementSimplex < TDim, TDim+1 > >(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            // Assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = it_elem->pGetGeometry();
            
            mp_distance_model_part->Elements().push_back(p_element);
        }

        // Next is for mpi (but mpi would also imply calling an mpi strategy)
        Communicator::Pointer pComm = base_model_part.GetCommunicator().Create();
        mp_distance_model_part->SetCommunicator(pComm);
        
        // Resize the arrays
        const auto n_nodes = mp_distance_model_part->NumberOfNodes();
        mv.resize(n_nodes);
        mvold.resize(n_nodes);
        mold_dist.resize(n_nodes);

        mdistance_part_is_initialized = true;

        KRATOS_CATCH("")
    }


    unsigned int EvaluateNumberOfSubsteps(){
        // First of all compute the cfl number 
        const auto n_elem = mp_distance_model_part->NumberOfElements();
        const double dt = mp_distance_model_part->GetProcessInfo()[DELTA_TIME];
        
		// Vector where each thread will store its maximum (VS does not support OpenMP reduce max)
		int NumThreads = OpenMPUtils::GetNumThreads();
		std::vector<double> list_of_max_local_cfl(NumThreads, 0.0);

        //TODO: Update this loop to avoid using thread id
        #pragma omp parallel shared(list_of_max_local_cfl)
        for(int i_elem = 0; i_elem < static_cast<int>(n_elem); i_elem++){
            const auto it_elem = mp_distance_model_part->ElementsBegin() + i_elem;
            Geometry< Node<3> >& r_geom = it_elem->GetGeometry();

            double vol;
            array_1d<double, TDim+1 > N;
            BoundedMatrix<double, TDim+1, TDim > DN_DX;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, vol);

			int k = OpenMPUtils::ThisThread();
			double& max_cfl = list_of_max_local_cfl[k];

            // Compute h
            double h=0.0;
            for(unsigned int i=0; i<TDim+1; i++){
                double h_inv = 0.0;
                for(unsigned int k=0; k<TDim; k++){
                    h_inv += DN_DX(i,k)*DN_DX(i,k);
                }
                h += 1.0/h_inv;
            }
            h = sqrt(h)/static_cast<double>(TDim+1);

            // Get average velocity at the nodes
            array_1d<double, 3 > vgauss = ZeroVector(3);
            for(unsigned int i=0; i<TDim+1; i++){
                vgauss += N[i]* r_geom[i].FastGetSolutionStepValue(VELOCITY);
            }

            double cfl_local = norm_2(vgauss) / h;
            if(cfl_local > max_cfl){
                max_cfl = cfl_local;
            }
        }

		// Now we get the maximum at each thread level
        double max_cfl_found = 0.0;
		for (int k=0; k < NumThreads;k++){
		    if (max_cfl_found < list_of_max_local_cfl[k]){
                max_cfl_found = list_of_max_local_cfl[k];
            }
        }
        max_cfl_found *= dt;

        // Synchronize maximum CFL between processes
        mp_distance_model_part->GetCommunicator().MaxAll(max_cfl_found);
        
        int n_steps = static_cast<unsigned int>(max_cfl_found / mmax_allowed_cfl); 
        if(n_steps < 1){
            n_steps = 1;
        }

		// Now we compare with the maximum set
		if (mMaxSubsteps > 0 && mMaxSubsteps < n_steps){
            n_steps = mMaxSubsteps;
        }
        
        return n_steps;
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
    LevelSetConvectionProcess& operator=(LevelSetConvectionProcess const& rOther);

    /// Copy constructor.
    //LevelSetConvectionProcess(LevelSetConvectionProcess const& rOther);

    ///@}
}; // Class LevelSetConvectionProcess

// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver > const Kratos::Flags LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::PERFORM_STEP1(Kratos::Flags::Create(0));
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver > const Kratos::Flags LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::DO_EXPENSIVE_CHECKS(Kratos::Flags::Create(1));

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// Input stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator >> (
    std::istream& rIStream,
    LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// Output stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis){

    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED  defined 


