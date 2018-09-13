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
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer plinear_solver,
        const double max_cfl = 1.0,
        const double cross_wind_stabilization_factor = 0.7,
        const unsigned int max_substeps = 0)
        : mrBaseModelPart(rBaseModelPart),
          mrLevelSetVar(rLevelSetVar),
          mMaxAllowedCFL(max_cfl),
          mMaxSubsteps(max_substeps)
    {
        KRATOS_TRY
        
        // Check that there is at least one element and node in the model
        const auto n_nodes = rBaseModelPart.NumberOfNodes();
        const auto n_elems = rBaseModelPart.NumberOfElements();

        KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;

        VariableUtils().CheckVariableExists< Variable< double > >(rLevelSetVar, rBaseModelPart.Nodes());
        VariableUtils().CheckVariableExists< Variable< array_1d < double, 3 > > >(VELOCITY, rBaseModelPart.Nodes());

        if(TDim == 2){
            KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle) << 
                "In 2D the element type is expected to be a triangle" << std::endl;
        } else if(TDim == 3) {
            KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra) << 
                "In 3D the element type is expected to be a tetrahedra" << std::endl;
        }
        
        // Allocate if needed the variable DYNAMIC_TAU of the process info, and if it does not exist, set it to zero
        if( rBaseModelPart.GetProcessInfo().Has(DYNAMIC_TAU) == false){
            rBaseModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU,0.0);
        }

        // Allocate if needed the variable CONVECTION_DIFFUSION_SETTINGS of the process info, and create it if it does not exist
        if( rBaseModelPart.GetProcessInfo().Has(CONVECTION_DIFFUSION_SETTINGS) == false){
            ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
            rBaseModelPart.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
            p_conv_diff_settings->SetUnknownVariable(rLevelSetVar);
            p_conv_diff_settings->SetConvectionVariable(VELOCITY);
        }

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        mDistancePartIsInitialized = false;
        ReGenerateConvectionModelPart(rBaseModelPart);

        // Generate a linear strategy
        typename SchemeType::Pointer pscheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();
        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        BuilderSolverTypePointer pBuilderSolver = Kratos::make_shared< ResidualBasedBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver > >(plinear_solver);
        mpSolvingStrategy = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            *mpDistanceModelPart,
            pscheme,
            plinear_solver,
            pBuilderSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        mpSolvingStrategy->SetEchoLevel(0);
        
        rBaseModelPart.GetProcessInfo().SetValue(CROSS_WIND_STABILIZATION_FACTOR, cross_wind_stabilization_factor);
        
        //TODO: check flag DO_EXPENSIVE_CHECKS
        mpSolvingStrategy->Check();

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

        if(mDistancePartIsInitialized == false){
            ReGenerateConvectionModelPart(mrBaseModelPart);
        }

        // Evaluate steps needed to achieve target max_cfl
        const auto n_substep = EvaluateNumberOfSubsteps();

        // Save the variables to be employed so that they can be restored after the solution
        ProcessInfo& rCurrentProcessInfo = mpDistanceModelPart->GetProcessInfo();
        const auto & r_previous_var = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->GetUnknownVariable();
        const double previous_delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);

        // Save current level set value and current and previous step velocity values
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
            const auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            mVelocity[i_node] = it_node->FastGetSolutionStepValue(VELOCITY);
            mVelocityOld[i_node] = it_node->FastGetSolutionStepValue(VELOCITY,1);
            mOldDistance[i_node] = it_node->FastGetSolutionStepValue(mrLevelSetVar,1);
        }

        const double dt = previous_delta_time / static_cast<double>(n_substep);
        rCurrentProcessInfo.SetValue(DELTA_TIME, dt);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(mrLevelSetVar);

        const int rank = mrBaseModelPart.GetCommunicator().MyPID();

        for(unsigned int step = 1; step <= n_substep; ++step){

            KRATOS_INFO_IF("LevelSetConvectionProcess", mpSolvingStrategy->GetEchoLevel() > 0 && rank == 0) << 
                "Doing step "<< step << " of " << n_substep << std::endl;

            // Compute shape functions of old and new step
            const double Nold = 1.0 - static_cast<double>(step) / static_cast<double>(n_substep);
            const double Nnew = 1.0 - Nold;
            
            const double Nold_before = 1.0 - static_cast<double>(step-1) / static_cast<double>(n_substep);
            const double Nnew_before = 1.0 - Nold_before;            
            
            // Emulate clone time step by copying the new distance onto the old one
            #pragma omp parallel for
            for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                const array_1d<double,3>& v = mVelocity[i_node];
                const array_1d<double,3>& v_old = mVelocityOld[i_node];

                it_node->FastGetSolutionStepValue(VELOCITY) = Nold * v_old + Nnew * v;
                it_node->FastGetSolutionStepValue(VELOCITY, 1) = Nold_before * v_old + Nnew_before * v;
                it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
            }
            
            mpSolvingStrategy->Solve();
        }

        // Reset the processinfo to the original settings 
        rCurrentProcessInfo.SetValue(DELTA_TIME, previous_delta_time);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(r_previous_var);
        
        // Reset the velocities and levelset values to the one saved before the solution process
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
            auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(VELOCITY) = mVelocity[i_node];
            it_node->FastGetSolutionStepValue(VELOCITY,1) = mVelocityOld[i_node];
            it_node->FastGetSolutionStepValue(mrLevelSetVar,1) = mOldDistance[i_node];
        }
        
        KRATOS_CATCH("")
    }

    virtual void Clear(){
        mpDistanceModelPart->Nodes().clear();
        mpDistanceModelPart->Conditions().clear();
        mpDistanceModelPart->Elements().clear();
        // mpDistanceModelPart->GetProcessInfo().clear();
        mDistancePartIsInitialized = false;

        mpSolvingStrategy->Clear();
        
        mVelocity.clear();
        mVelocityOld.clear();
        mOldDistance.clear();
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

    ModelPart& mrBaseModelPart;
    Kratos::unique_ptr<ModelPart> mpDistanceModelPart;

    Variable<double>& mrLevelSetVar;

    const double mMaxAllowedCFL;
    
    bool mDistancePartIsInitialized;

	const unsigned int mMaxSubsteps;

    std::vector< double > mOldDistance;
    std::vector< array_1d<double,3> > mVelocity, mVelocityOld;

    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /// Constructor without linear solver for derived classes
    LevelSetConvectionProcess(
        Variable<double> &rLevelSetVar,
        ModelPart &rBaseModelPart,
        const double MaxCFL = 1.0,
        const unsigned int MaxSubSteps = 0)
        : mrBaseModelPart(rBaseModelPart),
          mrLevelSetVar(rLevelSetVar),
          mMaxAllowedCFL(MaxCFL),
          mMaxSubsteps(MaxSubSteps)
    {
        mDistancePartIsInitialized = false;
    }

    virtual void ReGenerateConvectionModelPart(ModelPart& rBaseModelPart){

        KRATOS_TRY

        // Check buffer size
        const auto base_buffer_size = rBaseModelPart.GetBufferSize();
        KRATOS_ERROR_IF(base_buffer_size < 2) << 
            "Base model part buffer size is " << base_buffer_size << ". Set it to a minimum value of 2." << std::endl;

        // Generate
        Kratos::unique_ptr<ModelPart> p_aux_model_part = Kratos::make_unique<ModelPart>("DistancePart");
        mpDistanceModelPart.swap(p_aux_model_part);

        mpDistanceModelPart->Nodes().clear();
        mpDistanceModelPart->Conditions().clear();
        mpDistanceModelPart->Elements().clear();

        mpDistanceModelPart->SetProcessInfo(rBaseModelPart.pGetProcessInfo());
        mpDistanceModelPart->SetBufferSize(base_buffer_size);
        mpDistanceModelPart->SetProperties(rBaseModelPart.pProperties());
        mpDistanceModelPart->Tables() = rBaseModelPart.Tables();

        // Assigning the nodes to the new model part
        mpDistanceModelPart->Nodes() = rBaseModelPart.Nodes();

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof< Variable < double> >(mrLevelSetVar, rBaseModelPart);

        // Generating the elements
        mpDistanceModelPart->Elements().reserve(rBaseModelPart.NumberOfElements());
        for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem){
            Element::Pointer p_element = Kratos::make_shared< LevelSetConvectionElementSimplex < TDim, TDim+1 > >(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            // Assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = it_elem->pGetGeometry();
            
            mpDistanceModelPart->Elements().push_back(p_element);
        }

        // Next is for mpi (but mpi would also imply calling an mpi strategy)
        Communicator::Pointer pComm = rBaseModelPart.GetCommunicator().Create();
        mpDistanceModelPart->SetCommunicator(pComm);
        
        // Resize the arrays
        const auto n_nodes = mpDistanceModelPart->NumberOfNodes();
        mVelocity.resize(n_nodes);
        mVelocityOld.resize(n_nodes);
        mOldDistance.resize(n_nodes);

        mDistancePartIsInitialized = true;

        KRATOS_CATCH("")
    }


    unsigned int EvaluateNumberOfSubsteps(){
        // First of all compute the cfl number 
        const auto n_elem = mpDistanceModelPart->NumberOfElements();
        const double dt = mpDistanceModelPart->GetProcessInfo()[DELTA_TIME];
        
		// Vector where each thread will store its maximum (VS does not support OpenMP reduce max)
		int NumThreads = OpenMPUtils::GetNumThreads();
		std::vector<double> list_of_max_local_cfl(NumThreads, 0.0);

        //TODO: Update this loop to avoid using thread id
        #pragma omp parallel shared(list_of_max_local_cfl)
        for(int i_elem = 0; i_elem < static_cast<int>(n_elem); i_elem++){
            const auto it_elem = mpDistanceModelPart->ElementsBegin() + i_elem;
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
        mpDistanceModelPart->GetCommunicator().MaxAll(max_cfl_found);
        
        unsigned int n_steps = static_cast<unsigned int>(max_cfl_found / mMaxAllowedCFL); 
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


