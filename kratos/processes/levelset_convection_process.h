/*
Kratos Multi-Physics

Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
All rights reserved.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

		Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
		Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
		in the documentation and/or other materials provided with the distribution.
		All advertising materials mentioning features or use of this software must display the following acknowledgement:
			This product includes Kratos Multi-Physics technology.
		Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#if !defined(KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED )
#define  KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes
#include "includes/kratos_flags.h"



// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/levelset_convection_element_simplex.h"

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

template< unsigned int TDim >
class LevelSetConvectionProcess
    : public Process
{
public:
    
    KRATOS_DEFINE_LOCAL_FLAG(PERFORM_STEP1);
    KRATOS_DEFINE_LOCAL_FLAG(DO_EXPENSIVE_CHECKS);
    
    ///@name Type Definitions
    ///@{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType,  LocalSpaceType > SchemeType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > SolvingStrategyType;

///@}
    ///@name Pointer Definitions
    /// Pointer definition of LevelSetConvectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(LevelSetConvectionProcess);




    ///@}
    ///@name Life Cycle
    ///@{

    /**
     */
    LevelSetConvectionProcess(Variable<double>& rLevelSetVar,
                              ModelPart& base_model_part,
                                          typename LinearSolverType::Pointer plinear_solver,
                                          double max_cfl = 1.0
                                         )
        :mr_base_model_part(base_model_part), mrLevelSetVar(rLevelSetVar), mmax_allowed_cfl(max_cfl)
    {
        KRATOS_TRY

        
        //check that there is at least one element and node in the model
        if(base_model_part.Nodes().size() == 0) KRATOS_THROW_ERROR(std::logic_error, "the model has no Nodes","");
        if(base_model_part.Elements().size() == 0) KRATOS_THROW_ERROR(std::logic_error, "the model has no Elements","");
        if(base_model_part.NodesBegin()->SolutionStepsDataHas(rLevelSetVar) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing rLevelSetVar variable on solution step data","");
        if(base_model_part.NodesBegin()->SolutionStepsDataHas(VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data","");

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
        
        //allocate if needed the variable DYNAMIC_TAU of the process info, and if it is not zero, set it to zero
        if( base_model_part.GetProcessInfo().Has(DYNAMIC_TAU) == false)
            base_model_part.GetProcessInfo().SetValue(DYNAMIC_TAU,0.0);
        if( base_model_part.GetProcessInfo().Has(CONVECTION_DIFFUSION_SETTINGS) == false)
        {
            ConvectionDiffusionSettings::Pointer psettings( new ConvectionDiffusionSettings() );
            base_model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, psettings);
            psettings->SetUnknownVariable(rLevelSetVar);
            psettings->SetConvectionVariable(VELOCITY);
        }
        //generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        mdistance_part_is_initialized = false;
        ReGenerateConvectionModelPart(base_model_part);

        //generate a linear strategy
        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType,LocalSpaceType >() );
        typedef typename BuilderAndSolver<SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer BuilderSolverTypePointer;

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        BuilderSolverTypePointer pBuilderSolver = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<SparseSpaceType,LocalSpaceType,LinearSolverType>(plinear_solver) );
        mp_solving_strategy = typename SolvingStrategyType::Pointer( new ResidualBasedLinearStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType >(*mp_distance_model_part,pscheme,plinear_solver,pBuilderSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag) );

        mp_solving_strategy->SetEchoLevel(0);
        
        //TODO: check flag DO_EXPENSIVE_CHECKS
        mp_solving_strategy->Check();

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~LevelSetConvectionProcess()
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

    virtual void Execute()
    {
        KRATOS_TRY;


        if(mdistance_part_is_initialized == false)
            ReGenerateConvectionModelPart(mr_base_model_part);

        //evaluate steps needed to achieve target max_cfl
        unsigned int nsubstep = EvaluateNumberOfSubsteps();

        //save the variables to be employed so that they can be restored after the solution
        ProcessInfo& rCurrentProcessInfo = mp_distance_model_part->GetProcessInfo();
        const Variable<double>& previous_var = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->GetUnknownVariable();
        const double previous_delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);

        //save velocities and old value of LevelSetFunction
        int i = 0;
        for (ModelPart::NodesContainerType::iterator iii = mp_distance_model_part->NodesBegin(); iii != mp_distance_model_part->NodesEnd(); iii++)
        {
            mold_dist[i] = iii->FastGetSolutionStepValue(mrLevelSetVar,1);
            mv[i] = iii->FastGetSolutionStepValue(VELOCITY);
            mvold[i] = iii->FastGetSolutionStepValue(VELOCITY,1);
            i++;
        }

        double dt = previous_delta_time/static_cast<double>(nsubstep);
        rCurrentProcessInfo.SetValue(DELTA_TIME, dt);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(mrLevelSetVar);
        for(unsigned int step = 1; step<=nsubstep; step++)
        {
            std::cout << "doing step "<< step << " of " << nsubstep << std::endl;
            //compute shape functions of old and new step
            double Nold = 1.0-static_cast<double>(step)/static_cast<double>(nsubstep);
            double Nnew = 1.0-Nold;
            
            double Nold_before = 1.0-static_cast<double>(step-1)/static_cast<double>(nsubstep);
            double Nnew_before = 1.0-Nold_before;            
            
            //emulate clone time step by copying the new distance onto the old one
            i=0;
            for (ModelPart::NodesContainerType::iterator iii = mp_distance_model_part->NodesBegin(); iii != mp_distance_model_part->NodesEnd(); iii++)
            {
                iii->FastGetSolutionStepValue(mrLevelSetVar,1) = iii->FastGetSolutionStepValue(mrLevelSetVar);

                const array_1d<double,3>& vold = mvold[i];
                const array_1d<double,3>& v = mv[i];
                iii->FastGetSolutionStepValue(VELOCITY,1) = Nold_before*vold + Nnew_before*v;
                iii->FastGetSolutionStepValue(VELOCITY) = Nold*vold + Nnew*v;
                i++;
            }
            
            mp_solving_strategy->Solve();
        }

        
        
        //reset the processinfo to the original settings 
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(previous_var);
        rCurrentProcessInfo.SetValue(DELTA_TIME, previous_delta_time);
        
        //reset the velocities and levelset values to the one saved before the solution process
        i = 0;
        for (ModelPart::NodesContainerType::iterator iii = mp_distance_model_part->NodesBegin(); iii != mp_distance_model_part->NodesEnd(); iii++)
        {
            iii->FastGetSolutionStepValue(mrLevelSetVar,1) = mold_dist[i];
            
            const array_1d<double,3>& vold = mvold[i];
            const array_1d<double,3>& v = mv[i];
            iii->FastGetSolutionStepValue(VELOCITY,1) = vold;
            iii->FastGetSolutionStepValue(VELOCITY) = v;
            
            i++;
        }
        
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
        
        mold_dist.clear();
        mv.clear();
        mvold.clear();

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
    virtual std::string Info() const
    {
        return "LevelSetConvectionProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "LevelSetConvectionProcess";
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

    ModelPart& mr_base_model_part;
    Variable<double>& mrLevelSetVar;
    double mmax_allowed_cfl;
    
    bool mdistance_part_is_initialized;
    unsigned int mmax_iterations;
    ModelPart::Pointer mp_distance_model_part;

    std::vector< double > mold_dist;
    std::vector< array_1d<double,3> > mv, mvold;


    SolvingStrategyType::Pointer mp_solving_strategy;


    ///@}
    ///@name Protected Operators
    ///@{
    void ReGenerateConvectionModelPart(ModelPart& base_model_part)
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
            iii->AddDof(mrLevelSetVar);
        }

        //generating the elements
        mp_distance_model_part->Elements().reserve(base_model_part.Elements().size());
        for (ModelPart::ElementsContainerType::iterator iii = base_model_part.ElementsBegin(); iii != base_model_part.ElementsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();
            Element::Pointer p_element = Element::Pointer(new LevelSetConvectionElementSimplex<TDim, TDim+1>(
                                             iii->Id(),
                                             iii->pGetGeometry(),
                                             iii->pGetProperties() ) );

            //assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = iii->pGetGeometry();
            
            mp_distance_model_part->Elements().push_back(p_element);
        }


        //next is for mpi (but mpi would also imply calling an mpi strategy)
        Communicator::Pointer pComm = base_model_part.GetCommunicator().Create();
        mp_distance_model_part->SetCommunicator(pComm);
        
        //resize the arrays 
        mold_dist.resize(mp_distance_model_part->Nodes().size());
        mv.resize(mp_distance_model_part->Nodes().size());
        mvold.resize(mp_distance_model_part->Nodes().size());

        mdistance_part_is_initialized = true;

        KRATOS_CATCH("")
    }


    unsigned int EvaluateNumberOfSubsteps()
    {
        //first of all compute the cfl number 
        ModelPart::ElementsContainerType::iterator el_begin = mp_distance_model_part->ElementsBegin();
        const unsigned int nelem = mp_distance_model_part->Elements().size();
        const double dt = mp_distance_model_part->GetProcessInfo()[DELTA_TIME];
        
        double max_cfl_found = 0.0;
// //         #pragma omp parallel for reduce(min: min_cfl_found)
        for(unsigned int it=0; it<nelem; it++)
        {
            ModelPart::ElementsContainerType::iterator iii = el_begin+it;
            
            Geometry< Node<3> >& geom = iii->GetGeometry();

            boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim > DN_DX;
            array_1d<double, TDim+1 > N;
            double vol;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, vol);
            
            //compute h
            double h=0.0;
            for(unsigned int i=0; i<TDim+1; i++)
            {
                double h_inv = 0.0;
                for(unsigned int k=0; k<TDim; k++)
                {
                    h_inv += DN_DX(i,k)*DN_DX(i,k);
                }
                h += 1.0/h_inv;
            }
            h = sqrt(h)/static_cast<double>(TDim+1);


            //get avg velocity at the nodes
            array_1d<double, 3 > vgauss = ZeroVector(3);
            for(unsigned int i=0; i<TDim+1; i++)
            {
                vgauss += N[i]* geom[i].FastGetSolutionStepValue(VELOCITY);
            }

            const double vnorm = norm_2(vgauss);

            double cfl_local = vnorm/h;
            
            if(cfl_local > max_cfl_found) max_cfl_found = cfl_local;
        }
        
        max_cfl_found *= dt;
        
        int nsteps = static_cast<unsigned int>(max_cfl_found/mmax_allowed_cfl); 
        if(nsteps < 1) nsteps=1;
        
        KRATOS_WATCH(nsteps)
        
        
        return nsteps;
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
    LevelSetConvectionProcess& operator=(LevelSetConvectionProcess const& rOther);

    /// Copy constructor.
    //LevelSetConvectionProcess(LevelSetConvectionProcess const& rOther);


    ///@}

}; // Class LevelSetConvectionProcess

//avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim > const Kratos::Flags LevelSetConvectionProcess<TDim>::PERFORM_STEP1(Kratos::Flags::Create(0));
template< unsigned int TDim > const Kratos::Flags LevelSetConvectionProcess<TDim>::DO_EXPENSIVE_CHECKS(Kratos::Flags::Create(1));



///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim>
inline std::istream& operator >> (std::istream& rIStream,
                                  LevelSetConvectionProcess<TDim>& rThis);

/// output stream function
template< unsigned int TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const LevelSetConvectionProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED  defined 


