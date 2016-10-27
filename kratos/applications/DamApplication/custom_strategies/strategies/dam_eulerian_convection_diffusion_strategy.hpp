//   
//   Project Name:                 KratosDamApplication $
//   Last Modified by:    $Author:       Lorenzo Gracia $
//   Date:                $Date:           October 2016 $
//   Revision:            $Revision:                1.0 $
//

#if !defined(KRATOS_DAM_EULERIAN_CONVECTION_DIFUSSION_STRATEGY)
#define KRATOS_DAM_EULERIAN_CONVECTION_DIFUSSION_STRATEGY

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "convection_diffusion_application.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "includes/convection_diffusion_settings.h"
#include "custom_elements/eulerian_conv_diff.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class DamEulerianConvectionDiffusionStrategy : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(DamEulerianConvectionDiffusionStrategy);
    
    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    DamEulerianConvectionDiffusionStrategy(
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        Parameters& rParameters,
        bool ReformDofAtEachIteration = false,
        int dimension = 3
        ) : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false),
        mrParameters(rParameters)
        {            
           KRATOS_TRY

            this->GenerateMeshPart(dimension);
            KRATOS_WATCH(*mpConvectionModelPart);
            mdimension = dimension;
            mOldDt = 0.00;

            Check();
            
            //initializing fractional velocity solution step
            typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
            typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
                                                   ( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

            bool CalculateReactions = false;
            bool CalculateNormDxFlag = true;

            typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

            BuilderSolverTypePointer pBuilderSolver = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(pNewLinearSolver) );
            mstep1 = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(*mpConvectionModelPart,pscheme,pNewLinearSolver,pBuilderSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag) );

            mstep1->SetEchoLevel(2);

            KRATOS_CATCH("")
           
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~DamEulerianConvectionDiffusionStrategy() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double Solve()
        {
          KRATOS_TRY
                
          //SOLVING THE PROBLEM	  
          double Dp_norm = mstep1->Solve();
        
          return Dp_norm;
          KRATOS_CATCH("")
        }

    //------------------------------------------------------------------------------------
    
    virtual void SetEchoLevel(int Level)
    {
        mstep1->SetEchoLevel(Level);
    }

    //------------------------------------------------------------------------------------

    virtual void Clear()
    {
        mstep1->Clear();
    }
    
    //------------------------------------------------------------------------------------

    virtual int Check()
    {
        KRATOS_TRY
        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        if (rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)==false)
			KRATOS_THROW_ERROR(std::logic_error, "no CONVECTION_DIFFUSION_SETTINGS in model_part", "");
        
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
		
		//DENSITY VARIABLE
		if(my_settings->IsDefinedDensityVariable()==true) 
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetDensityVariable()) == false)
				KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: Density Variable defined but not contained in the model part", "");
		}
		else
			std::cout << "No density variable assigned for ConvDiff. Assuming density=1" << std::endl;

		//DIFFUSION VARIABLE
		if(my_settings->IsDefinedDiffusionVariable()==true) 
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetDiffusionVariable()) == false)
				KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: Diffusion Variable defined but not contained in the model part", "");
		}
		else
			std::cout << "No diffusion variable assigned for ConvDiff. Assuming diffusivity=0" << std::endl;

		//UNKNOWN VARIABLE
		if(my_settings->IsDefinedUnknownVariable()==true) 
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetUnknownVariable()) == false)
				KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: Unknown Variable defined but not contained in the model part", "");
		}
		else
			KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: Unknown Variable not defined!", "");
		
		//VOLUME SOURCE VARIABLE
		if(my_settings->IsDefinedVolumeSourceVariable()==true) 
			KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: VolumeSource not yet implemented", "");
		
		//SURFACE SOURCE VARIABLE
		if(my_settings->IsDefinedSurfaceSourceVariable()==true) 
			KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: SurfaceSource not yet implemented", "");

		//PROJECTION VARIABLE
		if(my_settings->IsDefinedProjectionVariable()==true) 
			KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: ProjectionVariable not useed. Remove it", "");

		//CONVECTION VELOCITY VARIABLE
		if(my_settings->IsDefinedConvectionVariable()==true) 
			KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: ConvectionVariable not used. Use VelocityVariable instead", "");

		//MESH VELOCITY VARIABLE
		if(my_settings->IsDefinedMeshVelocityVariable()==true) 
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetMeshVelocityVariable()) == false)
				KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: MeshVelocity Variable defined but not contained in the model part", "");
		}
		else
			std::cout << "No MeshVelocity variable assigned for ConvDiff. Assuming MeshVelocity=0" << std::endl;
		
		//VELOCITY VARIABLE	
		if(my_settings->IsDefinedVelocityVariable()==true) 
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetVelocityVariable()) == false)
				KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: Velocity Variable defined but not contained in the model part", "");
		}
		else
			std::cout << "No Velocity variable assigned for ConvDiff. Assuming Velocity=0" << std::endl;
		
		//TRANSFER COEFFICIENT VARIABLE
		if(my_settings->IsDefinedTransferCoefficientVariable()==true) 
			KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: TransferCoefficient not yet implemented", "");
		
		//SPECIFIC HEAT VARIABLE	
		if(my_settings->IsDefinedSpecificHeatVariable()==true) 
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetSpecificHeatVariable()) == false)
				KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: SpecificHeat Variable defined but not contained in the model part", "");
		}
		else
			std::cout << "No SpecificHeat variable assigned for ConvDiff. Assuming SpecificHeat=1" << std::endl;

        return 0;

        KRATOS_CATCH("")

    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    Parameters& mrParameters;
    ModelPart::Pointer mpConvectionModelPart;
    
    void GenerateMeshPart(int dimension)
    {        
        mpConvectionModelPart = ModelPart::Pointer( new ModelPart("ConvectionPart",1) );

        mpConvectionModelPart->SetProcessInfo(  BaseType::GetModelPart().pGetProcessInfo() );
        mpConvectionModelPart->SetBufferSize( BaseType::GetModelPart().GetBufferSize());

        //initializing mesh nodes
        mpConvectionModelPart->Nodes() = BaseType::GetModelPart().Nodes();

        //creating mesh elements
        ModelPart::ElementsContainerType& MeshElems = mpConvectionModelPart->Elements();
        Element::Pointer pElem;
        
        unsigned int NumBodySubModelParts = mrParameters["problem_domain_body_sub_model_part_list"].size();
    
        KRATOS_WATCH(NumBodySubModelParts)
    
        for(unsigned int m = 0; m < NumBodySubModelParts; m++)
        {
            ModelPart& SubModelPart = BaseType::GetModelPart().GetSubModelPart(mrParameters["problem_domain_body_sub_model_part_list"][m].GetString());
            
            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();
            
            //// Information about the geometry            
            Element::GeometryType& rGeom = el_begin->GetGeometry();
            const unsigned int& Dim  = rGeom.WorkingSpaceDimension();
            const unsigned int& NumNodes = rGeom.size();            
                         
            if (Dim == 2)
            {
                if (NumNodes == 3)
                {
                    for(int i = 0; i < NElems; i++)
                    {
                        ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
                        pElem = Element::Pointer(new EulerianConvectionDiffusionElement<2,3>( (*itElem).Id(), (*itElem).pGetGeometry(),(*itElem).pGetProperties() ) );
                        MeshElems.push_back(pElem);
                    }
                }
                else if(NumNodes == 4)
                {
                    for(int i = 0; i < NElems; i++)
                    {
                        ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
                        pElem = Element::Pointer(new EulerianConvectionDiffusionElement<2,4>( (*itElem).Id(), (*itElem).pGetGeometry(),(*itElem).pGetProperties() ) );
                        MeshElems.push_back(pElem);
                    }
                }
            }        
            else
            {
                if(NumNodes == 4)
                {
                    for(int i = 0; i < NElems; i++)
                    {
                        ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
                        pElem = Element::Pointer(new EulerianConvectionDiffusionElement<3,4>( (*itElem).Id(), (*itElem).pGetGeometry(),(*itElem).pGetProperties() ) );
                        MeshElems.push_back(pElem);
                    }
                }
                else if(NumNodes == 8)
                {
                    for(int i = 0; i < NElems; i++)
                    {
                        ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
                        pElem = Element::Pointer(new EulerianConvectionDiffusionElement<3,8>( (*itElem).Id(), (*itElem).pGetGeometry(),(*itElem).pGetProperties() ) );
                        MeshElems.push_back(pElem);
                    }
                }
            }
        }
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


private:

    typename BaseType::Pointer mstep1;
    double mOldDt;
    int mdimension;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    DamEulerianConvectionDiffusionStrategy(const DamEulerianConvectionDiffusionStrategy& Other);


}; // Class DamEulerianConvectionDiffusionStrategy

} // namespace Kratos

#endif // KRATOS_DAM_EULERIAN_CONVECTION_DIFUSSION_STRATEGY  defined
