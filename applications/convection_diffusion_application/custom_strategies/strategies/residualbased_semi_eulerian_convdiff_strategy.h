// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

#if !defined(KRATOS_RESIDUALBASED_SEMI_EULERIAN_CONVECTION_DIFFUSION_STRATEGY )
#define  KRATOS_RESIDUALBASED_SEMI_EULERIAN_CONVECTION_DIFFUSION_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/kernel.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "convection_diffusion_application.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "includes/convection_diffusion_settings.h"
//#include "custom_utilities/convection_diffusion_settings.h"



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

/// This strategy is used to solve convection-diffusion problem.
/**   Detail class definition.

Exactly the same as the ResidualBasedEulerianConvectionDiffusionStrategy, except for the fact that the assembled elements are only transient diffusion elements.
* This strategy must be combined with an auxiliary tool to convect the particles, which can be particles ( like PFEM-2 ) or BFECC.
* For updates, simply work on the Eulerian strategy and copypaste to this strategy, only changing the name and the type of element.

*/
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class ResidualBasedSemiEulerianConvectionDiffusionStrategy
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedSemiEulerianConvectionDiffusionStrategy );

    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::TDataType TDataType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;



    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /**
    @param model_part Reference to the ModelPart that contains the problem.
     @param pNewLinearSolver pointer to the solver for the temperature system.
     @paramReformDofAtEachIteration=true.
     @param time_order=2.
     @param prediction_order == 2.
    */


    ResidualBasedSemiEulerianConvectionDiffusionStrategy(
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        bool ReformDofAtEachIteration = false,
        int dimension = 3
    )
        : 
        mrReferenceModelPart(model_part),
        SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false)
    {
        KRATOS_TRY

        if(mrReferenceModelPart.GetOwnerModel().HasModelPart("ConvectionDiffusionPart"))
            KRATOS_ERROR << "ConvectionDiffusionPart already exists when constructing ResidualBasedSemiEulerianConvectionDiffusionStrategy" << std::endl;

		GenerateMeshPart(dimension);
		mdimension = dimension;
        mOldDt = 0.00;

        //ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        //ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        //(mpConvectionModelPart->GetProcessInfo()).GetValue(CONVECTION_DIFFUSION_SETTINGS) = my_settings;
		Check();


        //initializing fractional velocity solution step
        typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
        typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
                                               ( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

        bool CalculateReactions = false;
        bool CalculateNormDxFlag = true;

        //choosing the solving strategy
//			mstep1 = typename BaseType::Pointer(
//				new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver >
//				(model_part,pscheme,pNewLinearSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
//			mstep1->SetEchoLevel(2);


        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        //const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
        //BuilderSolverTypePointer componentwise_build = BuilderSolverTypePointer(new	ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace,TDenseSpace,TLinearSolver,Variable<double> > (pNewLinearSolver,rUnknownVar) );
        //mstep1 = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > 				(*mpConvectionModelPart,pscheme,pNewLinearSolver,componentwise_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );

        BuilderSolverTypePointer pBuilderSolver = BuilderSolverTypePointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(pNewLinearSolver) );
        mstep1 = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(*mpConvectionModelPart,pscheme,pNewLinearSolver,pBuilderSolver,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag) );


        mstep1->SetEchoLevel(2);

        KRATOS_CATCH("")
    }



    /** Destructor.
    */
    virtual ~ResidualBasedSemiEulerianConvectionDiffusionStrategy() 
    {
        mrReferenceModelPart.GetOwnerModel().DeleteModelPart("ConvectionDiffusionPart");
    }

    /** Destructor.
    */

    //*********************************************************************************
    //**********************************************************************
    double Solve() override
    {
      KRATOS_TRY

        //calculate the BDF coefficients
        //careful. Using information from the currentprocessinfo of the original model part (not the convection model part)
      //ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
      //double Dt = rCurrentProcessInfo[DELTA_TIME];
      //int stationary= rCurrentProcessInfo[STATIONARY];

	  //SOLVING THE PROBLEM
	  double Dp_norm = mstep1->Solve();

      return Dp_norm;
      KRATOS_CATCH("")
	}



    void SetEchoLevel(int Level) override
    {
        mstep1->SetEchoLevel(Level);
    }

    void Clear() override
    {
        mstep1->Clear();
    }

    int Check() override
    {
        KRATOS_TRY
        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
        if (rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)==false)
			KRATOS_THROW_ERROR(std::logic_error, "no CONVECTION_DIFFUSION_SETTINGS in model_part", "");
        //std::cout << "ConvDiff::Check(). If crashes, check CONVECTION_DIFFUSION_SETTINGS is defined" << std::endl;

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
		//if(my_settings->IsDefinedVolumeSourceVariable()==true)
		//{
		//	if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetVolumeSourceVariable()) == false)
		//		KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: VolumeSource Variable defined but not contained in the model part", "");
		//}
		//else
		//	std::cout << "No VolumeSource variable assigned for ConvDiff. Assuming VolumeSource=0" << std::endl;
		if(my_settings->IsDefinedVolumeSourceVariable()==true)
			KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: VolumeSource not yet implemented", "");

		//SURFACE SOURCE VARIABLE
		//if(my_settings->IsDefinedSurfaceSourceVariable()==true)
		//{
		//	if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetSurfaceSourceVariable()) == false)
		//		KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: SurfaceSource Variable defined but not contained in the model part", "");
		//}
		//else
		//	std::cout << "No SurfaceSource variable assigned for ConvDiff. Assuming SurfaceSource=0" << std::endl;
		if(my_settings->IsDefinedSurfaceSourceVariable()==true)
			KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: SurfaceSource not yet implemented", "");

		//PROJECTION VARIABLE
		//used as intermediate variable, is the variable at time n+1 but only accounting for the convective term.
		if(my_settings->IsDefinedProjectionVariable()==true)
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetProjectionVariable()) == false)
				KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: Projection Variable defined but not contained in the model part", "");
		}
		else
		{
			std::cout << "No Projection variable assigned for ConvDiff. Using PROJECTED_SCALAR1" << std::endl;
			my_settings->SetProjectionVariable(PROJECTED_SCALAR1);

		}
		//if(my_settings->IsDefinedProjectionVariable()==true)
		//	KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: ProjectionVariable not useed. Remove it", "");

		//CONVECTION VELOCITY VARIABLE
		//CURRENTLY WE ARE USING (VELOCITY -MESH_VELOCITY) TO CONVECT, so the ConvectionVariable must not be used:
		//if(my_settings->IsDefinedConvectionVariable()==true)
		//{
		//	if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetConvectionVariable()) == false)
		//		KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: Convection Variable defined but not contained in the model part", "");
		//}
		//else
		//	std::cout << "No Projection variable assigned for ConvDiff. Assuming Convection=0" << std::endl;
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
		//if(my_settings->IsDefinedTransferCoefficientVariable()==true)
		//{
		//	if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetTransferCoefficientVariable()) == false)
		//		KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: TransferCoefficient Variable defined but not contained in the model part", "");
		//}
		//else
		//	std::cout << "No TransferCoefficient variable assigned for ConvDiff. Assuming TransferCoefficient=0" << std::endl;
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
    ModelPart* mpConvectionModelPart;
    typename BaseType::Pointer mstep1;
    double mOldDt;
    int mdimension;





    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */

  void GenerateMeshPart(int dimension)
  {
    if(!mrReferenceModelPart.GetOwnerModel().HasModelPart("ConvectionDiffusionPart"))
        mrReferenceModelPart.GetOwnerModel().DeleteModelPart("ConvectionDiffusionPart");

    mpConvectionModelPart = &(mrReferenceModelPart.GetOwnerModel().CreateModelPart("ConvectionDiffusionPart"));

	mpConvectionModelPart->SetProcessInfo(  BaseType::GetModelPart().pGetProcessInfo() );
    mpConvectionModelPart->SetBufferSize( BaseType::GetModelPart().GetBufferSize());

    //initializing mesh nodes
    mpConvectionModelPart->Nodes() = BaseType::GetModelPart().Nodes();

    //creating mesh elements
    ModelPart::ElementsContainerType& MeshElems = mpConvectionModelPart->Elements();
    Element::Pointer pElem;

    if(dimension == 2)
    for(ModelPart::ElementsContainerType::iterator it
        = BaseType::GetModelPart().ElementsBegin();

        it != BaseType::GetModelPart().ElementsEnd(); ++it) {

      pElem = Element::Pointer(new EulerianDiffusionElement<2,3>(
              (*it).Id(),
              (*it).pGetGeometry(),
              (*it).pGetProperties() ) );
      MeshElems.push_back(pElem);
    }

    if(dimension == 3)
    for(ModelPart::ElementsContainerType::iterator it
        = BaseType::GetModelPart().ElementsBegin();

        it != BaseType::GetModelPart().ElementsEnd(); ++it) {

      pElem = Element::Pointer(new EulerianDiffusionElement<3,4>(
              (*it).Id(),
              (*it).pGetGeometry(),
              (*it).pGetProperties() ) );
      MeshElems.push_back(pElem);
    }
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
    ResidualBasedSemiEulerianConvectionDiffusionStrategy(const ResidualBasedSemiEulerianConvectionDiffusionStrategy& Other);


    /*@} */

}; /* Class ResidualBasedSemiEulerianConvectionDiffusionStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_SEMI_EULERIAN_CONVECTION_DIFFUSION_STRATEGY  defined */

