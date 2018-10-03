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

#if !defined(KRATOS_RESIDUALBASED_EULERIAN_CONVECTION_DIFFUSION_STRATEGY )
#define  KRATOS_RESIDUALBASED_EULERIAN_CONVECTION_DIFFUSION_STRATEGY


/* System includes */


/* External includes */

/* Project includes */
#include "includes/define.h"
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

The convection-diffusion problem
	\f$ \rho C M \frac{\partial T}{\partial t} + \rho C S T  = - \kappa L T \f$ (1)

      is completed with the standard boundary conditions of prescribed temperature and prescribed normal
      heat flux in the thermal problem. For surfaces exposed to fire conditions, energy losses
      due to radiation and convection must be taken into account, and the thermal boundary condition is

          \f$ \kappa \frac{\partial T}{\partial n} + \overline{q_n} = 0 \f$ (2)
      where

          \f$ \overline{q_n} = q_n - \varepsilon \sigma (T^4 - T_0^4) - \alpha_c (T - T_0)\f$ (3)

Then, this strategy is employed to solve the following equation


	Evaluates  \f$ L h s = \frac{\rho C}{\Delta t} (W, N) + (W, v. \nabla N) + \kappa (\nabla W, \nabla N) + 4 \epsilon \sigma T^3 \left\langle W, N \right\rangle + \alpha \left\langle
W, N \right\rangle \f$ and \f$ R h s = \rho (W, Q) + \frac{\rho C}{\Delta t} (W, T^n)- \left\langle W, q \right\rangle - \epsilon \sigma \left\langle W, T^4 - T_0^4 \right\rangle - \left\langle W, \alpha (T - T_0) \right\rangle - L h s \ast T \f$


*/
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class ResidualBasedEulerianConvectionDiffusionStrategy
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedEulerianConvectionDiffusionStrategy );

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


    ResidualBasedEulerianConvectionDiffusionStrategy(
        ModelPart& model_part,
        typename TLinearSolver::Pointer pNewLinearSolver,
        bool ReformDofAtEachIteration = false,
        int dimension = 3
    )
        : 
        mrModelPart(model_part),
        SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false)
    {
        KRATOS_TRY

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
    virtual ~ResidualBasedEulerianConvectionDiffusionStrategy() 
    {
        Model& current_model = mrModelPart.GetOwnerModel();
        current_model.DeleteModelPart("ConvDiffPart");
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
            KRATOS_ERROR << "no CONVECTION_DIFFUSION_SETTINGS in model_part" << std::endl;
        //std::cout << "ConvDiff::Check(). If crashes, check CONVECTION_DIFFUSION_SETTINGS is defined" << std::endl;

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

		//DENSITY VARIABLE
		if(my_settings->IsDefinedDensityVariable()==true)
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetDensityVariable()) == false)
                KRATOS_ERROR << "ConvDiffSettings: Density Variable defined but not contained in the model part" << std::endl;
		}
		else
			std::cout << "No density variable assigned for ConvDiff. Assuming density=1" << std::endl;

		//DIFFUSION VARIABLE
		if(my_settings->IsDefinedDiffusionVariable()==true)
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetDiffusionVariable()) == false)
                KRATOS_ERROR << "ConvDiffSettings: Diffusion Variable defined but not contained in the model part" << std::endl;
		}
		else
			std::cout << "No diffusion variable assigned for ConvDiff. Assuming diffusivity=0" << std::endl;

		//UNKNOWN VARIABLE
		if(my_settings->IsDefinedUnknownVariable()==true)
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetUnknownVariable()) == false)
				KRATOS_ERROR << "ConvDiffSettings: Unknown Variable defined but not contained in the model part" << std::endl;
		}
		else
			KRATOS_ERROR << "ConvDiffSettings: Unknown Variable not defined!" << std::endl;

		//VOLUME SOURCE VARIABLE
		if(my_settings->IsDefinedVolumeSourceVariable()==true)
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetVolumeSourceVariable()) == false)
				KRATOS_ERROR << "ConvDiffSettings: VolumeSource Variable defined but not contained in the model part" << std::endl;
		}
		else
			std::cout << "No VolumeSource variable assigned for ConvDiff. Assuming VolumeSource=0" << std::endl;


		//SURFACE SOURCE VARIABLE
		//if(my_settings->IsDefinedSurfaceSourceVariable()==true)
		//{
		//	if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetSurfaceSourceVariable()) == false)
		//		KRATOS_ERROR << "ConvDiffSettings: SurfaceSource Variable defined but not contained in the model part" << std::endl;
		//}
		//else
		//	std::cout << "No SurfaceSource variable assigned for ConvDiff. Assuming SurfaceSource=0" << std::endl;
        if(my_settings->IsDefinedSurfaceSourceVariable()==true)
			KRATOS_ERROR << "ConvDiffSettings: SurfaceSource not yet implemented" << std::endl;

		//PROJECTION VARIABLE
		//if(my_settings->IsDefinedProjectionVariable()==true)
		//{
		//	if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetProjectionVariable()) == false)
		//		KRATOS_THROW_ERROR(std::logic_error, "ConvDiffSettings: Projection Variable defined but not contained in the model part", "");
		//}
		//else
		//	std::cout << "No Projection variable assigned for ConvDiff. Assuming Projection=0" << std::endl;
		if(my_settings->IsDefinedProjectionVariable()==true)
			KRATOS_ERROR << "ConvDiffSettings: ProjectionVariable not useed. Remove it" << std::endl;

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
			KRATOS_ERROR << "ConvDiffSettings: ConvectionVariable not used. Use VelocityVariable instead" << std::endl;

		//MESH VELOCITY VARIABLE
		if(my_settings->IsDefinedMeshVelocityVariable()==true)
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetMeshVelocityVariable()) == false)
				KRATOS_ERROR << "ConvDiffSettings: MeshVelocity Variable defined but not contained in the model part" << std::endl;
		}
		else
			std::cout << "No MeshVelocity variable assigned for ConvDiff. Assuming MeshVelocity=0" << std::endl;

		//VELOCITY VARIABLE
		if(my_settings->IsDefinedVelocityVariable()==true)
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetVelocityVariable()) == false)
				KRATOS_ERROR << "ConvDiffSettings: Velocity Variable defined but not contained in the model part" << std::endl;
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
			KRATOS_ERROR << "ConvDiffSettings: TransferCoefficient not yet implemented" << std::endl;

		//SPECIFIC HEAT VARIABLE
		if(my_settings->IsDefinedSpecificHeatVariable()==true)
		{
			if (BaseType::GetModelPart().NodesBegin()->SolutionStepsDataHas(my_settings->GetSpecificHeatVariable()) == false)
				KRATOS_ERROR << "ConvDiffSettings: SpecificHeat Variable defined but not contained in the model part" << std::endl;
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

  virtual void GenerateMeshPart(int dimension)
  {
    Model& current_model = mrModelPart.GetOwnerModel();
    if(current_model.HasModelPart("ConvDiffPart"))
        current_model.DeleteModelPart("ConvDiffPart");
        
    // Generate
    mpConvectionModelPart = &(current_model.CreateModelPart("ConvDiffPart"));


	mpConvectionModelPart->SetProcessInfo(  BaseType::GetModelPart().pGetProcessInfo() );
    mpConvectionModelPart->SetBufferSize( BaseType::GetModelPart().GetBufferSize());

    //initializing mesh nodes
    mpConvectionModelPart->Nodes() = BaseType::GetModelPart().Nodes();

    //creating mesh elements
    ModelPart::ElementsContainerType& MeshElems = mpConvectionModelPart->Elements();
    Element::Pointer pElem;

    if(dimension == 2)
    {
        for(ModelPart::ElementsContainerType::iterator it= BaseType::GetModelPart().ElementsBegin(); it != BaseType::GetModelPart().ElementsEnd(); ++it)
         {
            Element::GeometryType& rGeom = it->GetGeometry();
            const unsigned int& NumNodes = rGeom.size();

            if (NumNodes == 3)
            {
              pElem = Element::Pointer(new EulerianConvectionDiffusionElement<2,3>(
                      (*it).Id(),
                      (*it).pGetGeometry(),
                      (*it).pGetProperties() ) );
              MeshElems.push_back(pElem);
            }
            else if(NumNodes == 4)
            {
                pElem = Element::Pointer(new EulerianConvectionDiffusionElement<2,4>(
                      (*it).Id(),
                      (*it).pGetGeometry(),
                      (*it).pGetProperties() ) );
              MeshElems.push_back(pElem);
            }
        }
    }
    else
    {
        for(ModelPart::ElementsContainerType::iterator it= BaseType::GetModelPart().ElementsBegin(); it != BaseType::GetModelPart().ElementsEnd(); ++it)
        {
            Element::GeometryType& rGeom = it->GetGeometry();
            const unsigned int& NumNodes = rGeom.size();

            if (NumNodes == 4)
            {
              pElem = Element::Pointer(new EulerianConvectionDiffusionElement<3,4>(
                      (*it).Id(),
                      (*it).pGetGeometry(),
                      (*it).pGetProperties() ) );
              MeshElems.push_back(pElem);
            }
            else if(NumNodes == 8)
            {
                pElem = Element::Pointer(new EulerianConvectionDiffusionElement<3,8>(
                      (*it).Id(),
                      (*it).pGetGeometry(),
                      (*it).pGetProperties() ) );
              MeshElems.push_back(pElem);
            }
         }
    }
  }

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
    ModelPart& mrModelPart;
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
    ResidualBasedEulerianConvectionDiffusionStrategy(const ResidualBasedEulerianConvectionDiffusionStrategy& Other);


    /*@} */

}; /* Class ResidualBasedEulerianConvectionDiffusionStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_EULERIAN_CONVECTION_DIFFUSION_STRATEGY  defined */
