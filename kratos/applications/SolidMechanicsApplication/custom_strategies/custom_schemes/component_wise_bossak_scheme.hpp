//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPONENT_BASED_BOSSAK_SCHEME )
#define  KRATOS_COMPONENT_BASED_BOSSAK_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "solid_mechanics_application.h"

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
/*@} */

template<class TSparseSpace,  class TDenseSpace >
class ComponentBasedBossakScheme: public ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>
{
public:

    struct LocalSystemContributions
    {
    private:
    
      LocalSystemMatrixType *mpLeftHandSideContribution;
      LocalSystemVectorType *mpRightHandSideContribution;

      std::vector<LocalSystemMatrixType> *mpLHS_Components;
      const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Variables;

      std::vector<LocalSystemVectorType> *mpRHS_Components;
      const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Variables;
      
    public:
      
      //setting pointer variables
      void SetLeftHandSideContribution ( LocalSystemMatrixType& rLeftHandSideContribution ) { mpLeftHandSideContribution = &rLeftHandSideContribution; };
      void SetRightHandSideContribution ( LocalSystemVectorType& rRightHandSideContribution ) { mpRightHandSideContribution = &rRightHandSideContribution; };
      void SetLHS_Components ( std::vector<LocalSystemMatrixType>& rLHS_Components ) { mpLHS_Components = &rLHS_Components; };
      void SetLHS_Variables     ( std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables ) { mpLHS_Variables = &rLHS_Variables; };
      void SetRHS_Components ( std::vector<LocalSystemVectorType>& rRHS_Components ) { mpRHS_Components = &rRHS_Components; };
      void SetRHS_Variables     ( std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables ) { mpRHS_Variables = &rRHS_Variables; };

      //getting pointer variables
      LocalSystemMatrixType& GetLeftHandSideContribution  () { return *mpLeftHandSideContribution; };
      LocalSystemVectorType& GetRightHandSideContribution () { return *mpRightHandSideContribution; };

      std::vector<LocalSystemMatrixType>& GetLHS_Components() { return *mpLHS_Components; };
      std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Variables() { return *mpLHS_Variables; };
      std::vector<LocalSystemVectorType>& GetRHS_Components() { return *mpRHS_Components; };
      std::vector< Variable< LocalSystemVectorType > >& GetRHS_Variables() { return *mpRHS_Variables; };

    }
      
    struct LHS_LocalSystemContributions
    {
    private:
    
      std::vector<LocalSystemMatrixType> *mpLHS_Components;
      const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Variables;

      
    public:
      
      LocalSystemMatrixType LeftHandSideContribution;
 
      //setting pointer variables
      void SetLHS_Components ( std::vector<LocalSystemMatrixType>& rLHS_Components ) { mpLHS_Components = &rLHS_Components; };
      void SetLHS_Variables     ( std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables ) { mpLHS_Variables = &rLHS_Variables; };
 
      //getting pointer variables
      std::vector<LocalSystemMatrixType>& GetLHS_Components() { return *mpLHS_Components; };
      std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Variables() { return *mpLHS_Variables; };

    }

    struct RHS_LocalSystemContributions
    {
    private:
    
      std::vector<LocalSystemVectorType> *mpRHS_Components;
      const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Variables;
      
    public:
      
      LocalSystemVectorType RightHandSideContribution;

      //setting pointer variables
      void SetRHS_Components ( std::vector<LocalSystemVectorType>& rRHS_Components ) { mpRHS_Components = &rRHS_Components; };
      void SetRHS_Variables     ( std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables ) { mpRHS_Variables = &rRHS_Variables; };

      //getting pointer variables
      std::vector<LocalSystemVectorType>& GetRHS_Components() { return *mpRHS_Components; };
      std::vector< Variable< LocalSystemVectorType > >& GetRHS_Variables() { return *mpRHS_Variables; };

    }


    /**@name Type Definitions */

    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION(ComponentBasedBossakScheme);

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;


    /*@} */

    /**
     * Constructor.
     * The bossak method
     */
    ComponentBasedBossakScheme(double rAlpham=0,double rDynamic=1)
      :ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>(rAlpham,rDynamic)
    {
    }

    /** Destructor.
     */
    virtual ~ComponentBasedBossakScheme
    () {}


    /*@} */
    /**@name Operators
     */
    /*@{ */


    //***************************************************************************
    //***************************************************************************

    /** this function is designed to be called in the builder and solver to introduce*/

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
	SystemContributions& rLocalSystem,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //(rCurrentElement) -> InitializeNonLinearIteration(rCurrentProcessInfo);

	LocalSystemMatrixType& rLeftHandSideContribution  = rLocalSystem.GetLeftHandSideContribution();
	LocalSystemVectorType& rRightHandSideContribution = rLocalSystem.GetRightHandSideContribution();

	std::vector<LocalSystemMatrixType>& rLHS_Components = rLocalSystem.GetLHS_Components();
	std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables = rLocalSystem.GetLHS_Variables();
	std::vector<LocalSystemVectorType>& rRHS_Components = rLocalSystem.GetRHS_Components();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = rLocalSystem.GetRHS_Variables();

        (rCurrentElement) -> CalculateLocalSystem(rLHS_Components, rLHS_Variables, rRHS_Components, rRHS_Variables, rCurrentProcessInfo);

	for( unsigned int i=0; i<rLHS_Components.size(); i++ )
	  {	    
	    rLeftHandSideContribution += rLHS_Components[i];
	  }
	
	for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	  {
	    rRightHandSideContribution += rRHS_Components[i];
	  }
	

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentElement) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        }

        (rCurrentElement) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS (rLeftHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

            AddDynamicsToRHS (rCurrentElement, rRightHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

        }

        //AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
	RHS_SystemContributions& rRHS_LocalSystem,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {

        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

	LocalSystemVectorType& rRightHandSideContribution = rLocalSystem.GetRightHandSideContribution();

	std::vector<LocalSystemVectorType>& rRHS_Components = rRHS_LocalSystem.GetRHS_Components();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = rRHS_LocalSystem.GetRHS_Variables();

        //basic operations for the element considered
        (rCurrentElement) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);

	
	for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	  {
	    rRightHandSideContribution += rRHS_Components[i];
	  }
	

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentElement) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);
        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToRHS (rCurrentElement, rRightHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
        }

        KRATOS_CATCH("")

    }

    //***************************************************************************
    //***************************************************************************

    /** functions totally analogous to the precedent but applied to
          the "condition" objects
    */
    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
	SystemContributions& rLocalSystem,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {


        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

	LocalSystemMatrixType& rLeftHandSideContribution  = rLocalSystem.GetLeftHandSideContribution();
	LocalSystemVectorType& rRightHandSideContribution = rLocalSystem.GetRightHandSideContribution();

	std::vector<LocalSystemMatrixType>& rLHS_Components = rLocalSystem.GetLHS_Components();
	std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables = rLocalSystem.GetLHS_Variables();
	std::vector<LocalSystemVectorType>& rRHS_Components = rLocalSystem.GetRHS_Components();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = rLocalSystem.GetRHS_Variables();

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateLocalSystem(rLHS_Components, rLHS_Variables, rRHS_Components, rRHS_Variables, rCurrentProcessInfo);


	for( unsigned int i=0; i<rLHS_Components.size(); i++ )
	  {	    
	    rLeftHandSideContribution += rLHS_Components[i];
	  }
	
	for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	  {
	    rRightHandSideContribution += rRHS_Components[i];
	  }

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentCondition) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentCondition) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS  (rLeftHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

            AddDynamicsToRHS  (rCurrentCondition, rRightHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
        }

        //AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
	RHS_SystemContributions& rRHS_LocalSystem,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current condition
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

	LocalSystemVectorType& rRightHandSideContribution = rLocalSystem.GetRightHandSideContribution();

	std::vector<LocalSystemVectorType>& rRHS_Components = rRHS_LocalSystem.GetRHS_Components();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = rRHS_LocalSystem.GetRHS_Variables();

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);
	
	for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	  {
	    rRightHandSideContribution += rRHS_Components[i];
	  }

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentCondition) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentCondition) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

        //adding the dynamic contributions (static is already included)

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToRHS  (rCurrentCondition, rRightHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

        }

        KRATOS_CATCH("")
    }


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

protected:

    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
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
private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */
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
    /**@name Unaccessible methods */
    /*@{ */
}; /* Class ComponentBasedBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_COMPONENT_BASED_BOSSAK_SCHEME defined */


