//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPONENT_WISE_BOSSAK_SCHEME )
#define KRATOS_COMPONENT_WISE_BOSSAK_SCHEME

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
class ComponentWiseBossakScheme: public ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>
{
public:

    /**@name Type Definitions */

    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ComponentWiseBossakScheme );

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
    ComponentWiseBossakScheme(double rAlpham=0,double rDynamic=1)
      :ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>(rAlpham,rDynamic)
    {
    }

    /** Destructor.
     */
    virtual ~ComponentWiseBossakScheme() 
    {}


    /*@} */
    /**@name Operators
     */
    /*@{ */

    LocalSystemComponents& GetLocalSystemComponents()
    {
      return mLocalSystem;
    }

    //***************************************************************************
    //***************************************************************************

    /** this function is designed to be called in the builder and solver to introduce*/

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //(rCurrentElement) -> InitializeNonLinearIteration(rCurrentProcessInfo);

	std::vector<LocalSystemMatrixType>& rLHS_Components = mLocalSystem.GetLHS_Element_Components();
	std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables = mLocalSystem.GetLHS_Element_Variables();
	std::vector<LocalSystemVectorType>& rRHS_Components = mLocalSystem.GetRHS_Element_Components();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = mLocalSystem.GetRHS_Element_Variables();

	if( rRHS_Variables.size() > 0 && rLHS_Variables.size() > 0)
	  {

	    //basic operations for the element considered
	    (rCurrentElement) -> CalculateLocalSystem(rLHS_Components, rLHS_Variables, rRHS_Components, rRHS_Variables, rCurrentProcessInfo);


	    for( unsigned int i=0; i<rLHS_Components.size(); i++ )
	      {	    
		rLHS_Contribution += rLHS_Components[i];
	      }
	
	    for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	      {
		rRHS_Contribution += rRHS_Components[i];
	      }

	  }
	else
	  {
	    if( rRHS_Variables.size() > 0 )
	      {
		(rCurrentElement) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);

		(rCurrentElement) -> CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
		
		for( unsigned int i=0; i<rRHS_Components.size(); i++ )
		  {
		    rRHS_Contribution += rRHS_Components[i];
		  }
		

	      }
	    else if ( rRHS_Variables.size() > 0 )
	      {
		
		(rCurrentElement) -> CalculateLeftHandSide(rLHS_Components, rLHS_Variables, rCurrentProcessInfo);

		(rCurrentElement) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
		
		for( unsigned int i=0; i<rLHS_Components.size(); i++ )
		  {
		    rLHS_Contribution += rLHS_Components[i];
		  }
		

	      }
	    else
	      {

		(rCurrentElement) -> CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

	      }

	  }
	

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentElement) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        }

        (rCurrentElement) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS (rLHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

            AddDynamicsToRHS (rCurrentElement, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

        }

        //AssembleTimeSpaceLHS(rCurrentElement, rLHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
	LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {

        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

	std::vector<LocalSystemVectorType>& rRHS_Components = mLocalSystem.GetRHS_Element_Components();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = mLocalSystem.GetRHS_Element_Variables();

	if( rRHS_Variables.size() > 0 )
	  {
	    //basic operations for the element considered
	    (rCurrentElement) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);
	
	    for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	      {
		rRHS_Contribution += rRHS_Components[i];
	      }
	  }
	else
	  {
	    //basic operations for the element considered
	    (rCurrentElement) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
	  }
	

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentElement) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);
        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToRHS (rCurrentElement, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
        }

        KRATOS_CATCH( "" )

    }

    //***************************************************************************
    //***************************************************************************

    /** functions totally analogous to the precedent but applied to
          the "condition" objects
    */
    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {


        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

	std::vector<LocalSystemMatrixType>& rLHS_Components = mLocalSystem.GetLHS_Condition_Components();
	std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables = mLocalSystem.GetLHS_Condition_Variables();
	std::vector<LocalSystemVectorType>& rRHS_Components = mLocalSystem.GetRHS_Condition_Components();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = mLocalSystem.GetRHS_Condition_Variables();

	if( rRHS_Variables.size() > 0 && rLHS_Variables.size() > 0)
	  {

	    //basic operations for the element considered
	    (rCurrentCondition) -> CalculateLocalSystem(rLHS_Components, rLHS_Variables, rRHS_Components, rRHS_Variables, rCurrentProcessInfo);


	    for( unsigned int i=0; i<rLHS_Components.size(); i++ )
	      {	    
		rLHS_Contribution += rLHS_Components[i];
	      }
	
	    for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	      {
		rRHS_Contribution += rRHS_Components[i];
	      }

	  }
	else
	  {
	    if( rRHS_Variables.size() > 0 )
	      {
		(rCurrentCondition) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);

		(rCurrentCondition) -> CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
		
		for( unsigned int i=0; i<rRHS_Components.size(); i++ )
		  {
		    rRHS_Contribution += rRHS_Components[i];
		  }
		

	      }
	    else if ( rRHS_Variables.size() > 0 )
	      {
		
		(rCurrentCondition) -> CalculateLeftHandSide(rLHS_Components, rLHS_Variables, rCurrentProcessInfo);

		(rCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
		
		for( unsigned int i=0; i<rLHS_Components.size(); i++ )
		  {
		    rLHS_Contribution += rLHS_Components[i];
		  }
		

	      }
	    else
	      {

		(rCurrentCondition) -> CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

	      }

	  }
	    	    
        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentCondition) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentCondition) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS  (rLHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

            AddDynamicsToRHS  (rCurrentCondition, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
        }

        //AssembleTimeSpaceLHS_Condition(rCurrentCondition, rLHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current condition
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

	std::vector<LocalSystemVectorType>& rRHS_Components = mLocalSystem.GetRHS_Condition_Components();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = mLocalSystem.GetRHS_Condition_Variables();

	if( rRHS_Variables.size() > 0 )
	  {
	    //basic operations for the element considered
	    (rCurrentCondition) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);
	
	    for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	      {
		rRHS_Contribution += rRHS_Components[i];
	      }
	  }
	else
	  {
	    //basic operations for the element considered
	    (rCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
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

            AddDynamicsToRHS  (rCurrentCondition, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

        }

        KRATOS_CATCH( "" )
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
    LocalSystemComponents mLocalSystem;
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
}; /* Class ComponentWiseBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_COMPONENT_WISE_BOSSAK_SCHEME defined */


