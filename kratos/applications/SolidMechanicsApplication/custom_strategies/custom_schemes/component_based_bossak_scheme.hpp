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
    
      std::vector<LocalSystemMatrixType> *mpLHS_Contributions;
      const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Variables;

      std::vector<LocalSystemVectorType> *mpRHS_Contributions;
      const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Variables;
      
    public:
      
      LocalSystemMatrixType LeftHandSideContribution;
      LocalSystemVectorType RightHandSideContribution;

      //setting pointer variables
      void SetLHS_Contributions ( std::vector<LocalSystemMatrixType>& rLHS_Contributions ) { mpLHS_Contributions = &rLHS_Contributions; };
      void SetLHS_Variables     ( std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables ) { mpLHS_Variables = &rLHS_Variables; };
      void SetRHS_Contributions ( std::vector<LocalSystemVectorType>& rRHS_Contributions ) { mpRHS_Contributions = &rRHS_Contributions; };
      void SetRHS_Variables     ( std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables ) { mpRHS_Variables = &rRHS_Variables; };

      //getting pointer variables
      std::vector<LocalSystemMatrixType>& GetLHS_Contributions() { return *mpLHS_Contributions; };
      std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Variables() { return *mpLHS_Variables; };
      std::vector<LocalSystemVectorType>& GetRHS_Contributions() { return *mpRHS_Contributions; };
      std::vector< Variable< LocalSystemVectorType > >& GetRHS_Variables() { return *mpRHS_Variables; };

    }
      
    struct LHS_LocalSystemContributions
    {
    private:
    
      std::vector<LocalSystemMatrixType> *mpLHS_Contributions;
      const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Variables;

      
    public:
      
      LocalSystemMatrixType LeftHandSideContribution;
 
      //setting pointer variables
      void SetLHS_Contributions ( std::vector<LocalSystemMatrixType>& rLHS_Contributions ) { mpLHS_Contributions = &rLHS_Contributions; };
      void SetLHS_Variables     ( std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables ) { mpLHS_Variables = &rLHS_Variables; };
 
      //getting pointer variables
      std::vector<LocalSystemMatrixType>& GetLHS_Contributions() { return *mpLHS_Contributions; };
      std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Variables() { return *mpLHS_Variables; };

    }

    struct RHS_LocalSystemContributions
    {
    private:
    
      std::vector<LocalSystemVectorType> *mpRHS_Contributions;
      const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Variables;
      
    public:
      
      LocalSystemVectorType RightHandSideContribution;

      //setting pointer variables
      void SetRHS_Contributions ( std::vector<LocalSystemVectorType>& rRHS_Contributions ) { mpRHS_Contributions = &rRHS_Contributions; };
      void SetRHS_Variables     ( std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables ) { mpRHS_Variables = &rRHS_Variables; };

      //getting pointer variables
      std::vector<LocalSystemVectorType>& GetRHS_Contributions() { return *mpRHS_Contributions; };
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
	SystemContributions& rLocalSystemContributions,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //(rCurrentElement) -> InitializeNonLinearIteration(rCurrentProcessInfo);

	std::vector<LocalSystemMatrixType>& rLHS_Contributions = rSystemContributions.GetLHS_Contributions();
	std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables = rSystemContributions.GetLHS_Variables();
	std::vector<LocalSystemVectorType>& rRHS_Contributions = rSystemContributions.GetRHS_Contributions();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = rSystemContributions.GetRHS_Variables();

        (rCurrentElement) -> CalculateLocalSystem(rLHS_Contribution, rLHS_Variables, rRHS_Contribution, rRHS_Variables, rCurrentProcessInfo);

	for( unsigned int i=0; i<rLHS_Contributions.size(); i++ )
	  {	    
	    rSystemContributions.LeftHandSideContribution += rLHS_Contributions[i];
	  }
	
	for( unsigned int i=0; i<rRHS_Contributions.size(); i++ )
	  {
	    rSystemContributions.RightHandSideContribution += rRHS_Contributions[i];
	  }
	

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentElement) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        }

        (rCurrentElement) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS (rSystemContributions.LeftHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

            AddDynamicsToRHS (rCurrentElement, rSystemContributions.RightHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

        }

        //AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
	RHS_SystemContributions& rRHS_LocalSystemContributions,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {

        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

	std::vector<LocalSystemVectorType>& rRHS_Contributions = rRHS_SystemContributions.GetRHS_Contributions();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = rRHS_SystemContributions.GetRHS_Variables();

        //basic operations for the element considered
        (rCurrentElement) -> CalculateRightHandSide(rRHS_Contributions, rRHS_Variables, rCurrentProcessInfo);

	
	for( unsigned int i=0; i<rRHS_Contributions.size(); i++ )
	  {
	    rSystemContributions.RightHandSideContribution += rRHS_Contributions[i];
	  }
	

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentElement) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);
        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToRHS (rCurrentElement, rSystemContributions.RightHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
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
	SystemContributions& rLocalSystemContributions,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {


        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

	std::vector<LocalSystemMatrixType>& rLHS_Contributions = rSystemContributions.GetLHS_Contributions();
	std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables = rSystemContributions.GetLHS_Variables();
	std::vector<LocalSystemVectorType>& rRHS_Contributions = rSystemContributions.GetRHS_Contributions();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = rSystemContributions.GetRHS_Variables();

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateLocalSystem(rLHS_Contributions, rLHS_Variables, rRHS_Contributions, rRHS_Variables, rCurrentProcessInfo);


	for( unsigned int i=0; i<rLHS_Contributions.size(); i++ )
	  {	    
	    rSystemContributions.LeftHandSideContribution += rLHS_Contributions[i];
	  }
	
	for( unsigned int i=0; i<rRHS_Contributions.size(); i++ )
	  {
	    rSystemContributions.RightHandSideContribution += rRHS_Contributions[i];
	  }

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentCondition) -> MassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentCondition) -> DampMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS  (rSystemContributions.LeftHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

            AddDynamicsToRHS  (rCurrentCondition, rSystemContributions.RightHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
        }

        //AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
	RHS_SystemContributions& rRHS_LocalSystemContributions,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current condition
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

	std::vector<LocalSystemVectorType>& rRHS_Contributions = rRHS_SystemContributions.GetRHS_Contributions();
	std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = rRHS_SystemContributions.GetRHS_Variables();

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateRightHandSide(rRHS_Contributions, rRHS_Variables, rCurrentProcessInfo);
	
	for( unsigned int i=0; i<rRHS_Contributions.size(); i++ )
	  {
	    rSystemContributions.RightHandSideContribution += rRHS_Contributions[i];
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

            AddDynamicsToRHS  (rCurrentCondition, rSystemContributions.RightHandSideContribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

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


