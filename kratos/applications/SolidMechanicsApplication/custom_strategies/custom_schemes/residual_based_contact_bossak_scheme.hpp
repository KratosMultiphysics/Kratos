//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_CONTACT_BOSSAK_SCHEME )
#define  KRATOS_RESIDUAL_BASED_CONTACT_BOSSAK_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/kratos_flags.h"
#include "custom_strategies/custom_schemes/residual_based_bossak_scheme.hpp"

#include "solid_mechanics_application_variables.h"

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
class ResidualBasedContactBossakScheme: public ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>
{
public:


    /**@name Type Definitions */

    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedContactBossakScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                                    BaseType;

    typedef typename BaseType::TDataType                                       TDataType;

    typedef typename BaseType::DofsArrayType                               DofsArrayType;

    typedef typename Element::DofsVectorType                              DofsVectorType;

    typedef typename BaseType::TSystemMatrixType                       TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                       TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType               LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType               LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType                           ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                       ConditionsArrayType;
    
    typedef typename BaseType::Pointer                                   BaseTypePointer;
  
    typedef ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>          DerivedBaseType;
    
    typedef typename BaseType::LocalSystemComponents           LocalSystemComponentsType;

    /*@} */

    /**
     * Constructor.
     * The bossak method
     */
    ResidualBasedContactBossakScheme(double rAlpham=0,double rDynamic=1)
      :DerivedBaseType(rAlpham, rDynamic)
    {
    }


    /** Copy Constructor.
     */
    ResidualBasedContactBossakScheme(ResidualBasedContactBossakScheme& rOther)
      :DerivedBaseType(rOther)
    {
    }


    /** Destructor.
     */
    virtual ~ResidualBasedContactBossakScheme
    () {}

   /*@} */
    /**@name Operators
     */
    /*@{ */


    /**
     * Clone 
     */
    virtual BaseTypePointer Clone()
    {
      return BaseTypePointer( new ResidualBasedContactBossakScheme(*this) );
    }



    //***************************************************************************
    //***************************************************************************

    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo)
    {

        KRATOS_TRY

	  if( (rCurrentCondition)->Is(CONTACT) ){

	    int thread = OpenMPUtils::ThisThread();

	    //Initializing the non linear iteration for the current element
	    //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

	    //Set local system to compute condition component wise: (to get contact forces vector)
	    std::vector< Variable< LocalSystemVectorType > > RHS_Variables;
	    RHS_Variables.push_back(CONTACT_FORCES_VECTOR);

	    std::vector<LocalSystemVectorType> RHS_Components;
	    RHS_Components.resize( RHS_Variables.size() );

	    for( unsigned int i=0; i<RHS_Components.size(); i++ )
	      {
		RHS_Components[i] = LocalSystemVectorType(0);
	      }

	    //basic operations for the element considered
	    (rCurrentCondition) -> CalculateRightHandSide(RHS_Components, RHS_Variables, rCurrentProcessInfo);
	    (rCurrentCondition) -> CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

	    if( rRHS_Contribution.size() != RHS_Components[0].size() )
	      rRHS_Contribution.resize(RHS_Components[0].size());
	
	    rRHS_Contribution.clear();
	
	    for( unsigned int i=0; i<RHS_Components.size(); i++ )
	      {
		rRHS_Contribution += RHS_Components[i];

		if( RHS_Variables[i] == CONTACT_FORCES_VECTOR ){
		  //add explicit components to nodes
		  (rCurrentCondition) -> AddExplicitContribution(RHS_Components.back(), RHS_Variables.back(), CONTACT_FORCE, rCurrentProcessInfo);
		}
	      }

	    if(this->mNewmark.static_dynamic !=0)
	      {

		(rCurrentCondition) -> CalculateMassMatrix(this->mMatrix.M[thread], rCurrentProcessInfo);

		(rCurrentCondition) -> CalculateDampingMatrix(this->mMatrix.D[thread], rCurrentProcessInfo);

	      }

	    (rCurrentCondition) -> EquationIdVector(rEquationId, rCurrentProcessInfo);
	
	    if(this->mNewmark.static_dynamic !=0)
	      {

		this->AddDynamicsToLHS  (rLHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);

		this->AddDynamicsToRHS  (rCurrentCondition, rRHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);
	      }

	    //AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);
	  }
	  else{

	    DerivedBaseType::Condition_CalculateSystemContributions(rCurrentCondition, rLHS_Contribution, rRHS_Contribution, rEquationId, rCurrentProcessInfo);

	  }

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

	  if( (rCurrentCondition)->Is(CONTACT) ){

	    int thread = OpenMPUtils::ThisThread();

	    //Initializing the non linear iteration for the current condition
	    //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

	    //Set local system to compute condition component wise: (to get contact forces vector)
	    std::vector< Variable< LocalSystemVectorType > > RHS_Variables;
	    RHS_Variables.push_back(CONTACT_FORCES_VECTOR);

	    std::vector<LocalSystemVectorType> RHS_Components;
	    RHS_Components.resize( RHS_Variables.size() );

	    for( unsigned int i=0; i<RHS_Components.size(); i++ )
	      {
		RHS_Components[i] = LocalSystemVectorType(0);
	      }

	    //basic operations for the element considered
	    (rCurrentCondition) -> CalculateRightHandSide(RHS_Components, RHS_Variables, rCurrentProcessInfo);

	    if( rRHS_Contribution.size() != RHS_Components[0].size() )
	      rRHS_Contribution.resize(RHS_Components[0].size());
	
	    rRHS_Contribution.clear();
	
	    for( unsigned int i=0; i<RHS_Components.size(); i++ )
	      {
		rRHS_Contribution += RHS_Components[i];

		if( RHS_Variables[i] == CONTACT_FORCES_VECTOR ){
		  //add explicit components to nodes
		  (rCurrentCondition) -> AddExplicitContribution(RHS_Components.back(), RHS_Variables.back(), CONTACT_FORCE, rCurrentProcessInfo);
		}

	      }


	    if(this->mNewmark.static_dynamic !=0)
	      {

		(rCurrentCondition) -> CalculateMassMatrix(this->mMatrix.M[thread], rCurrentProcessInfo);

		(rCurrentCondition) -> CalculateDampingMatrix(this->mMatrix.D[thread], rCurrentProcessInfo);

	      }

	    (rCurrentCondition) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

	    //adding the dynamic contributions (static is already included)

	    if(this->mNewmark.static_dynamic !=0)
	      {

		this->AddDynamicsToRHS  (rCurrentCondition, rRHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);

	      }
	
	  }
	  else{

	    DerivedBaseType::Condition_Calculate_RHS_Contribution(rCurrentCondition, rRHS_Contribution, rEquationId, rCurrentProcessInfo);

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
}; /* Class ResidualBasedContactBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_CONTACT_BOSSAK_SCHEME defined */


