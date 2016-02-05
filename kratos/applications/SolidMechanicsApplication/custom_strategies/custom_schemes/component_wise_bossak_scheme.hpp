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

/* Project includes */
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
class ComponentWiseBossakScheme: public ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>
{
public:

    /**@name Type Definitions */

    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( ComponentWiseBossakScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                            BaseType;

    typedef typename BaseType::TDataType                               TDataType;

    typedef typename BaseType::DofsArrayType                       DofsArrayType;

    typedef typename Element::DofsVectorType                      DofsVectorType;

    typedef typename BaseType::TSystemMatrixType               TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType               TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType       LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType       LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType                   ElementsArrayType;

    typedef ModelPart::ConditionsContainerType               ConditionsArrayType;

    typedef typename BaseType::LocalSystemComponents   LocalSystemComponentsType;

    typedef typename BaseType::Pointer                           BaseTypePointer;

    typedef ResidualBasedBossakScheme<TSparseSpace,TDenseSpace>  DerivedBaseType;

    /*@} */

    /**
     * Constructor.
     * The bossak method
     */
    ComponentWiseBossakScheme(double rAlpham=0,double rDynamic=1)
      :DerivedBaseType(rAlpham,rDynamic)
    {
      mLocalSystem.Initialize();
      //std::cout<<" component wise bossak scheme selected "<<std::endl;
    }


    /** Copy Constructor.
     */
    ComponentWiseBossakScheme(ComponentWiseBossakScheme& rOther)
      :DerivedBaseType(rOther)
    {
      mLocalSystem.Initialize();

      //set local component variables for the elements
      mLocalSystem.SetLHS_Element_Variables(rOther.mLocalSystem.GetLHS_Element_Variables());
      mLocalSystem.SetRHS_Element_Variables(rOther.mLocalSystem.GetRHS_Element_Variables());   

      //set local component variables for the conditions
      mLocalSystem.SetLHS_Condition_Variables(rOther.mLocalSystem.GetLHS_Condition_Variables());
      mLocalSystem.SetRHS_Condition_Variables(rOther.mLocalSystem.GetRHS_Condition_Variables());   
      
      //element components and condition components must be set by the builder_and_solver
        
    }

    /** Destructor.
     */
    virtual ~ComponentWiseBossakScheme() 
    {}


    /*@} */
    /**@name Operators
     */
    /*@{ */


    /**
     * Clone 
     */
    BaseTypePointer Clone()
    {
      return BaseTypePointer( new ComponentWiseBossakScheme(*this) );
    }


    LocalSystemComponentsType& GetLocalSystemComponents()
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

	bool LHS_Element_Components_Set = mLocalSystem.Are_LHS_Element_Components_Set();
	std::vector<LocalSystemMatrixType>& rLHS_Components = mLocalSystem.GetLHS_Element_Components();
	const std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables = mLocalSystem.GetLHS_Element_Variables();

	bool RHS_Element_Components_Set = mLocalSystem.Are_RHS_Element_Components_Set();
	std::vector<LocalSystemVectorType>& rRHS_Components = mLocalSystem.GetRHS_Element_Components();
	const std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = mLocalSystem.GetRHS_Element_Variables();


	if( LHS_Element_Components_Set && RHS_Element_Components_Set )
	  {

	    //basic operations for the element considered
	    (rCurrentElement) -> CalculateLocalSystem(rLHS_Components, rLHS_Variables, rRHS_Components, rRHS_Variables, rCurrentProcessInfo);

	    if( rLHS_Contribution.size1() != rLHS_Components[0].size1() )
	      rLHS_Contribution.resize(rLHS_Components[0].size1(), rLHS_Components[0].size2());	
	    
	    rLHS_Contribution.clear();

	    for( unsigned int i=0; i<rLHS_Components.size(); i++ )
	      {	    
		rLHS_Contribution += rLHS_Components[i];
	      }
	
	    if( rRHS_Contribution.size() != rRHS_Components[0].size() )
	      rRHS_Contribution.resize(rRHS_Components[0].size());

	    rRHS_Contribution.clear();

	    for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	      {
		rRHS_Contribution += rRHS_Components[i];
	      }

	  }
	else
	  {
	    if( !LHS_Element_Components_Set && RHS_Element_Components_Set )
	      {
		(rCurrentElement) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);

		(rCurrentElement) -> CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
		
		if( rRHS_Contribution.size() != rRHS_Components[0].size() )
		  rRHS_Contribution.resize(rRHS_Components[0].size());

		rRHS_Contribution.clear();

		for( unsigned int i=0; i<rRHS_Components.size(); i++ )
		  {
		    rRHS_Contribution += rRHS_Components[i];
		  }
		

	      }
	    else if ( LHS_Element_Components_Set && !RHS_Element_Components_Set )
	      {
	
		KRATOS_THROW_ERROR( std::logic_error, " scheme asks for a unusual Element LHS components not implemented ", "" )

		(rCurrentElement) -> CalculateLeftHandSide(rLHS_Components, rLHS_Variables, rCurrentProcessInfo);

		(rCurrentElement) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
		
		if( rLHS_Contribution.size1() != rLHS_Components[0].size1() )
		  rLHS_Contribution.resize(rLHS_Components[0].size1(), rLHS_Components[0].size2());	

		rLHS_Contribution.clear();

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
	

        if(this->mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> CalculateMassMatrix(this->mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentElement) -> CalculateDampingMatrix(this->mMatrix.D[thread], rCurrentProcessInfo);

        }

        (rCurrentElement) -> EquationIdVector(rEquationId, rCurrentProcessInfo);

        if(this->mNewmark.static_dynamic !=0)
        {

            this->AddDynamicsToLHS (rLHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);

            this->AddDynamicsToRHS (rCurrentElement, rRHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);

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

	bool RHS_Element_Components_Set = mLocalSystem.Are_RHS_Element_Components_Set();
	std::vector<LocalSystemVectorType>& rRHS_Components = mLocalSystem.GetRHS_Element_Components();
	const std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = mLocalSystem.GetRHS_Element_Variables();

	if( RHS_Element_Components_Set )
	  {
	    //basic operations for the element considered
	    (rCurrentElement) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);
	
	    if( rRHS_Contribution.size() != rRHS_Components[0].size() )
	      rRHS_Contribution.resize(rRHS_Components[0].size());

	    rRHS_Contribution.clear();

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
	

        if(this->mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> CalculateMassMatrix(this->mMatrix.M[thread], rCurrentProcessInfo);

            (rCurrentElement) -> CalculateDampingMatrix(this->mMatrix.D[thread], rCurrentProcessInfo);
        }

        (rCurrentElement) -> EquationIdVector(rEquationId,rCurrentProcessInfo);

        if(this->mNewmark.static_dynamic !=0)
        {

            this->AddDynamicsToRHS (rCurrentElement, rRHS_Contribution, this->mMatrix.D[thread], this->mMatrix.M[thread], rCurrentProcessInfo);
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

	bool LHS_Condition_Components_Set = mLocalSystem.Are_LHS_Condition_Components_Set();
	std::vector<LocalSystemMatrixType>& rLHS_Components = mLocalSystem.GetLHS_Condition_Components();
	const std::vector< Variable< LocalSystemMatrixType > >& rLHS_Variables = mLocalSystem.GetLHS_Condition_Variables();

	bool RHS_Condition_Components_Set = mLocalSystem.Are_RHS_Condition_Components_Set();
	std::vector<LocalSystemVectorType>& rRHS_Components = mLocalSystem.GetRHS_Condition_Components();
	const std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = mLocalSystem.GetRHS_Condition_Variables();

	if( LHS_Condition_Components_Set && RHS_Condition_Components_Set )
	  {

	    //basic operations for the element considered
	    (rCurrentCondition) -> CalculateLocalSystem(rLHS_Components, rLHS_Variables, rRHS_Components, rRHS_Variables, rCurrentProcessInfo);

	    if( rLHS_Contribution.size1() != rLHS_Components[0].size1() )
	      rLHS_Contribution.resize(rLHS_Components[0].size1(), rLHS_Components[0].size2());	

	    rLHS_Contribution.clear();

	    for( unsigned int i=0; i<rLHS_Components.size(); i++ )
	      {	    
		rLHS_Contribution += rLHS_Components[i];
	      }
	
	    if( rRHS_Contribution.size() != rRHS_Components[0].size() )
	      rRHS_Contribution.resize(rRHS_Components[0].size());

	    rRHS_Contribution.clear();
		
	    for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	      {
		rRHS_Contribution += rRHS_Components[i];

		if( rRHS_Variables[i] == CONTACT_FORCES_VECTOR ){
		  //add explicit components to nodes
		  (rCurrentCondition) -> AddExplicitContribution(rRHS_Components[i], rRHS_Variables[i], CONTACT_FORCE, rCurrentProcessInfo);
		}
	      }

	  }
	else
	  {
	    if( !LHS_Condition_Components_Set && RHS_Condition_Components_Set )
	      {
		(rCurrentCondition) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);

		(rCurrentCondition) -> CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
		
		if( rRHS_Contribution.size() != rRHS_Components[0].size() )
		  rRHS_Contribution.resize(rRHS_Components[0].size());

		rRHS_Contribution.clear();

		for( unsigned int i=0; i<rRHS_Components.size(); i++ )
		  {
		    rRHS_Contribution += rRHS_Components[i];

		    if( rRHS_Variables[i] == CONTACT_FORCES_VECTOR ){
		      //add explicit components to nodes
		      (rCurrentCondition) -> AddExplicitContribution(rRHS_Components[i], rRHS_Variables[i], CONTACT_FORCE, rCurrentProcessInfo);
		    }
		    
		    //std::cout<<" Condition ["<<(rCurrentCondition) -> Id()<<"] :"<<rRHS_Contribution<<" and "<<rRHS_Components[i]<<std::endl;
		  }

		

	      }
	    else if ( LHS_Condition_Components_Set && !RHS_Condition_Components_Set )
	      {
		
		KRATOS_THROW_ERROR( std::logic_error, " scheme asks for a unusual Condition LHS components not implemented ", "" )

		(rCurrentCondition) -> CalculateLeftHandSide(rLHS_Components, rLHS_Variables, rCurrentProcessInfo);

		(rCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
		
		if( rLHS_Contribution.size1() != rLHS_Components[0].size1() )
		  rLHS_Contribution.resize(rLHS_Components[0].size1(), rLHS_Components[0].size2());	
		
		rLHS_Contribution.clear();
	    
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

	bool RHS_Condition_Components_Set = mLocalSystem.Are_RHS_Condition_Components_Set();
	std::vector<LocalSystemVectorType>& rRHS_Components = mLocalSystem.GetRHS_Condition_Components();
	const std::vector< Variable< LocalSystemVectorType > >& rRHS_Variables = mLocalSystem.GetRHS_Condition_Variables();

	if( RHS_Condition_Components_Set )
	  {
	    //basic operations for the element considered
	    (rCurrentCondition) -> CalculateRightHandSide(rRHS_Components, rRHS_Variables, rCurrentProcessInfo);
	
	    if( rRHS_Contribution.size() != rRHS_Components[0].size() )
	      rRHS_Contribution.resize(rRHS_Components[0].size());

	    rRHS_Contribution.clear();

	    for( unsigned int i=0; i<rRHS_Components.size(); i++ )
	      {
		rRHS_Contribution += rRHS_Components[i];

		if( rRHS_Variables[i] == CONTACT_FORCES_VECTOR ){
		  //add explicit components to nodes
		  (rCurrentCondition) -> AddExplicitContribution(rRHS_Components[i], rRHS_Variables[i], CONTACT_FORCE, rCurrentProcessInfo);
		}

	      }
	  }
	else
	  {
	    //basic operations for the element considered
	    (rCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);
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
    LocalSystemComponentsType mLocalSystem;
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


