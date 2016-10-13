//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_SCHEME )
#define  KRATOS_SCHEME


/* System includes */
#include <set>

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

#include "includes/condition.h"
//#include "vectorial_spaces/vector.h"
//#include "vectorial_spaces/matrix.h"

//default dense space
//#include "vectorial_spaces/dense_space.h"

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

/** Short class definition.

  This class provides the implementation of the basic tasks that are needed by the solution strategy.
  It is intended to be the place for tailoring the solution strategies to problem specific tasks.

        Detail class definition.

          \URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

                \URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

                  \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

                        \URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


                                \URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

                                  \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

                                        \URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

                                          \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


 */
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class Scheme
{
public:

    /**@name Type Definitions */
    /*@{ */
    typedef typename TSparseSpace::DataType TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    typedef Dof<double> TDofType;

    typedef ModelPart::DofsArrayType DofsArrayType;
    typedef typename PointerVectorSet<TDofType, IndexedObject>::iterator DofIterator;
    typedef typename PointerVectorSet<TDofType, IndexedObject>::const_iterator DofConstantIterator;

    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    //typedef PointerVectorSet<TDofType, IndexedObject> DofsArrayType;
    //typedef Node::DofArrayType DofArrayType;



    /**
     * This struct is used in the component wise calculation only
     * is defined here and is used to declare a member variable in the component wise schemes
     * private pointers can only be accessed by means of set and get functions
     * this allows to set and not copy the Element_Variables and Condition_Variables
     * which will be asked and set by another strategy object
     */

    struct LocalSystemComponents
    {
    private:

      //elements
      std::vector<LocalSystemMatrixType> *mpLHS_Element_Components;
      const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Element_Variables;

      std::vector<LocalSystemVectorType> *mpRHS_Element_Components;
      const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Element_Variables;
      
      //conditions
      std::vector<LocalSystemMatrixType> *mpLHS_Condition_Components;
      const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Condition_Variables;

      std::vector<LocalSystemVectorType> *mpRHS_Condition_Components;
      const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Condition_Variables;
      
    public:
      
      void Initialize()
      {
	mpLHS_Element_Components = NULL;
	mpLHS_Element_Variables  = NULL;

	mpRHS_Element_Components = NULL;
	mpRHS_Element_Variables  = NULL;

	mpLHS_Condition_Components = NULL;
	mpLHS_Condition_Variables  = NULL;

	mpRHS_Condition_Components = NULL;
	mpRHS_Condition_Variables  = NULL;
      }

      //setting pointer variables

      //elements
      void SetLHS_Element_Components ( std::vector<LocalSystemMatrixType>& rLHS_Element_Components ) { mpLHS_Element_Components = &rLHS_Element_Components; };
      void SetLHS_Element_Variables     ( const std::vector< Variable< LocalSystemMatrixType > >& rLHS_Element_Variables ) { mpLHS_Element_Variables = &rLHS_Element_Variables; };
      void SetRHS_Element_Components ( std::vector<LocalSystemVectorType>& rRHS_Element_Components ) { mpRHS_Element_Components = &rRHS_Element_Components; };
      void SetRHS_Element_Variables     ( const std::vector< Variable< LocalSystemVectorType > >& rRHS_Element_Variables ) { mpRHS_Element_Variables = &rRHS_Element_Variables; };

      bool Are_LHS_Element_Components_Set() { if( mpLHS_Element_Variables == NULL ) return false; else return true; };
      bool Are_RHS_Element_Components_Set() { if( mpRHS_Element_Variables == NULL ) return false; else return true; };

      //conditions
      void SetLHS_Condition_Components ( std::vector<LocalSystemMatrixType>& rLHS_Condition_Components ) { mpLHS_Condition_Components = &rLHS_Condition_Components; };
      void SetLHS_Condition_Variables     ( const std::vector< Variable< LocalSystemMatrixType > >& rLHS_Condition_Variables ) { mpLHS_Condition_Variables = &rLHS_Condition_Variables; };
      void SetRHS_Condition_Components ( std::vector<LocalSystemVectorType>& rRHS_Condition_Components ) { mpRHS_Condition_Components = &rRHS_Condition_Components; };
      void SetRHS_Condition_Variables     ( const std::vector< Variable< LocalSystemVectorType > >& rRHS_Condition_Variables ) { mpRHS_Condition_Variables = &rRHS_Condition_Variables; };

      bool Are_LHS_Condition_Components_Set() { if( mpLHS_Condition_Variables == NULL ) return false; else return true; };
      bool Are_RHS_Condition_Components_Set() { if( mpRHS_Condition_Variables == NULL ) return false; else return true; };

      //getting pointer variables

      //elements
      std::vector<LocalSystemMatrixType>& GetLHS_Element_Components() { return *mpLHS_Element_Components; };
      const std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Element_Variables() { return *mpLHS_Element_Variables; };
      std::vector<LocalSystemVectorType>& GetRHS_Element_Components() { return *mpRHS_Element_Components; };
      const std::vector< Variable< LocalSystemVectorType > >& GetRHS_Element_Variables() { return *mpRHS_Element_Variables; };

      //conditions
      std::vector<LocalSystemMatrixType>& GetLHS_Condition_Components() { return *mpLHS_Condition_Components; };
      const std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Condition_Variables() { return *mpLHS_Condition_Variables; };
      std::vector<LocalSystemVectorType>& GetRHS_Condition_Components() { return *mpRHS_Condition_Components; };
      const std::vector< Variable< LocalSystemVectorType > >& GetRHS_Condition_Variables() { return *mpRHS_Condition_Variables; };

    };

    //pointer definition

    KRATOS_CLASS_POINTER_DEFINITION(Scheme);


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    Scheme()
    {
        mSchemeIsInitialized = false;
        mElementsAreInitialized = false;
	mConditionsAreInitialized = false;
    }


    /** Copy Constructor.
     */
    Scheme(Scheme& rOther)
      :mSchemeIsInitialized(rOther.mSchemeIsInitialized)
      ,mElementsAreInitialized(rOther.mElementsAreInitialized)
      ,mConditionsAreInitialized(rOther.mConditionsAreInitialized)
    {
    }


    /** Destructor.
     */
    virtual ~Scheme()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    /**
     * Clone 
     */
    virtual Pointer Clone()
    {
      return Pointer( new Scheme(*this) );
    }

    /**
     * Component wise components Get method
     */

    virtual LocalSystemComponents& GetLocalSystemComponents()
    {
      KRATOS_THROW_ERROR(std::logic_error, "Asking for Local Components to the SCHEME base class which is not component wise and not contains this member variable","")
    }

    /**
    this is the place to initialize the Scheme.
    This is intended to be called just once when the strategy is initialized
     */
    virtual void Initialize(
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY
        mSchemeIsInitialized = true;
        KRATOS_CATCH("")
    }

    bool SchemeIsInitialized()
    {
        return mSchemeIsInitialized;
    }

    bool ElementsAreInitialized()
    {
        return mElementsAreInitialized;
    }

    bool ConditionsAreInitialized()
    {
        return mConditionsAreInitialized;
    }

    void SetSchemeIsInitialized(bool SchemeIsInitializedFlag = true)
    {
        mSchemeIsInitialized = SchemeIsInitializedFlag;
    }

    void SetElementsAreInitialized(bool ElementsAreInitializedFlag = true)
    {
        mElementsAreInitialized = ElementsAreInitializedFlag;
    }

    void SetConditionsAreInitialized(bool ConditionsAreInitializedFlag = true)
    {
        mConditionsAreInitialized = ConditionsAreInitializedFlag;
    }

    /**
    this is the place to initialize the elements.
    This is intended to be called just once when the strategy is initialized
     */
    virtual void InitializeElements(
        ModelPart& rModelPart)
    {
        KRATOS_TRY
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++)
        {
            auto itElem = rModelPart.ElementsBegin() + i;
            itElem->Initialize(); 
        }

//         int NumThreads = OpenMPUtils::GetNumThreads();
//         OpenMPUtils::PartitionVector ElementPartition;
//         OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), NumThreads, ElementPartition);
// 
//         #pragma omp parallel
//         {
//             int k = OpenMPUtils::ThisThread();
//             ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin() + ElementPartition[k];
//             ElementsArrayType::iterator ElemEnd = rModelPart.Elements().begin() + ElementPartition[k + 1];
// 
//             for (ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)
//             {
//                 itElem->Initialize(); //function to initialize the element
//             }
// 
//         }

        mElementsAreInitialized = true;

        KRATOS_CATCH("")
    }


    /**
    this is the place to initialize the conditions.
    This is intended to be called just once when the strategy is initialized
    */
    virtual void InitializeConditions(
        ModelPart& rModelPart)
    {
        KRATOS_TRY

        if(mElementsAreInitialized==false)
	    KRATOS_THROW_ERROR(std::logic_error, "Before initilizing Conditions, initialize Elements FIRST","")
            
                #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++)
        {
            auto itElem = rModelPart.ConditionsBegin() + i;
            itElem->Initialize(); 
        }

//         int NumThreads = OpenMPUtils::GetNumThreads();
//         OpenMPUtils::PartitionVector ConditionPartition;
//         OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), NumThreads, ConditionPartition);
// 
//         #pragma omp parallel
//         {
//             int k = OpenMPUtils::ThisThread();
//             ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin() + ConditionPartition[k];
//             ConditionsArrayType::iterator CondEnd = rModelPart.Conditions().begin() + ConditionPartition[k + 1];
// 
//             for (ConditionsArrayType::iterator itCond = CondBegin; itCond != CondEnd; itCond++)
//             {
//                 itCond->Initialize(); //function to initialize the condition
//             }
// 
//         }

        mConditionsAreInitialized = true;
        KRATOS_CATCH("")
    }

    /**
    Function called once at the beginning of each solution step.
    The basic operations to be carried in there are the following:
    - managing variables to be kept constant over the time step
    (for example time-Scheme constants depending on the actual time step)
     */
    virtual void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        //initialize solution step for all of the elements
        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            (it) -> InitializeSolutionStep(CurrentProcessInfo);
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            (it) -> InitializeSolutionStep(CurrentProcessInfo);
        }
        KRATOS_CATCH("")
    }

    /**
    function called once at the end of a solution step, after convergence is reached if
    an iterative process is needed
     */
    virtual void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        //finalizes solution step for all of the elements
//         ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++)
        {
            auto itElem = rModelPart.ElementsBegin() + i;
            itElem->FinalizeSolutionStep(CurrentProcessInfo);
        }
        
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++)
        {
            auto itElem = rModelPart.ConditionsBegin() + i;
            itElem->FinalizeSolutionStep(CurrentProcessInfo);
        }        

//         int NumThreads = OpenMPUtils::GetNumThreads();
//         OpenMPUtils::PartitionVector ElementPartition;
//         OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);
// 
//         #pragma omp parallel
//         {
//             int k = OpenMPUtils::ThisThread();
// 
//             ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
//             ElementsArrayType::iterator ElementsEnd = rElements.begin() + ElementPartition[k + 1];
// 
//             for (ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
//             {
//                 itElem->FinalizeSolutionStep(CurrentProcessInfo);
//             }
//         }

//         ConditionsArrayType& rConditions = rModelPart.Conditions();
// 
//         OpenMPUtils::PartitionVector ConditionPartition;
//         OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
// 
//         #pragma omp parallel
//         {
//             int k = OpenMPUtils::ThisThread();
// 
//             ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
//             ConditionsArrayType::iterator ConditionsEnd = rConditions.begin() + ConditionPartition[k + 1];
// 
//             for (ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
//             {
//                 itCond->FinalizeSolutionStep(CurrentProcessInfo);
//             }
//         }
        KRATOS_CATCH("")
    }

    /**
    Completely analogous to the precedent,
    to be used when system is not explicitely defined, for example for fractional step
    strategies
     */
    /*		virtual void InitializeSolutionStep(
                            ModelPart& r_model_part
                            )
                    {
                            KRATOS_TRY
                            KRATOS_CATCH("")
                    }
     */
    /**
    Completely analogous to the precedent,
    to be used when system is not explicitely defined, for example for fractional step
    strategies
     */
    /*		virtual void FinalizeSolutionStep(
                            ModelPart& r_model_part
                            )
                    {
                            KRATOS_TRY
                            KRATOS_CATCH("")
                    }
     */
    /**
    executed before each fractional step
     */
    /*		virtual void InitializeFractionalSolutionStep(
                            ModelPart& r_model_part
                            )
                    {
                            KRATOS_TRY
                            KRATOS_CATCH("")
                    }
     */
    /**
    executed after each fractional step
     */
    /*		virtual void FinalizeFractionalSolutionStep(
                            ModelPart& r_model_part
                            )
                    {
                            KRATOS_TRY
                            KRATOS_CATCH("")
                    }
     */

    /**
    function to be called when it is needed to initialize an iteration.
    it is designed to be called at the beginning of each non linear iteration

      take care: the elemental function with the same name is NOT called here.
      The function is called in the builder for memory efficiency
     */
    virtual void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }
    
    virtual void InitializeNonLinearIteration(Condition::Pointer rCurrentCondition,
					      ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    virtual void InitializeNonLinearIteration(Element::Pointer rCurrentElement,
					      ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }


    /**
    function to be called when it is needed to finalize an iteration.
    it is designed to be called at the end of each non linear iteration
     */
    virtual void FinalizeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            (it) -> FinalizeNonLinearIteration(CurrentProcessInfo);
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            (it) -> FinalizeNonLinearIteration(CurrentProcessInfo);
        }
        KRATOS_CATCH("")
    }

    /**
    Performing the prediction of the solution.
     */
    virtual void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
    Performing the update of the solution.
     */
    virtual void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
    functions to be called to prepare the data needed for the output of results.
     */
    virtual void CalculateOutputData(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
    functions that cleans the results data.
     */
    virtual void CleanOutputData()
    {
    }

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed after the end of the solution step
     */
    virtual void Clean()
    {
    }

    /**
    function to clean up "element" scratch space after each element is built.
     */
    virtual void CleanMemory(Element::Pointer rCurrentElement)
    {
        rCurrentElement->CleanMemory();
    }

    /**
    function to clean up "condition" scratch space after each condition is built.
     */
    virtual void CleanMemory(Condition::Pointer rCurrentCondition)
    {
        rCurrentCondition->CleanMemory();
    }

    /**
     * Liberate internal storage.
     */
    virtual void Clear()
    {

    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    virtual int Check(ModelPart& r_model_part)
    {
        KRATOS_TRY

        for(ModelPart::ElementsContainerType::iterator it=r_model_part.ElementsBegin();
                it!=r_model_part.ElementsEnd(); it++)
        {
            it->Check(r_model_part.GetProcessInfo());
        }

        for(ModelPart::ConditionsContainerType::iterator it=r_model_part.ConditionsBegin();
                it!=r_model_part.ConditionsEnd(); it++)
        {
            it->Check(r_model_part.GetProcessInfo());
        }

        return 0;
        KRATOS_CATCH("");
    }

    /** this function is designed to be called in the builder and solver to introduce
    the selected time integration scheme. It "asks" the matrix needed to the element and
    performs the operations needed to introduce the seected time integration scheme.

      this function calculates at the same time the contribution to the LHS and to the RHS
      of the system
     */
    virtual void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
    }

    virtual void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
    }

    virtual void Calculate_LHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
    }

    virtual void EquationId(
        Element::Pointer rCurrentElement,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);
    }

    /** functions totally analogous to the precedent but applied to
    the "condition" objects
     */
    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
    }

    virtual void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
    }

    virtual void Condition_Calculate_LHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
    }

    virtual void Condition_EquationId(
        Condition::Pointer rCurrentCondition,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        (rCurrentCondition)->EquationIdVector(EquationId, CurrentProcessInfo);
    }

    /** Function that returns the list of Degrees of freedom to be
    assembled in the system for a Given Element
     */
    virtual void GetElementalDofList(
        Element::Pointer rCurrentElement,
        Element::DofsVectorType& ElementalDofList,
        ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement->GetDofList(ElementalDofList, CurrentProcessInfo);
    }

    /** Function that returns the list of Degrees of freedom to be
    assembled in the system for a Given Element
     */
    virtual void GetConditionDofList(
        Condition::Pointer rCurrentCondition,
        Element::DofsVectorType& ConditionDofList,
        ProcessInfo& CurrentProcessInfo)
    {
        rCurrentCondition->GetDofList(ConditionDofList, CurrentProcessInfo);
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


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /// flag to be used in controlling if the Scheme has been intialized or not
    bool mSchemeIsInitialized;

    /// flag taking in account if the elements were initialized correctly or not
    bool mElementsAreInitialized;

    /// flag taking in account if the conditions were initialized correctly or not
    bool mConditionsAreInitialized;

    /** Pointer to the Model.
     */

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


    /*@} */

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_SCHEME  defined */

