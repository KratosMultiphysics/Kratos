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
//

#if !defined(KRATOS_OR_CRITERIA_H)
#define  KRATOS_OR_CRITERIA_H


/* System includes */


/* External includes */


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

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
         class TDenseSpace
         >
class Or_Criteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    /**@name Type Definitions */
    /*@{ */

    /** Counted pointer of Or_Criteria */


    KRATOS_CLASS_POINTER_DEFINITION(Or_Criteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;


    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    Or_Criteria
    (
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer first_criterion,
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer second_criterion)
        :ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mpfirst_criterion   =  first_criterion;
        mpsecond_criterion  =  second_criterion;
    }

    /** Copy constructor.
    */
    Or_Criteria(Or_Criteria const& rOther)
      :BaseType(rOther)
     {
       mpfirst_criterion   =  rOther.mpfirst_criterion;
       mpsecond_criterion  =  rOther.mpsecond_criterion;      
     }

    /** Destructor.
    */
    virtual ~Or_Criteria () {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    /**level of echo for the convergence criterion
    0 -> mute... no echo at all
    1 -> print basic informations
    2 -> print extra informations
     */
    void SetEchoLevel(int Level)
    {
      BaseType::SetEchoLevel(Level);
      mpfirst_criterion->SetEchoLevel(Level);
      mpsecond_criterion->SetEchoLevel(Level);
    }


    /*Criteria that need to be called after getting the solution */
    bool PostCriteria(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
    )
    {
        bool first_criterion_result  = mpfirst_criterion ->PostCriteria(r_model_part,rDofSet,A,Dx,b);
        bool second_criterion_result = mpsecond_criterion ->PostCriteria(r_model_part,rDofSet,A,Dx,b);

        return (first_criterion_result || second_criterion_result);

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
    typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer mpfirst_criterion;
    typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer mpsecond_criterion;


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

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_AND_CRITERIA  defined */

