//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt ------------------------------->>> WHAT ABOUT THE LICENSE?¿
//
//  Main authors:    Ruben Zorrilla
//                    
//

#if !defined(KRATOS_CONVERGENCE_ACCELERATOR )
#define  KRATOS_CONVERGENCE_ACCELERATOR


/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "includes/define.h"
//~ #include "includes/model_part.h"


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
//~ template<class TSparseSpace,
         //~ class TDenseSpace,
         //~ class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         //~ >
template<class TSpace>
class ConvergenceAccelerator
{
public:
    /**@name Type Definitions */
    /*@{ */
    
    typedef typename TSpace::VectorType                             VectorType;
    typedef typename TSpace::MatrixType                             MatrixType;
    
    typedef typename TSpace::VectorPointerType               VectorPointerType;
    typedef typename TSpace::MatrixPointerType               MatrixPointerType;
    
    //~ typedef typename TSpace::MatrixType::size_type                    SizeType;
        
    //~ /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ConvergenceAccelerator );

    /*@} */

    /** Constructor.
     */

    /*@{ */
    ConvergenceAccelerator()
    {
    }
    /*@} */
    
    /** Copy constructor.
    */
    
    /*@{ */
    ConvergenceAccelerator(const ConvergenceAccelerator& Other);
    /*@{ */

    /** Destructor.
     */

    /*@{ */
    virtual ~ConvergenceAccelerator()
    {
    }
    /*@} */

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/
    /*@{ */

    /**
    Initialization of member variables and prior operations. Executed once.
     */
    virtual void Initialize()
    {
    }

    /**
    Performs all the required operations that should be done (for each step) before solving the solution step.
	A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep()
    {
    }

    /**
    Performs all the required operations that should be done at the beginning of each non-linear iteration.
     */
    virtual void InitializeNonLinearIteration()
    {
    }

    /**
    Computes the correction over the given iteration guess
     */
    virtual void UpdateSolution(
        const VectorType& rResidualVector,
        VectorType& rIterationGuess)
    {
    }
    
    /**
    Performs all the required operations that should be done at the end of each non-linear iteration.
     */
    virtual void FinalizeNonLinearIteration()
    {
    }

	/**
    Performs all the required operations that should be done (for each step) after solving the solution step.
	A member variable should be used as a flag to make sure this function is called only once per step.
    */
	virtual void FinalizeSolutionStep()
	{
	}

    /**
    Clears the internal storage
     */
    virtual void Clear()
    {
    }
	
    //*********************************************************************************

    /**level of echo for the convergence accelerator
    0 -> mute... no echo at all
    1 -> printing time and basic informations
    2 -> printing linear solver data
    3 -> Print of debug informations:
     */
    virtual void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    int GetEchoLevel()
    {
        return mEchoLevel;
    }

    //*********************************************************************************

    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    //level of echo for the convergence accelerator
    int mEchoLevel;

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

private:

    /*@} */
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    //~ ModelPart& mr_model_part;

    //~ bool mMoveMeshFlag;


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

}; /* Class ConvergenceAccelerator */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_CONVERGENCE_ACCELERATOR  defined */

