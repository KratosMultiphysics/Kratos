//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_OMP_SET_NUM_THREADS )
#define  KRATOS_OMP_SET_NUM_THREADS


/* System includes */
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif


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



class OMPSetNumThreads

{
public:
    /**@name Type Definitions */
    /*@{ */


    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    OMPSetNumThreads(){};
  
    /** Destructor.
     */
    ~OMPSetNumThreads(){};

    /** Operators.
     */



    //**************************************************************************
    //**************************************************************************


    void SetNumThreads( int num_threads = 1)
    {

 	int nthreads,tid, procs, maxt, inpar, dynamic, nested;
  
	/* Start parallel region */
  
	/* Set thread number */  
	omp_set_num_threads(num_threads);

    #pragma omp parallel private(nthreads, tid)
	{
	  /* Obtain thread number */
	  tid = omp_get_thread_num();
    
	  /* Only master thread does this */
	  if (tid == 0)
	    {
	      printf("Thread %d getting environment info...\n", tid);
	
	      /* Get environment information */
	      procs    = omp_get_num_procs();
	      nthreads = omp_get_num_threads();
	      maxt     = omp_get_max_threads();
	      inpar    = omp_in_parallel();
	      //omp_set_dynamic(true);
	      dynamic  = omp_get_dynamic();
	      //omp_set_nested(true);
	      nested   = omp_get_nested();
	
	      /* Print environment information */
	      printf( "  | ---------------OMP--------------- |\n");
	      printf( "  | Machine number of processors  = %d |\n", procs);
	      printf( "  | Number of threads set         = %d |\n", nthreads);
	      printf( "  | Max threads in use            = %d |\n", maxt);
	      printf( "  | In parallel?                  = %d |\n", inpar);
	      printf( "  | Dynamic threads enabled?      = %d |\n", dynamic);
	      printf( "  | Nested parallelism supported? = %d |\n", nested);
	      printf( "  | --------------------------------- |\n");
	
	    }
    
	}
  
  
      
      

       
    };


    //**************************************************************************
    //**************************************************************************


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

}; /* Class ComputeLineSearch */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_OMP_SET_NUM_THREADS defined */

