// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_RESIDUAL_BASED_STATIC_SIMP_SCHEME )
#define  KRATOS_RESIDUAL_BASED_STATIC_SIMP_SCHEME


// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

//#include "solid_mechanics_application.h"
#include "topology_optimization_application.h"

// Application includes
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


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

template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class ResidualBasedIncrementalUpdateStaticSIMPScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    //typedef Kratos::shared_ptr< ResidualBasedIncrementalUpdateStaticSIMPScheme<TSparseSpace,TDenseSpace> > Pointer;
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedIncrementalUpdateStaticSIMPScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ResidualBasedIncrementalUpdateStaticSIMPScheme()
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>()
    {}

    /** Destructor.
    */
    virtual ~ResidualBasedIncrementalUpdateStaticSIMPScheme() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- UPDATE LHS AND RHS ----------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// This function calculates a new Youngs Modulus based on the densities and multiplies it into the
    /// LHS and RHS contributions of the complete system
    virtual void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LSH,
        LocalSystemVectorType& RHS,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    )
    {
    	KRATOS_TRY
		//Initializing the non linear iteration for the current element
		(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

    	//basic operations for the element considered
    	(rCurrentElement)->CalculateLocalSystem(LSH,RHS,CurrentProcessInfo);


    	//Determine the new Youngs Modulus based on the assigned new density (X_PHYS)
    	double E_min     = (rCurrentElement)->GetValue(E_MIN);
    	double E_initial = (rCurrentElement)->GetValue(E_0);
    	double E_current = (rCurrentElement)->GetValue(YOUNG_MODULUS);
    	double penalty   = (rCurrentElement)->GetValue(PENAL);
    	double x_phys    = (rCurrentElement)->GetValue(X_PHYS);

    	double E_new     = (E_min + pow(x_phys, penalty) * (E_initial - E_min));

    	//Calculate the factor that needs to be multiplied on the RHS and LHS
    	double factor    = (1/E_current)*E_new;

    	// Factorize LHS and RHS according SIMP approach
    	// Note that when this function is called, all the contributions from the force conditions are missing.
    	// I.e. RHS = -K*u_init. Hence we can directly factorize LHS and RHS to obtained the modified stiffnesses
    	LSH *= factor;
    	RHS *= factor;

    	//Continuation of the basic operations
    	(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

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

}; /* Class Scheme */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_STATIC_SIMP_SCHEME  defined */

