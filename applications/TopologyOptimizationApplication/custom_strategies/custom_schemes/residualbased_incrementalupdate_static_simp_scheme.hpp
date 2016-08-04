//
//   Project Name:        $Project:     Topology_Optimization_Application $
//   Last modified by:    $Author:      Malfavón Farías, Baumgärtner      $
//   Date:                $Date:        May 2016                          $
//   Revision:            $Revision:    0.0                               $
//


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

    //typedef boost::shared_ptr< ResidualBasedIncrementalUpdateStaticSIMPScheme<TSparseSpace,TDenseSpace> > Pointer;
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
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    )
    {
        KRATOS_TRY
        //Initializing the non linear iteration for the current element
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentElement)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        //Determine the new Youngs Modulus based on the assigned new density (X_PHYS)
        double E_min     = (rCurrentElement)->GetValue(E_MIN);
        double E_initial = (rCurrentElement)->GetValue(E_0);
        double E_current = (rCurrentElement)->GetValue(YOUNG_MODULUS);
        double penalty   = (rCurrentElement)->GetValue(PENAL);
        double x_phys    = (rCurrentElement)->GetValue(X_PHYS);

        double E_new     = (E_min + pow(x_phys, penalty) * (E_initial - E_min));

        //Calculate the factor that needs to be multiplied on the RHS and LHS respectively
        double factor    = (1/E_current)*E_new;
        LHS_Contribution = LHS_Contribution * factor;
        RHS_Contribution = RHS_Contribution * factor;

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

