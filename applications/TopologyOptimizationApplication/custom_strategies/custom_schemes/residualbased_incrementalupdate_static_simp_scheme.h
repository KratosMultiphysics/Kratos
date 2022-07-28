//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
#if !defined(KRATOS_RESIDUALBASED_INCREMENTALUPDATE_STATIC_SIMP_SCHEME_H)
#define  KRATOS_RESIDUALBASED_INCREMENTALUPDATE_STATIC_SIMP_SCHEME_H


// External includes


// Project includes
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

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

template<class TSparseSpace,class TDenseSpace >
class ResidualBasedIncrementalUpdateStaticSIMPScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedIncrementalUpdateStaticSIMPScheme );

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;


    /// The definition of the vector containing the equation ids
    typedef Element::EquationIdVectorType                              EquationIdVectorType;

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
    ~ResidualBasedIncrementalUpdateStaticSIMPScheme() override {} 


    /*@} */
    /**@name Operators
    */
    /*@{ */

    // UPDATE LHS AND RHS
    /// This function calculates a new Youngs Modulus based on the densities and multiplies it into the
    /// LHS and RHS contributions of the complete system

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemVectorType& rRHSContribution,
        EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY
        //Initializing the non linear iteration for the current element
        rCurrentElement.InitializeNonLinearIteration(rCurrentProcessInfo);

        //basic operations for the element considered
        rCurrentElement.CalculateLocalSystem(rLHSContribution,rRHSContribution,rCurrentProcessInfo);
        
        //Determine the new Youngs Modulus based on the assigned new density (X_PHYS)
        const double E_min     = rCurrentElement.GetProperties()[YOUNGS_MODULUS_MIN];
        const double E_initial = rCurrentElement.GetProperties()[YOUNGS_MODULUS_0];
        //const std::string mat_interp = rCurrentElement.GetProperties()[MAT_INTERP];
        const std::string mat_interp = rCurrentElement.GetValue(MAT_INTERP);
        const double E_current = rCurrentElement.GetValue(YOUNG_MODULUS);
        const double penalty   = rCurrentElement.GetValue(PENAL);
        const double x_phys    = rCurrentElement.GetValue(X_PHYS);

        double E_new = 0;
        if (mat_interp == "simp")
        {
            E_new += E_initial*pow(x_phys, penalty);
        }
        else if (mat_interp == "simp_modified")
        {
            E_new += (E_min + pow(x_phys, penalty) * (E_initial - E_min));
        }
        else if (mat_interp == "ramp")
        {
            E_new += E_min + x_phys/(1+penalty*(1-x_phys)) * (E_initial - E_min);
        }
        else
        {
            KRATOS_ERROR << "Material interpolation option incorrectly chosen \nAvailable methods are: 'simp', 'simp_modified', 'ramp'." << std::endl;
        }

        //Calculate the factor that needs to be multiplied on the RHS and LHS
        const double factor    = E_new/E_current;

        // Factorize LHS and RHS according SIMP approach
        // Note that when this function is called, all the contributions from the force conditions are missing.
        // I.e. RHS = -K*u_init. Hence we can directly factorize LHS and RHS to obtained the modified stiffnesses
        rLHSContribution *= factor;
        rRHSContribution *= factor;
        //Continuation of the basic operations
        rCurrentElement.EquationIdVector(rEquationId,rCurrentProcessInfo);

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

