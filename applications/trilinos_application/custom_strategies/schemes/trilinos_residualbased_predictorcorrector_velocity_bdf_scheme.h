/*
==============================================================================
KratosTrilinosApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

#ifndef KRATOS_TRILINOS_RESIDUALBASED_PREDICTORCORRECTOR_VELOCITY_BDF2_SCHEME_H
#define KRATOS_TRILINOS_RESIDUALBASED_PREDICTORCORRECTOR_VELOCITY_BDF2_SCHEME_H


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "../../../FluidDynamicsApplication/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bdf_scheme_turbulent.h"
#include "includes/variables.h"
#include "containers/array_1d.h"

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

/// Trilinos version of ResidualBasedPredictorCorrectorBDFSchemeTurbulent.
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class TrilinosResidualBasedPredictorCorrectorBDFScheme
    : public ResidualBasedPredictorCorrectorBDFSchemeTurbulent<TSparseSpace, TDenseSpace>
{
public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION(TrilinosResidualBasedPredictorCorrectorBDFScheme );

    typedef ResidualBasedPredictorCorrectorBDFSchemeTurbulent<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename Element::DofsVectorType DofsVectorType;

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
    TrilinosResidualBasedPredictorCorrectorBDFScheme(unsigned int DomainSize,
                                                     Variable<double>& rSlipVar):
        ResidualBasedPredictorCorrectorBDFSchemeTurbulent<TSparseSpace,TDenseSpace>(DomainSize,rSlipVar)
    {
    }


    /** Destructor.
     */
    virtual ~TrilinosResidualBasedPredictorCorrectorBDFScheme()
    {
    }


    /*@} */
    /**@name Operators
     */
    /*@{ */

    virtual int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        int ErrorCode = BaseType::Check(rModelPart);
        if (ErrorCode != 0) return ErrorCode;

        // Check buffer size
        if (rModelPart.GetBufferSize() < 2)
            KRATOS_THROW_ERROR(std::logic_error, "BDF2 error: Insufficient buffer size for BDF2, should be at least 2, got ",rModelPart.GetBufferSize());

        // Check that all required variables were registered
        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if all applications were correctly registered.","");
        if(OSS_SWITCH.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"OSS_SWITCH Key is 0. Check if all applications were correctly registered.","");

        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT Key is 0. Check if all applications were correctly registered.","");
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if all applications were correctly registered.","");
        if(MESH_VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if all applications were correctly registered.","");

        // Checks on process info
//            const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

//            if(rCurrentProcessInfo.Has(DELTA_TIME) == 0)
//                KRATOS_THROW_ERROR(std::invalid_argument,"ProcessInfo does not contain a value for DELTA_TIME","");

        return 0;
        KRATOS_CATCH("");
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

} /* namespace Kratos.*/

#endif // KRATOS_TRILINOS_RESIDUALBASED_PREDICTORCORRECTOR_VELOCITY_BDF2_SCHEME_H
