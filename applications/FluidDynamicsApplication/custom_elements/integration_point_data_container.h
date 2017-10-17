//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_INTEGRATION_POINT_DATA_CONTAINER_H_INCLUDED)
#define KRATOS_INTEGRATION_POINT_DATA_CONTAINER_H_INCLUDED

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"

#include "fluid_dynamics_application_variables.h"
#include "fluid_element_data.h"

#define FLUID_ELEMENT_VARIABLES \
APPLY_TO_VARIABLE(mVelocityHandler,VELOCITY,NodalVectorType) \
APPLY_TO_VARIABLE(mPressureHandler,PRESSURE,NodalScalarType)

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

/**
 * @brief Auxiliary class to evaluate and hold data at integration points for elements
 * based on FluidElement.
 * Note that this class makes heavy use of the X-macro programming technique. 
 * (see for example https://en.wikipedia.org/wiki/X_Macro for reference).
 * The 'X macro' is called APPLY_TO_VARIABLE here.
 */
class IntegrationPointDataContainer
{
public:
    ///@name Type Definitions
    ///@{

    typedef NodalDataHandler< array_1d<double,3>, 3, boost::numeric::ublas::bounded_matrix<double,3,2> > NodalVectorType;

    typedef NodalDataHandler< double, 3, array_1d<double,3> > NodalScalarType;

    ///@}
    ///@name Life Cycle
    ///@{


    #define APPLY_TO_VARIABLE(Name,Variable,Handler) Name(Variable),
    
    IntegrationPointDataContainer():
    FLUID_ELEMENT_VARIABLES
    dummy(0)
    {

    }

    #undef APPLY_TO_VARIABLE

    ///@}
    ///@name Public members
    ///@{


    ///@}
    ///@name Public Operations
    ///@{

    void Initialize(Element& rElement, const ProcessInfo& rProcessInfo) {
        #define APPLY_TO_VARIABLE(Name,Variable,Handler) Name.Initialize(rElement,rProcessInfo);
        FLUID_ELEMENT_VARIABLES
        #undef APPLY_TO_VARIABLE
    }

    ///@}

    #define APPLY_TO_VARIABLE(Name,Variable,Handler) private: Handler Name;
    FLUID_ELEMENT_VARIABLES
    #undef APPLY_TO_VARIABLE

    #define APPLY_TO_VARIABLE(Name,Variable,Handler) public: \
    Handler& Get##Variable() \
    { \
        return Name; \
    }

    FLUID_ELEMENT_VARIABLES
    
    #undef APPLY_TO_VARIABLE

private:

    const int dummy;

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    IntegrationPointDataContainer& operator=(IntegrationPointDataContainer const& rOther);

    /// Copy constructor.
    IntegrationPointDataContainer(IntegrationPointDataContainer const& rOther);

    ///@}

}; // struct IntegrationPointDataContainer

///@}

///@} addtogroup block

} // namespace Kratos.

#undef FLUID_ELEMENT_VARIABLES

#endif // KRATOS_INTEGRATION_POINT_DATA_CONTAINER_H_INCLUDED  defined
