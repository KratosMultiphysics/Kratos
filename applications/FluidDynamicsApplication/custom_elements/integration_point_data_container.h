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

#define FLUID_ELEMENT_VARIABLES(MACRO_TO_APPLY) \
MACRO_TO_APPLY(VELOCITY,NodalVectorType) \
MACRO_TO_APPLY(PRESSURE,NodalScalarType)

#define MEMBER_NAME(Variable) m ## Variable ## Handler

#define DECLARE_CLASS_MEMBER_FOR_HANDLER(Variable,HandlerType) \
private: HandlerType MEMBER_NAME(Variable);

#define DEFINE_GET_FUNCTION_FOR_HANDLER(Variable,HandlerType) \
public: HandlerType& Get##Variable() \
{ \
    return MEMBER_NAME(Variable); \
}

#define CONSTRUCT_CLASS_MEMBER_FOR_HANDLER(Variable,HandlerType) \
MEMBER_NAME(Variable) (Variable), // mHandlerName(Variable),

#define INITIALIZE_HANDLER(Variable,HandlerType) \
MEMBER_NAME(Variable).Initialize(rElement,rProcessInfo);



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
    
    IntegrationPointDataContainer():
    FLUID_ELEMENT_VARIABLES(CONSTRUCT_CLASS_MEMBER_FOR_HANDLER)
    dummy(0)
    {

    }

    /// (deleted) assignment operator.
    IntegrationPointDataContainer& operator=(IntegrationPointDataContainer const& rOther) = delete;

    /// (deleted) copy constructor.
    IntegrationPointDataContainer(IntegrationPointDataContainer const& rOther) = delete;

    ///@}
    ///@name Public members
    ///@{


    ///@}
    ///@name Public Operations
    ///@{

    void Initialize(Element& rElement, const ProcessInfo& rProcessInfo) {
        FLUID_ELEMENT_VARIABLES(INITIALIZE_HANDLER)
    }

    ///@}

    FLUID_ELEMENT_VARIABLES(DECLARE_CLASS_MEMBER_FOR_HANDLER)

    FLUID_ELEMENT_VARIABLES(DEFINE_GET_FUNCTION_FOR_HANDLER)
    
private:

    const int dummy;

}; // struct IntegrationPointDataContainer

///@}

///@} addtogroup block

} // namespace Kratos.

#undef FLUID_ELEMENT_VARIABLES

#undef DECLARE_CLASS_MEMBER_FOR_HANDLER
#undef DEFINE_GET_FUNCTION_FOR_HANDLER
#undef CONSTRUCT_CLASS_MEMBER_FOR_HANDLER
#undef INITIALIZE_HANDLER

#endif // KRATOS_INTEGRATION_POINT_DATA_CONTAINER_H_INCLUDED  defined
