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
#include "custom_utilities/nodal_data_handler.h"

/**
 * Note that this file makes use of the X-macro programming technique.
 * Suggested reading:
 * https://www.embedded.com/design/programming-languages-and-tools/4403953/C-language-coding-errors-with-X-macros-Part-1
 * https://en.wikipedia.org/wiki/X_Macro
 */ 

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@brief Base class for data containers used within FluidElement and derived types.
template< size_t TDim, size_t TNumNodes >
class FluidElementDataContainer
{
public:
    ///@name Type Definitions
    ///@{

    using NodalScalar = NodalDataHandler< double, TNumNodes, array_1d<double,TNumNodes> >;

    // Note: the first template argument is always array_1d<double,3>, as it refers to the variable type, not dimension
    using NodalVector = NodalDataHandler< array_1d<double,3>, TNumNodes, boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim> >;

    constexpr static unsigned int Dim = TDim;

    constexpr static unsigned int NumNodes = TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    FluidElementDataContainer() {
    }

    /// Destructor
    ~FluidElementDataContainer() {
    }

    /// (deleted) assignment operator.
    FluidElementDataContainer& operator=(FluidElementDataContainer const& rOther) = delete;

    /// (deleted) copy constructor.
    FluidElementDataContainer(FluidElementDataContainer const& rOther) = delete;

    ///@}
    ///@name Public Operations
    ///@{

    void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) {
    }

    int Check(const Element& rElement) {
        return 0;
    }

    ///@}
};

}

// Variable -> mVariableHandler
#define MEMBER_NAME(Variable) m ## Variable ## Handler

#define DECLARE_CLASS_MEMBER_FOR_HANDLER(Variable,HandlerType) \
private: HandlerType MEMBER_NAME(Variable);

#define DEFINE_GET_FUNCTION_FOR_HANDLER(Variable,HandlerType) \
public: HandlerType& Get##Variable() \
{ \
    return MEMBER_NAME(Variable); \
}

// Adds ", mVariableHandler(Variable)" (used for constructor initialization list)
#define CONSTRUCT_CLASS_MEMBER_FOR_HANDLER(Variable,HandlerType) \
, MEMBER_NAME(Variable) (Variable) 

#define INITIALIZE_HANDLER(Variable,HandlerType) \
MEMBER_NAME(Variable).Initialize(rElement,rProcessInfo);

#define CHECK_HANDLER(Variable,HandlerType) \
out = MEMBER_NAME(Variable).Check(rElement);

#define MAKE_FLUID_ELEMENT_DATA_CONTAINER(ClassName, Dimension, NumNodes, HANDLER_LIST) \
class ClassName : public FluidElementDataContainer<Dimension, NumNodes>                 \
{                                                                                       \
public:                                                                                 \
    ClassName()                                                                         \
    : FluidElementDataContainer()                                                       \
    HANDLER_LIST(CONSTRUCT_CLASS_MEMBER_FOR_HANDLER) {                                  \
    }                                                                                   \
                                                                                        \
    ~ClassName() {                                                                      \
    }                                                                                   \
                                                                                        \
    void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo) {         \
        HANDLER_LIST(INITIALIZE_HANDLER)                                                \
    }                                                                                   \
                                                                                        \
    int Check(const Element& rElement) {                                                \
        int out = 0;                                                                    \
        HANDLER_LIST(CHECK_HANDLER)                                                     \
        return out;                                                                     \
    }                                                                                   \
                                                                                        \
    HANDLER_LIST(DECLARE_CLASS_MEMBER_FOR_HANDLER)                                      \
                                                                                        \
    HANDLER_LIST(DEFINE_GET_FUNCTION_FOR_HANDLER)                                       \
};

#endif // KRATOS_INTEGRATION_POINT_DATA_CONTAINER_H_INCLUDED  defined
    