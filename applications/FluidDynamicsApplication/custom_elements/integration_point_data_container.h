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

#define FLUID_ELEMENT_INTEGRATION_POINT_VARIABLES \
X(mVelocityHandler,VELOCITY,NodalVectorType) \
X(mPressureHandler,PRESSURE,NodalScalarType)

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

/// Auxiliary class to evaluate and hold data at integration points for elements
/// based on FluidElement
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

#define X(Name,Variable,Handler) Name(Variable),

    IntegrationPointDataContainer():
    FLUID_ELEMENT_INTEGRATION_POINT_VARIABLES
    dummy(0)
    {

    }

#undef X

    ///@}
    ///@name Public members
    ///@{


    ///@}
    ///@name Public Operations
    ///@{

    template< unsigned int N > void ConstructHandler();

    void Initialize(Element& rElement, const ProcessInfo& rProcessInfo) {
        #define X(Name,Variable,Handler) Name.Initialize(rElement,rProcessInfo);
        FLUID_ELEMENT_INTEGRATION_POINT_VARIABLES
        #undef X
    }

    ///@}

    #define X(Name,Variable,Handler) private: Handler Name;
    FLUID_ELEMENT_INTEGRATION_POINT_VARIABLES
    #undef X

public: 
    NodalDataHandler< array_1d<double,3>, 3, boost::numeric::ublas::bounded_matrix<double,3,2> >& GetVelocity()
    {
        return mVelocityHandler;
    }

    NodalDataHandler< double, 3, array_1d<double,3> >& GetPressure()
    {
        return mPressureHandler;
    }

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

#endif // KRATOS_INTEGRATION_POINT_DATA_CONTAINER_H_INCLUDED  defined
