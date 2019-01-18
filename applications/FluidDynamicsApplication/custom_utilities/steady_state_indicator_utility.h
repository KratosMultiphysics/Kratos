#if !defined(KRATOS_STEADY_STATE_INDICATOR_UTILITY_H)
#define KRATOS_STEADY_STATE_INDICATOR_UTILITY_H

// External includes 

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/find_nodal_h_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// A set of functions to indicate if steady state has been reached by anaylising
/// the change of quantities integrated within the domain
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) SteadyStateIndicatorUtility {

public:
    /// Pointer definition of SteadyStateIndicatorUtility
    KRATOS_CLASS_POINTER_DEFINITION(SteadyStateIndicatorUtility);

    /// Default constructor
    SteadyStateIndicatorUtility(ModelPart& rModelPart);

    /// Destructor.
    ~SteadyStateIndicatorUtility(){};

    ///@}
    ///@name Operators
    ///@{

    /// Calculating change of solution values in time over the domain
    void EstimateQuantityChangesInTime();

    /// Return memember variables
    double GetVelocityChange(){return mChangeInVelocity;}
    double GetPressureChange(){return mChangeInPressure;}

private:

    /// members
    ModelPart&                               mrModelPart;
    double                                   mChangeInVelocity;
    double                                   mChangeInPressure;
    bool                                     mIsSteady;

    double NodalVelocityChange(ModelPart::NodeType& inode);
    double NodalPressureChange(ModelPart::NodeType& inode);

};  // Class SteadyStateIndicatorUtility

}  // namespace Kratos.

#endif // KRATOS_STEADY_STATE_INDICATOR_UTILITY_H defined 


