//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:      BSD License
//                Kratos default license: kratos/license.txt
//
//  Main authors: Miguel Angel Celigueta
//                Aditya Ghantasala
//

#ifndef POST_PROCESS_UTILITIES_H
#define POST_PROCESS_UTILITIES_H

#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"

namespace Kratos
{
  ///@addtogroup FluidDynamicsApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{
  ///
  /** Auxiliary utility for computing different values as post processing.
   */
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidPostProcessUtilities {
  public:

    ///@name Type Definitions
    ///@{

    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                 GeometryType;
    typedef IntegrationPoint<3>                       IntegrationPointType;
    typedef std::vector<IntegrationPointType>   IntegrationPointsArrayType;
    typedef GeometryData::IntegrationMethod          IntegrationMethodType;

    /// Pointer definition of DragUtilities
    KRATOS_CLASS_POINTER_DEFINITION(FluidPostProcessUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    FluidPostProcessUtilities() {};

    /// Destructor.
    virtual ~FluidPostProcessUtilities() = default;


    /// Assignment operator.
    FluidPostProcessUtilities& operator=(FluidPostProcessUtilities const& rOther) = delete;

    /// Copy constructor.
    FluidPostProcessUtilities(FluidPostProcessUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     *  @brief This method calculates the flow throught the given modelpart (surface/line) in the normal
     *          direction.
     *  @param[in] rModelPart The model part instance where statistics are recorded.
     *  @return flow through the modelpart.
     */
    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please use the one from the \'FluidAuxiliaryUtilities\'.") double CalculateFlow(const ModelPart& rModelPart);

    ///@}

}; // Class FluidPostProcessUtilities

} // namespace Kratos.

#endif // POST_PROCESS_UTILITIES_H
