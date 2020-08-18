//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta


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
    class KRATOS_API(FLUID_DYNAMICS_APPLICATION) PostProcessUtilities {
  public:

    ///@name Type Definitions
    ///@{

    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                 GeometryType;
    typedef IntegrationPoint<3>                       IntegrationPointType;
    typedef std::vector<IntegrationPointType>   IntegrationPointsArrayType;

    /// Pointer definition of DragUtilities
    KRATOS_CLASS_POINTER_DEFINITION(PostProcessUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    PostProcessUtilities() {};

    /// Destructor.
    virtual ~PostProcessUtilities() {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    double ComputeFlow(const ModelPart& rModelPart);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    PostProcessUtilities& operator=(PostProcessUtilities const& rOther);

    /// Copy constructor.
    PostProcessUtilities(PostProcessUtilities const& rOther);

    ///@}


}; // Class PostProcessUtilities

} // namespace Kratos.

#endif // POST_PROCESS_UTILITIES_H
