//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_DRAG_UTILITIES_H_INCLUDED )
#define  KRATOS_DRAG_UTILITIES_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
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

  /// Auxiliary utility to compute the drag force.
  /** For embedded formulations, this utility iterates all the elements of a provided model part. In this iteration
   * calls the calculate method of each element to compute the value of the variable DRAG_FORCE. If the element is split,
   * this method computes the integration of the stress term over the interface. Otherwise, the value is just zero.
   * The obtained values are accumulated to get the total drag force in the model part.
   *
   * Note that if there is more than one embedded object, one just needs to save the surrounding elements to each embedded
   * object in different submodelparts and call this process for each one of that submodelparts.
   *
   * For the body fitted slip case, it integrates the pressure stress term over the given submodelpart conditions (the
   * shear stress term is assumed to be zero).
   */
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) DragUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef Geometry<Node<3>>                                 GeometryType;
    typedef IntegrationPoint<3>                       IntegrationPointType;
    typedef std::vector<IntegrationPointType>   IntegrationPointsArrayType;

    /// Pointer definition of DragUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DragUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    DragUtilities() {};

    /// Destructor.
    ~DragUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * Computes the integral of the pressure stress term normal projection over the conditions
    * of the given modelpart
    * @param rModelPart reference to the model part in where the drag is to be computed
    * @return An array containing the drag force value.
    */
    array_1d<double, 3> CalculateBodyFittedDrag(ModelPart &rModelPart);

    /**
    * Computes the integral of the Cauchy stress term normal projection in the given modelpart elements.
    * @param rModelPart reference to the model part in where the drag is to be computed
    * @return An array containing the drag force value.
    */
    array_1d<double, 3> CalculateEmbeddedDrag(ModelPart &rModelPart);

    /**
    * Calculates the drag force location in embedded formulations
    * @param rModelPart reference to the model part in where the drag force location is to be computed
    * @return An array containing the drag force location coordinates.
    */
    array_1d<double, 3> CalculateEmbeddedDragCenter(const ModelPart &rModelPart);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

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
    DragUtilities& operator=(DragUtilities const& rOther);

    /// Copy constructor.
    DragUtilities(DragUtilities const& rOther);

    ///@}

}; // Class DragUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DragUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DRAG_UTILITIES_H_INCLUDED  defined
