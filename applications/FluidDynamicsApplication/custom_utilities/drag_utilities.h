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

  /// Auxiliary utility to compute the drag force in embedded formulations.
  /** This utility iterates all the elements of a provided model part. In this iteration calls the calculate method of
   * each element to compute the value of the variable DRAG_FORCE. If the element is split, this method computes the 
   * integration of the stress term over the interface. Otherwise, the value is just zero. The obtained values are 
   * accumulated to get the total drag force in the model part.
   * 
   * Note that if there is more than one embedded object, one just needs to save the surrounding elements to each embedded
   * object in different submodelparts and call this process for each one of that submodelparts.
   */
  class DragUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

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

    array_1d<double, 3> CalculateEmbeddedDrag(ModelPart &rModelPart);

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
