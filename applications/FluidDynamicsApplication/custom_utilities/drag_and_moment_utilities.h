//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Anoop Kodakkal
//

#if !defined(KRATOS_DRAG_AND_MOMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_DRAG_AND_MOMENT_UTILITIES_H_INCLUDED

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

  /// Auxiliary utility to compute the drag force and the base moment.
  /**
   * Add description here.
   */
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) DragAndMomentUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef Geometry<Node<3>>                                 GeometryType;
    typedef IntegrationPoint<3>                       IntegrationPointType;
    typedef std::vector<IntegrationPointType>   IntegrationPointsArrayType;

    /// Pointer definition of DragAndMomentUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DragAndMomentUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    DragAndMomentUtilities() {};

    /// Destructor.
    ~DragAndMomentUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * Computes the integral of the pressure stress term normal projection over the conditions
    * of the given modelpart and the base moment.
    * @param rModelPart reference to the model part in where the drag is to be computed
    * @return An array containing the drag force value.
    */
    array_1d<double, 6> CalculateBodyFittedDragAndMoment(ModelPart &rModelPart, array_1d<double, 3> rReferencePoint);

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
    DragAndMomentUtilities& operator=(DragAndMomentUtilities const& rOther);

    /// Copy constructor.
    DragAndMomentUtilities(DragAndMomentUtilities const& rOther);

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
    const DragAndMomentUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DRAG_UTILITIES_H_INCLUDED  defined
