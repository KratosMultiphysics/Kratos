//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Simon Wenczowski
//

#if !defined(KRATOS_EXTRUSION_HEIGHT_UTILITIES_H_INCLUDED )
#define  KRATOS_EXTRUSION_HEIGHT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geodata_processing_application_variables.h"
#include "includes/checks.h"
#include "includes/model_part.h"


namespace Kratos
{
  ///@addtogroup GeodataProcessingApplication
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

  /// Auxiliary utility to compute the vertical extrusion height of the terrain.

  class KRATOS_API(GEODATA_PROCESSING_APPLICATION) ExtrusionHeightUtilities
  {

  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ExtrusionHeightUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ExtrusionHeightUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ExtrusionHeightUtilities( ModelPart& rModelPart ) : mrModelPart(rModelPart)
    {
        KRATOS_CHECK_VARIABLE_KEY( EXTRUSION_HEIGHT );
    };

    /// Destructor.
    ~ExtrusionHeightUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief The point cloud inside the model part is analysed and the each point is asigned an extrusion height.
     * Points can be deleted from the terrain point cloud
     *
     * @param radius Radius for the search of a minimal height that serves as a reference
     * @param max_etrusion_height The etrusion height over the found point with the minimal terrain height
     * @param free_board The minimal vertical distance between a terrain point and the top of the domain
     */
    void SetExtrusionHeight( const double radius, const double max_etrusion_height, const double free_board );

    /**
     * @brief The field of the variable extrusion height is smoothed to avoid regional refinement at the topper
     *
     * @param radius Radius to define a neighborhood of points in an unstructured point cloud
     * @param iterations Iterations in the smoothing mechanism
     * @param free_board Minimal distance between topper and the highest terrain points
     */
    void SmoothExtrusionHeight( const double radius, const int iterations, const double free_board );

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

    ModelPart& mrModelPart;

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
    ExtrusionHeightUtilities& operator=(ExtrusionHeightUtilities const& rOther);

    /// Copy constructor.
    ExtrusionHeightUtilities(ExtrusionHeightUtilities const& rOther);

    ///@}

}; // Class ExtrusionHeightUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ExtrusionHeightUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_EXTRUSION_HEIGHT_UTILITIES_H_INCLUDED  defined
