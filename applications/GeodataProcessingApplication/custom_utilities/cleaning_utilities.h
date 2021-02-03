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

#if !defined(KRATOS_CLEANING_UTILITIES_H_INCLUDED )
#define  KRATOS_CLEANING_UTILITIES_H_INCLUDED

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

  /// Auxiliary utility to maintain the quality of the model part

  class KRATOS_API(GEODATA_PROCESSING_APPLICATION) CleaningUtilities
  {

  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of CleaningUtilities
    KRATOS_CLASS_POINTER_DEFINITION(CleaningUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    CleaningUtilities( ModelPart& rModelPart ) : mrModelPart(rModelPart)
    { };

    /// Destructor.
    ~CleaningUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to erase all nodes that do not belong to an element
     *
     */
    void CleanIsolatedNodes();

    /**
     * @brief [NG] Function to erase all conditions with nodes no longer in model part
     *
     */
    void CleanConditions();

    /**
     * @brief [NG] Function to erase all conditions in BottomModelPart with nodes belong to SKIN_ISOSURFACE
     * nodes in the angles
     *
     */
    void CleanConditionsAngles();

    /**
     * @brief [NG] Function to fill the new conditions into bottom sub model part
     *
     */
    void FillBottom();

    /**
     * @brief Hard copy the content between model parts
     *
     */
    ModelPart& HardCopyBeforeSurfaceDiscretization( ModelPart& OriginalModelPart, ModelPart& NewModelPart );

    /**
     * @brief Hard copy the content between model parts
     *
     */
    ModelPart& HardCopyAfterSurfaceDiscretization( ModelPart& OriginalModelPart, ModelPart& NewModelPart );

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
    CleaningUtilities& operator=(CleaningUtilities const& rOther);

    /// Copy constructor.
    CleaningUtilities(CleaningUtilities const& rOther);

    ///@}

}; // Class CleaningUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CleaningUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CLEANING_UTILITIES_H_INCLUDED  defined
