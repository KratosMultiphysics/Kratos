//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Nicola Germano
//

#if !defined(KRATOS_FILL_MODELPART_UTILITIES_H_INCLUDED )
#define  KRATOS_FILL_MODELPART_UTILITIES_H_INCLUDED

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

  class KRATOS_API(GEODATA_PROCESSING_APPLICATION) FillModelpartUtilities
  {

  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of FillModelpartUtilities
    KRATOS_CLASS_POINTER_DEFINITION(FillModelpartUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    FillModelpartUtilities( ModelPart& rModelPart ) : mrModelPart(rModelPart)
    { };

    /// Destructor.
    ~FillModelpartUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to fill "Parts_Fluid" sub model part
     *
     */
    void FillFluid();

    /**
     * @brief Function to fill "Inlet" sub model part
     *
     */
    void FillInlet();

    /**
     * @brief Function to fill "Outlet" sub model part
     *
     */
    void FillSlip();

    /**
     * @brief Function to fill "Slip" sub model part
     *
     */
    void FillSlip();


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
    FillModelpartUtilities& operator=(FillModelpartUtilities const& rOther);

    /// Copy constructor.
    FillModelpartUtilities(FillModelpartUtilities const& rOther);

    ///@}

}; // Class FillModelpartUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const FillModelpartUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILL_MODELPART_UTILITIES_H_INCLUDED  defined
