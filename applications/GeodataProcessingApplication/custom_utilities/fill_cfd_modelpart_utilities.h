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

#if !defined(KRATOS_FILL_CFD_MODELPART_UTILITIES_H_INCLUDED)
#define  KRATOS_FILL_CFD_MODELPART_UTILITIES_H_INCLUDED

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

  class KRATOS_API(GEODATA_PROCESSING_APPLICATION) FillCfdModelpartUtilities
  {

  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of FillCfdModelpartUtilities
    KRATOS_CLASS_POINTER_DEFINITION(FillCfdModelpartUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    FillCfdModelpartUtilities( ModelPart& rModelPart ) : mrModelPart(rModelPart)
    { };

    /// Destructor.
    ~FillCfdModelpartUtilities() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to fill cfd model part
     * @param NewModelPart The cfd model part
     */
    void FillModelPart(ModelPart& NewModelPart);

    /**
     * @brief Function to fill "Parts_Fluid" sub model part
     *
     */
    void FillPartsFluid(ModelPart& OriginModelPart, const std::string& ElementModelName);
    // ModelPart& FillPartsFluid (ModelPart& OriginModelPart, const std::string& ElementModelName);

    /**
     * @brief Function to fill "Inlet" sub model part
     *
     */
    void FillInlet(ModelPart& OriginModelPart, const std::string& ConditionModelName);
    // ModelPart& FillInlet (ModelPart& OriginModelPart, const std::string& ConditionModelName);

    /**
     * @brief Function to fill "Outlet" sub model part
     *
     */
    void FillOutlet(ModelPart& OriginModelPart, const std::string& ConditionModelName);
    // ModelPart& FillOutlet (ModelPart& OriginModelPart, const std::string& ConditionModelName);

    /**
     * @brief Function to fill "Slip" sub model part
     *
     */
    void FillSlip(const std::string& ConditionModelName);

    /**
     * @brief Function to fill "NoSlip" sub model part
     *
     */
    void FillNoslip(const std::string& ConditionModelName);


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
    FillCfdModelpartUtilities& operator=(FillCfdModelpartUtilities const& rOther);

    /// Copy constructor.
    FillCfdModelpartUtilities(FillCfdModelpartUtilities const& rOther);

    ///@}

}; // Class FillCfdModelpartUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const FillCfdModelpartUtilities& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CLEANING_UTILITIES_H_INCLUDED  defined
