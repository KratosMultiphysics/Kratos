//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_NEGATIVE_HEIGHT_WETTING_MODEL
#define KRATOS_NEGATIVE_HEIGHT_WETTING_MODEL


// System includes


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
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

/** 
 * @ingroup ShallowWaterApplication
 * @class NegativeHeightWettingModel
 * @brief A simple wetting model. Heniche et al. "A two-dimensional ®nite element drying-wetting shallow water model for rivers and estuaries"
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) NegativeHeightWettingModel : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(NegativeHeightWettingModel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor with model and parameters
    NegativeHeightWettingModel(ModelPart& rModelPart, Parameters ThisParameters);

    /// @brief Constructor with model and doubles
    NegativeHeightWettingModel(ModelPart& rModelPart, double Beta);

    /// Destructor.
    virtual ~NegativeHeightWettingModel() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// @brief The scheme calls this method at the begining of the non-linear loop
    void ExecuteInitializeSolutionStep() override;

    ///@brief The scheme calls this method after convergence is reached in an iterative process
    void ExecuteFinalizeSolutionStep() override;

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
    std::string Info() const override
    {
        return "NegativeHeightWettingModel";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ModelPart& mrModelPart;
    double mBeta;
    double mDryHeight;

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
    NegativeHeightWettingModel& operator=(NegativeHeightWettingModel const& rOther) = delete;

    /// Copy constructor.
    NegativeHeightWettingModel(NegativeHeightWettingModel const& rOther) = delete;


    ///@}

}; // Class NegativeHeightWettingModel

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                NegativeHeightWettingModel& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const NegativeHeightWettingModel& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEGATIVE_HEIGHT_WETTING_MODEL  defined
