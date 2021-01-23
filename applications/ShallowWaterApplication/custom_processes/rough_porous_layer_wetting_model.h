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

#ifndef KRATOS_ROUGH_POROUS_LAYER_WETTING_MODEL
#define KRATOS_ROUGH_POROUS_LAYER_WETTING_MODEL


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
 * @brief A simple wetting model. Barros, Rosman, Telles, Azevedo "A simple wetting and drying method for shallow water flow with application in the Vitoria Bay estuary, Brazil"
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) RoughPorousLayerWettingModel : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RoughPorousLayerWettingModel
    KRATOS_CLASS_POINTER_DEFINITION(RoughPorousLayerWettingModel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor with model and parameters
    RoughPorousLayerWettingModel(ModelPart& rModelPart, Parameters ThisParameters);

    /// @brief Constructor with model and doubles
    RoughPorousLayerWettingModel(ModelPart& rModelPart, double LayerThickness, double RoughnessFactor);

    /// Destructor.
    virtual ~RoughPorousLayerWettingModel() = default;

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
        return "RoughPorousLayerWettingModel";
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
    double mLayerThickness;
    double mRoughnessFactor;

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
    RoughPorousLayerWettingModel& operator=(RoughPorousLayerWettingModel const& rOther) = delete;

    /// Copy constructor.
    RoughPorousLayerWettingModel(RoughPorousLayerWettingModel const& rOther) = delete;


    ///@}

}; // Class RoughPorousLayerWettingModel

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                RoughPorousLayerWettingModel& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const RoughPorousLayerWettingModel& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ROUGH_POROUS_LAYER_WETTING_MODEL  defined
