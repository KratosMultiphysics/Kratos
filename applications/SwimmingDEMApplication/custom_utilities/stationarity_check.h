//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas
//

#if !defined( KRATOS_SWIMMING_DEM_FLOW_STATIONARITY_CHECK )
#define  KRATOS_SWIMMING_DEM_FLOW_STATIONARITY_CHECK

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/define.h"

namespace Kratos {
///@addtogroup StructuralMechanicsApplication
///@{
///@name Kratos Classes
///@{

/**
 * @class FlowStationarityCheck
 * @ingroup SwimmingDEMApplication
 * @brief This defines a class to assess whether stationarity has been reached in the fluid
 * @details It compares a spacial and historically-averaged pressure time derivative to the current spacially averaged value
 * @author Guillermo Casas
*/
class KRATOS_API(SWIMMING_DEM_APPLICATION) FlowStationarityCheck
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FlowStationarityCheck
    KRATOS_CLASS_POINTER_DEFINITION(FlowStationarityCheck);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FlowStationarityCheck(ModelPart& rFluidModelPart, const double tolerance)
        : mrModelPart(rFluidModelPart),
          mAveragingStep(0),
          mTolerance(tolerance),
          mCharacteristicPressureRate(1.0),
          mCurrentPressureRate(2.0){}

    /// Copy constructor.
    FlowStationarityCheck(FlowStationarityCheck const& rOther) = delete;

    /// Destructor.
    virtual ~FlowStationarityCheck() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    FlowStationarityCheck &operator=(FlowStationarityCheck const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function calculates the current time derivative average and compares it to the accumulated historical average
     *  @details This function can be called every few fluid time steps an
     *  @return true if the stationary state has been reached according to the chosen tolerance
     */
    bool AssessStationarity();
    double GetCharacteristicPressureDerivative();
    double GetCurrentPressureDerivative();
    double GetTolerance();
    double GetTransienceMeasure();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "FlowStationarityCheck";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FlowStationarityCheck";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

  protected:
    ///@name Protected Operations
    ///@{
    ///@}

  private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; /// The model part where to apply the constraints
    int mAveragingStep;
    double mTolerance;
    double mCharacteristicPressureRate;
    double mCurrentPressureRate;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}

}; // Class FlowStationarityCheck

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream& rIStream,
                                FlowStationarityCheck& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream& rOStream,
                                const FlowStationarityCheck& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}
///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_SWIMMING_DEM_FLOW_STATIONARITY_CHECK  defined
