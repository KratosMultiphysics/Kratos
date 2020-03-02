//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_LOGARITHMIC_Y_PLUS_VELOCITY_SENSITIVITIES_PROCESS_H_INCLUDED)
#define KRATOS_RANS_LOGARITHMIC_Y_PLUS_VELOCITY_SENSITIVITIES_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "processes/process.h"

// Application includes


namespace Kratos
{
///@addtogroup RANSApplication
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
 * @brief Calculates $y^+$ sensitivities w.r.t. velocity
 *
 * This process calculates sensitivities of $y^+$ w.r.t. velocity, In here
 * The $y^+$ should be calculated based on logarithmic law.
 *
 * @see RansLogarithmicYPlusCalculationProcess
 */

class KRATOS_API(RANS_APPLICATION) RansLogarithmicYPlusVelocitySensitivitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::NodeType NodeType;

    typedef ModelPart::ElementType ElementType;

    typedef Geometry<NodeType> GeometryType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Pointer definition of RansLogarithmicYPlusVelocitySensitivitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansLogarithmicYPlusVelocitySensitivitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansLogarithmicYPlusVelocitySensitivitiesProcess(Model& rModel, Parameters& rParameters);
    /// Destructor.
    ~RansLogarithmicYPlusVelocitySensitivitiesProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

    void Execute() override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Model& mrModel;
    Parameters& mrParameters;
    std::string mModelPartName;

    int mEchoLevel;

    double mVonKarman;
    double mBeta;

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
    RansLogarithmicYPlusVelocitySensitivitiesProcess& operator=(
        RansLogarithmicYPlusVelocitySensitivitiesProcess const& rOther);

    /// Copy constructor.
    RansLogarithmicYPlusVelocitySensitivitiesProcess(
        RansLogarithmicYPlusVelocitySensitivitiesProcess const& rOther);

    ///@}

}; // namespace Kratos

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansLogarithmicYPlusVelocitySensitivitiesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_LOGARITHMIC_Y_PLUS_VELOCITY_SENSITIVITIES_PROCESS_H_INCLUDED defined
