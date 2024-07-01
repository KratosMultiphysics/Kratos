//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala Pascual
//
//

#ifndef KRATOS_MICROFLUIDIC_TUBE_PROCESSH
#define KRATOS_MICROFLUIDIC_TUBE_PROCESSH

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

#include "includes/node.h"
#include "utilities/math_utils.h"
#include "spatial_containers/spatial_containers.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
///@addtogroup SwimmingDEMApplication
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

class MicrofluidicTubeProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MicrofluidicTubeProcess
    KRATOS_CLASS_POINTER_DEFINITION(MicrofluidicTubeProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MicrofluidicTubeProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    MicrofluidicTubeProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    MicrofluidicTubeProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~MicrofluidicTubeProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Operations
    ///@{

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
        std::stringstream buffer;
        buffer << "MicrofluidicTubeProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "MicrofluidicTubeProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{

    ///@}

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart&                                       mrModelPart;
    ModelPart*                                       mpOutletModelPart = nullptr;
    std::string                                      mGeometryType;
    Vector                                           mDimensions;
    double                                           mPressureOutlet;
    double                                           mDensity;
    double                                           mSpecificHeat;
    double                                           mConductivity;
    double                                           mViscosity;
    double                                           mDiffusionCoefficient;

    ///@}
    ///@name Protected Operators
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    const Parameters GetDefaultParameters() const override;

    double GetVelocityValue(Node<3> node);

    void SetOutletVelocityAndPressure();

    void SetSystemProperties();

    ///@}
private:

    ///@name Private  Access
    ///@{

    Vector GetSurfaceCenter();

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    MicrofluidicTubeProcess() = delete;

    /// Assignment operator.
    MicrofluidicTubeProcess& operator=(MicrofluidicTubeProcess const& rOther) = delete;

    /// Copy constructor.
    MicrofluidicTubeProcess(MicrofluidicTubeProcess const& rOther) = delete;

    ///@}

}; // Class MicrofluidicTubeProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_MICROFLUIDIC_TUBE_PROCESSH