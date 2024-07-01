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

#ifndef KRATOS_MANUFACTURED_BODY_FORCE_PROCESSH
#define KRATOS_MANUFACTURED_BODY_FORCE_PROCESSH

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

class ManufacturedBodyForceProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ManufacturedBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(ManufacturedBodyForceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ManufacturedBodyForceProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    ManufacturedBodyForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    ManufacturedBodyForceProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~ManufacturedBodyForceProcess() override {}

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
        buffer << "ManufacturedBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "ManufacturedBodyForceProcess";}

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
    bool                                             mImposeConditionsInlet;
    bool                                             mWriteExactSolutionOutput;
    bool                                             mInletTemperature;
    bool                                             mCavity;
    ModelPart*                                       mpInletModelPart = nullptr;
    ModelPart*                                       mpOutletModelPart = nullptr;
    double                                           mDiffusionCoefficient;
    double                                           mDensity;
    double                                           mSpecificHeat;
    double                                           mConductivity;
    double                                           mViscosity;
    double                                           mv0;
    double                                           mRad = 0.5e-3;
    double                                           mLen = 2e-3;
    bool                                             mAddHeatFluxSource;
    bool                                             mAddBodyForceSource;
    double                                           mTGrad;
    ModelPart*                                       mpNoSlipModelPart;
    ModelPart*                                       mpOutetModelPart;

    ///@}
    ///@name Protected Operators
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    const Parameters GetDefaultParameters() const override;

    double GetFakeOutletPressure();

    void SetInitialBodyForceAndMassSource();

    void SetBodyForceAndMassSource();

    void SetSystemProperties();

    ///@}
private:

    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    ManufacturedBodyForceProcess() = delete;

    /// Assignment operator.
    ManufacturedBodyForceProcess& operator=(ManufacturedBodyForceProcess const& rOther) = delete;

    /// Copy constructor.
    ManufacturedBodyForceProcess(ManufacturedBodyForceProcess const& rOther) = delete;

    ///@}

}; // Class ManufacturedBodyForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_MANUFACTURED_BODY_FORCE_PROCESSH