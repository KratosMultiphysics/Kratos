//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//
//

#ifndef KRATOS_BUMP_TRANSIENT_POROSITY_SOLUTION_BODY_FORCE_PROCESS_H
#define KRATOS_BUMP_TRANSIENT_POROSITY_SOLUTION_BODY_FORCE_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

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

class KRATOS_API(SWIMMING_DEM_APPLICATION) BumpTransientPorositySolutionBodyForceProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BumpTransientPorositySolutionBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(BumpTransientPorositySolutionBodyForceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    BumpTransientPorositySolutionBodyForceProcess();
    /// Constructor.
    BumpTransientPorositySolutionBodyForceProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    BumpTransientPorositySolutionBodyForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    BumpTransientPorositySolutionBodyForceProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~BumpTransientPorositySolutionBodyForceProcess() override {}

    ///@}

    ModelPart&                                       mrModelPart;
    double                                              mDensity;
    double                                            mViscosity;
    double                                         mPermeability;
    double                                                mUchar;
    double                                           mDeltaAlpha;
    double                                               mLength;
    double                                   mMaxSqueezeFraction;
    double                                                mOmega;
    double                                     mSqueezeAmplitude;
    double                                              mNSafety;
    double                                             mX1Origin;
    double                                             mX2Origin;
    double                                       mReynoldsNumber;
    double                                      mDamKohlerNumber;
    bool                                      mInitialConditions;
    bool                                 mAlternativeFormulation;
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

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    const Parameters GetDefaultParameters() const override;

    void CalculateKinematicViscosity();

    void CalculatePermeability(double &dynamic_viscosity);

    void SetInitialBodyForceAndPorosityField();

    void SetBodyForceAndPorosityField();

    void SetFluidProperties();

    bool IsInsideEllipticalSupport(
        const double x1,
        const double x2,
        const double c,
        const double R);

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
        buffer << "BumpBumpTransientPorositySolutionBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "BumpTransientPorositySolutionBodyForceProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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

    ///@name Protected Operators
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

    /// Default constructor.

    /// Assignment operator.

    /// Copy constructor.
    ///@}

}; // Class BumpBumpTransientPorositySolutionBodyForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_BUMP_TRANSIENT_POROSITY_SOLUTION_BODY_FORCE_PROCESS_H