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

#ifndef KRATOS_CASAS_SOLUTION_BODY_FORCE_PROCESS_H
#define KRATOS_CASAS_SOLUTION_BODY_FORCE_PROCESS_H

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

class KRATOS_API(SWIMMING_DEM_APPLICATION) CasasSolutionBodyForceProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3>> > ComponentType;

    /// Pointer definition of CasasSolutionBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(CasasSolutionBodyForceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    CasasSolutionBodyForceProcess(
        ModelPart& rModelPart,
        const double Density,
        const double Viscosity,
        const double DeltaAlpha,
        const double Length,
        const double SqueezeAmplitude,
        const double X1Origin,
        const double X2Origin);

    /// Constructor with Kratos parameters.
    CasasSolutionBodyForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    CasasSolutionBodyForceProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~CasasSolutionBodyForceProcess() override {}

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
        buffer << "CasasSolutionBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "CasasSolutionBodyForceProcess";}

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

    ModelPart&                                       mrModelPart;
    double                                              mDensity;
    double                                            mViscosity;
    double                                           mDeltaAlpha;
    double                                               mLength;
    double                                   mMaxSqueezeFraction;
    double                                                mOmega;
    double                                     mSqueezeAmplitude;
    double                                             mX1Origin;
    double                                             mX2Origin;
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    void SetInitialBodyForceAndPorosityField();

    void SetBodyForceAndPorosityField();

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
    CasasSolutionBodyForceProcess() = delete;

    /// Assignment operator.
    CasasSolutionBodyForceProcess& operator=(CasasSolutionBodyForceProcess const& rOther) = delete;

    /// Copy constructor.
    CasasSolutionBodyForceProcess(CasasSolutionBodyForceProcess const& rOther) = delete;

    ///@}

}; // Class CasasSolutionBodyForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_CASAS_SOLUTION_BODY_FORCE_PROCESS_H