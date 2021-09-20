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

#ifndef KRATOS_SPATIAL_DEPENDANT_HYPERBOLIC_TANGENTIAL_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H
#define KRATOS_SPATIAL_DEPENDANT_HYPERBOLIC_TANGENTIAL_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H

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

class KRATOS_API(SWIMMING_DEM_APPLICATION) SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess();
    /// Constructor.
    SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess() override {}

    ///@}

    ModelPart&                                       mrModelPart;
    double                                              mDensity;
    double                                            mViscosity;
    double                                         mPermeability;
    double                                                mUchar;
    double                                               mLength;
    double                                            mMeanAlpha;
    double                                             mMinAlpha;
    double                                               mHeight;
    double                                       mReynoldsNumber;
    double                                      mDamKohlerNumber;
    double                                         mMaxGradAlpha;
    double                                       mFirstParameter;
    double                                      mSecondParameter;
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

    void CalculateKinematicViscosity(
        double &rReynoldsNumber,
        double &viscosity);

    void CalculatePermeability(
        double &rDamKohlerNumber,
        double &dynamic_viscosity,
        double &permeability);

    void CalculateFunctionParameters(
        double &mFirstParameter,
        double &mSecondParameter);

    void SetInitialBodyForceAndPorosityField();

    void SetFluidProperties();

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
        buffer << "SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess";}

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

}; // Class SpatialDependantHyperbolicTangentialPorositySolutionAndBodyForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_SPATIAL_DEPENDANT_HYPERBOLIC_TANGENTIAL_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H