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

#ifndef KRATOS_SPATIAL_DEPENDANT_POROSITY_SOLUTION_BODY_FORCE_PROCESS_H
#define KRATOS_SPATIAL_DEPENDANT_POROSITY_SOLUTION_BODY_FORCE_PROCESS_H

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

class KRATOS_API(SWIMMING_DEM_APPLICATION) SpatialDependantPorositySolutionBodyForceProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpatialDependantPorositySolutionBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(SpatialDependantPorositySolutionBodyForceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SpatialDependantPorositySolutionBodyForceProcess(
        ModelPart& rModelPart,
        const double Density,
        const double Viscosity,
        const double IndependentTerm,
        const double MaximumAlpha,
        const double Centerx1,
        const double Centerx2);

    /// Constructor with Kratos parameters.
    SpatialDependantPorositySolutionBodyForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    SpatialDependantPorositySolutionBodyForceProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~SpatialDependantPorositySolutionBodyForceProcess() override {}

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
        buffer << "SpatialDependantPorositySolutionBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "SpatialDependantPorositySolutionBodyForceProcess";}

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
    double                                      mIndependentTerm;
    double                                         mMaximumAlpha;
    double                                             mCenterx1;
    double                                             mCenterx2;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    const Parameters GetDefaultParameters() const override;

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
    SpatialDependantPorositySolutionBodyForceProcess() = delete;

    /// Assignment operator.
    SpatialDependantPorositySolutionBodyForceProcess& operator=(SpatialDependantPorositySolutionBodyForceProcess const& rOther) = delete;

    /// Copy constructor.
    SpatialDependantPorositySolutionBodyForceProcess(SpatialDependantPorositySolutionBodyForceProcess const& rOther) = delete;

    ///@}

}; // Class SpatialDependantPorositySolutionBodyForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_SPATIAL_DEPENDANT_POROSITY_SOLUTION_BODY_FORCE_PROCESS_H