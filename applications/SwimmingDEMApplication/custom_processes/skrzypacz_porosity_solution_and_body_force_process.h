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

#ifndef KRATOS_SKRZYPACZ_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H
#define KRATOS_SKRZYPACZ_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H

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

class KRATOS_API(SWIMMING_DEM_APPLICATION) SkrzypaczPorositySolutionAndBodyForceProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    typedef Node NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;
    /// Pointer definition of SkrzypaczPorositySolutionAndBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(SkrzypaczPorositySolutionAndBodyForceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    SkrzypaczPorositySolutionAndBodyForceProcess();
    /// Constructor.
    SkrzypaczPorositySolutionAndBodyForceProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    SkrzypaczPorositySolutionAndBodyForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    SkrzypaczPorositySolutionAndBodyForceProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~SkrzypaczPorositySolutionAndBodyForceProcess() override {}

    ///@}

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
        buffer << "SkrzypaczPorositySolutionAndBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "SkrzypaczPorositySolutionAndBodyForceProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ModelPart&                                       mrModelPart;
    double                                              mDensity;
    double                                            mViscosity;
    bool                                      mInitialConditions;
    bool                                 mAlternativeFormulation;


    ///@}
    ///@name Protected Operators
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    const Parameters GetDefaultParameters() const override;

    void SetInitialBodyForceAndPorosityField();

    void SetValuesOnIntegrationPoints();

    void SetFluidProperties();
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

}; // Class SkrzypaczPorositySolutionAndBodyForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_SKRZYPACZ_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H