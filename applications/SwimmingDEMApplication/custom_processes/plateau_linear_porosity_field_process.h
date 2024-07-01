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

#ifndef KRATOS_PLATEAU_LINEAR_POROSITY_FIELD_PROCESS_H
#define KRATOS_PLATEAU_LINEAR_POROSITY_FIELD_PROCESS_H

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

class KRATOS_API(SWIMMING_DEM_APPLICATION) PlateauLinearPorosityFieldProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Node type (default is: Node<3>)
    typedef Node NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Pointer definition of PlateauLinearPorosityFieldProcess
    KRATOS_CLASS_POINTER_DEFINITION(PlateauLinearPorosityFieldProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PlateauLinearPorosityFieldProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    PlateauLinearPorosityFieldProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    PlateauLinearPorosityFieldProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~PlateauLinearPorosityFieldProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    void SetValuesOnIntegrationPoints();

    void SetPorosityField();

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
        buffer << "PlateauLinearPorosityFieldProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "PlateauLinearPorosityFieldProcess";}

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

    ModelPart&                                        mrModelPart;
    double                                              mAlphaMin;
    double                                              mAlphaMax;
    double                                        mBeginningSlope;
    double                                      mBeginningPlateau;
    double                                         mEndingPlateau;
    double                                           mEndingSlope;
    double                                         mInletVelocity;


    ///@}
    ///@name Protected Operators
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    const Parameters GetDefaultParameters() const override;


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
    PlateauLinearPorosityFieldProcess() = delete;

    /// Assignment operator.
    PlateauLinearPorosityFieldProcess& operator=(PlateauLinearPorosityFieldProcess const& rOther) = delete;

    /// Copy constructor.
    PlateauLinearPorosityFieldProcess(PlateauLinearPorosityFieldProcess const& rOther) = delete;

    ///@}

}; // Class PlateauLinearPorosityFieldProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_PLATEAU_LINEAR_POROSITY_FIELD_PROCESS_H