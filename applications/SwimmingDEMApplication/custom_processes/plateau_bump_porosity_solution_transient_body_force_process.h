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

#ifndef KRATOS_PLATEAU_BUMP_POROSITY_SOLUTION_TRANSIENT_BODY_FORCE_PROCESS_H
#define KRATOS_PLATEAU_BUMP_POROSITY_SOLUTION_TRANSIENT_BODY_FORCE_PROCESS_H

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

class KRATOS_API(SWIMMING_DEM_APPLICATION) PlateauBumpPorositySolutionTransientBodyForceProcess : public Process
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

    /// Pointer definition of PlateauBumpPorositySolutionTransientBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(PlateauBumpPorositySolutionTransientBodyForceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    PlateauBumpPorositySolutionTransientBodyForceProcess();
    /// Constructor.
    PlateauBumpPorositySolutionTransientBodyForceProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    PlateauBumpPorositySolutionTransientBodyForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    PlateauBumpPorositySolutionTransientBodyForceProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~PlateauBumpPorositySolutionTransientBodyForceProcess() override {}

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
        buffer << "PlateauBumpPorositySolutionTransientBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "PlateauBumpPorositySolutionTransientBodyForceProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ModelPart&                                       mrModelPart;
    double                                              mDensity;
    double                                            mViscosity;
    double                                             mAlphaMin;
    double                                             mAlphaMax;
    double                                                mSigma;
    double                                                mUchar;
    double                                             mX1Origin;
    double                                             mX2Origin;
    double                                           mBumpRadius;
    double                                        mPlateauRadius;
    bool                                      mInitialConditions;
    bool                                 mAlternativeFormulation;

    ///@}
    ///@name Protected Operators
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    const Parameters GetDefaultParameters() const override;

    void SetBodyForceAndPorosityField();

    void SetValuesOnIntegrationPoints();

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

}; // Class PlateauBumpPorositySolutionTransientBodyForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_PLATEAU_BUMP_POROSITY_SOLUTION_TRANSIENT_BODY_FORCE_PROCESS_H