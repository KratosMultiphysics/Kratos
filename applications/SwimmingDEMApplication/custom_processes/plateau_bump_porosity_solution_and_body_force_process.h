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

//#pragma once
#ifndef KRATOS_PLATEAU_BUMP_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H
#define KRATOS_PLATEAU_BUMP_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H
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
    //extern DenseVector<std::vector<double>> mExactScalar;
    // extern DenseVector<std::vector<double>> mExactPorosity;
    // extern DenseVector<std::vector<double>> mExactPorosityRate;
    extern DenseVector<Matrix> mExactBodyForce;
    // extern DenseVector<Matrix> mExactPorosityGradient;
    // extern DenseVector<Matrix> mExactVector;
    //extern DenseVector<Matrix> mExactScalarGradient;
    //extern DenseVector<DenseVector<Matrix>> mExactVectorGradient;
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

class KRATOS_API(SWIMMING_DEM_APPLICATION) PlateauBumpPorositySolutionAndBodyForceProcess : public Process
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

    /// Pointer definition of PlateauBumpPorositySolutionAndBodyForceProcess
    KRATOS_CLASS_POINTER_DEFINITION(PlateauBumpPorositySolutionAndBodyForceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    PlateauBumpPorositySolutionAndBodyForceProcess();
    /// Constructor.
    PlateauBumpPorositySolutionAndBodyForceProcess(
        ModelPart& rModelPart);

    /// Constructor with Kratos parameters.
    PlateauBumpPorositySolutionAndBodyForceProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    PlateauBumpPorositySolutionAndBodyForceProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~PlateauBumpPorositySolutionAndBodyForceProcess() override {}

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
        buffer << "PlateauBumpPorositySolutionAndBodyForceProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "PlateauBumpPorositySolutionAndBodyForceProcess";}

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

    void SetInitialBodyForceAndPorosityField();

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

}; // Class PlateauBumpPorositySolutionAndBodyForceProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.
#endif // KRATOS_PLATEAU_BUMP_POROSITY_SOLUTION_AND_BODY_FORCE_PROCESS_H