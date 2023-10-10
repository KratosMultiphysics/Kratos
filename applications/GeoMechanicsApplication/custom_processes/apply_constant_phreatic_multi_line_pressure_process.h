// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//                   Jonathan Nuttall

#pragma once

#include <algorithm>
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyConstantPhreaticMultiLinePressureProcess : public Process
{

public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantPhreaticMultiLinePressureProcess);

    ApplyConstantPhreaticMultiLinePressureProcess(ModelPart& model_part, Parameters rParameters);
    ~ApplyConstantPhreaticMultiLinePressureProcess() override = default;
    ApplyConstantPhreaticMultiLinePressureProcess& operator=(ApplyConstantPhreaticMultiLinePressureProcess const&) = delete;
    ApplyConstantPhreaticMultiLinePressureProcess(ApplyConstantPhreaticMultiLinePressureProcess const&) = delete;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    std::string Info() const override;
    void PrintInfo(std::ostream& rOStream) const override;
    const std::string& VariableName() const;
    bool IsFixed() const;
    bool IsSeepage() const;
    unsigned int GravityDirection() const;
    double SpecificWeight() const;
    unsigned int OutOfPlaneDirection() const;
    unsigned int HorizontalDirection() const;
    const Vector& HorizontalDirectionCoordinates() const;
    const Vector& GravityDirectionCoordinates() const;
    const Vector& XCoordinates() const;
    const Vector& YCoordinates() const;
    const Vector& ZCoordinates() const;
    double PressureTensionCutOff() const;
    bool IsFixedProvided() const;

protected:
    ModelPart& mrModelPart;
    int findIndex(const Node& rNode) const;
    double CalculatePressure(const Node& rNode, std::vector<double> deltaH = {}) const;

private:
    std::string mVariableName;
    bool mIsFixed;
    bool mIsFixedProvided;
    bool mIsSeepage;
    unsigned int mGravityDirection;
    double mSpecificWeight;
    unsigned int mOutOfPlaneDirection;
    unsigned int mHorizontalDirection;
    Vector mHorizontalDirectionCoordinates;
    Vector mGravityDirectionCoordinates;
    Vector mXCoordinates;
    Vector mYCoordinates;
    Vector mZCoordinates;
    double mPressureTensionCutOff;

    void InitializeHorizontalDirection();
    void InitializeGravityDirection(const Parameters &rParameters);
    void InitializeCoordinates(const Parameters &rParameters);
    void ValidateCoordinates(const Parameters &rParameters) const;

    void InitializeParameters(Parameters &rParameters) const;
}; // Class ApplyConstantPhreaticMultiLinePressureProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantPhreaticMultiLinePressureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantPhreaticMultiLinePressureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.