//  KRATOS  ____       _               _
//         / ___|  ___(_)___ _ __ ___ (_) ___
//         \___ \ / _ \ / __| '_ ` _ \| |/ __|
//          ___) |  __/ \__ \ | | | | | | (__
//         |____/ \___|_|___/_| |_| |_|_|\___|
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_classes/fiber_beam_column_component.h"

namespace Kratos {

FiberBeamColumnComponent::FiberBeamColumnComponent(IndexType NewId)
    : mId(NewId) {}


void FiberBeamColumnComponent::Initialize()
{
    KRATOS_ERROR <<  "Calling virtual function for Initialize"<< std::endl;
}

void FiberBeamColumnComponent::GetGlobalFlexibilityMatrix(Matrix& rGlobalFlexibilityMatrix)
{
    KRATOS_ERROR <<  "Calling virtual function for GetGlobalFlexibilityMatrix"<< std::endl;
}

void FiberBeamColumnComponent::GetGlobalStiffnessMatrix(Matrix& rGlobalStiffnessMatrix)
{
    KRATOS_ERROR <<  "Calling virtual function for GetGlobalStiffnessMatrix"<< std::endl;
}

void FiberBeamColumnComponent::GetGlobalInternalForces(Vector& rGlobalForces)
{
    KRATOS_ERROR <<  "Calling virtual function for GetGlobalInternalForces"<< std::endl;
}

bool FiberBeamColumnComponent::StateDetermination(const Vector& rElementForceIncrements)
{
    KRATOS_ERROR <<  "Calling virtual function for StateDetermination"<< std::endl;
}

bool FiberBeamColumnComponent::StateDetermination(const Vector& rElementForceIncrements, const double Tolerance)
{
    KRATOS_ERROR <<  "Calling virtual function for StateDetermination"<< std::endl;
}

void FiberBeamColumnComponent::FinalizeSolutionStep()
{
    KRATOS_ERROR <<  "Calling virtual function for FinalizeSolutionStep"<< std::endl;
}

std::string FiberBeamColumnComponent::Info() const {
    std::stringstream buffer;
    buffer << "FiberBeamColumnComponent #" << mId;
    return buffer.str();
}

void FiberBeamColumnComponent::PrintInfo(std::ostream& rOStream) const {
    rOStream << Info();
}

void FiberBeamColumnComponent::PrintData(std::ostream& rOStream) const {
}

void FiberBeamColumnComponent::save(Serializer& rSerializer) const
{
    rSerializer.save("mId", mId);
}
void FiberBeamColumnComponent::load(Serializer& rSerializer)
{
    rSerializer.load("mId", mId);
}

} // namespace Kratos.
