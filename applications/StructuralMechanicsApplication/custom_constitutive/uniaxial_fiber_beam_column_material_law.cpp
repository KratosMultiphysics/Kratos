// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
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
#include "custom_constitutive/uniaxial_fiber_beam_column_material_law.hpp"

namespace Kratos {

UniaxialFiberBeamColumnMaterialLaw::UniaxialFiberBeamColumnMaterialLaw() {}

UniaxialFiberBeamColumnMaterialLaw::UniaxialFiberBeamColumnMaterialLaw(PropertiesType::Pointer pProperties)
    : mpProperties(pProperties)
{
}

void UniaxialFiberBeamColumnMaterialLaw::Initialize()
{
    KRATOS_ERROR << "Calling virtual function for Initialize" << std::endl;
}

double& UniaxialFiberBeamColumnMaterialLaw::GetTangentModulus()
{
    KRATOS_ERROR << "Calling virtual function for GetTangentModulus" << std::endl;
}

double& UniaxialFiberBeamColumnMaterialLaw::GetStress()
{
    KRATOS_ERROR << "Calling virtual function for GetStress" << std::endl;
}

void UniaxialFiberBeamColumnMaterialLaw::SetStrain(const double Strain)
{
    KRATOS_ERROR << "Calling virtual function for SetStrain" << std::endl;
}

void UniaxialFiberBeamColumnMaterialLaw::CalculateMaterialResponse()
{
    KRATOS_ERROR << "Calling virtual function for CalculateMaterialResponse" << std::endl;
}

void UniaxialFiberBeamColumnMaterialLaw::FinalizeMaterialResponse()
{
    KRATOS_ERROR << "Calling virtual function for FinalizeMaterialResponse" << std::endl;
}


std::string UniaxialFiberBeamColumnMaterialLaw::Info() const {
    std::stringstream buffer;
    buffer << "UniaxialFiberBeamColumnMaterialLaw";
    return buffer.str();
}

void UniaxialFiberBeamColumnMaterialLaw::PrintInfo(std::ostream& rOStream) const {
    rOStream << "UniaxialFiberBeamColumnMaterialLaw";
}

void UniaxialFiberBeamColumnMaterialLaw::PrintData(std::ostream& rOStream) const {
    rOStream << "UniaxialFiberBeamColumnMaterialLaw has no data" << std::endl;
}

void UniaxialFiberBeamColumnMaterialLaw::save(Serializer& rSerializer) const
{
    rSerializer.save("mpProperties", mpProperties);
}
void UniaxialFiberBeamColumnMaterialLaw::load(Serializer& rSerializer)
{
    rSerializer.load("mpProperties", mpProperties);
}

} // namespace Kratos.
