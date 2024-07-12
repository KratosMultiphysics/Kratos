// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes


// Application includes
#include "axisymmetric_thermal_face.h"

namespace Kratos
{

// Public Life Cycle //////////////////////////////////////////////////////////

AxisymmetricThermalFace::AxisymmetricThermalFace(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry)
    : ThermalFace(NewId, pGeometry)
{
}

AxisymmetricThermalFace::AxisymmetricThermalFace(
    IndexType NewId,
    typename GeometryType::Pointer pGeometry,
    typename PropertiesType::Pointer pProperties)
    : ThermalFace(NewId, pGeometry, pProperties)
{
}

// Public Operations //////////////////////////////////////////////////////////

Condition::Pointer AxisymmetricThermalFace::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    typename PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<AxisymmetricThermalFace>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer AxisymmetricThermalFace::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    typename PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<AxisymmetricThermalFace>(NewId, pGeom, pProperties);
}

// Protected Operations //////////////////////////////////////////////////////////

void AxisymmetricThermalFace::SetIntegrationWeight(
    const IndexType IntegrationPointIndex,
    const typename GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const Vector &rJacobianDeterminantsVector,
    ConditionDataStruct &rData)
{
    double radius = 0.0;
    const auto& r_geom = this->GetGeometry();
    for (IndexType i = 0; i < r_geom.PointsNumber(); ++i) {
        radius += rData.N[i] * r_geom[i].Y();
    }
    rData.Weight = 2.0 * Globals::Pi * radius * rJacobianDeterminantsVector[IntegrationPointIndex] * rIntegrationPoints[IntegrationPointIndex].Weight();
}

// Input and Output ///////////////////////////////////////////////////////////

std::string AxisymmetricThermalFace::Info() const
{
    std::stringstream buffer;
    buffer << "AxisymmetricThermalFace #" << Id();
    return buffer.str();
}

void AxisymmetricThermalFace::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "AxisymmetricThermalFace #" << Id();
}

void AxisymmetricThermalFace::PrintData(std::ostream& rOStream) const
{
    rOStream << "AxisymmetricThermalFace #" << Id() << std::endl;
    this->GetGeometry().PrintData(rOStream);
}

// Serialization //////////////////////////////////////////////////////////////

AxisymmetricThermalFace::AxisymmetricThermalFace():
    BaseType()
{
}

void AxisymmetricThermalFace::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

void AxisymmetricThermalFace::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

}
