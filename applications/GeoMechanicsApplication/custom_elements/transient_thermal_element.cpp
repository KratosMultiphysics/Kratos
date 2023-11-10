// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Gennady Markelov
//

#include "custom_elements/transient_thermal_element.h"
#include "custom_constitutive/thermal_dispersion_law.h"
#include "includes/condition.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(IndexType NewId)
    : Element(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(IndexType NewId,
                                                                  GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(
    IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::~TransientThermalElement() = default;

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientThermalElement(
        NewId, GetGeometry().Create(rThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new TransientThermalElement(NewId, pGeom, pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                          const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int N_DOF = TNumNodes;
    if (rElementalDofList.size() != N_DOF) {
        rElementalDofList.resize(N_DOF);
    }

    const GeometryType& rGeom = GetGeometry();
    for (unsigned int i = 0; i < N_DOF; ++i) {
        rElementalDofList[i] = rGeom[i].pGetDof(TEMPERATURE);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int N_DOF = TNumNodes;
    if (rResult.size() != N_DOF) {
        rResult.resize(N_DOF, false);
    }

    const GeometryType& rGeom = GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = rGeom[i].GetDof(TEMPERATURE).EquationId();
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
int TransientThermalElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    CheckDomainSize();
    const GeometryType& rGeom = GetGeometry();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        CheckSolutionStepsData(i, TEMPERATURE);
        CheckSolutionStepsData(i, DT_TEMPERATURE);
        if (!rGeom[i].HasDofFor(TEMPERATURE)) {
            KRATOS_ERROR << "missing degree of freedom for TEMPERATURE on node "
                         << rGeom[i].Id() << std::endl;
        }
    }

    VerifyProperty(DENSITY_WATER);
    VerifyProperty(POROSITY);
    VerifyProperty(SATURATION);
    VerifyProperty(DENSITY_SOLID);
    VerifyProperty(SPECIFIC_HEAT_CAPACITY_WATER);
    VerifyProperty(SPECIFIC_HEAT_CAPACITY_SOLID);
    VerifyProperty(THERMAL_CONDUCTIVITY_WATER);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XX);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_YY);
    VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XY);
    VerifyProperty(LONGITUDINAL_DISPERSIVITY);
    VerifyProperty(TRANSVERSE_DISPERSIVITY);
    VerifyProperty(SOLID_COMPRESSIBILITY);

    if (TDim == 2) {
        auto pos = std::find_if(rGeom.begin(), rGeom.end(),
                                [](const auto& node) { return node.Z() != 0.0; });
        if (pos != rGeom.end()) {
            KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << pos->Id()
                         << std::endl;
        }
    }

    if (TDim > 2) {
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_ZZ);
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_YZ);
        VerifyProperty(THERMAL_CONDUCTIVITY_SOLID_XZ);
    }

    KRATOS_CATCH("")

    return 0;
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes>
TransientThermalElement<TDim, TNumNodes>::CalculateCapacityMatrix(const Vector& rIntegrationCoefficients) const
{
    const auto& r_properties = GetProperties();

    const double cWater = r_properties[POROSITY] * r_properties[SATURATION] *
                          r_properties[DENSITY_WATER] * r_properties[SPECIFIC_HEAT_CAPACITY_WATER];
    const double cSolid = (1.0 - r_properties[POROSITY]) *
                          r_properties[DENSITY_SOLID] * r_properties[SPECIFIC_HEAT_CAPACITY_SOLID];

    const auto& IntegrationPoints = GetGeometry().IntegrationPoints(GetIntegrationMethod());

    const auto NContainer = Matrix{GetGeometry().ShapeFunctionsValues(GetIntegrationMethod())};

    auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
    for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        const auto N = Vector{row(NContainer, GPoint)};
        result += (cWater + cSolid) * outer_prod(N, N) * rIntegrationCoefficients[GPoint];
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
Vector TransientThermalElement<TDim, TNumNodes>::CalculateIntegrationCoefficients(const Vector& detJContainer) const
{
    const auto& IntegrationPoints = GetGeometry().IntegrationPoints(GetIntegrationMethod());

    Vector result{IntegrationPoints.size()};
    for (unsigned int GPoint = 0; GPoint < IntegrationPoints.size(); ++GPoint) {
        result[GPoint] = IntegrationPoints[GPoint].Weight() * detJContainer[GPoint];
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes>
TransientThermalElement<TDim, TNumNodes>::CalculateConductivityMatrix(const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
                                                                      const Vector& rIntegrationCoefficients,
                                                                      const ProcessInfo& rCurrentProcessInfo) const
{
    GeoThermalDispersionLaw geo(TDim);
    const auto constitutive_matrix =
        geo.CalculateThermalDispersionMatrix(GetProperties(), rCurrentProcessInfo, GetGeometry());

    auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
    for (unsigned int GPoint = 0; GPoint < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod()); ++GPoint) {
        BoundedMatrix<double, TDim, TNumNodes> Temp = prod(constitutive_matrix, trans(rShapeFunctionGradients[GPoint]));
        result += prod(rShapeFunctionGradients[GPoint], Temp) * rIntegrationCoefficients[GPoint];
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int N_DOF = TNumNodes;
    noalias(rLeftHandSideMatrix)  = ZeroMatrix(N_DOF, N_DOF);
    noalias(rRightHandSideVector) = ZeroVector(N_DOF);

    const auto& rGeom = GetGeometry();
    const unsigned int NumGPoints = rGeom.IntegrationPointsNumber(GetIntegrationMethod());

    GeometryType::ShapeFunctionsGradientsType DN_DXContainer;
    Vector detJContainer{NumGPoints};
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, detJContainer, GetIntegrationMethod());

    const auto integration_coefficients = CalculateIntegrationCoefficients(detJContainer);

    const auto conductivity_matrix =
        CalculateConductivityMatrix(DN_DXContainer, integration_coefficients, rCurrentProcessInfo);
    const auto capacity_matrix = CalculateCapacityMatrix(integration_coefficients);

    AddContributionsToLhsMatrix(rLeftHandSideMatrix, conductivity_matrix, capacity_matrix,
                                rCurrentProcessInfo[DT_TEMPERATURE_COEFFICIENT]);

    const auto capacity_vector =
        array_1d<double, TNumNodes>{-prod(capacity_matrix, GetNodalValuesOf(DT_TEMPERATURE))};
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, capacity_vector);

    const auto conductivity_vector =
        array_1d<double, TNumNodes>{-prod(conductivity_matrix, GetNodalValuesOf(TEMPERATURE))};
    GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, conductivity_vector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod TransientThermalElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    using Data = GeometryData::IntegrationMethod;

    switch (TNumNodes) {
    case 3:
    case 6:
        return Data::GI_GAUSS_2;
    case 10:
        return Data::GI_GAUSS_4;
    case 15:
        return Data::GI_GAUSS_5;
    default:
        return Data::GI_GAUSS_2;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::VerifyProperty(Kratos::Variable<double>& rVariable) const
{
    const PropertiesType& rProp = GetProperties();
    if (!rProp.Has(rVariable)) {
        KRATOS_ERROR << rVariable.Name()
                     << " does not exist in the material properties." << std::endl;
    }
    else if (rProp[rVariable] < 0.0) {
        KRATOS_ERROR << rVariable.Name() << " has an invalid value at element"
                     << Id() << "." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CheckDomainSize() const
{
    const GeometryType& rGeom = GetGeometry();
    if (rGeom.DomainSize() < 1.0e-15) {
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << Id() << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CheckSolutionStepsData(
    int rId, Kratos::Variable<double>& rVariable) const
{
    const GeometryType& rGeom = GetGeometry();
    if (rGeom[rId].SolutionStepsDataHas(rVariable)) {
        KRATOS_ERROR << "missing variable " << rVariable.Name() << " on node "
                     << rGeom[rId].Id() << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
}

template class TransientThermalElement<2, 3>;
template class TransientThermalElement<2, 4>;
template class TransientThermalElement<2, 6>;
template class TransientThermalElement<2, 8>;
template class TransientThermalElement<2, 9>;
template class TransientThermalElement<2, 10>;
template class TransientThermalElement<2, 15>;
template class TransientThermalElement<3, 4>;
template class TransientThermalElement<3, 8>;
template class TransientThermalElement<3, 10>;
template class TransientThermalElement<3, 20>;
template class TransientThermalElement<3, 27>;

} // Namespace Kratos
