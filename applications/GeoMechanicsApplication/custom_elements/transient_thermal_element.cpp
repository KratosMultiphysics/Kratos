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
#include "custom_constitutive/thermal_filter_law.h"
#include "custom_retention/retention_law_factory.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/variables_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(IndexType NewId) : Element(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(IndexType             NewId,
                                                                  GeometryType::Pointer pGeometry,
                                                                  std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : Element(NewId, pGeometry), mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
{
}

template <unsigned int TDim, unsigned int TNumNodes>
TransientThermalElement<TDim, TNumNodes>::TransientThermalElement(IndexType             NewId,
                                                                  GeometryType::Pointer pGeometry,
                                                                  PropertiesType::Pointer pProperties,
                                                                  std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : Element(NewId, pGeometry, pProperties),
      mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                  const NodesArrayType& rThisNodes,
                                                                  PropertiesType::Pointer pProperties) const
{
    return make_intrusive<TransientThermalElement>(NewId, GetGeometry().Create(rThisNodes), pProperties,
                                                   this->CloneIntegrationCoefficientModifier());
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                  GeometryType::Pointer pGeom,
                                                                  PropertiesType::Pointer pProperties) const
{
    return make_intrusive<TransientThermalElement>(NewId, pGeom, pProperties,
                                                   this->CloneIntegrationCoefficientModifier());
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                                const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                                    VectorType& rRightHandSideVector,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType::ShapeFunctionsGradientsType dN_dX_container;
    Vector                                    det_J_container;

    // ShapeFunctionsIntegrationsPointsGradients does not allow for the line element in 2D/3D
    // configuration and will produce errors. To circumvent this, the dN_dX_container is
    // separately computed with correct dimensions for the line element.
    if (GetGeometry().LocalSpaceDimension() == 1) {
        GetGeometry().DeterminantOfJacobian(det_J_container, this->GetIntegrationMethod());
        dN_dX_container = GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        std::transform(dN_dX_container.begin(), dN_dX_container.end(), det_J_container.begin(),
                       dN_dX_container.begin(), std::divides<>());
    } else {
        GetGeometry().ShapeFunctionsIntegrationPointsGradients(dN_dX_container, det_J_container,
                                                               GetIntegrationMethod());
    }

    const auto integration_coefficients = mIntegrationCoefficientsCalculator.Run<Vector>(
        GetGeometry().IntegrationPoints(GetIntegrationMethod()), det_J_container, this);
    const auto conductivity_matrix = CalculateConductivityMatrix(dN_dX_container, integration_coefficients);
    const auto capacity_matrix = CalculateCapacityMatrix(integration_coefficients);

    AddContributionsToLhsMatrix(rLeftHandSideMatrix, conductivity_matrix, capacity_matrix,
                                rCurrentProcessInfo[DT_TEMPERATURE_COEFFICIENT]);
    AddContributionsToRhsVector(rRightHandSideVector, conductivity_matrix, capacity_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod TransientThermalElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    switch (this->GetGeometry().GetGeometryOrderType()) {
        using enum GeometryData::KratosGeometryOrderType;
        using enum GeometryData::KratosGeometryType;
        using enum GeometryData::IntegrationMethod;
    case Kratos_Cubic_Order:
        return this->GetGeometry().GetGeometryType() == Kratos_Triangle2D10 ? GI_GAUSS_4 : GI_GAUSS_3;
    case Kratos_Quartic_Order:
        return GI_GAUSS_5;
    default:
        return GI_GAUSS_2;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int TransientThermalElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto element_Id = this->Id();
    CheckUtilities::CheckDomainSize(GetGeometry().DomainSize(), element_Id);

    CheckUtilities::CheckHasNodalSolutionStepData(
        GetGeometry(), {std::cref(TEMPERATURE), std::cref(DT_TEMPERATURE)});
    CheckUtilities::CheckHasDofs(GetGeometry(), {std::cref(TEMPERATURE)});

    const auto&           r_properties = this->GetProperties();
    const CheckProperties check_properties(r_properties, "properties", element_Id,
                                           CheckProperties::Bounds::AllInclusive);
    check_properties.Check(DENSITY_WATER);
    constexpr auto max_value = 1.0;
    check_properties.Check(POROSITY, max_value);
    check_properties.Check(DENSITY_SOLID);
    check_properties.Check(SPECIFIC_HEAT_CAPACITY_WATER);
    check_properties.Check(SPECIFIC_HEAT_CAPACITY_SOLID);
    check_properties.Check(THERMAL_CONDUCTIVITY_WATER);
    check_properties.Check(THERMAL_CONDUCTIVITY_SOLID_XX);
    check_properties.Check(THERMAL_CONDUCTIVITY_SOLID_YY);
    check_properties.Check(THERMAL_CONDUCTIVITY_SOLID_XY);

    if constexpr (TDim == 3) {
        check_properties.Check(THERMAL_CONDUCTIVITY_SOLID_ZZ);
        check_properties.Check(THERMAL_CONDUCTIVITY_SOLID_YZ);
        check_properties.Check(THERMAL_CONDUCTIVITY_SOLID_XZ);
    }

    if constexpr (TDim == 2) CheckUtilities::CheckForNonZeroZCoordinateIn2D(GetGeometry());

    check_properties.CheckAvailabilityAndEquality(RETENTION_LAW, "SaturatedLaw");
    return RetentionLawFactory::Clone(r_properties)->Check(r_properties, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
std::unique_ptr<IntegrationCoefficientModifier> TransientThermalElement<TDim, TNumNodes>::CloneIntegrationCoefficientModifier() const
{
    return mIntegrationCoefficientsCalculator.CloneModifier();
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::AddContributionsToLhsMatrix(
    MatrixType&                                        rLeftHandSideMatrix,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rConductivityMatrix,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rCapacityMatrix,
    double                                             DtTemperatureCoefficient)
{
    rLeftHandSideMatrix = rConductivityMatrix;
    rLeftHandSideMatrix += (DtTemperatureCoefficient * rCapacityMatrix);
}

template <unsigned int TDim, unsigned int TNumNodes>
void TransientThermalElement<TDim, TNumNodes>::AddContributionsToRhsVector(
    VectorType&                                        rRightHandSideVector,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rConductivityMatrix,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rCapacityMatrix) const
{
    const auto capacity_vector     = array_1d<double, TNumNodes>{-prod(
        rCapacityMatrix, VariablesUtilities::GetNodalValuesOf<TNumNodes>(DT_TEMPERATURE, this->GetGeometry()))};
    rRightHandSideVector           = capacity_vector;
    const auto conductivity_vector = array_1d<double, TNumNodes>{-prod(
        rConductivityMatrix, VariablesUtilities::GetNodalValuesOf<TNumNodes>(TEMPERATURE, this->GetGeometry()))};
    rRightHandSideVector += conductivity_vector;
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> TransientThermalElement<TDim, TNumNodes>::CalculateConductivityMatrix(
    const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients, const Vector& rIntegrationCoefficients) const
{
    const auto law                 = CreateThermalLaw();
    const auto constitutive_matrix = law->CalculateThermalDispersionMatrix(GetProperties());

    auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
    for (unsigned int integration_point_index = 0;
         integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
         ++integration_point_index) {
        BoundedMatrix<double, TDim, TNumNodes> Temp =
            prod(constitutive_matrix, trans(rShapeFunctionGradients[integration_point_index]));
        result += prod(rShapeFunctionGradients[integration_point_index], Temp) *
                  rIntegrationCoefficients[integration_point_index];
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> TransientThermalElement<TDim, TNumNodes>::CalculateCapacityMatrix(
    const Vector& rIntegrationCoefficients) const
{
    const auto&              r_properties = GetProperties();
    RetentionLaw::Parameters parameters(r_properties);
    auto                     retention_law = RetentionLawFactory::Clone(r_properties);
    const double             saturation    = retention_law->CalculateSaturation(parameters);
    const auto c_water = r_properties[POROSITY] * saturation * r_properties[DENSITY_WATER] *
                         r_properties[SPECIFIC_HEAT_CAPACITY_WATER];
    const auto c_solid = (1.0 - r_properties[POROSITY]) * r_properties[DENSITY_SOLID] *
                         r_properties[SPECIFIC_HEAT_CAPACITY_SOLID];

    auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
    const auto& r_N_container = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());
    for (unsigned int integration_point_index = 0;
         integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
         ++integration_point_index) {
        const auto N = Vector{row(r_N_container, integration_point_index)};
        result += (c_water + c_solid) * outer_prod(N, N) * rIntegrationCoefficients[integration_point_index];
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
std::unique_ptr<GeoThermalLaw> TransientThermalElement<TDim, TNumNodes>::CreateThermalLaw() const
{
    const std::size_t              number_of_dimensions = GetGeometry().LocalSpaceDimension();
    std::unique_ptr<GeoThermalLaw> law;

    if (GetProperties().Has(THERMAL_LAW_NAME)) {
        const std::string& rThermalLawName = GetProperties()[THERMAL_LAW_NAME];
        if (rThermalLawName == "GeoThermalDispersionLaw") {
            law = std::make_unique<GeoThermalDispersionLaw>(number_of_dimensions);
        } else if (rThermalLawName == "GeoThermalFilterLaw") {
            law = std::make_unique<GeoThermalFilterLaw>();
        } else {
            KRATOS_ERROR << "Undefined THERMAL_LAW_NAME! " << rThermalLawName << std::endl;
        }
    } else {
        law = std::make_unique<GeoThermalDispersionLaw>(number_of_dimensions);
    }
    return law;
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::DofsVectorType TransientThermalElement<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), TEMPERATURE);
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

template class TransientThermalElement<2, 2>;
template class TransientThermalElement<2, 3>;
template class TransientThermalElement<2, 4>;
template class TransientThermalElement<2, 5>;
template class TransientThermalElement<2, 6>;
template class TransientThermalElement<2, 8>;
template class TransientThermalElement<2, 9>;
template class TransientThermalElement<2, 10>;
template class TransientThermalElement<2, 15>;
template class TransientThermalElement<3, 2>;
template class TransientThermalElement<3, 3>;
template class TransientThermalElement<3, 4>;
template class TransientThermalElement<3, 8>;
template class TransientThermalElement<3, 10>;
template class TransientThermalElement<3, 20>;
template class TransientThermalElement<3, 27>;

} // Namespace Kratos
