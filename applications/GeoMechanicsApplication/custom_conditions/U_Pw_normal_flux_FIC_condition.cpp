// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// Application includes
#include "custom_conditions/U_Pw_normal_flux_FIC_condition.h"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"

#include <numeric>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
UPwNormalFluxFICCondition<TDim, TNumNodes>::UPwNormalFluxFICCondition()
    : UPwNormalFluxFICCondition(0, nullptr, nullptr)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwNormalFluxFICCondition<TDim, TNumNodes>::UPwNormalFluxFICCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : UPwNormalFluxFICCondition(NewId, pGeometry, nullptr)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
UPwNormalFluxFICCondition<TDim, TNumNodes>::UPwNormalFluxFICCondition(IndexType NewId,
                                                                      GeometryType::Pointer pGeometry,
                                                                      PropertiesType::Pointer pProperties)
    : UPwNormalFluxCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer UPwNormalFluxFICCondition<TDim, TNumNodes>::Create(IndexType NewId,
                                                                      NodesArrayType const& ThisNodes,
                                                                      PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(
        new UPwNormalFluxFICCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod UPwNormalFluxFICCondition<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxFICCondition<TDim, TNumNodes>::CalculateAll(Matrix& rLeftHandSideMatrix,
                                                              Vector& rRightHandSideVector,
                                                              const ProcessInfo& CurrentProcessInfo)
{
    const auto&                                     r_prop = this->GetProperties();
    const auto&                                     r_geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int number_of_integration_points = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& n_container = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType j_container(number_of_integration_points);
    for (auto& j : j_container)
        j.resize(TDim, r_geom.LocalSpaceDimension(), false);
    r_geom.Jacobian(j_container, this->GetIntegrationMethod());

    // Condition variables
    NormalFluxVariables    Variables;
    NormalFluxFICVariables FICVariables;
    FICVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];
    this->CalculateElementLength(FICVariables.ElementLength, r_geom);
    const double BulkModulusSolid = r_prop[BULK_MODULUS_SOLID];
    const double Porosity         = r_prop[POROSITY];
    const double BulkModulus = r_prop[YOUNG_MODULUS] / (3.0 * (1.0 - 2.0 * r_prop[POISSON_RATIO]));
    const double BiotCoefficient = 1.0 - BulkModulus / BulkModulusSolid;
    FICVariables.BiotModulusInverse =
        (BiotCoefficient - Porosity) / BulkModulusSolid + Porosity / r_prop[BULK_MODULUS_FLUID];

    auto normal_flux_vector = VariablesUtilities::GetNodalValues<TNumNodes>(r_geom, NORMAL_FLUID_FLUX);
    FICVariables.DtPressureVector = VariablesUtilities::GetNodalValues<TNumNodes>(r_geom, DT_WATER_PRESSURE);

    // Loop over integration points
    for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
         integration_point++) {
        noalias(Variables.Np) = row(n_container, integration_point);

        const auto normal_flux =
            std::inner_product(Variables.Np.begin(), Variables.Np.end(), normal_flux_vector.begin(), 0.0);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient = ConditionUtilities::CalculateIntegrationCoefficient(
            j_container[integration_point], IntegrationPoints[integration_point].Weight());

        // Contributions to the left hand side
        this->CalculateAndAddLHSStabilization(rLeftHandSideMatrix, Variables, FICVariables);

        // Contributions to the right hand side
        GeoElementUtilities::AssemblePBlockVector(
            rRightHandSideVector, -normal_flux * Variables.Np * Variables.IntegrationCoefficient);

        this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxFICCondition<TDim, TNumNodes>::CalculateRHS(Vector& rRightHandSideVector,
                                                              const ProcessInfo& CurrentProcessInfo)
{
    const auto&                                     r_prop = this->GetProperties();
    const auto&                                     r_geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int number_of_integration_points = IntegrationPoints.size();

    // Containers of variables at all integration points
    const Matrix& n_container = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType j_container(number_of_integration_points);
    for (auto& j : j_container)
        j.resize(TDim, r_geom.LocalSpaceDimension(), false);
    r_geom.Jacobian(j_container, this->GetIntegrationMethod());

    // Condition variables
    NormalFluxVariables    Variables;
    NormalFluxFICVariables FICVariables;
    FICVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];
    this->CalculateElementLength(FICVariables.ElementLength, r_geom);
    const double BulkModulusSolid = r_prop[BULK_MODULUS_SOLID];
    const double Porosity         = r_prop[POROSITY];
    const double BulkModulus = r_prop[YOUNG_MODULUS] / (3.0 * (1.0 - 2.0 * r_prop[POISSON_RATIO]));
    const double BiotCoefficient = 1.0 - BulkModulus / BulkModulusSolid;
    FICVariables.BiotModulusInverse =
        (BiotCoefficient - Porosity) / BulkModulusSolid + Porosity / r_prop[BULK_MODULUS_FLUID];
    const auto normal_flux_vector = VariablesUtilities::GetNodalValues<TNumNodes>(r_geom, NORMAL_FLUID_FLUX);
    FICVariables.DtPressureVector = VariablesUtilities::GetNodalValues<TNumNodes>(r_geom, DT_WATER_PRESSURE);

    for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
         integration_point++) {
        noalias(Variables.Np) = row(n_container, integration_point);

        const auto normal_flux =
            std::inner_product(Variables.Np.begin(), Variables.Np.end(), normal_flux_vector.begin(), 0.0);

        // Compute weighting coefficient for integration
        Variables.IntegrationCoefficient = ConditionUtilities::CalculateIntegrationCoefficient(
            j_container[integration_point], IntegrationPoints[integration_point].Weight());

        // Contributions to the right hand side
        GeoElementUtilities::AssemblePBlockVector(
            rRightHandSideVector, -normal_flux * Variables.Np * Variables.IntegrationCoefficient);

        this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
    }
}

template <>
void UPwNormalFluxFICCondition<2, 2>::CalculateElementLength(double& rElementLength, const GeometryType& rGeom)
{
    rElementLength = rGeom.Length();
}

template <>
void UPwNormalFluxFICCondition<3, 3>::CalculateElementLength(double& rElementLength, const GeometryType& rGeom)
{
    rElementLength = sqrt(4.0 * rGeom.Area() / Globals::Pi);
}

template <>
void UPwNormalFluxFICCondition<3, 4>::CalculateElementLength(double& rElementLength, const GeometryType& rGeom)
{
    rElementLength = sqrt(4.0 * rGeom.Area() / Globals::Pi);
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxFICCondition<TDim, TNumNodes>::CalculateAndAddLHSStabilization(
    Matrix& rLeftHandSideMatrix, NormalFluxVariables& rVariables, NormalFluxFICVariables& rFICVariables)
{
    this->CalculateAndAddBoundaryMassMatrix(rLeftHandSideMatrix, rVariables, rFICVariables);
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxFICCondition<TDim, TNumNodes>::CalculateAndAddBoundaryMassMatrix(
    Matrix& rLeftHandSideMatrix, const NormalFluxVariables& rVariables, const NormalFluxFICVariables& rFICVariables)
{
    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rFICVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    // Distribute boundary mass matrix into the elemental matrix
    // it seems the factor of 1/6 comes when Eq. 2.56 substituted into Eqs.2.69/2.70 in Pouplana's PhD thesis.
    GeoElementUtilities::AssemblePPBlockMatrix(
        rLeftHandSideMatrix, compressibility_matrix * rFICVariables.DtPressureCoefficient *
                                 rFICVariables.ElementLength / 6.0);
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxFICCondition<TDim, TNumNodes>::CalculateAndAddRHSStabilization(
    Vector& rRightHandSideVector, NormalFluxVariables& rVariables, NormalFluxFICVariables& rFICVariables)
{
    this->CalculateAndAddBoundaryMassFlow(rRightHandSideVector, rVariables, rFICVariables);
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxFICCondition<TDim, TNumNodes>::CalculateAndAddBoundaryMassFlow(
    Vector& rRightHandSideVector, NormalFluxVariables& rVariables, const NormalFluxFICVariables& rFICVariables)
{
    const auto compressibility_matrix = GeoTransportEquationUtilities::CalculateCompressibilityMatrix(
        rVariables.Np, rFICVariables.BiotModulusInverse, rVariables.IntegrationCoefficient);

    // it seems the factor of 1/6 comes when Eq. 2.56 substituted into Eqs.2.69/2.70 in Pouplana's
    // PhD thesis. Distribute boundary mass flow vector into elemental vector
    GeoElementUtilities::AssemblePBlockVector(
        rRightHandSideVector,
        prod(compressibility_matrix * rFICVariables.ElementLength / 6.0, rFICVariables.DtPressureVector));
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string UPwNormalFluxFICCondition<TDim, TNumNodes>::Info() const
{
    return "UPwNormalFluxFICCondition";
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxFICCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFluxFICCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}
template class UPwNormalFluxFICCondition<2, 2>;
template class UPwNormalFluxFICCondition<3, 3>;
template class UPwNormalFluxFICCondition<3, 4>;

} // Namespace Kratos.
