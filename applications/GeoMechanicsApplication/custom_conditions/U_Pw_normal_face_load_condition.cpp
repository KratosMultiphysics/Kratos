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
#include "custom_conditions/U_Pw_normal_face_load_condition.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"

#include <numeric>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer UPwNormalFaceLoadCondition<TDim, TNumNodes>::Create(IndexType NewId,
                                                                       NodesArrayType const& ThisNodes,
                                                                       PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(
        new UPwNormalFaceLoadCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFaceLoadCondition<TDim, TNumNodes>::CalculateRHS(Vector& rRightHandSideVector,
                                                               const ProcessInfo& CurrentProcessInfo)
{
    const auto& r_geometry           = this->GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const auto  number_of_integration_points = r_integration_points.size();

    // Containers of variables at all integration points
    const Matrix& r_n_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType j_container(number_of_integration_points);
    for (auto& j : j_container)
        j.resize(TDim, r_geometry.LocalSpaceDimension(), false);
    r_geometry.Jacobian(j_container, this->GetIntegrationMethod());

    // Condition variables
    NormalFaceLoadVariables variables;
    this->InitializeConditionVariables(variables, r_geometry);

    // Loop over integration points
    for (unsigned int integration_point = 0; integration_point < number_of_integration_points; ++integration_point) {
        // Compute traction vector
        array_1d<double, TDim> traction_vector;
        this->CalculateTractionVector(traction_vector, j_container[integration_point],
                                      r_n_container, variables, integration_point);

        // Compute Nu Matrix
        BoundedMatrix<double, TDim, TNumNodes * TDim> nu = ZeroMatrix(TDim, TNumNodes * TDim);
        ConditionUtilities::CalculateNuMatrix<TDim, TNumNodes>(nu, r_n_container, integration_point);

        // Contributions to the right hand side
        GeoElementUtilities::AssembleUBlockVector(
            rRightHandSideVector,
            prod(trans(nu), traction_vector) *
                this->CalculateIntegrationCoefficient(integration_point, r_integration_points));
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFaceLoadCondition<TDim, TNumNodes>::InitializeConditionVariables(NormalFaceLoadVariables& rVariables,
                                                                               const GeometryType& rGeom)
{
    VariablesUtilities::GetNodalValues(rGeom, NORMAL_CONTACT_STRESS, rVariables.NormalStressVector.begin());

    if constexpr (TDim == 2) {
        VariablesUtilities::GetNodalValues(rGeom, TANGENTIAL_CONTACT_STRESS,
                                           rVariables.TangentialStressVector.begin());
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFaceLoadCondition<TDim, TNumNodes>::CalculateTractionVector(array_1d<double, TDim>& rTractionVector,
                                                                          const Matrix& Jacobian,
                                                                          const Matrix& NContainer,
                                                                          const NormalFaceLoadVariables& Variables,
                                                                          const unsigned int& GPoint)
{
    const auto shape_function_values = row(NContainer, GPoint);
    const auto normal_stress =
        std::inner_product(shape_function_values.begin(), shape_function_values.end(),
                           Variables.NormalStressVector.begin(), 0.0);

    Vector normal_vector = ZeroVector(3);
    if constexpr (TDim == 2) {
        const auto tangential_stress =
            std::inner_product(shape_function_values.begin(), shape_function_values.end(),
                               Variables.TangentialStressVector.begin(), 0.0);

        Vector tangential_vector = ZeroVector(3);
        std::copy_n(column(Jacobian, 0).begin(), TDim, tangential_vector.begin());
        Vector out_of_plane_vector = ZeroVector(3);
        out_of_plane_vector[2]     = 1.0;
        MathUtils<>::CrossProduct(normal_vector, out_of_plane_vector, tangential_vector);
        auto traction_vector = tangential_stress * tangential_vector + normal_stress * normal_vector;
        std::copy_n(traction_vector.begin(), TDim, rTractionVector.begin());
    } else if constexpr (TDim == 3) {
        MathUtils<>::CrossProduct(normal_vector, column(Jacobian, 0), column(Jacobian, 1));
        rTractionVector = -normal_stress * normal_vector;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double UPwNormalFaceLoadCondition<TDim, TNumNodes>::CalculateIntegrationCoefficient(
    const IndexType PointNumber, const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const
{
    return IntegrationPoints[PointNumber].Weight();
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string UPwNormalFaceLoadCondition<TDim, TNumNodes>::Info() const
{
    return "UPwNormalFaceLoadCondition";
}

template class UPwNormalFaceLoadCondition<2, 2>;
template class UPwNormalFaceLoadCondition<2, 3>;
template class UPwNormalFaceLoadCondition<2, 4>;
template class UPwNormalFaceLoadCondition<2, 5>;
template class UPwNormalFaceLoadCondition<3, 3>;
template class UPwNormalFaceLoadCondition<3, 4>;
template class UPwNormalFaceLoadCondition<3, 6>;
template class UPwNormalFaceLoadCondition<3, 8>;

} // Namespace Kratos.
