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
#include "custom_conditions/U_Pw_face_load_condition.hpp"
#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer UPwFaceLoadCondition<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                 NodesArrayType const& ThisNodes,
                                                                 PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwFaceLoadCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwFaceLoadCondition<TDim, TNumNodes>::CalculateRHS(Vector&            rRightHandSideVector,
                                                         const ProcessInfo& CurrentProcessInfo)
{
    const GeometryType&                             r_geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        r_geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = r_integration_points.size();

    // Containers of variables at all integration points
    const Matrix& r_n_container = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType j_container(NumGPoints);
    for (auto& j : j_container)
        j.resize(TDim, r_geom.LocalSpaceDimension(), false);
    r_geom.Jacobian(j_container, this->GetIntegrationMethod());

    // Condition variables
    array_1d<double, TNumNodes * TDim> face_load_vector;
    ConditionUtilities::GetFaceLoadVector<TDim, TNumNodes>(face_load_vector, r_geom);

    for (unsigned int integration_point = 0; integration_point < NumGPoints; ++integration_point) {
        // Compute traction vector
        array_1d<double, TDim> traction_vector;
        ConditionUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            traction_vector, r_n_container, face_load_vector, integration_point);

        // Compute Nu Matrix
        BoundedMatrix<double, TDim, TNumNodes * TDim> nu = ZeroMatrix(TDim, TNumNodes * TDim);
        ConditionUtilities::CalculateNuMatrix<TDim, TNumNodes>(nu, r_n_container, integration_point);

        // Compute weighting coefficient for integration
        auto integration_coefficient = ConditionUtilities::CalculateIntegrationCoefficient<TDim, TNumNodes>(
            j_container[integration_point], r_integration_points[integration_point].Weight());

        // Contributions to the right hand side
        GeoElementUtilities::AssembleUBlockVector(
            rRightHandSideVector, prod(trans(nu), traction_vector) * integration_coefficient);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string UPwFaceLoadCondition<TDim, TNumNodes>::Info() const
{
    return "UPwFaceLoadCondition";
}

template class UPwFaceLoadCondition<2, 2>;
template class UPwFaceLoadCondition<2, 3>;
template class UPwFaceLoadCondition<2, 4>;
template class UPwFaceLoadCondition<2, 5>;
template class UPwFaceLoadCondition<3, 3>;
template class UPwFaceLoadCondition<3, 4>;

} // Namespace Kratos.
