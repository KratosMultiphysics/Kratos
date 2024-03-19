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

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwNormalFaceLoadCondition<TDim,TNumNodes>::
    Create( IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwNormalFaceLoadCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFaceLoadCondition<TDim,TNumNodes>::
    CalculateRHS(VectorType& rRightHandSideVector,
                 const ProcessInfo& CurrentProcessInfo)
{
    //Previous definitions
    const GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = rGeom.LocalSpaceDimension();

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType JContainer(NumGPoints);
    for (unsigned int i = 0; i<NumGPoints; ++i)
        (JContainer[i]).resize(TDim,LocalDim,false);
    rGeom.Jacobian(JContainer, this->GetIntegrationMethod());

    //Condition variables
    NormalFaceLoadVariables Variables;
    this->InitializeConditionVariables(Variables, rGeom);
    array_1d<double,TDim> TractionVector;
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TNumNodes*TDim> UVector;

    //Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute traction vector
        this->CalculateTractionVector(TractionVector, JContainer[GPoint], NContainer, Variables, GPoint);

        //Compute Nu Matrix
        ConditionUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nu,NContainer,GPoint);

        //Compute weighting coefficient for integration
        const double IntegrationCoefficient = 
            this->CalculateIntegrationCoefficient(GPoint, IntegrationPoints);

        //Contributions to the right hand side
        noalias(UVector) = prod(trans(Nu),TractionVector) * IntegrationCoefficient;
        GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, UVector);
    }
}


template<unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFaceLoadCondition<TDim, TNumNodes>::InitializeConditionVariables(
    NormalFaceLoadVariables& rVariables,
    const GeometryType& rGeom)
{
    std::transform(rGeom.begin(), rGeom.end(), rVariables.NormalStressVector.begin(), [](const auto& node) {
        return node.FastGetSolutionStepValue(NORMAL_CONTACT_STRESS); });

    if constexpr (TDim == 2) {
        std::transform(rGeom.begin(), rGeom.end(), rVariables.TangentialStressVector.begin(), [](const auto& node) {
            return node.FastGetSolutionStepValue(TANGENTIAL_CONTACT_STRESS); });
    }
}


template<unsigned int TDim, unsigned int TNumNodes>
void UPwNormalFaceLoadCondition<TDim, TNumNodes>::CalculateTractionVector(
    array_1d<double, TDim>& rTractionVector,
    const Matrix& Jacobian,
    const Matrix& NContainer,
    const NormalFaceLoadVariables& Variables,
    const unsigned int& GPoint)
{
    Vector NormalVector = ZeroVector(TDim);
    const double NormalStress = MathUtils<>::Dot(row(NContainer, GPoint), Variables.NormalStressVector);

    if constexpr (TDim == 2) {
        const double TangentialStress = MathUtils<>::Dot(row(NContainer, GPoint), Variables.TangentialStressVector);
        NormalVector = column(Jacobian, 0);
        rTractionVector[0] = TangentialStress * NormalVector[0] - NormalStress * NormalVector[1];
        rTractionVector[1] = NormalStress * NormalVector[0] + TangentialStress * NormalVector[1];
    }
    else if constexpr (TDim == 3) {
        MathUtils<double>::CrossProduct(NormalVector, column(Jacobian, 0), column(Jacobian, 1));
        rTractionVector = NormalStress * NormalVector;
    }
}


template< unsigned int TDim, unsigned int TNumNodes >
double UPwNormalFaceLoadCondition<TDim,TNumNodes>::CalculateIntegrationCoefficient(
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const
{
    return IntegrationPoints[PointNumber].Weight();
}


template class UPwNormalFaceLoadCondition<2,2>;
template class UPwNormalFaceLoadCondition<2,3>;
template class UPwNormalFaceLoadCondition<2,4>;
template class UPwNormalFaceLoadCondition<2,5>;
template class UPwNormalFaceLoadCondition<3,3>;
template class UPwNormalFaceLoadCondition<3,4>;

} // Namespace Kratos.
