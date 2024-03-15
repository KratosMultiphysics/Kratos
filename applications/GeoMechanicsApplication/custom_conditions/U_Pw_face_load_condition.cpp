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

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwFaceLoadCondition<TDim,TNumNodes>::
    Create( IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer(new UPwFaceLoadCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwFaceLoadCondition<TDim,TNumNodes>::
    CalculateRHS( VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{        
    //Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i<NumGPoints; ++i)
        (JContainer[i]).resize(TDim, LocalDim, false);
    Geom.Jacobian(JContainer, this->GetIntegrationMethod());

    //Condition variables
    array_1d<double,TNumNodes*TDim> FaceLoadVector;
    ConditionUtilities::GetFaceLoadVector<TNumNodes>(FaceLoadVector, Geom);
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TDim> TractionVector;
    array_1d<double,TNumNodes*TDim> UVector;

    //Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute traction vector 
        ConditionUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(TractionVector,
                                                                               NContainer,
                                                                               FaceLoadVector,
                                                                               GPoint);

        //Compute Nu Matrix
        ConditionUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nu,NContainer,GPoint);

        //Compute weighting coefficient for integration
        double integration_coefficient = 
            ConditionUtilities::CalculateIntegrationCoefficient<TDim, TNumNodes>(
            JContainer[GPoint], IntegrationPoints[GPoint].Weight());

        //Contributions to the right hand side
        noalias(UVector) = prod(trans(Nu),TractionVector) * integration_coefficient;
        GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector, UVector);
    }
}


template class UPwFaceLoadCondition<2,2>;
template class UPwFaceLoadCondition<2,3>;
template class UPwFaceLoadCondition<2,4>;
template class UPwFaceLoadCondition<2,5>;
template class UPwFaceLoadCondition<3,3>;
template class UPwFaceLoadCondition<3,4>;

} // Namespace Kratos.
