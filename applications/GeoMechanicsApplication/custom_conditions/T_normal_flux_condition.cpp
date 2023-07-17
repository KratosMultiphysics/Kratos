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
//  Main authors:    Mohamed Nabi
//                   
//                   
//


// Application includes
#include "utilities/math_utils.h"
#include "custom_conditions/T_normal_flux_condition.hpp"

namespace Kratos
{

// ============================================================================================
// ============================================================================================
// Default constructor
template<unsigned int TDim, unsigned int TNumNodes>
TNormalFluxCondition<TDim, TNumNodes>::TNormalFluxCondition() : TCondition<TDim, TNumNodes>() {}

// Constructor 1
template<unsigned int TDim, unsigned int TNumNodes>
TNormalFluxCondition<TDim, TNumNodes>::TNormalFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : TCondition<TDim, TNumNodes>(NewId, pGeometry) {}

// Constructor 2
template<unsigned int TDim, unsigned int TNumNodes>
TNormalFluxCondition<TDim, TNumNodes>::TNormalFluxCondition(IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties) : TCondition<TDim, TNumNodes>(NewId, pGeometry, pProperties) {}

// Destructor
template<unsigned int TDim, unsigned int TNumNodes>
TNormalFluxCondition<TDim, TNumNodes>::~TNormalFluxCondition() = default;

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer TNormalFluxCondition<TDim,TNumNodes>::Create(
    IndexType NewId,NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new TNormalFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim,TNumNodes>::CalculateRHS(
    VectorType& rRightHandSideVector,
    const ProcessInfo& CurrentProcessInfo)
{        
    //Previous definitions
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(this->GetIntegrationMethod());
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i < NumGPoints; ++i)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian( JContainer, this->GetIntegrationMethod());
    
    //Condition variables
    array_1d<double,TNumNodes> normal_flux_vector;
    for(unsigned int i = 0; i < TNumNodes; ++i)
    {
        normal_flux_vector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_HEAT_FLUX);
    }
    NormalFluxVariables Variables;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
        //Compute normal flux 
        Variables.NormalFlux = 0.0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            Variables.NormalFlux += NContainer(GPoint,i) * normal_flux_vector[i];
        }
        
        //Obtain Np
        noalias(Variables.Np) = row(NContainer,GPoint);
                
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient,
                                              JContainer[GPoint],
                                              IntegrationPoints[GPoint].Weight() );
                
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim,TNumNodes>::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    NormalFluxVariables& rVariables )
{
    noalias(rVariables.TVector) = rVariables.NormalFlux * rVariables.Np * rVariables.IntegrationCoefficient;

    GeoElementUtilities::
        AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.TVector);
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TNormalFluxCondition<TDim,TNumNodes>::CalculateIntegrationCoefficient(
    double& rIntegrationCoefficient,
    const Matrix& Jacobian,
    const double& Weight)
{
    Vector NormalVector = ZeroVector(TDim);

    if constexpr (TDim == 2)
    {
        NormalVector = column(Jacobian, 0);
    }
    else if constexpr (TDim == 3)
    {
        MathUtils<double>::CrossProduct(NormalVector, column(Jacobian, 0), column(Jacobian, 1));
    }
    rIntegrationCoefficient = Weight * MathUtils<double>::Norm(NormalVector);
}

// ============================================================================================
// ============================================================================================
template class TNormalFluxCondition<2,2>;
template class TNormalFluxCondition<2,3>;
template class TNormalFluxCondition<3,3>;
template class TNormalFluxCondition<3,4>;

} // Namespace Kratos.
