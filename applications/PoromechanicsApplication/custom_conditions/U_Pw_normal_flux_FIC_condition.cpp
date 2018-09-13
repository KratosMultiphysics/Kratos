//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


// Application includes
#include "custom_conditions/U_Pw_normal_flux_FIC_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwNormalFluxFICCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwNormalFluxFICCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFluxFICCondition<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{    
    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i<NumGPoints; i++)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    
    //Condition variables
    array_1d<double,TNumNodes> NormalFluxVector;
    NormalFluxVariables Variables;
    NormalFluxFICVariables FICVariables;
    FICVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];
    this->CalculateElementLength(FICVariables.ElementLength,Geom);
    const double& BulkModulusSolid = Prop[BULK_MODULUS_SOLID];
    const double& Porosity = Prop[POROSITY];
    const double BulkModulus = Prop[YOUNG_MODULUS]/(3.0*(1.0-2.0*Prop[POISSON_RATIO]));
    const double BiotCoefficient = 1.0-BulkModulus/BulkModulusSolid;
    FICVariables.BiotModulusInverse = (BiotCoefficient-Porosity)/BulkModulusSolid + Porosity/Prop[BULK_MODULUS_FLUID];
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        NormalFluxVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
        FICVariables.DtPressureVector[i] = Geom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute normal flux
        Variables.NormalFlux = 0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            Variables.NormalFlux += NContainer(GPoint,i)*NormalFluxVector[i];
        }

        //Obtain Np
        noalias(Variables.Np) = row(NContainer,GPoint);
                
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, JContainer[GPoint], integration_points[GPoint].Weight() );

        //Contributions to the left hand side        
        this->CalculateAndAddLHSStabilization(rLeftHandSideMatrix, Variables, FICVariables);
        
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        
        this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFluxFICCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{        
    //Previous definitions
    const PropertiesType& Prop = this->GetProperties();
    const GeometryType& Geom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = integration_points.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i<NumGPoints; i++)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    
    //Condition variables
    array_1d<double,TNumNodes> NormalFluxVector;
    NormalFluxVariables Variables;
    NormalFluxFICVariables FICVariables;
    FICVariables.DtPressureCoefficient = CurrentProcessInfo[DT_PRESSURE_COEFFICIENT];
    this->CalculateElementLength(FICVariables.ElementLength,Geom);
    const double& BulkModulusSolid = Prop[BULK_MODULUS_SOLID];
    const double& Porosity = Prop[POROSITY];
    const double BulkModulus = Prop[YOUNG_MODULUS]/(3.0*(1.0-2.0*Prop[POISSON_RATIO]));
    const double BiotCoefficient = 1.0-BulkModulus/BulkModulusSolid;
    FICVariables.BiotModulusInverse = (BiotCoefficient-Porosity)/BulkModulusSolid + Porosity/Prop[BULK_MODULUS_FLUID];
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        NormalFluxVector[i] = Geom[i].FastGetSolutionStepValue(NORMAL_FLUID_FLUX);
        FICVariables.DtPressureVector[i] = Geom[i].FastGetSolutionStepValue(DT_WATER_PRESSURE);
    }
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute normal flux
        Variables.NormalFlux = 0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            Variables.NormalFlux += NContainer(GPoint,i)*NormalFluxVector[i];
        }
        
        //Obtain Np
        noalias(Variables.Np) = row(NContainer,GPoint);
                
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(Variables.IntegrationCoefficient, JContainer[GPoint], integration_points[GPoint].Weight() );
                
        //Contributions to the right hand side
        this->CalculateAndAddRHS(rRightHandSideVector, Variables);
        
        this->CalculateAndAddRHSStabilization(rRightHandSideVector, Variables, FICVariables);
    }
}

//----------------------------------------------------------------------------------------

template<>
void UPwNormalFluxFICCondition<2,2>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = Geom.Length();
}

//----------------------------------------------------------------------------------------

template<>
void UPwNormalFluxFICCondition<3,3>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = sqrt(4.0*Geom.Area()/Globals::Pi);
}

//----------------------------------------------------------------------------------------

template<>
void UPwNormalFluxFICCondition<3,4>::CalculateElementLength(double& rElementLength, const GeometryType& Geom)
{
    rElementLength = sqrt(4.0*Geom.Area()/Globals::Pi);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFluxFICCondition<TDim,TNumNodes>::CalculateAndAddLHSStabilization(MatrixType& rLeftHandSideMatrix, NormalFluxVariables& rVariables, 
                                                                                    NormalFluxFICVariables& rFICVariables)
{
    this->CalculateAndAddBoundaryMassMatrix(rLeftHandSideMatrix, rVariables, rFICVariables);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFluxFICCondition<TDim,TNumNodes>::CalculateAndAddBoundaryMassMatrix(MatrixType& rLeftHandSideMatrix, NormalFluxVariables& rVariables, 
                                                                                    NormalFluxFICVariables& rFICVariables)
{
    noalias(rFICVariables.PMatrix) = -rFICVariables.DtPressureCoefficient*rFICVariables.ElementLength*rFICVariables.BiotModulusInverse/6.0*
                                        outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;
    
    //Distribute boundary mass matrix into the elemental matrix
    PoroElementUtilities::AssemblePBlockMatrix< BoundedMatrix<double,TNumNodes,TNumNodes> >(rLeftHandSideMatrix,rFICVariables.PMatrix,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFluxFICCondition<TDim,TNumNodes>::CalculateAndAddRHSStabilization(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables, 
                                                                                    NormalFluxFICVariables& rFICVariables)
{
    this->CalculateAndAddBoundaryMassFlow(rRightHandSideVector, rVariables, rFICVariables);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwNormalFluxFICCondition<TDim,TNumNodes>::CalculateAndAddBoundaryMassFlow(VectorType& rRightHandSideVector, NormalFluxVariables& rVariables, 
                                                                                    NormalFluxFICVariables& rFICVariables)
{
    noalias(rFICVariables.PMatrix) = rFICVariables.ElementLength*rFICVariables.BiotModulusInverse/6.0*
                                        outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;
    

    noalias(rVariables.PVector) = prod(rFICVariables.PMatrix,rFICVariables.DtPressureVector);
    
    //Distribute boundary mass flow vector into elemental vector
    PoroElementUtilities::AssemblePBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwNormalFluxFICCondition<2,2>;
template class UPwNormalFluxFICCondition<3,3>;
template class UPwNormalFluxFICCondition<3,4>;

} // Namespace Kratos.
