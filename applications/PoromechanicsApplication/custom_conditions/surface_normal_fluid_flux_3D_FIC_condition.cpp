//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_conditions/surface_normal_fluid_flux_3D_FIC_condition.hpp"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
SurfaceNormalFluidFlux3DFICCondition::SurfaceNormalFluidFlux3DFICCondition() : SurfaceNormalFluidFlux3DCondition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SurfaceNormalFluidFlux3DFICCondition::SurfaceNormalFluidFlux3DFICCondition(IndexType NewId, GeometryType::Pointer pGeometry) : SurfaceNormalFluidFlux3DCondition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SurfaceNormalFluidFlux3DFICCondition::SurfaceNormalFluidFlux3DFICCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : SurfaceNormalFluidFlux3DCondition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
SurfaceNormalFluidFlux3DFICCondition::~SurfaceNormalFluidFlux3DFICCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer SurfaceNormalFluidFlux3DFICCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceNormalFluidFlux3DFICCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SurfaceNormalFluidFlux3DFICCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables)
{
    this->CalculateAndAddBoundaryMassMatrix(rLeftHandSideMatrix, rVariables);
}

//----------------------------------------------------------------------------------------

void SurfaceNormalFluidFlux3DFICCondition::CalculateAndAddBoundaryMassMatrix(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables)
{
    Matrix BoundaryMassMatrix = rVariables.NewmarkCoefficient*rVariables.ElementLength*rVariables.BiotModulusInverse/6*outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;
    
    //Distribute boundary mass block matrix into the elemental matrix
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int Global_i, Global_j;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;
        
        for(unsigned int j = 0; j < number_of_nodes; j++)
        {
            Global_j = j * (dimension + 1) + dimension;
            
            rLeftHandSideMatrix(Global_i,Global_j) -= BoundaryMassMatrix(i,j);
        }
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void SurfaceNormalFluidFlux3DFICCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    this->CalculateAndAddConditionForce(rRightHandSideVector, rVariables);
    
    this->CalculateAndAddBoundaryMassFlow(rRightHandSideVector, rVariables);
}

//----------------------------------------------------------------------------------------

void SurfaceNormalFluidFlux3DFICCondition::CalculateAndAddBoundaryMassFlow(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{    
    Matrix BoundaryMassMatrix = rVariables.ElementLength*rVariables.BiotModulusInverse/6.0*outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    Vector PressureDtVector = ZeroVector(number_of_nodes);
    
    for(unsigned int i=0; i<number_of_nodes; i++)
    {
        PressureDtVector[i] = rGeom[i].FastGetSolutionStepValue(DERIVATIVE_WATER_PRESSURE);
    }

    Vector BoundaryMassFlow = prod(BoundaryMassMatrix,PressureDtVector);

    //Distribute boundary mass block vector into the elemental vector
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int Global_i;
    
    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;
        
        rRightHandSideVector[Global_i] += BoundaryMassFlow[i];
    }
}

} // Namespace Kratos.
