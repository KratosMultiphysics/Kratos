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
//  Main authors:    Aron Noordam
//

#include <math.h>
// Application includes
#include "custom_conditions/U_Pw_normal_lysmer_absorbing_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwLysmerAbsorbingCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwLysmerAbsorbingCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/// <summary>
/// Calculates LHS stiffness part of absorbing boundary
/// </summary>
/// <param name="rLhsMatrix"></param>
/// <param name="rRightHandSideVector"></param>
/// <param name="CurrentProcessInfo"></param>
template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLhsMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo)
{
    int n_dof = TNumNodes * TDim;
    BoundedMatrix<double, TNumNodes* TDim, TNumNodes* TDim> rAbsKMatrix = ZeroMatrix(n_dof, n_dof);

    //Previous definitions
    GeometryType& Geom = this->GetGeometry();
    PropertiesType& prop = this->GetProperties();

    //GeometryData::IntegrationMethod rIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
    GeometryData::IntegrationMethod rIntegrationMethod = this->mThisIntegrationMethod;
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(rIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(rIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    for (unsigned int i = 0; i < NumGPoints; ++i)
        (JContainer[i]).resize(TDim, LocalDim, false);
    Geom.Jacobian(JContainer, rIntegrationMethod);

    //Condition variables
    BoundedMatrix<double, TDim, TNumNodes* TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);

    NormalLysmerAbsorbingVariables rVariables;

    this->GetVariables(rVariables, CurrentProcessInfo);
    this->CalculateNodalStiffnessMatrix(rVariables, CurrentProcessInfo, Geom);

    BoundedMatrix<double, TDim, TDim>             KMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TNumNodes* TDim> AuxAbsKMatrix = ZeroMatrix(TDim, TNumNodes * TDim);

    //Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {

        // calculate displacement shape function matrix
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nu, NContainer, GPoint);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(rVariables.IntegrationCoefficient,
            JContainer[GPoint],
            IntegrationPoints[GPoint].Weight());


        // set stiffness part of absorbing matrix
        AuxAbsKMatrix = prod(rVariables.KAbsMatrix, Nu);
        rAbsKMatrix += prod(trans(Nu), AuxAbsKMatrix) * rVariables.IntegrationCoefficient;

    }
    rVariables.UMatrix = rAbsKMatrix;

    // assemble left hand side vector
    const SizeType ConditionSize = TNumNodes * TDim + TNumNodes;
    noalias(rLhsMatrix) = ZeroMatrix(ConditionSize);
    if (rLhsMatrix.size1() != ConditionSize)
        rLhsMatrix.resize(ConditionSize, ConditionSize, false);

    noalias(rLhsMatrix) = ZeroMatrix(ConditionSize, ConditionSize);
    this->CalculateAndAddLHS(rLhsMatrix, rVariables);

    // no righthand side contribution
    rRightHandSideVector = ZeroVector(ConditionSize);
}

/// <summary>
/// Calculates LHS Damping part of abosrbing boundary
/// </summary>
/// <param name="rDampingMatrix"></param>
/// <param name="CurrentProcessInfo"></param>
template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& CurrentProcessInfo)
{
    int n_dof = TNumNodes * TDim;
    BoundedMatrix<double, TNumNodes* TDim, TNumNodes* TDim> rAbsMatrix = ZeroMatrix(n_dof, n_dof);

    //Previous definitions
    GeometryType& Geom = this->GetGeometry();
    PropertiesType& prop = this->GetProperties();

    GeometryData::IntegrationMethod rIntegrationMethod = this->mThisIntegrationMethod;
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(rIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();

    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(rIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    for (unsigned int i = 0; i < NumGPoints; ++i)
        (JContainer[i]).resize(TDim, LocalDim, false);
    Geom.Jacobian(JContainer, rIntegrationMethod);

    //Condition variables
    BoundedMatrix<double, TDim, TNumNodes* TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);

    NormalLysmerAbsorbingVariables rVariables;
    this->GetVariables(rVariables, CurrentProcessInfo);
    this->CalculateNodalDampingMatrix(rVariables, CurrentProcessInfo, Geom);

    BoundedMatrix<double, TDim, TNumNodes* TDim> AuxAbsMatrix = ZeroMatrix(TDim, TNumNodes * TDim);

    //Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {

        // calculate displacement shape function matrix
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nu, NContainer, GPoint);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(rVariables.IntegrationCoefficient,
            JContainer[GPoint],
            IntegrationPoints[GPoint].Weight());

        // set damping part of absorbing matrix
        AuxAbsMatrix = prod(rVariables.CAbsMatrix, Nu);
        rAbsMatrix += prod(trans(Nu), AuxAbsMatrix) * rVariables.IntegrationCoefficient;
    }
    rVariables.UMatrix = rAbsMatrix;


    // assemble right hand side vector
    const SizeType ConditionSize = TNumNodes * TDim + TNumNodes;
    noalias(rDampingMatrix) = ZeroMatrix(ConditionSize);
    if (rDampingMatrix.size1() != ConditionSize)
        rDampingMatrix.resize(ConditionSize, ConditionSize, false);
    rDampingMatrix = ZeroMatrix(ConditionSize, ConditionSize);

    this->CalculateAndAddLHS(rDampingMatrix, rVariables);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

/// <summary>
/// Gets the velocity vector of the absorbing boundary
/// </summary>
/// <param name="rValues"></param>
/// <param name="Step"></param>
template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int N_DOF = TNumNodes * TDim + TNumNodes;

    if (rValues.size() != N_DOF)
        rValues.resize(N_DOF, false);

    if (TDim == 2) {
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rValues[index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
            rValues[index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
            rValues[index++] = 0.0;
        }
    }
    else if (TDim == 3) {
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rValues[index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
            rValues[index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
            rValues[index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z, Step);
            rValues[index++] = 0.0;
        }
    }
    else {
        KRATOS_ERROR << "undefined dimension in GetFirstDerivativesVector... illegal operation!!" << this->Id() << std::endl;
    }

    KRATOS_CATCH("")
}

/// <summary>
/// Calculates the damping constant in all directions 
/// </summary>
/// <param name="rVariables"></param>
/// <param name="CurrentProcessInfo"></param>
/// <param name="Geom"></param>
template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
CalculateNodalDampingMatrix(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, Element::GeometryType& Geom)
{
    array_1d<double, 2> rDampingConstants;

    // calculate rotation matrix
    BoundedMatrix<double, TDim, TDim> rotationMatrix;
    CalculateRotationMatrix(rotationMatrix, Geom);

    const int local_perpendicular_direction = TDim - 1;


    // calculate constant traction vector part
    rDampingConstants[0] = rVariables.vs * rVariables.rho * rVariables.s_factor;
    rDampingConstants[1] = rVariables.vp * rVariables.rho * rVariables.p_factor;

    BoundedMatrix<double, TDim, TDim>             localCMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TDim>             auxLocalCMatrix = ZeroMatrix(TDim, TDim);

    rVariables.CAbsMatrix = ZeroMatrix(TDim, TDim);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        localCMatrix(idim, idim) = rDampingConstants[0];
    }
    //localCMatrix(local_perpendicular_direction, local_perpendicular_direction) = 0;
    localCMatrix(local_perpendicular_direction, local_perpendicular_direction) = rDampingConstants[1];

    auxLocalCMatrix = prod(localCMatrix, rotationMatrix);
    rVariables.CAbsMatrix = prod(trans(rotationMatrix), auxLocalCMatrix);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
            rVariables.CAbsMatrix(idim, idim) = abs(rVariables.CAbsMatrix(idim, idim));
    }
}

/// <summary>
/// Calculates the stiffness constants in all directions
/// </summary>
/// <param name="rVariables"></param>
/// <param name="CurrentProcessInfo"></param>
/// <param name="Geom"></param>
template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
CalculateNodalStiffnessMatrix(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, Element::GeometryType& Geom)
{
    array_1d<double, 2> rStiffnessConstants;

    // calculate rotation matrix
    BoundedMatrix<double, TDim, TDim> rotationMatrix;
    CalculateRotationMatrix(rotationMatrix, Geom);

    const int local_perpendicular_direction = TDim - 1;


    // calculate constant traction vector part
    rStiffnessConstants[0] = rVariables.G / rVariables.virtual_thickness;
    rStiffnessConstants[1] = rVariables.Ec / rVariables.virtual_thickness;


    BoundedMatrix<double, TDim, TDim>             localKMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TDim>             auxLocalKMatrix = ZeroMatrix(TDim, TDim);

    rVariables.KAbsMatrix = ZeroMatrix(TDim, TDim);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        localKMatrix(idim, idim) = rStiffnessConstants[0];
    }
    localKMatrix(local_perpendicular_direction, local_perpendicular_direction) = rStiffnessConstants[1];

    auxLocalKMatrix = prod(localKMatrix, rotationMatrix);
    rVariables.KAbsMatrix = prod(trans(rotationMatrix), auxLocalKMatrix);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        rVariables.KAbsMatrix(idim, idim) = abs(rVariables.KAbsMatrix(idim, idim));
    }
}

/// <summary>
/// This method gets the average of the variables of all the neighbour elements of the condition. 
/// </summary>
/// <param name="rVariables"></param>
/// <param name="rCurrentProcessInfo"></param>
template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
GetNeighbourElementVariables(
    NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    
    double rho_s = 0;
    double rho_w = 0;
    double E = 0;
    double n = 0;
    double nu = 0;
    double rMeanDegreeOfSaturation = 0;

    // get mean degree of saturation of all integration points in all neighbour elements
    std::vector<double> SaturationVector;
    
    // get neighbour elements
    auto neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    
    int nElements = neighbours.size();
    
    // loop over all neighbour elements
    int nValues = 0;
    for (unsigned int i = 0; i < neighbours.size(); ++i) {
        
        auto const rPropNeighbour = neighbours[i].GetProperties();

        rho_s = rho_s + rPropNeighbour[DENSITY_SOLID];
        rho_w = rho_w + rPropNeighbour[DENSITY_WATER];
        n = n + rPropNeighbour[POROSITY];
        E = E + rPropNeighbour[YOUNG_MODULUS];
        nu = nu + rPropNeighbour[POISSON_RATIO];

        neighbours[i].CalculateOnIntegrationPoints(DEGREE_OF_SATURATION, SaturationVector, rCurrentProcessInfo);

        // loop over all integration points
        for (unsigned int j = 0; j < SaturationVector.size(); ++j)
        {
            rMeanDegreeOfSaturation = rMeanDegreeOfSaturation + SaturationVector[j];
            nValues = nValues + 1;
        }
    }

    // calculate mean of neighbour element variables
    rVariables.E = E / nElements;
    rVariables.nu = nu / nElements;

    n = n / nElements;
    rho_s = rho_s / nElements;
    rho_w = rho_w / nElements;
    rMeanDegreeOfSaturation = rMeanDegreeOfSaturation / nValues;

    // calculate density of mixture
    rVariables.rho = (rMeanDegreeOfSaturation * n * rho_w) + (1.0 - n) * rho_s;
}


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
GetVariables(
    NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    // gets average of variables as stored in neighbour elements
    this->GetNeighbourElementVariables(rVariables, rCurrentProcessInfo);

    // calculate P wave modulus and shear wave modulus
    rVariables.Ec = rVariables.E * (1 - rVariables.nu) / ((1 + rVariables.nu) * (1 - 2 * rVariables.nu));
    rVariables.G = rVariables.E / (2 * (1 + rVariables.nu));

    // calculate pressure wave and shear wave velocity
    rVariables.vp = sqrt(rVariables.Ec / rVariables.rho);
    rVariables.vs = sqrt(rVariables.G / rVariables.rho);

    // get condition specific variables
    Vector absorbing_factors = this->GetValue(ABSORBING_FACTORS);
    rVariables.p_factor = absorbing_factors(0);
    rVariables.s_factor = absorbing_factors(1);
    rVariables.virtual_thickness = this->GetValue(VIRTUAL_THICKNESS);
}



template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
    NormalLysmerAbsorbingVariables& rVariables)
{
    //Adding contribution to left hand side

    GeoElementUtilities::
        AssembleUBlockMatrix< TDim, TNumNodes >(rLeftHandSideMatrix,
            rVariables.UMatrix);
}

template< >
void UPwLysmerAbsorbingCondition<2, 2>::CalculateRotationMatrix( BoundedMatrix<double, 2, 2>& rRotationMatrix, const Element::GeometryType& Geom)
{
    //Line_2d_2
    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = Geom.GetPoint(1) - Geom.GetPoint(0);
    double norm_x = norm_2(Vx);

    Vx[0] *= 1.0 / norm_x;
    Vx[1] *= 1.0 / norm_x;

    //Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];

    // We need to determine the unitary vector in local y direction pointing towards the TOP face of the joint

    // Unitary vector in local x direction (3D)
    array_1d<double, 3> Vx3D;
    Vx3D[0] = Vx[0];
    Vx3D[1] = Vx[1];
    Vx3D[2] = 0.0;

    // Unitary vector in local y direction (first option)
    array_1d<double, 3> Vy3D;
    Vy3D[0] = -Vx[1];
    Vy3D[1] = Vx[0];
    Vy3D[2] = 0.0;

    // Vector in global z direction (first option)
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx3D, Vy3D);

    // Vz must have the same sign as vector (0,0,1)
    if (Vz[2] > 0.0)
    {
        rRotationMatrix(1, 0) = -Vx[1];
        rRotationMatrix(1, 1) = Vx[0];
    }
    else
    {
        rRotationMatrix(1, 0) = Vx[1];
        rRotationMatrix(1, 1) = -Vx[0];
    }

}


template< >
void UPwLysmerAbsorbingCondition<2, 3>::CalculateRotationMatrix(BoundedMatrix<double, 2, 2>& rRotationMatrix, const Element::GeometryType& Geom)
{
    //Line_2d_3
    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = Geom.GetPoint(1) - Geom.GetPoint(0);
    double norm_x = norm_2(Vx);

    Vx[0] *= 1.0 / norm_x;
    Vx[1] *= 1.0 / norm_x;

    //Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];

    // We need to determine the unitary vector in local y direction pointing towards the TOP face of the joint

    // Unitary vector in local x direction (3D)
    array_1d<double, 3> Vx3D;
    Vx3D[0] = Vx[0];
    Vx3D[1] = Vx[1];
    Vx3D[2] = 0.0;

    // Unitary vector in local y direction (first option)
    array_1d<double, 3> Vy3D;
    Vy3D[0] = -Vx[1];
    Vy3D[1] = Vx[0];
    Vy3D[2] = 0.0;

    // Vector in global z direction (first option)
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx3D, Vy3D);

    // Vz must have the same sign as vector (0,0,1)
    if (Vz[2] > 0.0)
    {
        rRotationMatrix(1, 0) = -Vx[1];
        rRotationMatrix(1, 1) = Vx[0];
    }
    else
    {
        rRotationMatrix(1, 0) = Vx[1];
        rRotationMatrix(1, 1) = -Vx[0];
    }
}

template< >
void UPwLysmerAbsorbingCondition<3, 3>::CalculateRotationMatrix(BoundedMatrix<double, 3, 3>& rRotationMatrix, const Element::GeometryType& Geom)
{

    ////triangle_3d_3
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    noalias(pmid0) = 0.5 * (Geom.GetPoint(0) + Geom.GetPoint(1));
    noalias(pmid1) = 0.5 * (Geom.GetPoint(0) + Geom.GetPoint(2));

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = Geom.GetPoint(1) - Geom.GetPoint(0);
    double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;
    Vx[2] *= inv_norm_x;

    //Unitary vector in local z direction
    array_1d<double, 3> Vy;
    noalias(Vy) = Geom.GetPoint(2) - Geom.GetPoint(0);

    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    double norm_z = norm_2(Vz);

    Vz[0] *= 1.0 / norm_z;
    Vz[1] *= 1.0 / norm_z;
    Vz[2] *= 1.0 / norm_z;

    //Unitary vector in local y direction
    MathUtils<double>::CrossProduct(Vy, Vz, Vx);

    //Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];
    rRotationMatrix(0, 2) = Vx[2];

    rRotationMatrix(1, 0) = Vy[0];
    rRotationMatrix(1, 1) = Vy[1];
    rRotationMatrix(1, 2) = Vy[2];

    rRotationMatrix(2, 0) = Vz[0];
    rRotationMatrix(2, 1) = Vz[1];
    rRotationMatrix(2, 2) = Vz[2];
}
//----------------------------------------------------------------------------------------

template< >
void UPwLysmerAbsorbingCondition<3, 4>::CalculateRotationMatrix( BoundedMatrix<double, 3, 3>& rRotationMatrix, const Element::GeometryType& Geom)
{
    //Quadrilateral_3d_4
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    array_1d<double, 3> P2 = Geom.GetPoint(2);
    noalias(pmid0) = 0.5 * (Geom.GetPoint(0) + Geom.GetPoint(3));
    noalias(pmid1) = 0.5 * (Geom.GetPoint(1) + P2);

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = pmid1 - pmid0;
    double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;
    Vx[2] *= inv_norm_x;

    //Unitary vector in local z direction
    array_1d<double, 3> Vy;
    noalias(Vy) = P2 - pmid0;
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    double norm_z = norm_2(Vz);

    Vz[0] *= 1.0 / norm_z;
    Vz[1] *= 1.0 / norm_z;
    Vz[2] *= 1.0 / norm_z;

    //Unitary vector in local y direction
    MathUtils<double>::CrossProduct(Vy, Vz, Vx);

    //Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];
    rRotationMatrix(0, 2) = Vx[2];

    rRotationMatrix(1, 0) = Vy[0];
    rRotationMatrix(1, 1) = Vy[1];
    rRotationMatrix(1, 2) = Vy[2];

    rRotationMatrix(2, 0) = Vz[0];
    rRotationMatrix(2, 1) = Vz[1];
    rRotationMatrix(2, 2) = Vz[2];
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwLysmerAbsorbingCondition<2,2>;
template class UPwLysmerAbsorbingCondition<2,3>;
template class UPwLysmerAbsorbingCondition<3,3>;
template class UPwLysmerAbsorbingCondition<3,4>;

} // Namespace Kratos.
