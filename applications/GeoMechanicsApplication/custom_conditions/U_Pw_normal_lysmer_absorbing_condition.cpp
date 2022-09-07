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

    //Previous definitions
    GeometryType& rGeom = this->GetGeometry();

    GeometryData::IntegrationMethod rIntegrationMethod = this->mThisIntegrationMethod;
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(rIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = rGeom.LocalSpaceDimension();

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(rIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    for (unsigned int i = 0; i < NumGPoints; ++i)
        (JContainer[i]).resize(TDim, LocalDim, false);
    rGeom.Jacobian(JContainer, rIntegrationMethod);

    //Condition variables
    BoundedMatrix<double, TDim, N_DOF> Nu = ZeroMatrix(TDim, N_DOF);

    NormalLysmerAbsorbingVariables rVariables;

    this->GetVariables(rVariables, CurrentProcessInfo);
    

    BoundedMatrix<double, TDim, N_DOF> AuxAbsKMatrix;
    BoundedMatrix<double, N_DOF, N_DOF> rAbsKMatrix = ZeroMatrix(N_DOF, N_DOF);

    //Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {

        // calculate
        rVariables.Ec = 0.0;
        rVariables.G = 0.0;
        for (unsigned int node = 0; node < rGeom.size(); ++node)
        {
            rVariables.Ec += NContainer(GPoint, node) * rVariables.EcNodes[node];
            rVariables.G += NContainer(GPoint, node) * rVariables.GNodes[node];
        }

        this->CalculateNodalStiffnessMatrix(rVariables, CurrentProcessInfo, rGeom);

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
    if (rLhsMatrix.size1() != CONDITION_SIZE)
        rLhsMatrix.resize(CONDITION_SIZE, CONDITION_SIZE, false);

    noalias(rLhsMatrix) = ZeroMatrix(CONDITION_SIZE, CONDITION_SIZE);
    this->CalculateAndAddLHS(rLhsMatrix, rVariables);

    // no righthand side contribution
    if (rRightHandSideVector.size() != CONDITION_SIZE)
        rRightHandSideVector.resize(CONDITION_SIZE, false);
    noalias(rRightHandSideVector) = ZeroVector(CONDITION_SIZE);
}

/// <summary>
/// Calculates LHS Damping part of abosrbing boundary
/// </summary>
/// <param name="rDampingMatrix"></param>
/// <param name="CurrentProcessInfo"></param>
template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& CurrentProcessInfo)
{

    //Previous definitions
    GeometryType& rGeom = this->GetGeometry();

    GeometryData::IntegrationMethod rIntegrationMethod = this->mThisIntegrationMethod;
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(rIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = rGeom.LocalSpaceDimension();

    //Containers of variables at all integration points
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(rIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    for (unsigned int i = 0; i < NumGPoints; ++i)
        (JContainer[i]).resize(TDim, LocalDim, false);
    rGeom.Jacobian(JContainer, rIntegrationMethod);

    //Condition variables
    BoundedMatrix<double, TDim, N_DOF> Nu = ZeroMatrix(TDim, N_DOF);

    NormalLysmerAbsorbingVariables rVariables;
    this->GetVariables(rVariables, CurrentProcessInfo);

    BoundedMatrix<double, TDim, N_DOF> AuxAbsMatrix;
    BoundedMatrix<double, N_DOF, N_DOF> rAbsMatrix = ZeroMatrix(N_DOF, N_DOF);

    //Loop over integration points
    for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {

        // calculate
        rVariables.rho = 0.0;
        rVariables.Ec = 0.0;
        rVariables.G = 0.0;
        for (unsigned int node = 0; node < rGeom.size(); ++node)
        {
            rVariables.rho += NContainer(GPoint, node) * rVariables.rhoNodes[node];
            rVariables.Ec += NContainer(GPoint, node) * rVariables.EcNodes[node];
            rVariables.G += NContainer(GPoint, node) * rVariables.GNodes[node];
        }
        rVariables.vp = sqrt(rVariables.Ec / rVariables.rho);
        rVariables.vs = sqrt(rVariables.G / rVariables.rho);

        this->CalculateNodalDampingMatrix(rVariables, CurrentProcessInfo, rGeom);

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


    // assemble left hand side vector
    if (rDampingMatrix.size1() != CONDITION_SIZE)
        rDampingMatrix.resize(CONDITION_SIZE, CONDITION_SIZE, false);
    noalias(rDampingMatrix) = ZeroMatrix(CONDITION_SIZE, CONDITION_SIZE);

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

    if (rValues.size() != CONDITION_SIZE)
        rValues.resize(CONDITION_SIZE, false);

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
CalculateNodalDampingMatrix(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, const Element::GeometryType& Geom)
{
    array_1d<double, 2> dampingConstants;

    // calculate rotation matrix
    BoundedMatrix<double, TDim, TDim> rotationMatrix;
    CalculateRotationMatrix(rotationMatrix, Geom);

    const int local_perpendicular_direction = TDim - 1;

    // calculate constant traction vector part
    dampingConstants[0] = rVariables.vs * rVariables.rho * rVariables.s_factor;
    dampingConstants[1] = rVariables.vp * rVariables.rho * rVariables.p_factor;

    BoundedMatrix<double, TDim, TDim> localCMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TDim> auxLocalCMatrix = ZeroMatrix(TDim, TDim);

    rVariables.CAbsMatrix = ZeroMatrix(TDim, TDim);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        localCMatrix(idim, idim) = dampingConstants[0];
    }
    localCMatrix(local_perpendicular_direction, local_perpendicular_direction) = dampingConstants[1];

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
CalculateNodalStiffnessMatrix(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, const Element::GeometryType& Geom)
{
    array_1d<double, 2> stiffnessConstants;

    // calculate rotation matrix
    BoundedMatrix<double, TDim, TDim> rotationMatrix;
    CalculateRotationMatrix(rotationMatrix, Geom);

    const int local_perpendicular_direction = TDim - 1;

    // calculate constant traction vector part
    stiffnessConstants[0] = rVariables.G / rVariables.virtual_thickness;
    stiffnessConstants[1] = rVariables.Ec / rVariables.virtual_thickness;

    BoundedMatrix<double, TDim, TDim> localKMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TDim> auxLocalKMatrix = ZeroMatrix(TDim, TDim);

    rVariables.KAbsMatrix = ZeroMatrix(TDim, TDim);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        localKMatrix(idim, idim) = stiffnessConstants[0];
    }
    localKMatrix(local_perpendicular_direction, local_perpendicular_direction) = stiffnessConstants[1];

    auxLocalKMatrix = prod(localKMatrix, rotationMatrix);
    rVariables.KAbsMatrix = prod(trans(rotationMatrix), auxLocalKMatrix);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        rVariables.KAbsMatrix(idim, idim) = abs(rVariables.KAbsMatrix(idim, idim));
    }
}

/// <summary>
/// Calculates the extrapolation matrix for neighbour elements. Values from integration points are extrapolated to the nodes
/// </summary>
/// <param name="NeighbourElement"></param>
/// <returns></returns>
template< unsigned int TDim, unsigned int TNumNodes >
Matrix UPwLysmerAbsorbingCondition<TDim, TNumNodes >::CalculateExtrapolationMatrixNeighbour(const Element& NeighbourElement)
{
    const GeometryData::IntegrationMethod rIntegrationMethodNeighbour = NeighbourElement.GetIntegrationMethod();
    const GeometryType& rNeighbourGeom = NeighbourElement.GetGeometry();
    const IndexType rNumNodesNeighbour = rNeighbourGeom.size();
    const IndexType NumGPointsNeighbour = rNeighbourGeom.IntegrationPointsNumber(rIntegrationMethodNeighbour);

    Matrix rExtrapolationMatrix = ZeroMatrix(rNumNodesNeighbour, NumGPointsNeighbour);

    // Calculate extrapolation matrix for 2d elements
    if (TDim == 2)
    {
        if (rNumNodesNeighbour == 3)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixTriangle(rExtrapolationMatrix, rIntegrationMethodNeighbour);
            return rExtrapolationMatrix;
        }
        if (rNumNodesNeighbour == 4)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixQuad(rExtrapolationMatrix, rIntegrationMethodNeighbour);
            return rExtrapolationMatrix;
        }
       
    }
    // Calculate extrapolation matrix for 3d elements
    if (TDim == 3)
    {
        if (rNumNodesNeighbour == 4)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixTetra(rExtrapolationMatrix, rIntegrationMethodNeighbour);
            return rExtrapolationMatrix;
        }
        if (rNumNodesNeighbour == 8)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixHexa(rExtrapolationMatrix, rIntegrationMethodNeighbour);
            return rExtrapolationMatrix;
        }

    }

    // if no extrapolation matrix is implemented, take average values at gauss points
    double averaging_factor = 1 / NumGPointsNeighbour;
    for (unsigned int node = 0; node < rNumNodesNeighbour; ++node)
    {
        for (unsigned int GPoint = 0; GPoint < NumGPointsNeighbour; ++GPoint)
        {
            rExtrapolationMatrix(node, GPoint) = averaging_factor;
        }
    }
    return rExtrapolationMatrix;
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
    
    // get neighbour elements
    auto neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    GeometryType& rGeom = this->GetGeometry();

    if (neighbours.size() == 0)
    {
        KRATOS_ERROR << "Condition: " <<this->Id()<<" does not have neighbour elements" <<  std::endl;
    }
 
    rVariables.EcNodes.resize(TNumNodes);
    rVariables.GNodes.resize(TNumNodes);
    rVariables.rhoNodes.resize(TNumNodes);

    // only get values from first neighbour
    Element& rNeighbour = neighbours[0];
    auto rPropNeighbour = rNeighbour.GetProperties();

    const GeometryType& rNeighbourGeom = rNeighbour.GetGeometry();
    const GeometryData::IntegrationMethod rIntegrationMethodNeighbour = rNeighbour.GetIntegrationMethod();
    const IndexType rNumNodesNeighbour = rNeighbourGeom.size();
    const IndexType NumGPointsNeighbour = rNeighbourGeom.IntegrationPointsNumber(rIntegrationMethodNeighbour);

    // get density and porosity from element
    std::vector<double> SaturationVector;
    vector<double> DensityVector(NumGPointsNeighbour);
    vector<double> ConfinedStiffness(NumGPointsNeighbour);
    vector<double> ShearStiffness(NumGPointsNeighbour);

    std::vector<double> stdConfinedStiffness(NumGPointsNeighbour);
    std::vector<double> stdShearStiffness(NumGPointsNeighbour);

    // get parameters at neighbour element integration points
    rNeighbour.CalculateOnIntegrationPoints(DEGREE_OF_SATURATION, SaturationVector, rCurrentProcessInfo);
    rNeighbour.CalculateOnIntegrationPoints(CONFINED_STIFFNESS, stdConfinedStiffness, rCurrentProcessInfo);
    rNeighbour.CalculateOnIntegrationPoints(SHEAR_STIFFNESS, stdShearStiffness, rCurrentProcessInfo);

    for (unsigned int GPoint = 0; GPoint < NumGPointsNeighbour; ++GPoint) {

        // transform std vectors to boost vectors
        ConfinedStiffness[GPoint] = stdConfinedStiffness[GPoint];
        ShearStiffness[GPoint] = stdShearStiffness[GPoint];

        // calculate density mixture
        DensityVector[GPoint] = (SaturationVector[GPoint] * rPropNeighbour[POROSITY] * rPropNeighbour[DENSITY_WATER]) + (1.0 - rPropNeighbour[POROSITY]) * rPropNeighbour[DENSITY_SOLID];
    }

    Matrix rExtrapolationMatrix = CalculateExtrapolationMatrixNeighbour(rNeighbour);
   
    // project parameters on neighbour nodes
    vector<double> EcNodesNeighbour = prod(rExtrapolationMatrix, ConfinedStiffness);
    vector<double> GNodesNeighbour = prod(rExtrapolationMatrix, ShearStiffness);
    vector<double> rhoNodesNeighbour = prod(rExtrapolationMatrix, DensityVector);

    // add parameters to condition nodes
    for (unsigned int k = 0; k < rNumNodesNeighbour; ++k)
    {
        for (unsigned int node = 0; node < TNumNodes; ++node)
        {
            // add parameter if neighbour node id and condition node id are the same
            if (rNeighbourGeom[k].Id() == rGeom[node].Id())
            {
                rVariables.EcNodes[node] = EcNodesNeighbour(k);
                rVariables.GNodes[node] = GNodesNeighbour(k);
                rVariables.rhoNodes[node] = rhoNodesNeighbour(k);
            }
        }
    }   
}

/// <summary>
/// Gets condition variables
/// </summary>
/// <param name="rVariables"></param>
/// <param name="rCurrentProcessInfo"></param>
template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
GetVariables(
    NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    // gets variables from neighbour elements
    this->GetNeighbourElementVariables(rVariables, rCurrentProcessInfo);

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
void UPwLysmerAbsorbingCondition<2, 2>::CalculateRotationMatrix( BoundedMatrix<double, 2, 2>& rRotationMatrix, const Element::GeometryType& rGeom)
{
    //Line_2d_2
    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = rGeom.GetPoint(1) - rGeom.GetPoint(0);
    const double norm_x = norm_2(Vx);

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
void UPwLysmerAbsorbingCondition<2, 3>::CalculateRotationMatrix(BoundedMatrix<double, 2, 2>& rRotationMatrix, const Element::GeometryType& rGeom)
{
    //Line_2d_3
    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = rGeom.GetPoint(1) - rGeom.GetPoint(0);
    const double norm_x = norm_2(Vx);

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
void UPwLysmerAbsorbingCondition<3, 3>::CalculateRotationMatrix(BoundedMatrix<double, 3, 3>& rRotationMatrix, const Element::GeometryType& rGeom)
{

    ////triangle_3d_3
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    noalias(pmid0) = 0.5 * (rGeom.GetPoint(0) + rGeom.GetPoint(1));
    noalias(pmid1) = 0.5 * (rGeom.GetPoint(0) + rGeom.GetPoint(2));

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = rGeom.GetPoint(1) - rGeom.GetPoint(0);
    const double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;
    Vx[2] *= inv_norm_x;

    //Unitary vector in local z direction
    array_1d<double, 3> Vy;
    noalias(Vy) = rGeom.GetPoint(2) - rGeom.GetPoint(0);

    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    const double norm_z = norm_2(Vz);

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
void UPwLysmerAbsorbingCondition<3, 4>::CalculateRotationMatrix( BoundedMatrix<double, 3, 3>& rRotationMatrix, const Element::GeometryType& rGeom)
{
    //Quadrilateral_3d_4
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    const array_1d<double, 3> P2 = array_1d<double, 3>(rGeom.GetPoint(2));
    noalias(pmid0) = 0.5 * (rGeom.GetPoint(0) + rGeom.GetPoint(3));
    noalias(pmid1) = 0.5 * (rGeom.GetPoint(1) + P2);

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = pmid1 - pmid0;
    const double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;
    Vx[2] *= inv_norm_x;

    //Unitary vector in local z direction
    array_1d<double, 3> Vy;
    noalias(Vy) = P2 - pmid0;
    array_1d<double, 3> Vz;
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    const double norm_z = norm_2(Vz);

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
//

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwLysmerAbsorbingCondition<2,2>;
template class UPwLysmerAbsorbingCondition<2,3>;
template class UPwLysmerAbsorbingCondition<3,3>;
template class UPwLysmerAbsorbingCondition<3,4>;

} // Namespace Kratos.
