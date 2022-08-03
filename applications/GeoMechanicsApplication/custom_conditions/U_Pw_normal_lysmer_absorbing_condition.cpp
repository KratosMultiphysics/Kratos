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
        for (unsigned int node = 0; node < Geom.size(); ++node)
        {
            rVariables.Ec += NContainer(GPoint, node) * rVariables.EcNodes[node];
            rVariables.G += NContainer(GPoint, node) * rVariables.GNodes[node];
        }


        this->CalculateNodalStiffnessMatrix(rVariables, CurrentProcessInfo, Geom);

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
    noalias(rLhsMatrix) = ZeroMatrix(CONDITION_SIZE);
    if (rLhsMatrix.size1() != CONDITION_SIZE)
        rLhsMatrix.resize(CONDITION_SIZE, CONDITION_SIZE, false);

    noalias(rLhsMatrix) = ZeroMatrix(CONDITION_SIZE, CONDITION_SIZE);
    this->CalculateAndAddLHS(rLhsMatrix, rVariables);

    // no righthand side contribution
    rRightHandSideVector = ZeroVector(CONDITION_SIZE);
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
    GeometryType& Geom = this->GetGeometry();
    PropertiesType& prop = this->GetProperties();

    //GeometryData::IntegrationMethod rIntegrationMethod = this->mThisIntegrationMethod;
    GeometryData::IntegrationMethod rIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
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
        for (unsigned int node = 0; node < Geom.size(); ++node)
        {
            rVariables.rho += NContainer(GPoint, node) * rVariables.rhoNodes[node];
            rVariables.Ec += NContainer(GPoint, node) * rVariables.EcNodes[node];
            rVariables.G += NContainer(GPoint, node) * rVariables.GNodes[node];
        }
        rVariables.vp = sqrt(rVariables.Ec / rVariables.rho);
        rVariables.vs = sqrt(rVariables.G / rVariables.rho);


        this->CalculateNodalDampingMatrix(rVariables, CurrentProcessInfo, Geom);

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
    noalias(rDampingMatrix) = ZeroMatrix(CONDITION_SIZE);
    if (rDampingMatrix.size1() != CONDITION_SIZE)
        rDampingMatrix.resize(CONDITION_SIZE, CONDITION_SIZE, false);
    rDampingMatrix = ZeroMatrix(CONDITION_SIZE, CONDITION_SIZE);

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
CalculateNodalStiffnessMatrix(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, const Element::GeometryType& Geom)
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
template< unsigned int TDim, unsigned int TNumNodes >
Matrix UPwLysmerAbsorbingCondition<TDim, TNumNodes >::CalculateExtrapolationMatrixNeighbour(const Element& NeighbourElement)
{
    const GeometryData::IntegrationMethod rIntegrationMethodNeighbour = NeighbourElement.GetIntegrationMethod();
    const GeometryType& rNeighbourGeom = NeighbourElement.GetGeometry();
    const IndexType rNumNodesNeighbour = rNeighbourGeom.size();
    const IndexType NumGPointsNeighbour = rNeighbourGeom.IntegrationPointsNumber(rIntegrationMethodNeighbour);

    Matrix rExtrapolationMatrix = ZeroMatrix(rNumNodesNeighbour, NumGPointsNeighbour);

    if (TDim == 2)
    {
        if (rNumNodesNeighbour == 3)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixTriangle(rExtrapolationMatrix, rIntegrationMethodNeighbour);
        }
        else if (rNumNodesNeighbour == 4)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixQuad(rExtrapolationMatrix, rIntegrationMethodNeighbour);
        }
    }
    if (TDim == 3)
    {
        if (rNumNodesNeighbour == 4)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixTetra(rExtrapolationMatrix, rIntegrationMethodNeighbour);
        }
        else if (rNumNodesNeighbour == 8)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixHexa(rExtrapolationMatrix, rIntegrationMethodNeighbour);
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
    
    //double rho_s = 0;
    //double rho_w = 0;
    double Ec = 0;
    double G = 0;
    //double n = 0;
    double rMeanDegreeOfSaturation = 0;

    // get mean degree of saturation of all integration points in all neighbour elements
    
    
    // get neighbour elements
    auto neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    
    int nElements = neighbours.size();
    
    GeometryType& rGeom = this->GetGeometry();
    GeometryData::IntegrationMethod rIntegrationMethod = this->GetIntegrationMethod();
    const Matrix& NContainer = rGeom.ShapeFunctionsValues(rIntegrationMethod);

    std::vector<double> rEcNodes(TNumNodes, 0);
    //std::vector<double> rGNodes;
    //rEcNodes.resize(rNumNodesNeighbour);
    std::vector<double> rGNodes(TNumNodes, 0);
    std::vector<double> rSaturationNodes(TNumNodes, 0);
    std::vector<double> rRhoNodes(TNumNodes, 0);

    rVariables.EcNodes.resize(TNumNodes);
    rVariables.GNodes.resize(TNumNodes);
    rVariables.rhoNodes.resize(TNumNodes);

    // only get values from first neighbour
    int nValues = 0;

    Element rNeighbour = neighbours[0];

    auto rPropNeighbour = rNeighbour.GetProperties();

    const GeometryType& rNeighbourGeom = rNeighbour.GetGeometry();
    const GeometryData::IntegrationMethod rIntegrationMethodNeighbour = rNeighbour.GetIntegrationMethod();
    const Matrix& NContainerNeighbour = rNeighbourGeom.ShapeFunctionsValues(rIntegrationMethodNeighbour);

    const IndexType rNumNodesNeighbour = rNeighbourGeom.size();
    const IndexType NumGPointsNeighbour = rNeighbourGeom.IntegrationPointsNumber(rIntegrationMethodNeighbour);


    // get density and porosity from element
    double rho_s = rPropNeighbour[DENSITY_SOLID];
    double rho_w = rPropNeighbour[DENSITY_WATER];
    double n = rPropNeighbour[POROSITY];

    std::vector<double> SaturationVector;
    vector<double> rDensityVector(NumGPointsNeighbour);
    vector<double> rConfinedStiffness(NumGPointsNeighbour);
    vector<double> rShearStiffness(NumGPointsNeighbour);
        

    rNeighbour.CalculateOnIntegrationPoints(DEGREE_OF_SATURATION, SaturationVector, rCurrentProcessInfo);

    // get Ec and G from constitutive matrix
    Matrix ConstitutiveMatrix;
        
    const ConstitutiveLaw::Pointer pConstitutiveLaw = rPropNeighbour[CONSTITUTIVE_LAW];
    std::vector<ConstitutiveLaw::Pointer> constitutiveLawVector;
    rNeighbour.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutiveLawVector, rCurrentProcessInfo);
        

    for (unsigned int GPoint = 0; GPoint < constitutiveLawVector.size(); ++GPoint) {


        ConstitutiveLaw::Features rLawFeatures;
        constitutiveLawVector[GPoint]->GetLawFeatures(rLawFeatures);
        Flags lawOptions = rLawFeatures.GetOptions();

        ConstitutiveLaw::Parameters ConstitutiveParameters(rNeighbourGeom, rPropNeighbour, rCurrentProcessInfo);

        constitutiveLawVector[GPoint]->CalculateValue(ConstitutiveParameters, CONSTITUTIVE_MATRIX, ConstitutiveMatrix);

        if (TDim == 2)
        {
            int rEcIdx;
            int rGIdx;
            //ConstitutiveLaw::AXISYMMETRIC_LAW;
            bool tmp11 = lawOptions.Is(ConstitutiveLaw::INFINITESIMAL_STRAINS);

            if (lawOptions.Is(ConstitutiveLaw::PLANE_STRAIN_LAW))
            {
                //ConstitutiveLaw::INFINITESIMAL_STRAINS;
                rEcIdx = INDEX_2D_PLANE_STRAIN_XX;
            }
            else if (lawOptions.Is(ConstitutiveLaw::PLANE_STRESS_LAW))
            {
                rEcIdx = INDEX_2D_PLANE_STRESS_XX;
            }
            else
            {
                KRATOS_ERROR << "In 2D, lysmer absorbing boundary is only valid in plane strain" << std::endl;
            }
            //const int EcIdx = ConstitutiveLaw::PLANE_STRAIN_LAW ? INDEX_2D_PLANE_STRAIN_XX : ConstitutiveLaw::PLANE_STRESS_LAW ? INDEX_2D_PLANE_STRESS_XX : "i is over 0"

            Ec = ConstitutiveMatrix(rEcIdx, rEcIdx);
            G = ConstitutiveMatrix(INDEX_2D_PLANE_STRAIN_XY, INDEX_2D_PLANE_STRAIN_XY);

        }
        else if (TDim == 3)
        {
            Ec = ConstitutiveMatrix(INDEX_3D_XX, INDEX_3D_XX);
            G = ConstitutiveMatrix(INDEX_3D_XZ, INDEX_3D_XZ);
        }
        else
        {
            KRATOS_ERROR << "Lysmer absorbing boundary is only valid in a 2D or 3D space" << std::endl;
        }

        // calculate density of mixture
        double rho = (SaturationVector[GPoint] * rPropNeighbour[POROSITY] * rPropNeighbour[DENSITY_WATER]) + (1.0 - rPropNeighbour[POROSITY]) * rPropNeighbour[DENSITY_SOLID];

        rDensityVector[GPoint] = rho;
        rConfinedStiffness[GPoint] = Ec;
        rShearStiffness[GPoint] = G;

    }

    Matrix rExtrapolationMatrix = CalculateExtrapolationMatrixNeighbour(rNeighbour);
    /*Matrix rExtrapolationMatrix = ZeroMatrix(rNumNodesNeighbour, NumGPointsNeighbour);

    if (TDim == 2)
    {
        if (rNumNodesNeighbour == 3)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixTriangle(rExtrapolationMatrix, rIntegrationMethodNeighbour);
        }
        else if (rNumNodesNeighbour == 4)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixQuad(rExtrapolationMatrix, rIntegrationMethodNeighbour);
        }
    }
    if (TDim == 3)
    {
        if (rNumNodesNeighbour == 4)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixTetra(rExtrapolationMatrix, rIntegrationMethodNeighbour);
        }
        else if (rNumNodesNeighbour == 8)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixHexa(rExtrapolationMatrix, rIntegrationMethodNeighbour);
        }
    }*/

    auto EcNodesNeighbour = prod(rExtrapolationMatrix, rConfinedStiffness);
    auto GNodesNeighbour = prod(rExtrapolationMatrix, rShearStiffness);
    auto rhoNodesNeighbour = prod(rExtrapolationMatrix, rDensityVector);

    for (unsigned int k = 0; k < rNumNodesNeighbour; ++k)
    {
        for (unsigned int node = 0; node < TNumNodes; ++node)
        {
            if (rNeighbourGeom[k].Id() == rGeom[node].Id())
            {
                rVariables.EcNodes[node] = EcNodesNeighbour(k);
                rVariables.GNodes[node] = GNodesNeighbour(k);
                rVariables.rhoNodes[node] = rhoNodesNeighbour(k);
            }
        }
    }   
}


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
GetVariables(
    NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    // gets average of variables as stored in neighbour elements
    this->GetNeighbourElementVariables(rVariables, rCurrentProcessInfo);

    // calculate pressure wave and shear wave velocity
    //rVariables.vp = sqrt(rVariables.Ec / rVariables.rho);
    //rVariables.vs = sqrt(rVariables.G / rVariables.rho);

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
//

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwLysmerAbsorbingCondition<2,2>;
template class UPwLysmerAbsorbingCondition<2,3>;
template class UPwLysmerAbsorbingCondition<3,3>;
template class UPwLysmerAbsorbingCondition<3,4>;

} // Namespace Kratos.
