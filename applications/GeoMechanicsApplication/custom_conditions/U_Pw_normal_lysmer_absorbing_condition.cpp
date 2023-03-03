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

template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLhsMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{

    BoundedMatrix<double, N_DOF, N_DOF> stiffness_matrix;

    this->CalculateConditionStiffnessMatrix(stiffness_matrix, rCurrentProcessInfo);

    this->AddLHS(rLhsMatrix, stiffness_matrix);

    this->CalculateAndAddRHS(rRightHandSideVector, rLhsMatrix);
}

template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::CalculateConditionStiffnessMatrix(BoundedMatrix<double, N_DOF, N_DOF>& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    //Previous definitions
    GeometryType& r_geom = this->GetGeometry();

    GeometryData::IntegrationMethod integration_method = this->mThisIntegrationMethod;
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(integration_method);
    const unsigned int num_g_points = r_integration_points.size();
    const unsigned int local_dim = r_geom.LocalSpaceDimension();

    //Containers of variables at all integration points
    const Matrix& r_n_container = r_geom.ShapeFunctionsValues(integration_method);
    GeometryType::JacobiansType j_container(num_g_points);
    for (unsigned int i = 0; i < num_g_points; ++i)
        j_container[i].resize(TDim, local_dim, false);
    r_geom.Jacobian(j_container, integration_method);

    //Condition variables
    BoundedMatrix<double, TDim, N_DOF> nu_matrix = ZeroMatrix(TDim, N_DOF);

    NormalLysmerAbsorbingVariables r_variables;

    this->GetVariables(r_variables, rCurrentProcessInfo);

    BoundedMatrix<double, TDim, N_DOF> aux_abs_k_matrix;
    rStiffnessMatrix = ZeroMatrix(N_DOF, N_DOF);

    //Loop over integration points
    for (unsigned int g_point = 0; g_point < num_g_points; ++g_point) {

        // calculate
        r_variables.Ec = 0.0;
        r_variables.G = 0.0;
        for (unsigned int node = 0; node < r_geom.size(); ++node)
        {
            r_variables.Ec += r_n_container(g_point, node) * r_variables.EcNodes[node];
            r_variables.G += r_n_container(g_point, node) * r_variables.GNodes[node];
        }

        this->CalculateNodalStiffnessMatrix(r_variables, r_geom);

        // calculate displacement shape function matrix
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(nu_matrix, r_n_container, g_point);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(r_variables.IntegrationCoefficient,
            j_container[g_point],
            r_integration_points[g_point].Weight());

        // set stiffness part of absorbing matrix
        aux_abs_k_matrix = prod(r_variables.KAbsMatrix, nu_matrix);
        rStiffnessMatrix += prod(trans(nu_matrix), aux_abs_k_matrix) * r_variables.IntegrationCoefficient;
    }
}


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    BoundedMatrix<double, N_DOF, N_DOF> stiffness_matrix;

    this->CalculateConditionStiffnessMatrix(stiffness_matrix, rCurrentProcessInfo);

    MatrixType global_stiffness_matrix = ZeroMatrix(CONDITION_SIZE, CONDITION_SIZE);
    GeoElementUtilities::AssembleUBlockMatrix< TDim, TNumNodes >(global_stiffness_matrix, stiffness_matrix);

    this->CalculateAndAddRHS(rRightHandSideVector, global_stiffness_matrix);
}


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{

    //Previous definitions
    GeometryType& r_geom = this->GetGeometry();

    GeometryData::IntegrationMethod r_integration_method = this->mThisIntegrationMethod;
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geom.IntegrationPoints(r_integration_method);
    const unsigned int num_g_points = r_integration_points.size();
    const unsigned int local_dim = r_geom.LocalSpaceDimension();

    //Containers of variables at all integration points
    const Matrix& r_n_container = r_geom.ShapeFunctionsValues(r_integration_method);
    GeometryType::JacobiansType j_container(num_g_points);
    for (unsigned int i = 0; i < num_g_points; ++i)
        (j_container[i]).resize(TDim, local_dim, false);
    r_geom.Jacobian(j_container, r_integration_method);

    //Condition variables
    BoundedMatrix<double, TDim, N_DOF> nu_matrix = ZeroMatrix(TDim, N_DOF);

    NormalLysmerAbsorbingVariables r_variables;
    this->GetVariables(r_variables, rCurrentProcessInfo);

    BoundedMatrix<double, TDim, N_DOF> aux_abs_matrix;
    BoundedMatrix<double, N_DOF, N_DOF> abs_matrix = ZeroMatrix(N_DOF, N_DOF);

    //Loop over integration points
    for (unsigned int g_point = 0; g_point < num_g_points; ++g_point) {

        // calculate
        r_variables.rho = 0.0;
        r_variables.Ec = 0.0;
        r_variables.G = 0.0;
        for (unsigned int node = 0; node < r_geom.size(); ++node)
        {
            r_variables.rho += r_n_container(g_point, node) * r_variables.rhoNodes[node];
            r_variables.Ec += r_n_container(g_point, node) * r_variables.EcNodes[node];
            r_variables.G += r_n_container(g_point, node) * r_variables.GNodes[node];
        }
        r_variables.vp = sqrt(r_variables.Ec / r_variables.rho);
        r_variables.vs = sqrt(r_variables.G / r_variables.rho);

        this->CalculateNodalDampingMatrix(r_variables, r_geom);

        // calculate displacement shape function matrix
        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(nu_matrix, r_n_container, g_point);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(r_variables.IntegrationCoefficient,
            j_container[g_point],
            r_integration_points[g_point].Weight());

        // set damping part of absorbing matrix
        aux_abs_matrix = prod(r_variables.CAbsMatrix, nu_matrix);
        abs_matrix += prod(trans(nu_matrix), aux_abs_matrix) * r_variables.IntegrationCoefficient;
    }

    // assemble left hand side vector
    this->AddLHS(rDampingMatrix, abs_matrix);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{

    KRATOS_TRY

        const GeometryType& r_geom = this->GetGeometry();

        if (rValues.size() != CONDITION_SIZE) {
            rValues.resize(CONDITION_SIZE, false);
        }

        if constexpr (TDim == 2) {
            unsigned int index = 0;
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                rValues[index++] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
                rValues[index++] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
                rValues[index++] = 0.0;
            }
        }
        else if constexpr (TDim == 3) {
            unsigned int index = 0;
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                rValues[index++] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
                rValues[index++] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
                rValues[index++] = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT_Z, Step);
                rValues[index++] = 0.0;
            }
        }
        else {
            KRATOS_ERROR << "undefined dimension in GetValuesVector... illegal operation!!" << this->Id() << std::endl;
        }
    KRATOS_CATCH("")
}


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY

    const GeometryType& r_geom = this->GetGeometry();

    if (rValues.size() != CONDITION_SIZE)
        rValues.resize(CONDITION_SIZE, false);

    if constexpr(TDim == 2) {
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rValues[index++] = r_geom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
            rValues[index++] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
            rValues[index++] = 0.0;
        }
    }
    else if constexpr (TDim == 3) {
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rValues[index++] = r_geom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
            rValues[index++] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
            rValues[index++] = r_geom[i].FastGetSolutionStepValue(VELOCITY_Z, Step);
            rValues[index++] = 0.0;
        }
    }
    else {
        KRATOS_ERROR << "undefined dimension in GetFirstDerivativesVector... illegal operation!!" << this->Id() << std::endl;
    }

    KRATOS_CATCH("")
}


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
CalculateNodalDampingMatrix(NormalLysmerAbsorbingVariables& rVariables, const Element::GeometryType& rGeom)
{
    array_1d<double, 2> damping_constants;

    // calculate rotation matrix
    BoundedMatrix<double, TDim, TDim> rotation_matrix;
    CalculateRotationMatrix(rotation_matrix, rGeom);

    const int local_perpendicular_direction = TDim - 1;

    // calculate constant traction vector part
    damping_constants[0] = rVariables.vs * rVariables.rho * rVariables.s_factor;
    damping_constants[1] = rVariables.vp * rVariables.rho * rVariables.p_factor;

    BoundedMatrix<double, TDim, TDim> local_c_matrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TDim> aux_local_c_matrix = ZeroMatrix(TDim, TDim);

    rVariables.CAbsMatrix = ZeroMatrix(TDim, TDim);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        local_c_matrix(idim, idim) = damping_constants[0];
    }
    local_c_matrix(local_perpendicular_direction, local_perpendicular_direction) = damping_constants[1];

    aux_local_c_matrix = prod(local_c_matrix, rotation_matrix);
    rVariables.CAbsMatrix = prod(trans(rotation_matrix), aux_local_c_matrix);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
            rVariables.CAbsMatrix(idim, idim) = abs(rVariables.CAbsMatrix(idim, idim));
    }
}


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
CalculateNodalStiffnessMatrix(NormalLysmerAbsorbingVariables& rVariables, const Element::GeometryType& rGeom)
{
    array_1d<double, 2> stiffness_constants;

    // calculate rotation matrix
    BoundedMatrix<double, TDim, TDim> rotation_matrix;
    CalculateRotationMatrix(rotation_matrix, rGeom);

    const int local_perpendicular_direction = TDim - 1;

    // calculate constant traction vector part
    stiffness_constants[0] = rVariables.G / rVariables.virtual_thickness;
    stiffness_constants[1] = rVariables.Ec / rVariables.virtual_thickness;

    BoundedMatrix<double, TDim, TDim> local_k_matrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TDim> aux_local_k_matrix = ZeroMatrix(TDim, TDim);

    rVariables.KAbsMatrix = ZeroMatrix(TDim, TDim);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        local_k_matrix(idim, idim) = stiffness_constants[0];
    }
    local_k_matrix(local_perpendicular_direction, local_perpendicular_direction) = stiffness_constants[1];

    aux_local_k_matrix = prod(local_k_matrix, rotation_matrix);
    rVariables.KAbsMatrix = prod(trans(rotation_matrix), aux_local_k_matrix);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        rVariables.KAbsMatrix(idim, idim) = abs(rVariables.KAbsMatrix(idim, idim));
    }
}

/// <summary>
/// Calculates the extrapolation matrix for neighbour elements. Values from integration points are extrapolated to the nodes
/// </summary>
/// <param name="rNeighbourElement"></param>
/// <returns></returns>
template< unsigned int TDim, unsigned int TNumNodes >
Matrix UPwLysmerAbsorbingCondition<TDim, TNumNodes >::CalculateExtrapolationMatrixNeighbour(const Element& rNeighbourElement)
{
    const GeometryData::IntegrationMethod integration_method_neighbour = rNeighbourElement.GetIntegrationMethod();
    const GeometryType& r_neighbour_geom = rNeighbourElement.GetGeometry();
    const IndexType r_num_nodes_neighbour = r_neighbour_geom.size();
    const IndexType num_g_points_neighbour = r_neighbour_geom.IntegrationPointsNumber(integration_method_neighbour);

    Matrix extrapolation_matrix = ZeroMatrix(r_num_nodes_neighbour, num_g_points_neighbour);

    // Calculate extrapolation matrix for 2d elements
    if constexpr (TDim == 2)
    {
        if (r_num_nodes_neighbour == 3)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixTriangle(extrapolation_matrix, integration_method_neighbour);
            return extrapolation_matrix;
        }
        if (r_num_nodes_neighbour == 4)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixQuad(extrapolation_matrix, integration_method_neighbour);
            return extrapolation_matrix;
        }
       
    }
    // Calculate extrapolation matrix for 3d elements
    if constexpr (TDim == 3)
    {
        if (r_num_nodes_neighbour == 4)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixTetra(extrapolation_matrix, integration_method_neighbour);
            return extrapolation_matrix;
        }
        if (r_num_nodes_neighbour == 8)
        {
            GeoElementUtilities::CalculateExtrapolationMatrixHexa(extrapolation_matrix, integration_method_neighbour);
            return extrapolation_matrix;
        }

    }

    // if no extrapolation matrix is implemented, take average values at gauss points
    const double averaging_factor = 1 / num_g_points_neighbour;
    for (unsigned int node = 0; node < r_num_nodes_neighbour; ++node)
    {
        for (unsigned int g_point = 0; g_point < num_g_points_neighbour; ++g_point)
        {
            extrapolation_matrix(node, g_point) = averaging_factor;
        }
    }
    return extrapolation_matrix;
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
    std::vector <double> SaturationVector;
    Vector DensityVector(NumGPointsNeighbour);
    Vector ConfinedStiffness(NumGPointsNeighbour);
    Vector ShearStiffness(NumGPointsNeighbour);

    std::vector <double> stdConfinedStiffness(NumGPointsNeighbour);
    std::vector <double> stdShearStiffness(NumGPointsNeighbour);

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
    Vector EcNodesNeighbour = prod(rExtrapolationMatrix, ConfinedStiffness);
    Vector GNodesNeighbour = prod(rExtrapolationMatrix, ShearStiffness);
    Vector rhoNodesNeighbour = prod(rExtrapolationMatrix, DensityVector);

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
CalculateAndAddRHS(VectorType& rRightHandSideVector, const MatrixType& rStiffnessMatrix)
{

    if (rRightHandSideVector.size() != CONDITION_SIZE)
        rRightHandSideVector.resize(CONDITION_SIZE, false);
    noalias(rRightHandSideVector) = ZeroVector(CONDITION_SIZE);

    Vector displacements_vector = ZeroVector(CONDITION_SIZE);
    this->GetValuesVector(displacements_vector, 0);

    rRightHandSideVector -= prod(rStiffnessMatrix, displacements_vector);
}


template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
AddLHS(MatrixType& rLeftHandSideMatrix, const BoundedMatrix<double, N_DOF, N_DOF>& rUMatrix)
{
	// assemble left hand side vector
    if (rLeftHandSideMatrix.size1() != CONDITION_SIZE)
        rLeftHandSideMatrix.resize(CONDITION_SIZE, CONDITION_SIZE, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(CONDITION_SIZE, CONDITION_SIZE);

    //Adding contribution to left hand side
    GeoElementUtilities::
        AssembleUBlockMatrix< TDim, TNumNodes >(rLeftHandSideMatrix,
            rUMatrix);
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
