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
//#include "custom_processes/find_neighbour_elements_of_conditions_process.cpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwLysmerAbsorbingCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwLysmerAbsorbingCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    // reset bool which indicates if an abosrbing force is present on current node 
    GeometryType& Geom = this->GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i) {
    
        Geom[i].SetValue(IS_ABSORBING, false);
    }
}

template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim,TNumNodes>::
    CalculateRHS(VectorType& rRightHandSideVector,
                 const ProcessInfo& CurrentProcessInfo)
{        

    int n_dof = TNumNodes * TDim;
    BoundedMatrix<double, TNumNodes* TDim, TNumNodes* TDim> rAbsMatrix = ZeroMatrix(n_dof, n_dof);
    BoundedMatrix<double, TNumNodes* TDim, TNumNodes* TDim> rAbsKMatrix = ZeroMatrix(n_dof, n_dof);

    //Previous definitions
    GeometryType& Geom = this->GetGeometry();
    PropertiesType& prop = this->GetProperties();

    GeometryData::IntegrationMethod rIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(rIntegrationMethod);
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues(rIntegrationMethod);
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i<NumGPoints; ++i)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian( JContainer, rIntegrationMethod);
    
    //Condition variables
    BoundedMatrix<double, TDim, TNumNodes* TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
    array_1d<double, TNumNodes* TDim> rTractionVector;

    NormalLysmerAbsorbingVariables rVariables;

    this->GetVariables(rVariables, CurrentProcessInfo);
    this->CalculateTractionVector(rVariables, CurrentProcessInfo, Geom, rTractionVector);

    BoundedMatrix<double, TDim, TDim>             VpMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TNumNodes* TDim> AuxAbsMatrix = ZeroMatrix(TDim, TNumNodes * TDim);

    BoundedMatrix<double, TDim, TDim>             KMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TNumNodes* TDim> AuxAbsKMatrix = ZeroMatrix(TDim, TNumNodes * TDim);

    //for (unsigned int idim = 0; idim < TDim; ++idim) {
    //    for (unsigned int jdim = 0; jdim < TDim; ++jdim) {
    //        VpMatrix(idim, jdim) = abs(rTractionVector[0]);
    //        KMatrix(idim, jdim) = rVariables.Ec / thickness_virtual_layer;
    //        //VpMatrix(idim, jdim) = 0;
    //    }
    //}
    //VpMatrix(0, 0) = abs(rTractionVector[0]);
    //KMatrix(0, 0) = rVariables.Ec/ thickness_virtual_layer;
    rVariables.UVector = ZeroVector(TNumNodes * TDim);

    //const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ConditionSize = TNumNodes * TDim + TNumNodes;
    noalias(rRightHandSideVector) = ZeroVector(ConditionSize);

    array_1d<double, TNumNodes* TDim> rGlobalVelocityVector;
    array_1d<double, TNumNodes* TDim> rGlobalVelocityVector_prev;
    array_1d<double, TNumNodes* TDim> rDeltaGlobalVelocity;
    array_1d<double, TNumNodes* TDim> rGlobalDisplacementVector;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rGlobalVelocityVector, Geom, VELOCITY, 0);
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rGlobalVelocityVector_prev, Geom, VELOCITY, 1);


    // set velocity to zero, if the corresponding node is already set as an absorbing boundary
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        const unsigned int Global_i = i * (TDim);
        if (Geom[i].GetValue(IS_ABSORBING)) {
            for (unsigned int idim = 0; idim < TDim; ++idim)
            {
                rGlobalVelocityVector_prev[Global_i + idim] = 0;
            }
        }
    }
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {

        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nu, NContainer, GPoint);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(rVariables.IntegrationCoefficient,
                                              JContainer[GPoint],
                                              IntegrationPoints[GPoint].Weight() );


        AuxAbsKMatrix = prod(rVariables.KAbsMatrix, Nu);
        rAbsKMatrix += prod(trans(Nu), AuxAbsKMatrix) * rVariables.IntegrationCoefficient;

        AuxAbsMatrix = prod(rVariables.CAbsMatrix, Nu);
        rAbsMatrix = prod(trans(Nu), AuxAbsMatrix)* rVariables.IntegrationCoefficient;


        rVariables.UVector += prod(trans(rAbsMatrix), -rGlobalVelocityVector_prev);
        //rAbsMatrix += prod(trans(Nu), Nu) * rVariables.IntegrationCoefficient * - rVariables.vp* rVariables.rho/2;
    }
    for (unsigned int i = 0; i < TNumNodes; ++i) {
    
        Geom[i].SetValue(IS_ABSORBING, true);
    }
    
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rGlobalDisplacementVector, Geom, DISPLACEMENT, 1);

    //double tmp = Geom[0].FastGetSolutionStepValue(VELOCITY_X, 1);
    //rVariables.UVector = prod(trans(rAbsMatrix), -rGlobalVelocityVector) + prod(trans(rAbsKMatrix), -rGlobalDisplacementVector);
    //rVariables.UVector =  prod(trans(rAbsKMatrix), -rGlobalDisplacementVector);
    this->CalculateAndAddRHS(rRightHandSideVector, rVariables);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
CalculateTractionVector(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, Element::GeometryType& Geom, 
    array_1d<double, TNumNodes* TDim>& rTractionVector)
{
    double alpha1 = 0.986;
    double alpha2 = 0.747;


    array_1d<double, 2> rDampingConstants;
    array_1d<double, 2> rStiffnessConstants;
    array_1d<double, TDim> rNodalVelocityVector;
    array_1d<double, TDim> rLocalVelocityVector;


    // calculate rotation matrix
    BoundedMatrix<double, TDim, TDim> rotationMatrix;
    CalculateRotationMatrix(rotationMatrix, Geom);

    const int local_perpendicular_direction = TDim-1;


    // calculate constant traction vector part


    rDampingConstants[0] = rVariables.vs * rVariables.rho * rVariables.s_factor;
    rDampingConstants[1] = rVariables.vp * rVariables.rho * rVariables.p_factor;


    rStiffnessConstants[0] = rVariables.G / rVariables.virtual_thickness;
    rStiffnessConstants[1] = rVariables.Ec / rVariables.virtual_thickness;


    BoundedMatrix<double, TDim, TDim>             localCMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TDim>             localKMatrix = ZeroMatrix(TDim, TDim);

    BoundedMatrix<double, TDim, TDim>             auxLocalCMatrix = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, TDim>             auxLocalKMatrix = ZeroMatrix(TDim, TDim);

    rVariables.CAbsMatrix = ZeroMatrix(TDim, TDim);
    rVariables.KAbsMatrix = ZeroMatrix(TDim, TDim);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        for (unsigned int jdim = 0; jdim < TDim; ++jdim) {
            localCMatrix(idim, idim) = rDampingConstants[0];
            localKMatrix(idim, idim) = rStiffnessConstants[0];
        }
    }
    localCMatrix(local_perpendicular_direction, local_perpendicular_direction) = rDampingConstants[1];
    localKMatrix(local_perpendicular_direction, local_perpendicular_direction) = rStiffnessConstants[1];
  
    auxLocalCMatrix = prod(localCMatrix, rotationMatrix);
    auxLocalKMatrix = prod(localKMatrix, rotationMatrix);
    rVariables.CAbsMatrix = prod(trans(rotationMatrix), auxLocalCMatrix);
    rVariables.KAbsMatrix = prod(trans(rotationMatrix), auxLocalKMatrix);

    for (unsigned int idim = 0; idim < TDim; ++idim) {
        for (unsigned int jdim = 0; jdim < TDim; ++jdim) {
            rVariables.CAbsMatrix(idim, jdim) = abs(rVariables.CAbsMatrix(idim, jdim));
            rVariables.KAbsMatrix(idim, jdim) = abs(rVariables.KAbsMatrix(idim, jdim));
        }
    }
}


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
void UPwLysmerAbsorbingCondition<TDim,TNumNodes>::
    CalculateAndAddRHS( VectorType& rRightHandSideVector,
        NormalLysmerAbsorbingVariables& rVariables )
{
    ////Adding contribution to right hand side

    GeoElementUtilities::
        AssembleUBlockVector< TDim, TNumNodes >(rRightHandSideVector,
                                                rVariables.UVector);
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
    //todo
    // 
    // 
    ////triangle_3d_3
    //array_1d<double, 3> pmid0;
    //array_1d<double, 3> pmid1;
    //array_1d<double, 3> P2 = Geom.GetPoint(2);
    //noalias(pmid0) = 0.5 * (Geom.GetPoint(0) + Geom.GetPoint(3));
    //noalias(pmid1) = 0.5 * (Geom.GetPoint(1) + P2);

    ////Unitary vector in local x direction
    //array_1d<double, 3> Vx;
    //noalias(Vx) = pmid1 - pmid0;
    //double inv_norm_x = 1.0 / norm_2(Vx);
    //Vx[0] *= inv_norm_x;
    //Vx[1] *= inv_norm_x;
    //Vx[2] *= inv_norm_x;

    ////Unitary vector in local z direction
    //array_1d<double, 3> Vy;
    //noalias(Vy) = P2 - pmid0;
    //array_1d<double, 3> Vz;
    //MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    //double norm_z = norm_2(Vz);

    //Vz[0] *= 1.0 / norm_z;
    //Vz[1] *= 1.0 / norm_z;
    //Vz[2] *= 1.0 / norm_z;

    ////Unitary vector in local y direction
    //MathUtils<double>::CrossProduct(Vy, Vz, Vx);

    ////Rotation Matrix
    //rRotationMatrix(0, 0) = Vx[0];
    //rRotationMatrix(0, 1) = Vx[1];
    //rRotationMatrix(0, 2) = Vx[2];

    //rRotationMatrix(1, 0) = Vy[0];
    //rRotationMatrix(1, 1) = Vy[1];
    //rRotationMatrix(1, 2) = Vy[2];

    //rRotationMatrix(2, 0) = Vz[0];
    //rRotationMatrix(2, 1) = Vz[1];
    //rRotationMatrix(2, 2) = Vz[2];
}
//----------------------------------------------------------------------------------------

template< >
void UPwLysmerAbsorbingCondition<3, 4>::CalculateRotationMatrix( BoundedMatrix<double, 3, 3>& rRotationMatrix, const Element::GeometryType& Geom)
{
    //Quadrilateral_interface_3d_4
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
