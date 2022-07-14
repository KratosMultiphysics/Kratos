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

    //Previous definitions
    GeometryType& Geom = this->GetGeometry();
    const PropertiesType& prop = this->GetProperties();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints( mThisIntegrationMethod );
    const unsigned int NumGPoints = IntegrationPoints.size();
    const unsigned int LocalDim = Geom.LocalSpaceDimension();
    
    //Containers of variables at all integration points
    const Matrix& NContainer = Geom.ShapeFunctionsValues( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer(NumGPoints);
    for(unsigned int i = 0; i<NumGPoints; ++i)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian( JContainer, mThisIntegrationMethod );
    
    //Condition variables
    BoundedMatrix<double, TDim, TNumNodes* TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
    array_1d<double, TNumNodes* TDim> rTractionVector;

    NormalLysmerAbsorbingVariables rVariables;

    this->GetVariables(rVariables, CurrentProcessInfo);
    this->CalculateTractionVector(rVariables, CurrentProcessInfo, Geom, rTractionVector);

    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {

        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nu, NContainer, GPoint);

        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(rVariables.IntegrationCoefficient,
                                              JContainer[GPoint],
                                              IntegrationPoints[GPoint].Weight() );

        rAbsMatrix += prod(trans(Nu), Nu) * rVariables.IntegrationCoefficient;
    }

    rVariables.UVector = prod(rAbsMatrix, rTractionVector);
    this->CalculateAndAddRHS(rRightHandSideVector, rVariables);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim, TNumNodes>::
CalculateTractionVector(NormalLysmerAbsorbingVariables& rVariables, const ProcessInfo& CurrentProcessInfo, Element::GeometryType& Geom, 
    array_1d<double, TNumNodes* TDim>& rTractionVector)
{
    double alpha1 = 1;
    double alpha2 = 1;


    array_1d<double, TDim> rTractionVectorConstants;
    array_1d<double, TDim> rNodalVelocityVector;
    array_1d<double, TDim> rLocalVelocityVector;

    array_1d<double, TNumNodes* TDim> rGlobalVelocityVector;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(rGlobalVelocityVector, Geom, VELOCITY);


    // calculate rotation matrix
    BoundedMatrix<double, TDim, TDim> rotationMatrix;
    CalculateRotationMatrix(rotationMatrix, Geom);


    // calculate constant traction vector part
    rTractionVectorConstants[0] = rVariables.vs * rVariables.rho * alpha2;
    if (TDim == 2)
    {
        rTractionVectorConstants[1] = rVariables.vp * rVariables.rho * alpha1;
    }
    else if(TDim == 3)
    {
        rTractionVectorConstants[1] = rVariables.vs * rVariables.rho * alpha2;
        rTractionVectorConstants[2] = rVariables.vp * rVariables.rho * alpha1;
    }

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        const unsigned int Local_i = i * TDim;

        // check if absorbing forces are already set at current node, if so, do not fill traction vector
        if (Geom[i].GetValue(IS_ABSORBING))
        {
            for (unsigned int dim = 0; dim < TDim; ++dim) {
                rTractionVector[Local_i + dim] = 0;
            }
        }
        else 
        {
            // get velocity vector at current node
            for (unsigned int dim = 0; dim < TDim; ++dim) {
                rNodalVelocityVector[dim] = rGlobalVelocityVector[Local_i + dim];
            }

            // rotate nodal velocity vector to local system
            rLocalVelocityVector = prod(rNodalVelocityVector, rotationMatrix);


            // calculate traction vector
            for (unsigned int dim = 0; dim < TDim; ++dim) {
                rTractionVector[Local_i + dim] = -rLocalVelocityVector[dim] * rTractionVectorConstants[dim];
            }
            Geom[i].SetValue(IS_ABSORBING, true);
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

    //auto absorbing_factors = this->GetValue(ABSORBING_FACTORS);

    //rVariables.alpha1 = absorbing_factors[0];
    //rVariables.alpha2 = absorbing_factors[1];
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
void UPwLysmerAbsorbingCondition<2, 2>::
CalculateLocalVelocityVector(array_1d<double, 2>& rTractionVector,
    const Matrix& Jacobian,
    const Matrix& NContainer,
    const NormalLysmerAbsorbingVariables& Variables,
    const unsigned int& GPoint)
{
    double NormalStress = 0.0;
    double TangentialStress = 0.0;


    double dx_dxi = Jacobian(0, 0);
    double dy_dxi = Jacobian(1, 0);

    rTractionVector[0] = TangentialStress * dx_dxi - NormalStress * dy_dxi;
    rTractionVector[1] = NormalStress * dx_dxi + TangentialStress * dy_dxi;
}

//----------------------------------------------------------------------------------------
template< >
void UPwLysmerAbsorbingCondition<3, 3>::
CalculateLocalVelocityVector(array_1d<double, 3>& rTractionVector,
    const Matrix& Jacobian,
    const Matrix& NContainer,
    const NormalLysmerAbsorbingVariables& Variables,
    const unsigned int& GPoint)
{
    double NormalStress = 0.0;

  /*  for (unsigned int i = 0; i < 3; ++i) {
        NormalStress += NContainer(GPoint, i) * Variables.NormalStressVector[i];
    }*/

    double NormalVector[3];

    NormalVector[0] = Jacobian(1, 0) * Jacobian(2, 1) - Jacobian(2, 0) * Jacobian(1, 1);

    NormalVector[1] = Jacobian(2, 0) * Jacobian(0, 1) - Jacobian(0, 0) * Jacobian(2, 1);

    NormalVector[2] = Jacobian(0, 0) * Jacobian(1, 1) - Jacobian(1, 0) * Jacobian(0, 1);

    rTractionVector[0] = NormalStress * NormalVector[0];
    rTractionVector[1] = NormalStress * NormalVector[1];
    rTractionVector[2] = NormalStress * NormalVector[2];
}

//----------------------------------------------------------------------------------------

template< >
void UPwLysmerAbsorbingCondition<3, 4>::
CalculateLocalVelocityVector(array_1d<double, 3>& rTractionVector,
    const Matrix& Jacobian,
    const Matrix& NContainer,
    const NormalLysmerAbsorbingVariables& Variables,
    const unsigned int& GPoint)
{
    double NormalStress = 0.0;

    //for (unsigned int i = 0; i < 4; ++i) {
    //    NormalStress += NContainer(GPoint, i) * Variables.NormalStressVector[i];
    //}

    double NormalVector[3];

    NormalVector[0] = Jacobian(1, 0) * Jacobian(2, 1) - Jacobian(2, 0) * Jacobian(1, 1);

    NormalVector[1] = Jacobian(2, 0) * Jacobian(0, 1) - Jacobian(0, 0) * Jacobian(2, 1);

    NormalVector[2] = Jacobian(0, 0) * Jacobian(1, 1) - Jacobian(1, 0) * Jacobian(0, 1);

    rTractionVector[0] = NormalStress * NormalVector[0];
    rTractionVector[1] = NormalStress * NormalVector[1];
    rTractionVector[2] = NormalStress * NormalVector[2];
}


template< >
void UPwLysmerAbsorbingCondition<2, 2>::CalculateRotationMatrix( BoundedMatrix<double, 2, 2>& rRotationMatrix, const Element::GeometryType& Geom)
{
    //Line_interface_2d_2
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
template class UPwLysmerAbsorbingCondition<3,3>;
template class UPwLysmerAbsorbingCondition<3,4>;

} // Namespace Kratos.
