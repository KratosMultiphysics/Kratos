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

template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim,TNumNodes>::
    CalculateRHS(VectorType& rRightHandSideVector,
                 const ProcessInfo& CurrentProcessInfo)
{        
    //Previous definitions
    const GeometryType& Geom = this->GetGeometry();
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
    array_1d<double, TNumNodes* TDim> VelocityVector;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(VelocityVector, Geom, VELOCITY);
    BoundedMatrix<double, TDim, TDim> rotationMatrix;

    // calculate rotation matrix
    CalculateRotationMatrix(rotationMatrix, Geom);
    BoundedMatrix<double, TDim, TNumNodes* TDim> Nu = ZeroMatrix(TDim, TNumNodes * TDim);
    array_1d<double, TDim> relVelVector;
    array_1d<double, TDim> localRelVelVector;
    //std::vector< Vector > mVelocityVector;
    //array_1d<double, TDim> nodalVelocityVector;
    //for(unsigned int i=0; i<TNumNodes; ++i)
    //{
    //    nodalVelocityVector = Geom[i].FastGetSolutionStepValue(VELOCITY);
    //    for (unsigned int j = 0; j < TDim; ++j)
    //    {
    //        VelocityVector(i,j) = nodalVelocityVector[j];
    //    }
    //}

    NormalLysmerAbsorbingVariables Variables;


    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {

        GeoElementUtilities::CalculateNuMatrix<TDim, TNumNodes>(Nu, NContainer, GPoint);

        noalias(relVelVector) = prod(Nu, VelocityVector);

        noalias(localRelVelVector) = prod(rotationMatrix, relVelVector);

        // normal velocity????
        //Compute normal flux 
       /* Variables.NormalFlux = 0.0;
        for (unsigned int i=0; i<TNumNodes; ++i) {
            Variables.NormalFlux += NContainer(GPoint,i)*NormalFluxVector[i];
        }*/
        
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

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwLysmerAbsorbingCondition<TDim,TNumNodes>::
    CalculateAndAddRHS( VectorType& rRightHandSideVector,
        NormalLysmerAbsorbingVariables& rVariables )
{


    ////Adding contribution to Mass matrix
    //noalias(rMassMatrix) += prod(trans(Nut), AuxDensityMatrix)
    //    * Variables.IntegrationCoefficientInitialConfiguration;

    double alpha = 1;
    double rho = 2000;
    double E = 1e9;
    double nu = 0.2;

    double Ec = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    //double vp = math::sqrt(Ec / rho);

    //noalias(rVariables.PVector) = rVariables.NormalFlux * rVariables.Np * rVariables.IntegrationCoefficient;




    GeoElementUtilities::
        AssemblePBlockVector< TDim, TNumNodes >(rRightHandSideVector,
                                                rVariables.PVector);
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

  /*  for (unsigned int i = 0; i < 2; ++i) {
        NormalStress += NContainer(GPoint, i) * Variables.NormalStressVector[i];
        TangentialStress += NContainer(GPoint, i) * Variables.TangentialStressVector[i];
    }*/

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
