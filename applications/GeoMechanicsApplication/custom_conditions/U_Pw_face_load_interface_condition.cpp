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
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// Application includes
#include "custom_conditions/U_Pw_face_load_interface_condition.hpp"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwFaceLoadInterfaceCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwFaceLoadInterfaceCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwFaceLoadInterfaceCondition<TDim,TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    const GeometryType& Geom = this->GetGeometry();

    //Compute initial gap of the joint
    this->CalculateInitialGap(Geom);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< >
void UPwFaceLoadInterfaceCondition<2,2>::CalculateInitialGap(const GeometryType& Geom)
{
    mInitialGap.resize(1);
    
    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 1 ) - Geom.GetPoint( 0 );
    mInitialGap[0] = norm_2(Vx);
}

//----------------------------------------------------------------------------------------

template< >
void UPwFaceLoadInterfaceCondition<3,4>::CalculateInitialGap(const GeometryType& Geom)
{
    mInitialGap.resize(2);
    
    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 3 ) - Geom.GetPoint( 0 );
    mInitialGap[0] = norm_2(Vx);

    noalias(Vx) = Geom.GetPoint( 2 ) - Geom.GetPoint( 1 );
    mInitialGap[1] = norm_2(Vx);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwFaceLoadInterfaceCondition<TDim,TNumNodes>::
    CalculateRHS(VectorType& rRightHandSideVector,
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
    for (unsigned int i = 0; i<NumGPoints; i++)
        (JContainer[i]).resize(TDim,LocalDim,false);
    Geom.Jacobian(JContainer, this->GetIntegrationMethod());
    
    //Condition variables
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    ConditionUtilities::GetDisplacementsVector(DisplacementVector,Geom);
    array_1d<double,TNumNodes*TDim> FaceLoadVector;
    ConditionUtilities::GetFaceLoadVector<TNumNodes>(FaceLoadVector,Geom);
    BoundedMatrix<double,TDim,TDim> RotationMatrix;
    const double& MinimumJointWidth = this->GetProperties()[MINIMUM_JOINT_WIDTH];
    bool ComputeJointWidth;
    double JointWidth;
    this->CheckJointWidth(JointWidth,ComputeJointWidth,RotationMatrix,MinimumJointWidth,Geom);
    array_1d<double,TDim> LocalRelDispVector;
    array_1d<double,TDim> RelDispVector;
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TDim> TractionVector;
    array_1d<double,TNumNodes*TDim> UVector;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute traction vector 
        ConditionUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(TractionVector,NContainer,FaceLoadVector,GPoint);
        
        //Compute Nu Matrix
        InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);
        
        if(ComputeJointWidth==true)
            this->CalculateJointWidth(JointWidth, Nu, DisplacementVector, RelDispVector, RotationMatrix, LocalRelDispVector, MinimumJointWidth,GPoint);
        
        //Compute weighting coefficient for integration
        double integration_coefficient = this->CalculateIntegrationCoefficient(JContainer[GPoint], IntegrationPoints[GPoint].Weight(), JointWidth);
                
        //Contributions to the right hand side
        noalias(UVector) = prod(trans(Nu),TractionVector) * integration_coefficient;
        GeoElementUtilities::AssembleUBlockVector(rRightHandSideVector,UVector);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwFaceLoadInterfaceCondition<2,2>::CheckJointWidth(double& rJointWidth, bool& rComputeJointWidth, BoundedMatrix<double,2,2>& rRotationMatrix,
                                                                const double& MinimumJointWidth, const Element::GeometryType& Geom)
{
    //Line_interface_2d_2
    //Unitary vector in local x direction
    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 1 ) - Geom.GetPoint( 0 );
    double norm_x = norm_2(Vx);
    if(norm_x > 1.0e-8)
    {
        Vx[0] *= 1.0/norm_x;
        Vx[1] *= 1.0/norm_x;
        
        //Rotation Matrix
        rRotationMatrix(0,0) = Vx[0];
        rRotationMatrix(0,1) = Vx[1];
        
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
        if(Vz[2] > 0.0)
        {
            rRotationMatrix(1,0) = -Vx[1];
            rRotationMatrix(1,1) = Vx[0];
        }
        else
        {
            rRotationMatrix(1,0) = Vx[1];
            rRotationMatrix(1,1) = -Vx[0];
        }

        rComputeJointWidth = true;
    }
    else
    {
        rJointWidth = MinimumJointWidth;
        rComputeJointWidth = false;
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPwFaceLoadInterfaceCondition<3,4>::CheckJointWidth(double& rJointWidth, bool& rComputeJointWidth, BoundedMatrix<double,3,3>& rRotationMatrix,
                                                            const double& MinimumJointWidth, const Element::GeometryType& Geom)
{
    //Quadrilateral_interface_3d_4
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    const array_1d<double, 3>& P2 = Geom.GetPoint( 2 );
    noalias(pmid0) = 0.5 * (Geom.GetPoint( 0 ) + Geom.GetPoint( 3 ));
    noalias(pmid1) = 0.5 * (Geom.GetPoint( 1 ) + P2);
    
    //Unitary vector in local x direction
    array_1d<double,3> Vx;
    noalias(Vx) = pmid1 - pmid0;
    double inv_norm_x = 1.0/norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;
    Vx[2] *= inv_norm_x;
    
    //Unitary vector in local z direction
    array_1d<double,3> Vy;
    noalias(Vy) = P2 - pmid0;
    array_1d<double,3> Vz;
    MathUtils<double>::CrossProduct(Vz,Vx,Vy);
    double norm_z = norm_2(Vz);
    if(norm_z > 1.0e-8)
    {
        Vz[0] *= 1.0/norm_z;
        Vz[1] *= 1.0/norm_z;
        Vz[2] *= 1.0/norm_z;
        
        //Unitary vector in local y direction
        MathUtils<double>::CrossProduct(Vy, Vz, Vx);
        
        //Rotation Matrix
        rRotationMatrix(0,0) = Vx[0];
        rRotationMatrix(0,1) = Vx[1];
        rRotationMatrix(0,2) = Vx[2];
        
        rRotationMatrix(1,0) = Vy[0];
        rRotationMatrix(1,1) = Vy[1];
        rRotationMatrix(1,2) = Vy[2];
        
        rRotationMatrix(2,0) = Vz[0];
        rRotationMatrix(2,1) = Vz[1];
        rRotationMatrix(2,2) = Vz[2];
        
        rComputeJointWidth = true;
    }
    else
    {
        rJointWidth = MinimumJointWidth;
        rComputeJointWidth = false;
    }
}

//----------------------------------------------------------------------------------------
    
template< >
void UPwFaceLoadInterfaceCondition<2,2>::CalculateJointWidth( double& rJointWidth, const BoundedMatrix<double,2,4>& Nu,
                                                                        const array_1d<double,4>& DisplacementVector, array_1d<double,2>& rRelDispVector,
                                                                        const BoundedMatrix<double,2,2>& RotationMatrix,
                                                                        array_1d<double,2>& rLocalRelDispVector, const double& MinimumJointWidth,
                                                                        const unsigned int& GPoint )
{
    //Line_interface_2d_2
    noalias(rRelDispVector) = prod(Nu,DisplacementVector);
    
    noalias(rLocalRelDispVector) = prod(RotationMatrix,rRelDispVector);

    rJointWidth = mInitialGap[GPoint] + rLocalRelDispVector[0]; //The joint width is obtained in the local x direction
        
    if(rJointWidth < MinimumJointWidth)
    {
        rJointWidth = MinimumJointWidth;
    }

}

//----------------------------------------------------------------------------------------

template< >
void UPwFaceLoadInterfaceCondition<3,4>::CalculateJointWidth( double& rJointWidth, const BoundedMatrix<double,3,12>& Nu,
                                                                        const array_1d<double,12>& DisplacementVector, array_1d<double,3>& rRelDispVector,
                                                                        const BoundedMatrix<double,3,3>& RotationMatrix,
                                                                        array_1d<double,3>& rLocalRelDispVector, const double& MinimumJointWidth,
                                                                        const unsigned int& GPoint )
{
    //Quadrilateral_interface_3d_4
    noalias(rRelDispVector) = prod(Nu,DisplacementVector);
    
    noalias(rLocalRelDispVector) = prod(RotationMatrix,rRelDispVector);

    rJointWidth = mInitialGap[GPoint] + rLocalRelDispVector[1]; //The joint width is obtained in the local y direction

    if(rJointWidth < MinimumJointWidth)
    {
        rJointWidth = MinimumJointWidth;
    }
}

//----------------------------------------------------------------------------------------

template< >
double UPwFaceLoadInterfaceCondition<2,2>::CalculateIntegrationCoefficient(const Matrix& Jacobian, const double& Weight, const double& JointWidth)
{
    // Note: since we cannot include the determinant of the Jacobian here (because the Jacobian could be singular), we do not multiply the JointWidth by the Weight 
    //       In a normal line load we would have |J| * w = L/2.0 * 2.0 = L
    return JointWidth;
}

//----------------------------------------------------------------------------------------

template< >
double UPwFaceLoadInterfaceCondition<3, 4>::CalculateIntegrationCoefficient(const Matrix& Jacobian, const double& Weight, const double& JointWidth)
{
    KRATOS_TRY

    double ds = MathUtils<double>::Norm(column(Jacobian, 0));
    return Weight * ds * JointWidth;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwFaceLoadInterfaceCondition<2,2>;
template class UPwFaceLoadInterfaceCondition<3,4>;

} // Namespace Kratos.
