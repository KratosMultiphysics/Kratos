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
#include "custom_conditions/one-phase_flow/U_Pl_face_load_interface_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlFaceLoadInterfaceCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlFaceLoadInterfaceCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlFaceLoadInterfaceCondition<TDim,TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    const GeometryType& Geom = this->GetGeometry();

    //Compute initial gap of the joint
    this->CalculateInitialGap(Geom);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< >
void UPlFaceLoadInterfaceCondition<2,2>::CalculateInitialGap(const GeometryType& Geom)
{
    const double& InitialJointWidth = this->GetProperties()[INITIAL_JOINT_WIDTH];
    const double Tolerance = std::numeric_limits<double>::epsilon();
    mInitialGap.resize(1);
    
    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 1 ) - Geom.GetPoint( 0 );
    mInitialGap[0] = norm_2(Vx);
    if (mInitialGap[0] <= (InitialJointWidth+Tolerance)) {
        mInitialGap[0] = InitialJointWidth;
    } else {
        KRATOS_THROW_ERROR( std::invalid_argument, "The value of INITIAL_JOINT_WIDTH is smaller than the geometrical width.", "" )
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadInterfaceCondition<3,4>::CalculateInitialGap(const GeometryType& Geom)
{
    const double& InitialJointWidth = this->GetProperties()[INITIAL_JOINT_WIDTH];
    const double Tolerance = std::numeric_limits<double>::epsilon();
    mInitialGap.resize(2);
    
    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 3 ) - Geom.GetPoint( 0 );
    mInitialGap[0] = norm_2(Vx);
    if (mInitialGap[0] <= (InitialJointWidth+Tolerance)) {
        mInitialGap[0] = InitialJointWidth;
    } else {
        KRATOS_THROW_ERROR( std::invalid_argument, "The value of INITIAL_JOINT_WIDTH is smaller than the geometrical width.", "" )
    }

    noalias(Vx) = Geom.GetPoint( 2 ) - Geom.GetPoint( 1 );
    mInitialGap[1] = norm_2(Vx);
    if (mInitialGap[1] <= (InitialJointWidth+Tolerance)) {
        mInitialGap[1] = InitialJointWidth;
    } else {
        KRATOS_THROW_ERROR( std::invalid_argument, "The value of INITIAL_JOINT_WIDTH is smaller than the geometrical width.", "" )
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlFaceLoadInterfaceCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{        
    //Previous definitions
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
    array_1d<double,TNumNodes*TDim> DisplacementVector;
    PoroConditionUtilities::GetNodalVariableVector(DisplacementVector,Geom,DISPLACEMENT);
    array_1d<double,TNumNodes*TDim> FaceLoadVector;
    PoroConditionUtilities::GetNodalVariableVector(FaceLoadVector,Geom,FACE_LOAD);
    BoundedMatrix<double,TDim,TDim> RotationMatrix;
    const double& InitialJointWidth = this->GetProperties()[INITIAL_JOINT_WIDTH];
    bool ComputeJointWidth;
    double JointWidth;
    this->CheckJointWidth(JointWidth,ComputeJointWidth,RotationMatrix,InitialJointWidth,Geom);
    array_1d<double,TDim> LocalRelDispVector;
    array_1d<double,TDim> RelDispVector;
    BoundedMatrix<double,TDim, TNumNodes*TDim> Nu = ZeroMatrix(TDim, TNumNodes*TDim);
    array_1d<double,TDim> TractionVector;
    array_1d<double,TNumNodes*TDim> UVector;
    double IntegrationCoefficient;
    
    //Loop over integration points
    for(unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++)
    {
        //Compute traction vector 
        PoroConditionUtilities::InterpolateVariableWithComponents(TractionVector,NContainer,FaceLoadVector,GPoint);
        
        //Compute Nu Matrix
        InterfaceElementUtilities::CalculateNuMatrix(Nu,NContainer,GPoint);
        
        if(ComputeJointWidth==true)
            this->CalculateJointWidth(JointWidth, Nu, DisplacementVector, RelDispVector, RotationMatrix, LocalRelDispVector, InitialJointWidth,GPoint);
        
        //Compute weighting coefficient for integration
        this->CalculateIntegrationCoefficient(IntegrationCoefficient, JContainer[GPoint], integration_points[GPoint].Weight(), JointWidth);
                
        //Contributions to the right hand side
        noalias(UVector) = prod(trans(Nu),TractionVector) * IntegrationCoefficient;
        PoroConditionUtilities::AssembleUBlockVector(rRightHandSideVector,UVector);
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadInterfaceCondition<2,2>::CheckJointWidth(double& rJointWidth, bool& rComputeJointWidth, BoundedMatrix<double,2,2>& rRotationMatrix,
                                                                const double& InitialJointWidth, const Element::GeometryType& Geom)
{
    //Line_interface_2d_2
    //Unitary vector in local x direction
    array_1d<double,3> Vx;
    noalias(Vx) = Geom.GetPoint( 1 ) - Geom.GetPoint( 0 );
    double norm_x = norm_2(Vx);
    if(norm_x > InitialJointWidth) {
        Vx[0] *= 1.0/norm_x;
        Vx[1] *= 1.0/norm_x;

        //Rotation Matrix
        rRotationMatrix(0,0) = Vx[0];
        rRotationMatrix(0,1) = Vx[1];

        rRotationMatrix(1,0) = -Vx[1];
        rRotationMatrix(1,1) = Vx[0];

        rComputeJointWidth = true;
    } else {
        rJointWidth = InitialJointWidth;
        rComputeJointWidth = false;
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadInterfaceCondition<3,4>::CheckJointWidth(double& rJointWidth, bool& rComputeJointWidth, BoundedMatrix<double,3,3>& rRotationMatrix,
                                                            const double& InitialJointWidth, const Element::GeometryType& Geom)
{
    //Quadrilateral_interface_3d_4
    array_1d<double, 3> pmid0;
    array_1d<double, 3> pmid1;
    array_1d<double,3> P2 = Geom.GetPoint( 2 );
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
        rJointWidth = InitialJointWidth;
        rComputeJointWidth = false;
    }
}

//----------------------------------------------------------------------------------------
    
template< >
void UPlFaceLoadInterfaceCondition<2,2>::CalculateJointWidth( double& rJointWidth, const BoundedMatrix<double,2,4>& Nu,
                                                                        const array_1d<double,4>& DisplacementVector, array_1d<double,2>& rRelDispVector,
                                                                        const BoundedMatrix<double,2,2>& RotationMatrix,
                                                                        array_1d<double,2>& rLocalRelDispVector, const double& InitialJointWidth,
                                                                        const unsigned int& GPoint )
{
    //Line_interface_2d_2
    noalias(rRelDispVector) = prod(Nu,DisplacementVector);
    
    noalias(rLocalRelDispVector) = prod(RotationMatrix,rRelDispVector);

    rJointWidth = mInitialGap[GPoint] + rLocalRelDispVector[0]; //The joint width is obtained in the local x direction
    
    // JointWidth can't be negative
    if (rJointWidth < 0.0) {
        rJointWidth = 0.0;
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadInterfaceCondition<3,4>::CalculateJointWidth( double& rJointWidth, const BoundedMatrix<double,3,12>& Nu,
                                                                        const array_1d<double,12>& DisplacementVector, array_1d<double,3>& rRelDispVector,
                                                                        const BoundedMatrix<double,3,3>& RotationMatrix,
                                                                        array_1d<double,3>& rLocalRelDispVector, const double& InitialJointWidth,
                                                                        const unsigned int& GPoint )
{
    //Quadrilateral_interface_3d_4
    noalias(rRelDispVector) = prod(Nu,DisplacementVector);
    
    noalias(rLocalRelDispVector) = prod(RotationMatrix,rRelDispVector);

    rJointWidth = mInitialGap[GPoint] + rLocalRelDispVector[1]; //The joint width is obtained in the local y direction

    // JointWidth can't be negative
    if (rJointWidth < 0.0) {
        rJointWidth = 0.0;
    }
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadInterfaceCondition<2,2>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight, const double& JointWidth)
{
    // Note: since we cannot include the determinant of the Jacobian here (because the Jacobian could be singular), we do not multiply the JointWidth by the Weight 
    //       In a normal line load we would have |J| * w = L/2.0 * 2.0 = L
    rIntegrationCoefficient = JointWidth;
}

//----------------------------------------------------------------------------------------

template< >
void UPlFaceLoadInterfaceCondition<3,4>::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, const Matrix& Jacobian, const double& Weight, const double& JointWidth)
{
    KRATOS_TRY
    
    double dx_dxi = Jacobian(0,0);
    
    double dy_dxi = Jacobian(1,0);
    
    double dz_dxi = Jacobian(2,0);
    
    double ds = sqrt(dx_dxi*dx_dxi + dy_dxi*dy_dxi + dz_dxi*dz_dxi);
    
    rIntegrationCoefficient = Weight * ds * JointWidth;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlFaceLoadInterfaceCondition<2,2>;
template class UPlFaceLoadInterfaceCondition<3,4>;

} // Namespace Kratos.
