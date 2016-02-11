//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// System includes
#include <math.h>

/* Project includes */
#include "custom_elements/small_strain_U_Pw_FIC_element.hpp"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
SmallStrainUPwFICElement::SmallStrainUPwFICElement() : SmallStrainUPwElement() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SmallStrainUPwFICElement::SmallStrainUPwFICElement( IndexType NewId, GeometryType::Pointer pGeometry ) : SmallStrainUPwElement( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SmallStrainUPwFICElement::SmallStrainUPwFICElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : SmallStrainUPwElement( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
SmallStrainUPwFICElement::~SmallStrainUPwFICElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Element::Pointer SmallStrainUPwFICElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallStrainUPwFICElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::Initialize()
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const SizeType NumNodes = rGeom.size();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != NumGPoints )
        mConstitutiveLawVector.resize( NumGPoints );

    if ( GetProperties()[CONSTITUTIVE_LAW_POINTER] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW_POINTER]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), rGeom,row( rGeom.ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element with ID ", this->Id() )

    
    //FIC member variables
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();

    if( (dimension==2 && number_of_nodes==4) || (dimension==3 && number_of_nodes==8)) // Quadrilaterals2D4N or Hexahedra3D8N
    {
        this->CalculateExtrapolationMatrix(mExtrapolationMatrix,dimension); 
    }

    mDN_DXContainer.resize(NumGPoints,false);
    for(SizeType i = 0; i<NumGPoints; i++)
        (mDN_DXContainer[i]).resize(NumNodes,Dim,false);
    mdetJContainer.resize(NumGPoints,false);
    rGeom.ShapeFunctionsIntegrationPointsGradients(mDN_DXContainer,mdetJContainer,mThisIntegrationMethod);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    rGeom.ShapeFunctionsIntegrationPointsGradients(mDN_DXContainer,mdetJContainer,mThisIntegrationMethod);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    rGeom.ShapeFunctionsIntegrationPointsGradients(mDN_DXContainer,mdetJContainer,mThisIntegrationMethod);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
void SmallStrainUPwFICElement::CalculateExtrapolationMatrix(Matrix& rExtrapolationMatrix, const unsigned int dimension)
{
    const GeometryType& rGeom = GetGeometry();
    std::vector<Vector> NodeLocalCoordinates;
    
    if(dimension==2)
    {
        rExtrapolationMatrix.resize(4,4,false);
        NodeLocalCoordinates.resize(4);
        
        (NodeLocalCoordinates[0]).resize(3,false);
        NodeLocalCoordinates[0][0] = -sqrt(3);
        NodeLocalCoordinates[0][1] = -sqrt(3);
        NodeLocalCoordinates[0][2] = 0.0;
        
        (NodeLocalCoordinates[1]).resize(3,false);
        NodeLocalCoordinates[1][0] =  sqrt(3);
        NodeLocalCoordinates[1][1] = -sqrt(3);
        NodeLocalCoordinates[1][2] = 0.0;
        
        (NodeLocalCoordinates[2]).resize(3,false);
        NodeLocalCoordinates[2][0] =  sqrt(3);
        NodeLocalCoordinates[2][1] =  sqrt(3);
        NodeLocalCoordinates[2][2] =  0.0;
        
        (NodeLocalCoordinates[3]).resize(3,false);
        NodeLocalCoordinates[3][0] = -sqrt(3);
        NodeLocalCoordinates[3][1] =  sqrt(3);
        NodeLocalCoordinates[3][2] =  0.0;
        
        //Element nodes loop
        for(unsigned int i=0; i<4; i++)
        {
            const Element::GeometryType::CoordinatesArrayType NodeLocalCoord = NodeLocalCoordinates[i];
            
            rExtrapolationMatrix(i,0) = rGeom.ShapeFunctionValue(0,NodeLocalCoord);
            rExtrapolationMatrix(i,1) = rGeom.ShapeFunctionValue(1,NodeLocalCoord);
            rExtrapolationMatrix(i,2) = rGeom.ShapeFunctionValue(2,NodeLocalCoord);
            rExtrapolationMatrix(i,3) = rGeom.ShapeFunctionValue(3,NodeLocalCoord);
        }
    }
    else
    {
        rExtrapolationMatrix.resize(8,8,false);
        NodeLocalCoordinates.resize(8);
        
        (NodeLocalCoordinates[0]).resize(3,false);
        NodeLocalCoordinates[0][0] = -sqrt(3);
        NodeLocalCoordinates[0][1] = -sqrt(3);
        NodeLocalCoordinates[0][2] = -sqrt(3);
        
        (NodeLocalCoordinates[1]).resize(3,false);
        NodeLocalCoordinates[1][0] =  sqrt(3);
        NodeLocalCoordinates[1][1] = -sqrt(3);
        NodeLocalCoordinates[1][2] = -sqrt(3);
        
        (NodeLocalCoordinates[2]).resize(3,false);
        NodeLocalCoordinates[2][0] =  sqrt(3);
        NodeLocalCoordinates[2][1] =  sqrt(3);
        NodeLocalCoordinates[2][2] = -sqrt(3);
        
        (NodeLocalCoordinates[3]).resize(3,false);
        NodeLocalCoordinates[3][0] = -sqrt(3);
        NodeLocalCoordinates[3][1] =  sqrt(3);
        NodeLocalCoordinates[3][2] = -sqrt(3);
        
        (NodeLocalCoordinates[4]).resize(3,false);
        NodeLocalCoordinates[4][0] = -sqrt(3);
        NodeLocalCoordinates[4][1] = -sqrt(3);
        NodeLocalCoordinates[4][2] =  sqrt(3);
        
        (NodeLocalCoordinates[5]).resize(3,false);
        NodeLocalCoordinates[5][0] =  sqrt(3);
        NodeLocalCoordinates[5][1] = -sqrt(3);
        NodeLocalCoordinates[5][2] =  sqrt(3);
        
        (NodeLocalCoordinates[6]).resize(3,false);
        NodeLocalCoordinates[6][0] =  sqrt(3);
        NodeLocalCoordinates[6][1] =  sqrt(3);
        NodeLocalCoordinates[6][2] =  sqrt(3);
        
        (NodeLocalCoordinates[7]).resize(3,false);
        NodeLocalCoordinates[7][0] = -sqrt(3);
        NodeLocalCoordinates[7][1] =  sqrt(3);
        NodeLocalCoordinates[7][2] =  sqrt(3);
        
        //Element nodes loop
        for(unsigned int i=0; i<8; i++)
        {
            const Element::GeometryType::CoordinatesArrayType NodeLocalCoord = NodeLocalCoordinates[i];
            
            rExtrapolationMatrix(i,0) = rGeom.ShapeFunctionValue(0,NodeLocalCoord);
            rExtrapolationMatrix(i,1) = rGeom.ShapeFunctionValue(1,NodeLocalCoord);
            rExtrapolationMatrix(i,2) = rGeom.ShapeFunctionValue(2,NodeLocalCoord);
            rExtrapolationMatrix(i,3) = rGeom.ShapeFunctionValue(3,NodeLocalCoord);
            rExtrapolationMatrix(i,4) = rGeom.ShapeFunctionValue(4,NodeLocalCoord);
            rExtrapolationMatrix(i,5) = rGeom.ShapeFunctionValue(5,NodeLocalCoord);
            rExtrapolationMatrix(i,6) = rGeom.ShapeFunctionValue(6,NodeLocalCoord);
            rExtrapolationMatrix(i,7) = rGeom.ShapeFunctionValue(7,NodeLocalCoord);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::InitializeElementalVariables (ElementalVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType NumNodes = rGeom.size();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    
    //Variables at all integration points
    (rVariables.NContainer).resize(NumGPoints,NumNodes,false);
    rVariables.NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );

    (rVariables.Np).resize(NumNodes,false);
    (rVariables.GradNpT).resize(NumNodes,Dim,false);
    
    (rVariables.detJContainer).resize(NumGPoints,false);
    noalias(rVariables.detJContainer) = mdetJContainer;

    //Variables computed at each integration point
    unsigned int voigtsize  = 3;
    if( Dim == 3 ) voigtsize  = 6;
    (rVariables.B).resize(voigtsize, NumNodes * Dim, false);
    noalias(rVariables.B) = ZeroMatrix( voigtsize, NumNodes * Dim );
    (rVariables.StrainVector).resize(voigtsize,false);
    (rVariables.ConstitutiveMatrix).resize(voigtsize, voigtsize, false);
    (rVariables.StressVector).resize(voigtsize,false);

    //Needed parameters for consistency with the general constitutive law
    rVariables.detF  = 1.0;
    (rVariables.F).resize(Dim, Dim, false);
    noalias(rVariables.F) = identity_matrix<double>(Dim);
    
    //Nodal variables
    this->InitializeNodalVariables(rVariables);
    
    //Properties variables
    this->InitializeProperties(rVariables);
    
    //ProcessInfo variables
    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
    rVariables.NewmarkCoefficient1 = rCurrentProcessInfo[GAMMA_NEWMARK]/(rCurrentProcessInfo[BETA_NEWMARK]*DeltaTime);
    rVariables.NewmarkCoefficient2 = 1.0/(rCurrentProcessInfo[THETA_NEWMARK]*DeltaTime);
    
    //FIC variables
    rVariables.ShearModulus = GetProperties()[YOUNG_MODULUS]/(2.0*(1.0+GetProperties()[POISSON_RATIO]));
    rVariables.SecondOrderStrainDerivatives = false;
    if(Dim==2)
    {
        rVariables.ElementLength = sqrt(4.0*rGeom.Area()/KRATOS_M_PI);
        
        if(NumNodes==4) // Quadrilaterals2D4N
        {
            rVariables.SecondOrderStrainDerivatives = true;
            this->ExtrapolateShapeFunctionsDerivatives(rVariables);
            rVariables.ShapeFunctionsSecondOrderDerivatives.resize(4); //NumNodes
            for(SizeType i=0; i<4;i++)
                ((rVariables.ShapeFunctionsSecondOrderDerivatives)[i]).resize(3); //voigtsize
            (rVariables.StrainDerivativeTerm).resize(2,8, false);
        }
    }
    else
    {
        rVariables.ElementLength = pow((6.0*rGeom.Volume()/KRATOS_M_PI),(1.0/3.0));
        
        if(NumNodes==8) // Hexahedra3D8N
        {
            rVariables.SecondOrderStrainDerivatives = true;
            this->ExtrapolateShapeFunctionsDerivatives(rVariables);
            rVariables.ShapeFunctionsSecondOrderDerivatives.resize(8); //NumNodes
            for(SizeType i=0; i<8;i++)
                ((rVariables.ShapeFunctionsSecondOrderDerivatives)[i]).resize(6); //voigtsize
            (rVariables.StrainDerivativeTerm).resize(3,24, false);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::ExtrapolateShapeFunctionsDerivatives(ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const unsigned int  number_of_nodes = rGeom.size();
    const unsigned int  dimension       = rGeom.WorkingSpaceDimension();
    const unsigned int& integration_points_number = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    unsigned int node;
    
    Matrix GPShapeFunctionsDerivatives = ZeroMatrix(integration_points_number, dimension*number_of_nodes);
    
    for(unsigned int i=0; i<integration_points_number; i++)
    {
        for(unsigned int j=0; j<number_of_nodes; j++)
        {
            node = j*dimension;
            
            GPShapeFunctionsDerivatives(i,node)   = mDN_DXContainer[i](j,0);
            GPShapeFunctionsDerivatives(i,node+1) = mDN_DXContainer[i](j,1);
            if(dimension==3)
                GPShapeFunctionsDerivatives(i,node+2) = mDN_DXContainer[i](j,2);
        }
    }
    
    Matrix NodalShapeFunctionsDerivativesMatrix = prod(mExtrapolationMatrix, GPShapeFunctionsDerivatives);
    
    rVariables.NodalShapeFunctionsDerivatives.resize(number_of_nodes);
    
    if(dimension==2)
    {
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            (rVariables.NodalShapeFunctionsDerivatives[i]).resize(8,false);
            node = i*dimension;
            
            rVariables.NodalShapeFunctionsDerivatives[i][0] = NodalShapeFunctionsDerivativesMatrix(0,node);
            rVariables.NodalShapeFunctionsDerivatives[i][1] = NodalShapeFunctionsDerivativesMatrix(0,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][2] = NodalShapeFunctionsDerivativesMatrix(1,node);
            rVariables.NodalShapeFunctionsDerivatives[i][3] = NodalShapeFunctionsDerivativesMatrix(1,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][4] = NodalShapeFunctionsDerivativesMatrix(2,node);
            rVariables.NodalShapeFunctionsDerivatives[i][5] = NodalShapeFunctionsDerivativesMatrix(2,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][6] = NodalShapeFunctionsDerivativesMatrix(3,node);
            rVariables.NodalShapeFunctionsDerivatives[i][7] = NodalShapeFunctionsDerivativesMatrix(3,node+1);
        }
    }
    else
    {
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            (rVariables.NodalShapeFunctionsDerivatives[i]).resize(24,false);
            node = i*dimension;
            
            rVariables.NodalShapeFunctionsDerivatives[i][0] = NodalShapeFunctionsDerivativesMatrix(0,node);
            rVariables.NodalShapeFunctionsDerivatives[i][1] = NodalShapeFunctionsDerivativesMatrix(0,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][2] = NodalShapeFunctionsDerivativesMatrix(0,node+2);
            rVariables.NodalShapeFunctionsDerivatives[i][3] = NodalShapeFunctionsDerivativesMatrix(1,node);
            rVariables.NodalShapeFunctionsDerivatives[i][4] = NodalShapeFunctionsDerivativesMatrix(1,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][5] = NodalShapeFunctionsDerivativesMatrix(1,node+2);
            rVariables.NodalShapeFunctionsDerivatives[i][6] = NodalShapeFunctionsDerivativesMatrix(2,node);
            rVariables.NodalShapeFunctionsDerivatives[i][7] = NodalShapeFunctionsDerivativesMatrix(2,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][8] = NodalShapeFunctionsDerivativesMatrix(2,node+2);
            rVariables.NodalShapeFunctionsDerivatives[i][9] = NodalShapeFunctionsDerivativesMatrix(3,node);
            rVariables.NodalShapeFunctionsDerivatives[i][10] = NodalShapeFunctionsDerivativesMatrix(3,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][11] = NodalShapeFunctionsDerivativesMatrix(3,node+2);
            rVariables.NodalShapeFunctionsDerivatives[i][12] = NodalShapeFunctionsDerivativesMatrix(4,node);
            rVariables.NodalShapeFunctionsDerivatives[i][13] = NodalShapeFunctionsDerivativesMatrix(4,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][14] = NodalShapeFunctionsDerivativesMatrix(4,node+2);
            rVariables.NodalShapeFunctionsDerivatives[i][15] = NodalShapeFunctionsDerivativesMatrix(5,node);
            rVariables.NodalShapeFunctionsDerivatives[i][16] = NodalShapeFunctionsDerivativesMatrix(5,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][17] = NodalShapeFunctionsDerivativesMatrix(5,node+2);
            rVariables.NodalShapeFunctionsDerivatives[i][18] = NodalShapeFunctionsDerivativesMatrix(6,node);
            rVariables.NodalShapeFunctionsDerivatives[i][19] = NodalShapeFunctionsDerivativesMatrix(6,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][20] = NodalShapeFunctionsDerivativesMatrix(6,node+2);
            rVariables.NodalShapeFunctionsDerivatives[i][21] = NodalShapeFunctionsDerivativesMatrix(7,node);
            rVariables.NodalShapeFunctionsDerivatives[i][22] = NodalShapeFunctionsDerivativesMatrix(7,node+1);
            rVariables.NodalShapeFunctionsDerivatives[i][23] = NodalShapeFunctionsDerivativesMatrix(7,node+2);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateKinematics(ElementalVariables& rVariables, unsigned int PointNumber)

{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int node;
    
    //Setting the vector of shape functions and the matrix of the shape functions global gradients
    noalias(rVariables.Np) = row( rVariables.NContainer, PointNumber);
    noalias(rVariables.GradNpT) = mDN_DXContainer[PointNumber];
    
    //Compute the deformation matrix B
    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 2 * i;

            rVariables.B( 0, node + 0 ) = rVariables.GradNpT( i, 0 );
            rVariables.B( 1, node + 1 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 2, node + 0 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 2, node + 1 ) = rVariables.GradNpT( i, 0 );
        }
    }
    else
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 3 * i;

            rVariables.B( 0, node + 0 ) = rVariables.GradNpT( i, 0 );
            rVariables.B( 1, node + 1 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 2, node + 2 ) = rVariables.GradNpT( i, 2 );
            rVariables.B( 3, node + 0 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 3, node + 1 ) = rVariables.GradNpT( i, 0 );
            rVariables.B( 4, node + 1 ) = rVariables.GradNpT( i, 2 );
            rVariables.B( 4, node + 2 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 5, node + 0 ) = rVariables.GradNpT( i, 2 );
            rVariables.B( 5, node + 2 ) = rVariables.GradNpT( i, 0 );
        }
    }

    //Compute infinitessimal strain
    rVariables.StrainVector = prod(rVariables.B,rVariables.DisplacementVector);
    
    //Compute strain derivative terms
    if(rVariables.SecondOrderStrainDerivatives==true)
    {
        this->CalculateStrainDerivativeTerm(rVariables);
    }    

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateStrainDerivativeTerm(ElementalVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int node;

    //Compute the voigt identity matrix Iv
    unsigned int voigtsize  = 3;
    if( dimension == 3 ) voigtsize  = 6;
    Matrix Iv = ZeroMatrix(voigtsize,voigtsize);
    Iv(0,0) = 1.0;
    Iv(1,1) = 1.0;
    Iv(2,2) = 0.5;
    if(dimension==3)
    {
        Iv(2,2) = 1.0;
        Iv(3,3) = 0.5;
        Iv(4,4) = 0.5;
        Iv(5,5) = 0.5;
    }

    //Compute the shape function second order derivatives
    for(unsigned int i=0; i<number_of_nodes; i++)
    {
        rVariables.ShapeFunctionsSecondOrderDerivatives[i] = prod(Matrix(prod(Iv,rVariables.B)),rVariables.NodalShapeFunctionsDerivatives[i]);
    }
    
    //Compute the StrainDerivative term
    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 2 * i;
            
            rVariables.StrainDerivativeTerm(0, node)   = rVariables.ShapeFunctionsSecondOrderDerivatives[i][0]+0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][1];
            rVariables.StrainDerivativeTerm(1, node+1) = 0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][0]+rVariables.ShapeFunctionsSecondOrderDerivatives[i][1];
            rVariables.StrainDerivativeTerm(0, node+1) = 0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][2];
            rVariables.StrainDerivativeTerm(1, node)   = rVariables.StrainDerivativeTerm(0, node+1);
        }
    }
    else
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
        node = 3 * i;
            
        rVariables.StrainDerivativeTerm(0, node)   = rVariables.ShapeFunctionsSecondOrderDerivatives[i][0]+0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][1]+0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][2];
        rVariables.StrainDerivativeTerm(1, node+1) = 0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][0]+rVariables.ShapeFunctionsSecondOrderDerivatives[i][1]+0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][2];
        rVariables.StrainDerivativeTerm(2, node+2) = 0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][0]+0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][1]+rVariables.ShapeFunctionsSecondOrderDerivatives[i][2];
        rVariables.StrainDerivativeTerm(0, node+1) = 0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][3];
        rVariables.StrainDerivativeTerm(1, node)   = rVariables.StrainDerivativeTerm(0, node+1);
        rVariables.StrainDerivativeTerm(1, node+2) = 0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][4];
        rVariables.StrainDerivativeTerm(2, node+1) = rVariables.StrainDerivativeTerm(1, node+2);
        rVariables.StrainDerivativeTerm(0, node+2) = 0.5*rVariables.ShapeFunctionsSecondOrderDerivatives[i][5];
        rVariables.StrainDerivativeTerm(2, node)   = rVariables.StrainDerivativeTerm(0, node+2);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{    
    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);
    
    this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);
    
    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);
    
    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
    
    //FIC terms
    this->CalculateAndAddStrainDerivativeMatrix(rLeftHandSideMatrix,rVariables);
    
    this->CalculateAndAddStressDerivativeMatrix(rLeftHandSideMatrix,rVariables);
    
    this->CalculateAndAddPressureDerivativeMatrix(rLeftHandSideMatrix,rVariables);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddStrainDerivativeMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    if(rVariables.SecondOrderStrainDerivatives==true)
    {
        const GeometryType& rGeom = GetGeometry();
        
        Matrix StrainDerivativeMatrix = rVariables.NewmarkCoefficient1*0.25*rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient*
                                        prod(rVariables.GradNpT,rVariables.StrainDerivativeTerm)*rVariables.IntegrationCoefficient;

        //Distribute strain derivative block matrix into the elemental matrix
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        const unsigned int number_of_nodes = rGeom.size();
        unsigned int Global_i, Global_j, Local_j;
        
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            Global_i = i * (dimension + 1) + dimension;
            
            for(unsigned int j = 0; j < number_of_nodes; j++)
            {
                Global_j = j * (dimension + 1);
                Local_j = j * dimension;
                
                rLeftHandSideMatrix(Global_i,Global_j)   -= StrainDerivativeMatrix(i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1) -= StrainDerivativeMatrix(i,Local_j+1);
                if(dimension == 3)
                    rLeftHandSideMatrix(Global_i,Global_j+2) -= StrainDerivativeMatrix(i,Local_j+2);
            }
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddStressDerivativeMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    if(rVariables.SecondOrderStrainDerivatives==true)
    {   
        const GeometryType& rGeom = GetGeometry();
        
        Matrix StressDerivativeTerm;
        this->CalculateStressDerivativeTerm(StressDerivativeTerm,rVariables);
        
        double StabilizationParameter = rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient/(8.0*rVariables.ShearModulus);
        
        Matrix StressDerivativeMatrix = rVariables.NewmarkCoefficient1*StabilizationParameter/3.0*prod(rVariables.GradNpT,StressDerivativeTerm)*rVariables.IntegrationCoefficient;

        //Distribute stress derivative block matrix into the elemental matrix
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        const unsigned int number_of_nodes = rGeom.size();
        unsigned int Global_i, Global_j, Local_j;
        
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            Global_i = i * (dimension + 1) + dimension;
            
            for(unsigned int j = 0; j < number_of_nodes; j++)
            {
                Global_j = j * (dimension + 1);
                Local_j = j * dimension;
                
                rLeftHandSideMatrix(Global_i,Global_j)   -= StressDerivativeMatrix(i,Local_j);
                rLeftHandSideMatrix(Global_i,Global_j+1) -= StressDerivativeMatrix(i,Local_j+1);
                if(dimension == 3)
                    rLeftHandSideMatrix(Global_i,Global_j+2) -= StressDerivativeMatrix(i,Local_j+2);
            }
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateStressDerivativeTerm(Matrix& rStressDerivativeTerm, const ElementalVariables& rVariables)
{
    //TODO: for the moment this is implemented only for linear elasticity
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int node;
    double D1 = rVariables.ConstitutiveMatrix(0,0), D2 = rVariables.ConstitutiveMatrix(0,1);
    
    rStressDerivativeTerm.resize(dimension,dimension*number_of_nodes,false);
    noalias(rStressDerivativeTerm) = ZeroMatrix(dimension,dimension*number_of_nodes);
    
    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 2 * i;
            
            rStressDerivativeTerm(0, node)   = (D1+D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][0];
            rStressDerivativeTerm(1, node+1) = (D1+D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][1];
            rStressDerivativeTerm(0, node+1) = (D1+D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][2];
            rStressDerivativeTerm(1, node)   = rStressDerivativeTerm(0, node+1);
        }
    }
    else
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 3 * i;
                
            rStressDerivativeTerm(0, node)   = (D1+2.0*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][0];
            rStressDerivativeTerm(1, node+1) = (D1+2.0*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][1];
            rStressDerivativeTerm(2, node+2) = (D1+2.0*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][2];
            rStressDerivativeTerm(0, node+1) = (D1+2.0*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][3];
            rStressDerivativeTerm(1, node)   = rStressDerivativeTerm(0, node+1);
            rStressDerivativeTerm(1, node+2) = (D1+2.0*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][4];
            rStressDerivativeTerm(2, node+1) = rStressDerivativeTerm(1, node+2);
            rStressDerivativeTerm(0, node+2) = (D1+2.0*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][5];
            rStressDerivativeTerm(2, node)   = rStressDerivativeTerm(0, node+2);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddPressureDerivativeMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    double StabilizationParameter = rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient/(8.0*rVariables.ShearModulus);
    
    Matrix PressureDerivativeMatrix = rVariables.NewmarkCoefficient2*StabilizationParameter*(rVariables.BiotCoefficient-2.0*rVariables.ShearModulus*rVariables.BiotModulusInverse/(3.0*rVariables.BiotCoefficient))*
                                      prod(rVariables.GradNpT,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;
    
    //Distribute pressure derivative block matrix into the elemental matrix
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
            
            rLeftHandSideMatrix(Global_i,Global_j) += PressureDerivativeMatrix(i,j);
        }
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);
    
    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);
    
    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);
    
    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);
    
    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);
    
    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
    
    //FIC terms
    this->CalculateAndAddStrainDerivativeFlow(rRightHandSideVector,rVariables);
    
    this->CalculateAndAddStressDerivativeFlow(rRightHandSideVector,rVariables);
    
    this->CalculateAndAddPressureDerivativeFlow(rRightHandSideVector,rVariables);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddStrainDerivativeFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    if(rVariables.SecondOrderStrainDerivatives==true)
    {
        Matrix StrainDerivativeMatrix = 0.25*rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient*
                                        prod(rVariables.GradNpT,rVariables.StrainDerivativeTerm)*rVariables.IntegrationCoefficient;
        
        Vector StrainDerivativeFlow = prod(StrainDerivativeMatrix,rVariables.VelocityVector);
                
        //Distribute strain derivative block vector into the elemental vector
        const GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.size();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        unsigned int Global_i;
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            Global_i = i * (dimension + 1) + dimension;
            
            rRightHandSideVector[Global_i] += StrainDerivativeFlow[i];
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddStressDerivativeFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    if(rVariables.SecondOrderStrainDerivatives==true)
    {
        Matrix StressDerivativeTerm;
        this->CalculateStressDerivativeTerm(StressDerivativeTerm, rVariables);
        
        double StabilizationParameter = rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient/(8.0*rVariables.ShearModulus);
        
        Matrix StressDerivativeMatrix = StabilizationParameter/3.0*prod(rVariables.GradNpT,StressDerivativeTerm)*rVariables.IntegrationCoefficient;

        Vector StressDerivativeFlow = prod(StressDerivativeMatrix,rVariables.VelocityVector);
                
        //Distribute stress derivative block vector into the elemental vector
        const GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.size();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();
        unsigned int Global_i;
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            Global_i = i * (dimension + 1) + dimension;
            
            rRightHandSideVector[Global_i] += StressDerivativeFlow[i];
        }

    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddPressureDerivativeFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    double StabilizationParameter = rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient/(8.0*rVariables.ShearModulus);
    
    Matrix PressureDerivativeMatrix = StabilizationParameter*(rVariables.BiotCoefficient-2.0*rVariables.ShearModulus*rVariables.BiotModulusInverse/(3.0*rVariables.BiotCoefficient))*
                                      prod(rVariables.GradNpT,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;
    
    Vector PressureDerivativeFlow = prod(PressureDerivativeMatrix,rVariables.PressureDtVector);

    //Distribute pressure derivative block vector into the elemental vector
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int Global_i;
    
    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;
        
        rRightHandSideVector[Global_i] -= PressureDerivativeFlow[i];
    }
}

} // Namespace Kratos
