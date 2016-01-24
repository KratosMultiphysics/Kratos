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

	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    if ( GetProperties()[CONSTITUTIVE_LAW_POINTER] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW_POINTER]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element with ID ", this->Id() )

    
    //FIC member variables
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( (dimension==2 && number_of_nodes==4) || (dimension==3 && number_of_nodes==8)) // Quadrilaterals2D4N or Hexahedra3D8N
    {
        mExtrapolationMatrix = this->CalculateExtrapolationMatrix(dimension); 
    }

    mShapeFunctionsDerivativesContainer.resize(integration_points.size());
    
    GeometryType::ShapeFunctionsGradientsType DN_DlocalContainer;
    DN_DlocalContainer = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer;
    GetGeometry().Jacobian( JContainer, mThisIntegrationMethod );
    Matrix InvJ, GradNpT;
    double detJ;
    
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        MathUtils<double>::InvertMatrix(JContainer[PointNumber], InvJ, detJ);
        if(detJ<0)
            KRATOS_THROW_ERROR( std::invalid_argument," |J|<0 , |J| = ", detJ )
        
        //Calculate shape functions global gradients
        GradNpT = prod(DN_DlocalContainer[PointNumber] , InvJ);
        
        mShapeFunctionsDerivativesContainer[PointNumber] = GradNpT;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    GeometryType::ShapeFunctionsGradientsType DN_DlocalContainer;
    DN_DlocalContainer = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer;
    GetGeometry().Jacobian( JContainer, mThisIntegrationMethod );
    Matrix InvJ, GradNpT;
    double detJ;
    
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        MathUtils<double>::InvertMatrix(JContainer[PointNumber], InvJ, detJ);
        if(detJ<0)
            KRATOS_THROW_ERROR( std::invalid_argument," |J|<0 , |J| = ", detJ )
        
        //Calculate shape functions global gradients
        GradNpT = prod(DN_DlocalContainer[PointNumber] , InvJ);
        
        mShapeFunctionsDerivativesContainer[PointNumber] = GradNpT;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    GeometryType::ShapeFunctionsGradientsType DN_DlocalContainer;
    DN_DlocalContainer = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    GeometryType::JacobiansType JContainer;
    GetGeometry().Jacobian( JContainer, mThisIntegrationMethod );
    Matrix InvJ, GradNpT;
    double detJ;
    
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        MathUtils<double>::InvertMatrix(JContainer[PointNumber], InvJ, detJ);
        if(detJ<0)
            KRATOS_THROW_ERROR( std::invalid_argument," |J|<0 , |J| = ", detJ )
        
        //Calculate shape functions global gradients
        GradNpT = prod(DN_DlocalContainer[PointNumber] , InvJ);
        
        mShapeFunctionsDerivativesContainer[PointNumber] = GradNpT;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
Matrix SmallStrainUPwFICElement::CalculateExtrapolationMatrix(const unsigned int dimension)
{
    Matrix ExtrapolationMatrix;
    std::vector<Vector> NodeLocalCoordinates;
    
    if(dimension==2)
    {
        ExtrapolationMatrix = ZeroMatrix(4,4);
        NodeLocalCoordinates.resize(4);
        
        NodeLocalCoordinates[0] = ZeroVector(3);
        NodeLocalCoordinates[0][0] = -sqrt(3);
        NodeLocalCoordinates[0][1] = -sqrt(3);
        
        NodeLocalCoordinates[1] = ZeroVector(3);
        NodeLocalCoordinates[1][0] =  sqrt(3);
        NodeLocalCoordinates[1][1] = -sqrt(3);
        
        NodeLocalCoordinates[2] = ZeroVector(3);
        NodeLocalCoordinates[2][0] =  sqrt(3);
        NodeLocalCoordinates[2][1] =  sqrt(3);
        
        NodeLocalCoordinates[3] = ZeroVector(3);
        NodeLocalCoordinates[3][0] = -sqrt(3);
        NodeLocalCoordinates[3][1] =  sqrt(3);
        
        //Element nodes loop
        for(unsigned int i=0; i<4; i++)
        {
            const Element::GeometryType::CoordinatesArrayType NodeLocalCoord = NodeLocalCoordinates[i];
            
            ExtrapolationMatrix(i,0) = GetGeometry().ShapeFunctionValue(0,NodeLocalCoord);
            ExtrapolationMatrix(i,1) = GetGeometry().ShapeFunctionValue(1,NodeLocalCoord);
            ExtrapolationMatrix(i,2) = GetGeometry().ShapeFunctionValue(2,NodeLocalCoord);
            ExtrapolationMatrix(i,3) = GetGeometry().ShapeFunctionValue(3,NodeLocalCoord);
        }
    }
    else
    {
        ExtrapolationMatrix = ZeroMatrix(8,8);
        NodeLocalCoordinates.resize(8);
        
        NodeLocalCoordinates[0] = ZeroVector(3);
        NodeLocalCoordinates[0][0] = -sqrt(3);
        NodeLocalCoordinates[0][1] = -sqrt(3);
        NodeLocalCoordinates[0][2] = -sqrt(3);
        
        NodeLocalCoordinates[1] = ZeroVector(3);
        NodeLocalCoordinates[1][0] =  sqrt(3);
        NodeLocalCoordinates[1][1] = -sqrt(3);
        NodeLocalCoordinates[1][2] = -sqrt(3);
        
        NodeLocalCoordinates[2] = ZeroVector(3);
        NodeLocalCoordinates[2][0] =  sqrt(3);
        NodeLocalCoordinates[2][1] =  sqrt(3);
        NodeLocalCoordinates[2][2] = -sqrt(3);
        
        NodeLocalCoordinates[3] = ZeroVector(3);
        NodeLocalCoordinates[3][0] = -sqrt(3);
        NodeLocalCoordinates[3][1] =  sqrt(3);
        NodeLocalCoordinates[3][2] = -sqrt(3);
        
        NodeLocalCoordinates[4] = ZeroVector(3);
        NodeLocalCoordinates[4][0] = -sqrt(3);
        NodeLocalCoordinates[4][1] = -sqrt(3);
        NodeLocalCoordinates[4][2] =  sqrt(3);
        
        NodeLocalCoordinates[5] = ZeroVector(3);
        NodeLocalCoordinates[5][0] =  sqrt(3);
        NodeLocalCoordinates[5][1] = -sqrt(3);
        NodeLocalCoordinates[5][2] =  sqrt(3);
        
        NodeLocalCoordinates[6] = ZeroVector(3);
        NodeLocalCoordinates[6][0] =  sqrt(3);
        NodeLocalCoordinates[6][1] =  sqrt(3);
        NodeLocalCoordinates[6][2] =  sqrt(3);
        
        NodeLocalCoordinates[7] = ZeroVector(3);
        NodeLocalCoordinates[7][0] = -sqrt(3);
        NodeLocalCoordinates[7][1] =  sqrt(3);
        NodeLocalCoordinates[7][2] =  sqrt(3);
        
        //Element nodes loop
        for(unsigned int i=0; i<8; i++)
        {
            const Element::GeometryType::CoordinatesArrayType NodeLocalCoord = NodeLocalCoordinates[i];
            
            ExtrapolationMatrix(i,0) = GetGeometry().ShapeFunctionValue(0,NodeLocalCoord);
            ExtrapolationMatrix(i,1) = GetGeometry().ShapeFunctionValue(1,NodeLocalCoord);
            ExtrapolationMatrix(i,2) = GetGeometry().ShapeFunctionValue(2,NodeLocalCoord);
            ExtrapolationMatrix(i,3) = GetGeometry().ShapeFunctionValue(3,NodeLocalCoord);
            ExtrapolationMatrix(i,4) = GetGeometry().ShapeFunctionValue(4,NodeLocalCoord);
            ExtrapolationMatrix(i,5) = GetGeometry().ShapeFunctionValue(5,NodeLocalCoord);
            ExtrapolationMatrix(i,6) = GetGeometry().ShapeFunctionValue(6,NodeLocalCoord);
            ExtrapolationMatrix(i,7) = GetGeometry().ShapeFunctionValue(7,NodeLocalCoord);
        }
    }
    
    return ExtrapolationMatrix;
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::InitializeElementalVariables (ElementalVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    //Variables at all integration points
    rVariables.NContainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    GetGeometry().DeterminantOfJacobian(rVariables.detJContainer,mThisIntegrationMethod);

    //Variables computed at each integration point
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int voigtsize  = 3;
    if( dimension == 3 ) voigtsize  = 6;
    rVariables.B = ZeroMatrix( voigtsize, number_of_nodes * dimension );
    rVariables.StrainVector = ZeroVector( voigtsize );
    rVariables.ConstitutiveMatrix = ZeroMatrix( voigtsize, voigtsize );
    rVariables.StressVector = ZeroVector( voigtsize );

    //Needed parameters for consistency with the general constitutive law
    rVariables.detF  = 1;
    rVariables.F     = identity_matrix<double>(dimension);
    
    //Nodal variables
    this->InitializeNodalVariables(rVariables);
    
    //Properties variables
    this->InitializeProperties(rVariables);
    
    //ProcessInfo variables
    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
    rVariables.NewmarkCoefficient1 = rCurrentProcessInfo[GAMMA_NEWMARK]/(rCurrentProcessInfo[BETA_NEWMARK]*DeltaTime);
    rVariables.NewmarkCoefficient2 = 1/(rCurrentProcessInfo[THETA_NEWMARK]*DeltaTime);
    
    //FIC variables
    rVariables.ShearModulus = GetProperties()[YOUNG_MODULUS]/(2*(1+GetProperties()[POISSON_RATIO]));
    rVariables.SecondOrderStrainDerivatives = false;
    if(dimension==2)
    {
        rVariables.ElementLength = sqrt(4*GetGeometry().Area()/M_PI);
        
        if(number_of_nodes==4) // Quadrilaterals2D4N
        {
            rVariables.SecondOrderStrainDerivatives = true;
            this->ExtrapolateShapeFunctionsDerivatives(rVariables);
            rVariables.ShapeFunctionsSecondOrderDerivatives.resize(4);
            rVariables.StrainDerivativeTerm = ZeroMatrix(2,8);
        }
    }
    else
    {
        rVariables.ElementLength = pow((6*GetGeometry().Volume()/M_PI),(1.0/3.0));
        
        if(number_of_nodes==8) // Hexahedra3D8N
        {
            rVariables.SecondOrderStrainDerivatives = true;
            this->ExtrapolateShapeFunctionsDerivatives(rVariables);
            rVariables.ShapeFunctionsSecondOrderDerivatives.resize(8);
            rVariables.StrainDerivativeTerm = ZeroMatrix(3,24);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::ExtrapolateShapeFunctionsDerivatives(ElementalVariables& rVariables)
{
    const unsigned int  number_of_nodes = GetGeometry().size();
    const unsigned int  dimension       = GetGeometry().WorkingSpaceDimension();
    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    unsigned int node;
    
    Matrix GPShapeFunctionsDerivatives = ZeroMatrix(integration_points_number, dimension*number_of_nodes);
    
    for(unsigned int i=0; i<integration_points_number; i++)
    {
        for(unsigned int j=0; j<number_of_nodes; j++)
        {
            node = j*dimension;
            
            GPShapeFunctionsDerivatives(i,node)   = mShapeFunctionsDerivativesContainer[i](j,0);
            GPShapeFunctionsDerivatives(i,node+1) = mShapeFunctionsDerivativesContainer[i](j,1);
            if(dimension==3)
                GPShapeFunctionsDerivatives(i,node+2) = mShapeFunctionsDerivativesContainer[i](j,2);
        }
    }
    
    Matrix NodalShapeFunctionsDerivativesMatrix = prod(mExtrapolationMatrix, GPShapeFunctionsDerivatives);
    
    rVariables.NodalShapeFunctionsDerivatives.resize(number_of_nodes);
    
    if(dimension==2)
    {
        for(unsigned int i = 0; i < number_of_nodes; i++)
        {
            rVariables.NodalShapeFunctionsDerivatives[i] = ZeroVector(8);
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
            rVariables.NodalShapeFunctionsDerivatives[i] = ZeroVector(24);
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

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int node;
    
    //Setting the shape function vector
    rVariables.Np = row( rVariables.NContainer, PointNumber);
            
    //Setting the detJ
    rVariables.detJ = rVariables.detJContainer[PointNumber];
    
    //Setting the shape functions global gradients
    rVariables.GradNpT = mShapeFunctionsDerivativesContainer[PointNumber];
    
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
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int node;

    //Compute the voigt identity matrix Iv
    unsigned int voigtsize  = 3;
    if( dimension == 3 ) voigtsize  = 6;
    Matrix Iv = ZeroMatrix(voigtsize,voigtsize);
    Iv(0,0) = 1;
    Iv(1,1) = 1;
    Iv(2,2) = 0.5;
    if(dimension==3)
    {
        Iv(2,2) = 1;
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
        Matrix StrainDerivativeMatrix = rVariables.NewmarkCoefficient1*0.25*rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient*
                                        prod(rVariables.GradNpT,rVariables.StrainDerivativeTerm)*rVariables.IntegrationCoefficient;

        //Distribute strain derivative block matrix into the elemental matrix
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_nodes = GetGeometry().size();
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
        Matrix StressDerivativeTerm = this->CalculateStressDerivativeTerm(rVariables);
        
        double StabilizationParameter = rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient/(8*rVariables.ShearModulus);
        
        Matrix StressDerivativeMatrix = rVariables.NewmarkCoefficient1*StabilizationParameter/3*prod(rVariables.GradNpT,StressDerivativeTerm)*rVariables.IntegrationCoefficient;

        //Distribute stress derivative block matrix into the elemental matrix
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_nodes = GetGeometry().size();
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

Matrix SmallStrainUPwFICElement::CalculateStressDerivativeTerm(ElementalVariables& rVariables)
{
    //TODO: for the moment this is implemented only for linear elasticity
    
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int node;
    double D1 = rVariables.ConstitutiveMatrix(0,0), D2 = rVariables.ConstitutiveMatrix(0,1);
    
    Matrix StressDerivativeTerm = ZeroMatrix(dimension,dimension*number_of_nodes);
    
    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 2 * i;
            
            StressDerivativeTerm(0, node)   = (D1+D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][0];
            StressDerivativeTerm(1, node+1) = (D1+D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][1];
            StressDerivativeTerm(0, node+1) = (D1+D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][2];
            StressDerivativeTerm(1, node)   = StressDerivativeTerm(0, node+1);
        }
    }
    else
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 3 * i;
                
            StressDerivativeTerm(0, node)   = (D1+2*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][0];
            StressDerivativeTerm(1, node+1) = (D1+2*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][1];
            StressDerivativeTerm(2, node+2) = (D1+2*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][2];
            StressDerivativeTerm(0, node+1) = (D1+2*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][3];
            StressDerivativeTerm(1, node)   = StressDerivativeTerm(0, node+1);
            StressDerivativeTerm(1, node+2) = (D1+2*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][4];
            StressDerivativeTerm(2, node+1) = StressDerivativeTerm(1, node+2);
            StressDerivativeTerm(0, node+2) = (D1+2*D2)*rVariables.ShapeFunctionsSecondOrderDerivatives[i][5];
            StressDerivativeTerm(2, node)   = StressDerivativeTerm(0, node+2);
        }
    }
    
    return StressDerivativeTerm;
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwFICElement::CalculateAndAddPressureDerivativeMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    double StabilizationParameter = rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient/(8*rVariables.ShearModulus);
    
    Matrix PressureDerivativeMatrix = rVariables.NewmarkCoefficient2*StabilizationParameter*(rVariables.BiotCoefficient-2*rVariables.ShearModulus*rVariables.BiotModulusInverse/(3*rVariables.BiotCoefficient))*
                                      prod(rVariables.GradNpT,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;
    
    //Distribute pressure derivative block matrix into the elemental matrix
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
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
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
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
        Matrix StressDerivativeTerm = this->CalculateStressDerivativeTerm(rVariables);
        
        double StabilizationParameter = rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient/(8*rVariables.ShearModulus);
        
        Matrix StressDerivativeMatrix = StabilizationParameter/3*prod(rVariables.GradNpT,StressDerivativeTerm)*rVariables.IntegrationCoefficient;

        Vector StressDerivativeFlow = prod(StressDerivativeMatrix,rVariables.VelocityVector);
                
        //Distribute stress derivative block vector into the elemental vector
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
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
    double StabilizationParameter = rVariables.ElementLength*rVariables.ElementLength*rVariables.BiotCoefficient/(8*rVariables.ShearModulus);
    
    Matrix PressureDerivativeMatrix = StabilizationParameter*(rVariables.BiotCoefficient-2*rVariables.ShearModulus*rVariables.BiotModulusInverse/(3*rVariables.BiotCoefficient))*
                                      prod(rVariables.GradNpT,trans(rVariables.GradNpT))*rVariables.IntegrationCoefficient;
    
    Vector PressureDerivativeFlow = prod(PressureDerivativeMatrix,rVariables.PressureDtVector);

    //Distribute pressure derivative block vector into the elemental vector
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int Global_i;
    
    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;
        
        rRightHandSideVector[Global_i] -= PressureDerivativeFlow[i];
    }
}

} // Namespace Kratos
