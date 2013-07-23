//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/large_displacement_U_P_3D_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LargeDisplacementUP3DElement::LargeDisplacementUP3DElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacement3DElement( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  LargeDisplacementUP3DElement::LargeDisplacementUP3DElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacement3DElement( NewId, pGeometry, pProperties )
  {
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  LargeDisplacementUP3DElement::LargeDisplacementUP3DElement( LargeDisplacementUP3DElement const& rOther)
    :LargeDisplacement3DElement(rOther)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  LargeDisplacementUP3DElement&  LargeDisplacementUP3DElement::operator=(LargeDisplacementUP3DElement const& rOther)
  {
    LargeDisplacement3DElement::operator=(rOther);

    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer LargeDisplacementUP3DElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Element::Pointer( new LargeDisplacementUP3DElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
  }


  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  LargeDisplacementUP3DElement::~LargeDisplacementUP3DElement()
  {
  }


  //************* GETTING METHODS
  //************************************************************************************
  //************************************************************************************



  void LargeDisplacementUP3DElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
  {
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
	rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ));
      }
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
  {
    int number_of_nodes = GetGeometry().size();
    unsigned int element_size = number_of_nodes * 3 + number_of_nodes;

    if ( rResult.size() != element_size )
      rResult.resize( element_size, false );

    for ( int i = 0; i < number_of_nodes; i++ )
      {
	int index = i * 3;
	rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
	rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
	rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
	rResult[index + 3] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
      }

  }

  //*********************************DISPLACEMENT***************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::GetValuesVector( Vector& rValues, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    unsigned int j = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = i * dimension;
	array_1d<double, 3>& Displacement = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT, Step );
	for ( j = 0; j < Displacement.size(); j++ )
	  {
	    rValues[ index + j ] = Displacement[j];
	  }

       	rValues[ index + j ] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
	
      }
  }


  //************************************VELOCITY****************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::GetFirstDerivativesVector( Vector& rValues, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    unsigned int j = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = i * dimension;
	array_1d<double, 3>& Velocity = GetGeometry()[i].GetSolutionStepValue( VELOCITY, Step );
	for ( j = 0; j < Velocity.size(); j++ )
	  {
	    rValues[ index + j ] = Velocity[j];
	  }

	rValues[ index + j ] = 0;
      }
  }

  //*********************************ACCELERATION***************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::GetSecondDerivativesVector( Vector& rValues, int Step )
  {
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    unsigned int j = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int index = i * dimension;
	array_1d<double, 3>& Acceleration = GetGeometry()[i].GetSolutionStepValue( ACCELERATION, Step );
	for ( j = 0; j < Acceleration.size(); j++ )
	  {
	    rValues[ index + j ] = Acceleration[j];
	  }

	rValues[ index + j ] = 0;
      }

  }


  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************


  void LargeDisplacementUP3DElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							Flags& rCalculationFlags)

  {

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;
      
    if ( rCalculationFlags.Is(LargeDisplacement3DElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
	if ( rLeftHandSideMatrix.size1() != MatSize )
	  rLeftHandSideMatrix.resize( MatSize, MatSize, false );
	  
	noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }
      
      
    //resizing as needed the RHS
    if ( rCalculationFlags.Is(LargeDisplacement3DElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
      {
	if ( rRightHandSideVector.size() != MatSize )
	  rRightHandSideVector.resize( MatSize, false );
	  
	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
      }
  }



  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
  {

    //contributions to stiffness matrix calculated on the reference config
    
    // operation performed: add Km to the rLefsHandSideMatrix

    //respect to the current configuration n+1   
    CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    
    // operation performed: add Kg to the rLefsHandSideMatrix
    CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    
    // operation performed: add Kup to the rLefsHandSideMatrix
    CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    
    // operation performed: add Kpu to the rLefsHandSideMatrix
    CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    
    // operation performed: add Kpp to the rLefsHandSideMatrix
    CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    
    // operation performed: add Kpp Stab to the rLefsHandSideMatrix
    CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //KRATOS_WATCH(rLeftHandSideMatrix)
  }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
  {

    //contribution to external forces

    // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
    CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);
	    
    // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
    CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rIntegrationWeight);
	    
    //KRATOS_WATCH(rRightHandSideVector)
  }


  //************************************************************************************
  //************************************************************************************

   void LargeDisplacementUP3DElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
									GeneralVariables& rVariables,
									Vector& rVolumeForce,
									double& rIntegrationWeight)
								    
  {
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //VectorType Fh=rRightHandSideVector;

    double Fext=0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	int indexup = dimension * i + i;

	array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);
	
	Fext = 0;
	for ( unsigned int j = 0; j < dimension; j++ )
	  {
	    Fext = rIntegrationWeight * rVariables.N[i] * rVolumeForce[j];
	    rRightHandSideVector[indexup + j] += Fext;
	    ExternalForce[j] +=Fext;
	  }
      }


    // std::cout<<std::endl;
    // std::cout<<" Fext "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

   void LargeDisplacementUP3DElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
									GeneralVariables & rVariables,
									double& rIntegrationWeight
									)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //VectorType Fh=rRightHandSideVector;

    Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexup = dimension * i + i;
        unsigned int indexu  = dimension * i;

	array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_INTERNAL);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
	    InternalForce[j] -= InternalForces[indexu + j];
        }
    }

    // std::cout<<std::endl;
    // //std::cout<<" Stress "<<rVariables.StressVector<<std::endl;
    // std::cout<<" Fint "<<rRightHandSideVector-Fh<<std::endl;
    KRATOS_CATCH( "" )
      }
  

  //************************************************************************************
  //************************************************************************************
  
   void LargeDisplacementUP3DElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
									     GeneralVariables & rVariables,
									     double& rIntegrationWeight)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    VectorType Fh=rRightHandSideVector;

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    double consistent=1;
    
    //double auxiliar = 0.5*(rVariables.detF0*rVariables.detF0-1)/rVariables.detF0); //(J²-1)
    double auxiliar = (std::log(rVariables.detF0)/rVariables.detF0);  //(ln(J))

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

            consistent=1;
            if(i==j)
	      consistent=2;

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            //rRightHandSideVector[indexp] += (1.0/BulkModulus) * rVariables.N[i] * rVariables.N[j] * Pressure * rIntegrationWeight;
	    rRightHandSideVector[indexp] += consistent * (1.0/BulkModulus) * (1.0/12.0) * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //2D
	    //std::cout<<" Pressure ["<<j<<"] : "<<Pressure<<" rhs "<<std::endl;

	  }
	
        //rRightHandSideVector[indexp] -=  auxiliar * rVariables.N[i] * rIntegrationWeight / rVariables.detF;

        rRightHandSideVector[indexp] -=  auxiliar * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

	//std::cout<< " Mpres "<<rRightHandSideVector[indexp]<<" Ppres "<<auxiliar * rVariables.N[i] * rIntegrationWeight / rVariables.detF <<std::endl;

        indexp += (dimension + 1);
      }


    // std::cout<<std::endl;
    // std::cout<<" auxiliar " <<auxiliar<<" F0 "<<rVariables.detF0<<std::endl;
    // std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************

   void LargeDisplacementUP3DElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
									     GeneralVariables & rVariables,
										 double& rIntegrationWeight)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    VectorType Fh=rRightHandSideVector;

    double AlphaStabilization = 4.0; //GetProperties()[STABILIZATION];

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);    

    //use of this variable for the complete parameter:
    AlphaStabilization=(AlphaStabilization/(18.0*LameMu));

    
    double consistent = 1;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

            consistent=(-1)*AlphaStabilization;
            if(i==j)
	      consistent=2*AlphaStabilization;

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
	  }


        indexp += (dimension + 1);
      }


    // std::cout<<std::endl;
    // std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************
 //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
							    GeneralVariables& rVariables,
							    double& rIntegrationWeight)                                                       
  {
    KRATOS_TRY

      //contributions to stiffness matrix calculated on the reference config
      Matrix Kuu = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) ); //to be optimized to remove the temporary
    
    //assemble into rk the material uu contribution:
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType Kh=rLeftHandSideMatrix;

    unsigned int indexi = 0;
    unsigned int indexj  = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
	  {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
	      {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
		  {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuu(indexi,indexj);
                    indexj++;
		  }
	      }
            indexi++;
	  }
      }

    // std::cout<<std::endl;
    // std::cout<<" Kmat "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
      }




  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
							GeneralVariables& rVariables,
							double& rIntegrationWeight)
    
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    int size = number_of_nodes * dimension;

    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix ReducedKg = prod( rVariables.DN_DX,  rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized

    Matrix Kuu = zero_matrix<double> (size);
    MathUtils<double>::ExpandAndAddReducedMatrix( Kuu, ReducedKg, dimension );

    MatrixType Kh=rLeftHandSideMatrix;

    //assemble into rLeftHandSideMatrix the geometric uu contribution:
    unsigned int indexi = 0;
    unsigned int indexj = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
	  {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
	      {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
		  {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuu(indexi,indexj);
                    indexj++;
		  }
	      }
            indexi++;
	  }
      }

    // std::cout<<std::endl;
    // std::cout<<" Kgeo "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
							  GeneralVariables& rVariables,
							  double& rIntegrationWeight)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType Kh=rLeftHandSideMatrix;
    //contributions to stiffness matrix calculated on the reference configuration
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int indexp  = dimension;
        unsigned int indexup = dimension * i + i;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

            for ( unsigned int k = 0; k < dimension; k++ )
	      {
                rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight * rVariables.detF;
	      }
            indexp += (dimension + 1);
	  }
      }

    // std::cout<<std::endl;
    // std::cout<<" Kup "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
							GeneralVariables& rVariables,
							double& rIntegrationWeight)
    
  {
    KRATOS_TRY

      //repasar

      const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexp = dimension;

    //double auxiliar = (rVariables.detF0*rVariables.detF0 + 1)/(rVariables.detF0*rVariables.detF0); //(J²-1)
    double auxiliar = (1.0-std::log(rVariables.detF0))/(rVariables.detF0*rVariables.detF0);   //(ln(J))


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {
            int indexup= dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
	      {
		rLeftHandSideMatrix(indexp,indexup+k) +=  auxiliar  * rVariables.N[i] * rVariables.DN_DX ( j , k ) * rIntegrationWeight * rVariables.detF;

		//std::cout<<" value ("<<indexp<<","<<indexup+k<<") "<<(2*detF) * rN[i] * rDN_DX ( j , k ) * rIntegrationWeight<<std::endl;
	      }
	  }
        indexp += (dimension + 1);
      }


    // std::cout<<std::endl;
    // std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;


    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
							  GeneralVariables& rVariables,
							  double& rIntegrationWeight)
  {
    KRATOS_TRY

      
      const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;
    double consistent = 1.0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {
            consistent=1;
            if(indexpi==indexpj)
	      consistent=2;

            //rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(BulkModulus)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
            rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/12.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D

            indexpj += (dimension + 1);
	  }

        indexpi += (dimension + 1);
      }

    // std::cout<<std::endl;
    // std::cout<<" Kpp "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
      }



  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
							      GeneralVariables & rVariables,
							      double& rIntegrationWeight)
  {
    KRATOS_TRY

      //repasar

      const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;
    double consistent = 1.0;

    double AlphaStabilization = 4.0; //GetProperties()[STABILIZATION];

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);    


    //use of this variable for the complete parameter:
    AlphaStabilization=(AlphaStabilization/(18.0*LameMu));

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	unsigned int indexpj = dimension;
	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {
	    consistent=(-1)*AlphaStabilization;
	    if(indexpi==indexpj)
	      consistent=2*AlphaStabilization;

	    rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF);     //2D

	    indexpj += (dimension + 1);
	  }

	indexpi += (dimension + 1);
      }

    // std::cout<<std::endl;
    // std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
      }



  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = 0;
    this->CalculateTotalMass( TotalMass );

    if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

    Vector LumpFact;

    LumpFact = GetGeometry().LumpingFactors( LumpFact );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = LumpFact[i] * TotalMass;

        unsigned int indexup = dimension * i + i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rMassMatrix( indexup+j , indexup+j ) = temp;
        }
    }

    //std::cout<<std::endl;
    //std::cout<<" Mass Matrix "<<rMassMatrix<<std::endl;

    KRATOS_CATCH( "" )
      }

  //************************************************************************************
  //************************************************************************************

  void LargeDisplacementUP3DElement::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY
      unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rDampMatrix.size1() != MatSize )
      rDampMatrix.resize( MatSize, MatSize, false );

    noalias( rDampMatrix ) = ZeroMatrix( MatSize, MatSize );

    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************

  int  LargeDisplacementUP3DElement::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY
      
    int correct = 0;
      
    correct = LargeDisplacement3DElement::Check(rCurrentProcessInfo);

    //verify that the variables are correctly initialized

    if ( PRESSURE.Key() == 0 )
      KRATOS_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" );

    return correct;

    KRATOS_CATCH( "" );
  }

  //************************************************************************************
  //************************************************************************************


  void LargeDisplacementUP3DElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacement3DElement );
   }

  void LargeDisplacementUP3DElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacement3DElement );
  }


} // Namespace Kratos


