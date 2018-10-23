//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                 -0.1 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_updated_lagrangian_U_J_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"


namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   AxisymUpdatedLagrangianUJElement::AxisymUpdatedLagrangianUJElement()
      : UpdatedLagrangianUJElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJElement::AxisymUpdatedLagrangianUJElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUJElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJElement::AxisymUpdatedLagrangianUJElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUJElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJElement::AxisymUpdatedLagrangianUJElement( AxisymUpdatedLagrangianUJElement const& rOther)
      :UpdatedLagrangianUJElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJElement&  AxisymUpdatedLagrangianUJElement::operator=(AxisymUpdatedLagrangianUJElement const& rOther)
   {
      UpdatedLagrangianUJElement::operator=(rOther);

      mDeformationGradientF0.clear();
      mDeformationGradientF0.resize(rOther.mDeformationGradientF0.size());

      for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
         mDeformationGradientF0[i] = rOther.mDeformationGradientF0[i];
      }

      mDeterminantF0 = rOther.mDeterminantF0;

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUJElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new AxisymUpdatedLagrangianUJElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUJElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      AxisymUpdatedLagrangianUJElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

      //-----------//

      NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


      if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
         NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

         if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }

      for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
         NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

      //-----------//

      if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
         NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

      for(unsigned int i=0; i<mDeformationGradientF0.size(); i++)
      {
         NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
      }

      NewElement.mDeterminantF0 = mDeterminantF0;

      NewElement.SetData(this->GetData());
      NewElement.SetFlags(this->GetFlags());

      return Element::Pointer( new AxisymUpdatedLagrangianUJElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUJElement::~AxisymUpdatedLagrangianUJElement()
   {
   }


   //************* GETTING METHODS
   //************************************************************************************
   //************************************************************************************



   //************************************************************************************
   //************************************************************************************

   int  AxisymUpdatedLagrangianUJElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = UpdatedLagrangianUJElement::Check(rCurrentProcessInfo);


      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-J element type ", " AxisymUpdatedLagrangianUJElement" )

            //verify that the variables are correctly initialized

            if ( JACOBIAN.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "JACOBIAN has Key zero! (check if the application is correctly registered", "" )


      if ( this->GetProperties().Has(THICKNESS) ) {
	      double thickness = this->GetProperties()[THICKNESS];
	      if ( thickness <= 0.0) {
		      this->GetProperties()[THICKNESS] = 1.0;
	      }
      } else {
	     this->GetProperties()[THICKNESS] = 1.0;
      } 

      return correct;

      KRATOS_CATCH( "" );
   }


   //************* STARTING - ENDING  METHODS
   //************************************************************************************
   //************************************************************************************
   void AxisymUpdatedLagrangianUJElement::Initialize()
   {
      KRATOS_TRY

      UpdatedLagrangianUJElement::Initialize();

      SizeType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

      //Resize historic deformation gradient
      if ( mDeformationGradientF0.size() != integration_points_number )
         mDeformationGradientF0.resize( integration_points_number );

      if ( mDeterminantF0.size() != integration_points_number )
         mDeterminantF0.resize( integration_points_number, false );

      for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {
         mDeterminantF0[PointNumber] = 1;
         mDeformationGradientF0[PointNumber] = identity_matrix<double> (3);
      }

      const unsigned int number_of_nodes = GetGeometry().size();
      for (unsigned int Node = 0; Node < number_of_nodes; Node++) {
         double& DetFNodal = GetGeometry()[Node].GetSolutionStepValue(JACOBIAN );
         DetFNodal = 1.0;
         double& DetFNodalPrev = GetGeometry()[Node].GetSolutionStepValue(JACOBIAN , 1);
         DetFNodalPrev = 1.0;
      }

      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::InitializeElementData (ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {
      // copy from a non-derived class
      const unsigned int number_of_nodes = GetGeometry().size();

      rVariables.detF  = 1;

      rVariables.detF0 = 1;

      rVariables.detH = 1;

      rVariables.B.resize( 4 , number_of_nodes * 2 );

      rVariables.F.resize( 3, 3 );

      rVariables.F0.resize( 3, 3 );

      rVariables.H.resize( 3, 3 );

      rVariables.ConstitutiveMatrix.resize( 4, 4 );

      rVariables.StrainVector.resize( 4 );

      rVariables.StressVector.resize( 4 );

      rVariables.DN_DX.resize( number_of_nodes, 2 );

      //set variables including all integration points values

      //reading shape functions
      rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

      //reading shape functions local gradients
      rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

      //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
      rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );


      //Calculate Delta Position
      rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

      //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
      rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    //stabilization factor
    double StabilizationFactor = 0.20;
    if( GetProperties().Has(STABILIZATION_FACTOR_J) ) {
        StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_J];
    }
    else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR_J) ) {
        StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR_J];
    }
    GetProperties().SetValue(STABILIZATION_FACTOR_J, StabilizationFactor);


   }

   //************************************************************************************
   //************************************************************************************


   void AxisymUpdatedLagrangianUJElement::CalculateDeformationMatrix(Matrix& rB,
         Matrix& rDN_DX,
         Vector& rN,
         double & rCurrentRadius)
   {
      KRATOS_TRY


      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      rB.clear(); //set all components to zero

      if( dimension == 2 )
      {

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rN[i]/rCurrentRadius;
            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

         }

      }
      else if( dimension == 3 )
      {

         std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

      }
      else
      {

         KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

      }

      KRATOS_CATCH( "" )
   }



   //*************************COMPUTE DEFORMATION GRADIENT*******************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::CalculateDeformationGradient(const Matrix& rDN_DX,
         Matrix& rF,
         Matrix& rDeltaPosition, 
         double & rCurrentRadius, 
         double & rReferenceRadius)
   {

      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

      rF = identity_matrix<double> ( 3 );

      if( dimension == 2 )
      {

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
         }

         rF ( 2 , 2 ) = rCurrentRadius/rReferenceRadius;
      }
      else if( dimension == 3)
      {

         std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
      }
      else
      {

         KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

      }


      KRATOS_CATCH( "" )
   }
//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianUJElement::CalculateGreenLagrangeStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    //Right Cauchy-Green Calculation
    Matrix C ( 3, 3 );
    noalias( C ) = prod( trans( rF ), rF );

    if( dimension == 2 )
    {

        //Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 4 ) rStrainVector.resize( 4, false );

        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );

        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );

        rStrainVector[2] = 0.5 * ( C( 2, 2 ) - 1.00 );

        rStrainVector[3] = C( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void AxisymUpdatedLagrangianUJElement::CalculateAlmansiStrain(const Matrix& rF,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen(rF.size1(),rF.size1());
    noalias(LeftCauchyGreen) = prod( rF, trans( rF ) );

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen(rF.size1(),rF.size1());
    double det_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    if( dimension == 2 )
    {

        //Almansi Strain Calculation
        if ( rStrainVector.size() != 4 ) rStrainVector.resize( 4, false );

        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );

        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );

        rStrainVector[2] = 0.5 * ( 1.00 - InverseLeftCauchyGreen( 2, 2 ) );

        rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy

    }
    else if( dimension == 3 )
    {

        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }


    KRATOS_CATCH( "" )
}



   //************* COMPUTING  METHODS
   //************************************************************************************
   //************************************************************************************


   //*********************************COMPUTE KINEMATICS*********************************
   //************************************************************************************
   void AxisymUpdatedLagrangianUJElement::CalculateKinematics(ElementDataType& rVariables,
         const double& rPointNumber)

   {
      KRATOS_TRY

      // copy of a non-derived class

      //Get the parent coodinates derivative [dN/d£]
      const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

      //Get the shape functions for the order of the integration method [N]
      const Matrix& Ncontainer = rVariables.GetShapeFunctions();

      //Parent to reference configuration
      rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

      //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
      Matrix InvJ;
      MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

      //Compute cartesian derivatives [dN/dx_n]
      noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

      //Set Shape Functions Values for this integration point
      noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

      //Calculate IntegrationPoint radius
      CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius, rVariables.N);

      //Current Deformation Gradient [dx_n+1/dx_n]
      CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition, rVariables.CurrentRadius, rVariables.ReferenceRadius);

      //Determinant of the deformation gradient F
      rVariables.detF  = MathUtils<double>::Det(rVariables.F);

      //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
      Matrix Invj;
      MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

      //Compute cartesian derivatives [dN/dx_n+1]
      noalias( rVariables.DN_DX ) =  prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

      //Determinant of the Deformation Gradient F0
      rVariables.detF0 = mDeterminantF0[rPointNumber];
      rVariables.F0    = mDeformationGradientF0[rPointNumber];

      //Compute the deformation matrix B
      CalculateDeformationMatrix(rVariables.B, rVariables.DN_DX, rVariables.N, rVariables.CurrentRadius);


      KRATOS_CATCH( "" )
   }

   //*************************COMPUTE AXYSIMMETRIC RADIUS********************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::CalculateRadius(double & rCurrentRadius,
         double & rReferenceRadius,
         const Vector& rN)


   {

      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      rCurrentRadius=0;
      rReferenceRadius=0;

      if ( dimension == 2 )
      {
         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
            array_1d<double, 3 > & CurrentPosition      = GetGeometry()[i].Coordinates();
            array_1d<double, 3 > ReferencePosition      = CurrentPosition - DeltaDisplacement;

            rCurrentRadius   += CurrentPosition[0]*rN[i];
            rReferenceRadius += ReferencePosition[0]*rN[i];
            //std::cout<<" node "<<i<<" -> DeltaDisplacement : "<<DeltaDisplacement<<std::endl;
         }
      }


      if ( dimension == 3 )
      {
         std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
      }

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];


      UpdatedLagrangianUJElement::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight);

   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {

      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      UpdatedLagrangianUJElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, IntegrationWeight);

   }


   //******************************** JACOBIAN FORCES  **********************************
   //************************************************************************************
   void AxisymUpdatedLagrangianUJElement::CalculateAndAddJacobianForces(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            const double& NodalJac = (GetGeometry()[j].GetSolutionStepValue( JACOBIAN) );

            rRightHandSideVector[indexp] +=   rVariables.N[i] * rVariables.N[j]  * NodalJac * rIntegrationWeight / rVariables.detH;

         }

         rRightHandSideVector[indexp] -= rVariables.N[i] * rVariables.detH * rIntegrationWeight / rVariables.detH;



         indexp += (dimension + 1);

      }

      KRATOS_CATCH( "" )

   }



   //****************** STABILIZATION *********************************************************
   //************************* defined in the Stab element ************************************

   void AxisymUpdatedLagrangianUJElement::CalculateAndAddStabilizedJacobian(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      // VectorType Fh=rRightHandSideVector;
      // std::cout<<" Element "<<this->Id()<<" "<<std::endl;

      //use of this variable for the complete parameter: (deffault: 4)
      double AlphaStabilization  = 4.0; 
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_J];
      AlphaStabilization *= StabilizationFactor; 

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus= YoungModulus/(3*(1-2*PoissonCoefficient));

      AlphaStabilization=(AlphaStabilization/(18.0*LameMu));

      AlphaStabilization *= BulkModulus;  // TIMES THE BULK MODULUS BECAUSE I HAVE ALL THE EQUATION MULTIPLIED BY THE BULK MODULUS

      if (YoungModulus < 0.00001)
      {
         AlphaStabilization = 4.0 * StabilizationFactor / 18.0;
         AlphaStabilization *= mElementStabilizationNumber; 

      }

      double consistent = 1;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            consistent=(-1.0)*AlphaStabilization;
            if(i==j)
               consistent=2.0*AlphaStabilization;

            double& Jacobian= GetGeometry()[j].FastGetSolutionStepValue(JACOBIAN);
            rRightHandSideVector[indexp] += consistent * Jacobian * rIntegrationWeight / (rVariables.detH);

            // std::cout<<" Pressure "<<Pressure<<std::endl;
         }


         indexp += (dimension + 1);
      }


      // std::cout<<std::endl;
      // std::cout<<" IntegrationWeight "<<rIntegrationWeight<<" detF "<<rVariables.detF0<<std::endl;
      // std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;

      KRATOS_CATCH( "" )

   }


   //******** Kuu Material************************************************************
   //***************** It includes the pw geometric stiffness ************************

   void AxisymUpdatedLagrangianUJElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //assemble into rk the material uu contribution:
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      double dimension_double = 3.0;

      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;

      unsigned int voigtsize = 4;

      Matrix DeviatoricTensor = ZeroMatrix(voigtsize,voigtsize);
      Vector Identity = ZeroVector(voigtsize);

      for (unsigned int i = 0; i < voigtsize ; ++i) {
         DeviatoricTensor(i,i) = 1.0;
      }
      for (unsigned int i = 0; i < 3; i++) {
         Identity(i) = 1.0;
         for (unsigned int j = 0; j < 3; j++) {
            DeviatoricTensor( i,j) -= 1.0/dimension_double;
         }
      }


      ConstitutiveMatrix = prod( ConstitutiveMatrix, DeviatoricTensor);

      if ( this->Id() == 0) 
      {
         std::cout << " CONS 0 " << rVariables.ConstitutiveMatrix << std::endl;
         std::cout << " CONS 1 " << ConstitutiveMatrix << std::endl;
      }
      Matrix AuxMatrix = ZeroMatrix(voigtsize,voigtsize);

      for (unsigned int i = 0; i < voigtsize; i++) {
         for (unsigned int j = 0; j < voigtsize; j++) {
            ConstitutiveMatrix(i,j)  += (1 -  2/dimension_double) * rVariables.StressVector(i) * Identity(j);
            AuxMatrix(i,j) += rVariables.StressVector(i) * Identity(j);
         }
      }

      if ( this->Id() == 0) 
      {
         std::cout << " CONS 2 " << ConstitutiveMatrix << std::endl;
         std::cout << " AUX MATRIX " << AuxMatrix << std::endl;
         std::cout << " stressVector " << rVariables.StressVector << std::endl;
         std::cout << std::endl;
      }


      Matrix Kuu = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( ConstitutiveMatrix, rVariables.B ) ) ); 

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
      //if ( this->Id() < 10) 
      //   std::cout<<" Kmat "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::CalculateAndAddKuJ (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         double& rIntegrationWeight)
   {

      
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      double dimension_double = 3.0;

      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      unsigned int voigtsize = 4;

    // Trying to do it new
    Vector Identity = ZeroVector(voigtsize);
    for (unsigned int i = 0; i < 3; i++)
        Identity(i) = 1.0;

    Vector ConstVector = prod( ConstitutiveMatrix, Identity);
    ConstVector /= dimension_double;

    ConstVector += ( 2.0/dimension_double-1.0) * rVariables.StressVector;

    double ElementJacobian = 0.0;

      for ( unsigned int i = 0; i <  number_of_nodes ; i++)
         ElementJacobian += GetGeometry()[i].GetSolutionStepValue( JACOBIAN ) * rVariables.N[i] ;

    ConstVector /= ElementJacobian;

    Vector KuJ = prod( trans( rVariables.B), (ConstVector) );

          Matrix SecondMatrix = ZeroMatrix( dimension*number_of_nodes, number_of_nodes);

    for (unsigned int i = 0; i < dimension*number_of_nodes; i++) {
        for (unsigned int j = 0; j < number_of_nodes; j++) {
            SecondMatrix(i,j) =KuJ(i) * rVariables.N[j];
         }
      }
      SecondMatrix *= rIntegrationWeight;


      // Add the matrix in its place
      MatrixType Kh=rLeftHandSideMatrix;
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        for (unsigned int idim = 0; idim < dimension; idim++) {
            for (unsigned int j = 0; j < number_of_nodes; j++) {
                rLeftHandSideMatrix(i*(dimension+1) + idim, (dimension+1)*(j+1) -1 ) += SecondMatrix( i*(dimension) + idim, j);
               }
         }
      }

      //std::cout << std::endl;
      //std::cout << " TRUE MATRIX " << rLeftHandSideMatrix - Kh << std::endl;
      //std::cout << std::endl;

      KRATOS_CATCH( "" )
   }



   //******************* Kuug ********************************************************
   //*********************************************************************************

   void AxisymUpdatedLagrangianUJElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      int size = number_of_nodes * dimension;

      Matrix Kuu = zero_matrix<double>(size,size);

      // axisymmetric geometric matrix

      double alpha1 = 0;
      double alpha2 = 0;
      double alpha3 = 0;

      unsigned int indexi = 0;
      unsigned int indexj = 0;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         indexj =0;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            alpha1 = rVariables.DN_DX(j,0) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[0] + rVariables.DN_DX(i,1) * rVariables.StressVector[3] );
            alpha2 = rVariables.DN_DX(j,1) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[3] + rVariables.DN_DX(i,1) * rVariables.StressVector[1] );
            alpha3 = rVariables.N[i] * rVariables.N[j] * rVariables.StressVector[2] * (1.0/rVariables.CurrentRadius*rVariables.CurrentRadius);

            Kuu(indexi,indexj)     = alpha1 + alpha2 + alpha3 ;
            Kuu(indexi+1,indexj+1) = alpha1 + alpha2 ;

            indexj+=2;
         }

         indexi+=2;

      }

      Kuu *= rIntegrationWeight;

      //std::cout<<std::endl;
      //std::cout<<" Kuu "<<Kuu<<std::endl;


      MatrixType Kh=rLeftHandSideMatrix;

      //assemble into rk the geometric uu contribution:
      indexi = 0;
      indexj = 0;
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

      //std::cout<<std::endl;
      //std::cout<<" Kgeo "<<rK-Kh<<std::endl;

      KRATOS_CATCH( "" )

   
   }



   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::CalculateAndAddKJu (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY



      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexp = dimension;


      for (unsigned int i = 0; i < number_of_nodes; i++)
      {
         for (unsigned int j = 0; j < number_of_nodes; j++)
         {
            int indexup = dimension*j + j;
            for (unsigned int k = 0; k < dimension; k++)
            {
               rLeftHandSideMatrix(indexp, indexup + k) += rVariables.N[i] * rVariables.DN_DX( j, k) * rIntegrationWeight ;
               if ( k == 0)
                  rLeftHandSideMatrix(indexp, indexup + k) += rVariables.N[i] * rVariables.N[j] * ( 1.0/ rVariables.CurrentRadius) * rIntegrationWeight; 
            }

         }
         indexp += (dimension + 1);
      }


      //std::cout<<std::endl;
      //std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;


      KRATOS_CATCH( "" )
   }



   // ^^^^^^^^^^^^^^^^^^^^^ KJJ ***************************************************
   // ********************************************************************************
   void AxisymUpdatedLagrangianUJElement::CalculateAndAddKJJ (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      Matrix TotalF = prod( rVariables.F, rVariables.F0);

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            rLeftHandSideMatrix(indexpi,indexpj)  -= rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / rVariables.detH;
            indexpj += (dimension+1);
         }

         indexpi += (dimension + 1);
      }

      KRATOS_CATCH( "" )
   }





   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::CalculateAndAddKJJStab (MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         double& rIntegrationWeight)
   {

      KRATOS_TRY

      //repasar

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;
      double consistent = 1.0;

      //use of this variable for the complete parameter: (deffault: 4)
      double AlphaStabilization  = 4.0; 
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_J];
      AlphaStabilization *= StabilizationFactor; 

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus= YoungModulus/(3*(1-2*PoissonCoefficient));

      AlphaStabilization=(AlphaStabilization/(18.0*LameMu));

      AlphaStabilization *= BulkModulus;  // TIMES THE BULK MODULUS BECAUSE I HAVE ALL THE EQUATION MULTIPLIED BY THE BULK MODULUS

      if (YoungModulus < 0.00001)
      {
         AlphaStabilization = 4.0 * StabilizationFactor / 18.0;
         AlphaStabilization *= mElementStabilizationNumber; 
      }

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent=(-1.0)*AlphaStabilization;
            if(indexpi==indexpj)
               consistent=2.0*AlphaStabilization;

            rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF);     //2D

            indexpj += (dimension + 1);
         }

         indexpi += (dimension + 1);
      }

      // std::cout<<std::endl;
      // std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )

   }

   //************************************************************************************
   //************************************************************************************
   // Calculate mass matrix. copied from AxisymUpdatedLagrangianUPElement
   void AxisymUpdatedLagrangianUJElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      //lumped
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      const unsigned int number_of_nodes = GetGeometry().size();
      unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

      if ( rMassMatrix.size1() != MatSize )
         rMassMatrix.resize( MatSize, MatSize, false );

      rMassMatrix = ZeroMatrix( MatSize, MatSize );

      // Not Lumped Mass Matrix (numerical integration):

      //reading integration points
      IntegrationMethod CurrentIntegrationMethod = mThisIntegrationMethod; //GeometryData::GI_GAUSS_2; //GeometryData::GI_GAUSS_1;

      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( CurrentIntegrationMethod  );

      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);


      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {

         //compute element kinematics
         this->CalculateKinematics( Variables, PointNumber );

         //getting informations for integration
         double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ * 2.0 * 3.141592654 * Variables.CurrentRadius;


         //compute point volume change
         double PointVolumeChange = 0;
         PointVolumeChange = this->CalculateVolumeChange( PointVolumeChange, Variables );

         double CurrentDensity = PointVolumeChange * GetProperties()[DENSITY];

         for ( unsigned int i = 0; i < number_of_nodes; i++ )
         {
            unsigned int indexupi = dimension * i + i;

            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
               unsigned int indexupj = dimension * j + j;

               for ( unsigned int k = 0; k < dimension; k++ )
               {
                  rMassMatrix( indexupi+k , indexupj+k ) += Variables.N[i] * Variables.N[j] * CurrentDensity * IntegrationWeight;
               }
            }
         }

      }

      // Lumped Mass Matrix:

      // double TotalMass = 0;

      // this->CalculateTotalMass( TotalMass, rCurrentProcessInfo );

      // Vector LumpFact = ZeroVector(number_of_nodes);

      // LumpFact = GetGeometry().LumpingFactors( LumpFact );

      // for ( unsigned int i = 0; i < number_of_nodes; i++ )
      // 	{
      // 	  double temp = LumpFact[i] * TotalMass;

      // 	  unsigned int indexup = dimension * i + i;

      // 	  for ( unsigned int j = 0; j < dimension; j++ )
      // 	    {
      // 	      rMassMatrix( indexup+j , indexup+j ) += temp;
      // 	    }
      // 	}

      // std::cout<<std::endl;
      // std::cout<<" Mass Matrix "<<rMassMatrix<<std::endl;

      KRATOS_CATCH( "" )
   }

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::GetHistoricalVariables( ElementDataType& rVariables, const double& rPointNumber )
   {
      UpdatedLagrangianUJElement::GetHistoricalVariables(rVariables,rPointNumber);

      //Deformation Gradient F0
      rVariables.detF0 = mDeterminantF0[rPointNumber];
      rVariables.F0    = mDeformationGradientF0[rPointNumber];
   }



   void AxisymUpdatedLagrangianUJElement::ComputeConstitutiveVariables(  ElementDataType& rVariables, Matrix& rFT, double& rDetFT)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().size();

      double dimension_double = 3.0; 

      rDetFT = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++) 
         rDetFT += GetGeometry()[i].GetSolutionStepValue( JACOBIAN ) * rVariables.N[i];

      rFT = rVariables.H;

      rFT *=  pow( rDetFT/ rVariables.detH, 1.0/ dimension_double );

      // COMPUTE THE EFFECT OF THE INTERPOLATION, LETS SEE
      std::vector< Matrix > EECCInverseDefGrad;
      ProcessInfo SomeProcessInfo;
      this->GetValueOnIntegrationPoints( INVERSE_DEFORMATION_GRADIENT, EECCInverseDefGrad, SomeProcessInfo);
    if(  EECCInverseDefGrad[0].size1() == 0)
        return;
      Matrix EECCInverseBig = EECCInverseDefGrad[0];
      Matrix EECCDefGradInverse = ZeroMatrix(3,3);

      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            EECCDefGradInverse(i,j) = EECCInverseBig(i,j);
         }
      }

      double det;
      Matrix EECCDefGrad;
      MathUtils<double>::InvertMatrix( EECCDefGradInverse, EECCDefGrad, det);

      double detF0 = 0;
      unsigned int step = 1;
      if ( this->Is(SolidElement::FINALIZED_STEP) ) 
         step = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++)
         detF0 += GetGeometry()[i].GetSolutionStepValue( JACOBIAN, step ) * rVariables.N[i];

      Matrix F0 = rVariables.F0;
      F0 *= pow( detF0 / rVariables.detF0, 1.0/ dimension_double );

      Matrix F0Inverse;
      MathUtils<double>::InvertMatrix( F0, F0Inverse, det);

      Matrix Update = prod( F0Inverse, EECCDefGrad);

      if ( this->Id() == 0)
      {
         std::cout << " TRY TO SEE WHAT I DID " << std::endl;
         std::cout << " CONSTITUTIVE INVERSE " << EECCInverseDefGrad[0] << std::endl;
         std::cout << "  CONSTITUTIVE " << EECCDefGrad << std::endl;
         std::cout << std::endl;
         std::cout << " FINALIZED ?: " << this->Is(SolidElement::FINALIZED_STEP) << std::endl;
         std::cout << " NODAL " << detF0 << std::endl;
         std::cout << " PREVIOUS DISPL F " << rVariables.F0 << std::endl;
         std::cout << "   MAYBE " << rVariables.H << std::endl;
         std::cout << " SO FINALLY F0 displ theta " << F0 << std::endl;
         std::cout << " and the Inverse is: " << F0Inverse << std::endl;
         std::cout << " UPDATE ?" << Update << std::endl;
         std::cout << " AND THE LAST ONE " << prod( rFT, Update) << std::endl;
         std::cout << std::endl;
      }

      // SO FINALLY I DO THAT
      rFT = prod( rFT, Update);


      KRATOS_CATCH( " " )
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUJElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianUJElement )
   }

   void AxisymUpdatedLagrangianUJElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianUJElement )
   }




}
