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
//#include "custom_elements/large_displacement_element.hpp"
//#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"
//#include "solid_mechanics_application.h"
#include "custom_elements/axisym_spatial_lagrangian_U_wP_element.hpp"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   AxisymSpatialLagrangianUwPElement::AxisymSpatialLagrangianUwPElement()
      : SpatialLagrangianUwPElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymSpatialLagrangianUwPElement::AxisymSpatialLagrangianUwPElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : SpatialLagrangianUwPElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymSpatialLagrangianUwPElement::AxisymSpatialLagrangianUwPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : SpatialLagrangianUwPElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   AxisymSpatialLagrangianUwPElement::AxisymSpatialLagrangianUwPElement( AxisymSpatialLagrangianUwPElement const& rOther)
      :SpatialLagrangianUwPElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   AxisymSpatialLagrangianUwPElement&  AxisymSpatialLagrangianUwPElement::operator=(AxisymSpatialLagrangianUwPElement const& rOther)
   {
      SpatialLagrangianUwPElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer AxisymSpatialLagrangianUwPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new AxisymSpatialLagrangianUwPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer AxisymSpatialLagrangianUwPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      AxisymSpatialLagrangianUwPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

      //-----------//

      NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


      if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
         NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

         if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
            KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
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

      return Element::Pointer( new AxisymSpatialLagrangianUwPElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   AxisymSpatialLagrangianUwPElement::~AxisymSpatialLagrangianUwPElement()
   {
   }



   //************************************************************************************
   //************************************************************************************

   void AxisymSpatialLagrangianUwPElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
   {

      const unsigned int number_of_nodes = GetGeometry().size();

      rVariables.detF  = 1;

      rVariables.detF0 = 1;

      rVariables.DomainSize = 1;

      rVariables.B.resize( 4 , number_of_nodes * 2 );

      rVariables.F.resize( 3, 3 );

      rVariables.F0.resize( 3, 3 );

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

      // SAVE THE TIME STEP, THAT WILL BE USED; BUT IS A BAD IDEA TO DO IT THIS WAY.
      mTimeStep = rCurrentProcessInfo[DELTA_TIME];

      //mCompressibleWater = false;

   }

   //************************************************************************************
   //************************************************************************************


   void AxisymSpatialLagrangianUwPElement::CalculateDeformationGradient(const Matrix& rDN_DX,
         Matrix&  rF,
         Matrix&  rDeltaPosition,
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

         KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

      }


      KRATOS_CATCH( "" )
   }



   //*************************COMPUTE DEFORMATION GRADIENT*******************************
   //************************************************************************************
   void AxisymSpatialLagrangianUwPElement::CalculateDeformationMatrix(Matrix& rB,
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

         KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

      }

      KRATOS_CATCH( "" )
   }


   //************* COMPUTING  METHODS
   //************************************************************************************
   //************************************************************************************
   //************************************************************************************
   //************************************************************************************

   //************************************************************************************
   void AxisymSpatialLagrangianUwPElement::SetGeneralVariables(GeneralVariables& rVariables,
         ConstitutiveLaw::Parameters& rValues,
         const int & rPointNumber)
   {
      SpatialLagrangianUwPElement::SetGeneralVariables(rVariables,rValues,rPointNumber);

      //Set extra options for the contitutive law
      Flags &ConstitutiveLawOptions=rValues.GetOptions();
      ConstitutiveLawOptions.Set(ConstitutiveLaw::FINAL_CONFIGURATION);

   }

   //*********************************COMPUTE KINEMATICS*********************************
   //************************************************************************************
   void AxisymSpatialLagrangianUwPElement::CalculateKinematics(GeneralVariables& rVariables,
         const double& rPointNumber)

   {
      KRATOS_TRY

         //Get the parent coodinates derivative [dN/d£]
         const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
      //Get the shape functions for the order of the integration method [N]
      const Matrix& Ncontainer = rVariables.GetShapeFunctions();

      //Parent to reference configuration
      rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

      //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
      Matrix InvJ;
      MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

      //Step domain size
      rVariables.DomainSize = rVariables.detJ;

      //Compute cartesian derivatives [dN/dx_n]
      noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

      //Set Shape Functions Values for this integration point
      rVariables.N=row( Ncontainer, rPointNumber);

      //Calculate IntegrationPoint radius
      CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius, rVariables.N);

      //Current Deformation Gradient [dx_n+1/dx_n]
      CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition, rVariables.CurrentRadius, rVariables.ReferenceRadius);

      //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
      Matrix Invj;
      MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

      //Compute cartesian derivatives [dN/dx_n+1]
      rVariables.DN_DX = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

      //Determinant of the Deformation Gradient F0
      rVariables.detF0 = mDeterminantF0[rPointNumber];
      rVariables.F0    = mDeformationGradientF0[rPointNumber];

      //Compute the deformation matrix B
      CalculateDeformationMatrix(rVariables.B, rVariables.DN_DX, rVariables.N, rVariables.CurrentRadius);


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
   {


      /* std::cout << " I HAVE IT ALL WRONG: detF0 " << 1.0 - rVariables.detF0 << " detF " << 1.0 - rVariables.detF << std::endl;
         std::cout << " PART 2 F  " << rVariables.F << std::endl;
         std::cout << " PART 3 F0 " << rVariables.F0 << std::endl;
         std::cout << std::endl;
       */
      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      // SOMETHING OF DET
      SpatialLagrangianUwPElement::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight );

      //KRATOS_WATCH( rLeftHandSideMatrix )
   }

   //************************************************************************************
   //************************************************************************************

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {
      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      // SOMETHING OF DET
      SpatialLagrangianUwPElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, IntegrationWeight );
      //KRATOS_WATCH( rRightHandSideVector )
   }


   //************************************************************************************
   //************************************************************************************
   void AxisymSpatialLagrangianUwPElement::CalculateRadius(double & rCurrentRadius,
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
         }
      }


      if ( dimension == 3 )
      {
         std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
      }

      KRATOS_CATCH( "" )
   }

   //********** INTERNAL FORCES**************************************************************************
   //********************** with the total stress tensor ************************************************

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight
         )
   {
      KRATOS_TRY

         const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // VectorType Fh=rRightHandSideVector;

      Vector TotalStressVector = rVariables.StressVector;

      double Pressure = 0.0;
      for (unsigned int i = 0; i < number_of_nodes; ++i)
         Pressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];


      for (unsigned int i = 0; i < 3; ++i) 
         TotalStressVector(i) += Pressure;

      Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), TotalStressVector );

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexup = dimension * i + i;
         unsigned int indexu  = dimension * i;

         for ( unsigned int j = 0; j < dimension; j++ )
         {
            rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
         }
      }


      KRATOS_CATCH( "" )
   }


   //*************** PRESSURE FORCES *********************************************************************
   //************************* aka: mass balance equation ************************************************

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      unsigned int indexp = dimension;

      VectorType Fh=rRightHandSideVector;


      Matrix K = ZeroMatrix(dimension);
      Vector b = ZeroVector(dimension);
      double WaterDensity = 0.0;
      b(dimension-1) = -WaterDensity;
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = Permeability;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            const double& CurrentPressure   = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE );
            const double& PreviousPressure = GetGeometry()[j].FastGetSolutionStepValue( WATER_PRESSURE , 1);
            double DeltaPressure = CurrentPressure - PreviousPressure;

            const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[j].FastGetSolutionStepValue( DISPLACEMENT );
            const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[j].FastGetSolutionStepValue( DISPLACEMENT , 1 );
            array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;

            rRightHandSideVector[indexp] += (1.0/WaterBulk) * rVariables.N[i] * rVariables.N[j] * DeltaPressure * rIntegrationWeight * ScalingConstant / rVariables.detF0;

            for ( unsigned int p = 0; p < dimension; ++p )
            {

               rRightHandSideVector[indexp] -= rVariables.N[i]*DeltaDisplacement[p] * rVariables.DN_DX(j,p) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
               if (p == 0)
                  rRightHandSideVector[indexp] -= rVariables.N[i]*DeltaDisplacement[p] * rVariables.N[j]*(1.0/ rVariables.CurrentRadius) * rIntegrationWeight * ScalingConstant/ rVariables.detF0;

               for ( unsigned int q = 0; q < dimension; ++q )
               {

                  rRightHandSideVector[indexp] += DeltaTime * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q)* CurrentPressure *rIntegrationWeight * ScalingConstant / rVariables.detF0;

                  if (j == 0) {
                     // Gravity term of the DARCY FLOW.
                     rRightHandSideVector[indexp] += DeltaTime * rVariables.DN_DX(i,p) * K(p,q)*b(q)*rIntegrationWeight* ScalingConstant / rVariables.detF0;
                  }

               }


            }

         }

         indexp += (dimension + 1);

      }


      KRATOS_CATCH( "" )

   }



   //*********** STABILIZATION *************************************************************************
   //****************** in the Stab element ************************************************************
   void AxisymSpatialLagrangianUwPElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
   }




   //******************* Geometric Stiffness matrix *************************************
   //************************************************************************************

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddKuug(MatrixType& rK,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY
         return;
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



      MatrixType Kh=rK;

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
                  rK(indexi+i,indexj+j)+=Kuu(indexi,indexj);
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


   //********** MATERIAL STIFFNESS MATRIX ***********************************************
   //*********************** has de Pw term *********************************************
   void AxisymSpatialLagrangianUwPElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         //assemble into rk the material uu contribution:
         const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // FIRST: Ctotal = Cmaterial + pJ(1*1 - 2I4)
      Matrix ConstitutiveMatrix = rVariables.ConstitutiveMatrix;

      double Pressure = 0.0;
      for (unsigned int i = 0; i < number_of_nodes; ++i)
         Pressure += GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) * rVariables.N[i];


      Matrix FourthOrderIdentity = ZeroMatrix(4);
      Matrix MatrixProduct = ZeroMatrix(4);
      for (unsigned int i = 0; i < 3 ; ++i) {
         FourthOrderIdentity(i,i) = 1.0;
         for (unsigned int j = 0; j < 3; ++j) {
            MatrixProduct(i,j) = 1.0;
         }
      }
      for (unsigned int i = 3; i < 4; ++i)
         FourthOrderIdentity(i,i) = 0.5;


      ConstitutiveMatrix += Pressure * rVariables.detF0 * (MatrixProduct - 2.0*FourthOrderIdentity); 

      Matrix Kuu = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( ConstitutiveMatrix, rVariables.B ) ) ); 

      // MatrixType Kh=rLeftHandSideMatrix;

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

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //MatrixType Kh=rLeftHandSideMatrix;
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

               if (k==0)
                  rLeftHandSideMatrix(indexup+k, indexp) +=rVariables.N(i) * rVariables.N[j] * (1.0/ rVariables.CurrentRadius)* rIntegrationWeight* rVariables.detF;
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

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY
         double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      Matrix K = ZeroMatrix(dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = Permeability;

      MatrixType Kh=rLeftHandSideMatrix;

      unsigned int indexp = dimension;

      for ( unsigned int i = 0; i < number_of_nodes; ++i )
      {
         for ( unsigned int j = 0; j < number_of_nodes; ++j )
         {
            int indexup = dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexp, indexup+k) += rVariables.N[i] * rVariables.DN_DX( j , k ) * rIntegrationWeight * ScalingConstant / rVariables.detF0;
               if (k==0)
                  rLeftHandSideMatrix(indexp, indexup+k) +=rVariables.N[i]*rVariables.N[j] * rIntegrationWeight*ScalingConstant* (1.0/ rVariables.CurrentRadius) / rVariables.detF0;

               // TRY TO PROGRAM LD TERMS...
               for (unsigned int m = 0; m < number_of_nodes; ++m ) 
               {
                  const array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[m].FastGetSolutionStepValue( DISPLACEMENT );
                  const array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[m].FastGetSolutionStepValue( DISPLACEMENT , 1 );
                  array_1d<double, 3 > DeltaDisplacement            = CurrentDisplacement-PreviousDisplacement;
                  const double & WaterPressure = GetGeometry()[m].FastGetSolutionStepValue( WATER_PRESSURE );
                  for (unsigned n = 0; n < dimension; n++) {

                     rLeftHandSideMatrix(indexp, indexup+k) -= rVariables.N[i] * DeltaDisplacement[n] * rVariables.DN_DX(m,k) * rVariables.DN_DX(j,n) * rIntegrationWeight * ScalingConstant / rVariables.detF0;

                     for (unsigned t = 0; t < dimension; ++t) {
                        rLeftHandSideMatrix(indexp, indexup+k) += (rVariables.DN_DX(i,k)*rVariables.DN_DX(m,n) + rVariables.DN_DX(i,n)*rVariables.DN_DX(m,k) ) * K(t,n) * rVariables.DN_DX(j,t) * WaterPressure * DeltaTime * rIntegrationWeight * ScalingConstant / rVariables.detF0;
                     }
                  }
               }

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

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

         const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double ScalingConstant;
      double Permeability; double WaterBulk; double DeltaTime;

      GetConstants(ScalingConstant, WaterBulk, DeltaTime, Permeability);
      Matrix K = ZeroMatrix(dimension);
      for (unsigned int i = 0; i < dimension; ++i)
         K(i,i) = Permeability;

      MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(WaterBulk)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * ScalingConstant / rVariables.detF0;
            for ( unsigned int p = 0; p < dimension; ++p ) {
               for ( unsigned int q = 0; q < dimension; ++q ) {
                  rLeftHandSideMatrix(indexpi, indexpj) -= DeltaTime * rVariables.DN_DX(i, p) * K( p, q) * rVariables.DN_DX(j,q) * rIntegrationWeight* ScalingConstant / rVariables.detF0;
               }
            }

            indexpj += (dimension + 1);
         }

         indexpi += (dimension + 1);
      }

      // std::cout<<std::endl;
      // std::cout<<" Kpp "<< (rLeftHandSideMatrix-Kh) <<std::endl;
      // std::cout<<" Kpp "<< (rLeftHandSideMatrix-Kh) / ScalingConstant / rIntegrationWeight* rVariables.detF0  * WaterBulk <<std::endl;


      KRATOS_CATCH( "" )
   }



   //************ STABILIZATION TANGENT MATRIX *******************************************
   //****************************** defined in the Stab element **************************

   void AxisymSpatialLagrangianUwPElement::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         double& rIntegrationWeight)
   {
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymSpatialLagrangianUwPElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.save("DeterminantF0",mDeterminantF0);
   }

   void AxisymSpatialLagrangianUwPElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
         rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
      rSerializer.load("DeterminantF0",mDeterminantF0);
   }







}  // END KRATOS NAMESPACE
