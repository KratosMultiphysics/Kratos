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
#include "custom_elements/axisym_updated_lagrangian_U_Pressure_element.hpp"
#include "pfem_solid_mechanics_application_variables.h"
//#include "includes/define.h"
//#include "utilities/math_utils.h"
//#include "includes/constitutive_law.h"


// axisym UP Element ( with coupled volumetric and deviatoric behaviour)
// (some part of the code is copied from AxisymUpdatedLagrangianUPElement (since it is not derived, I copy code)

namespace Kratos
{


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************
   // Aquest a l'altre no hi és....
   AxisymUpdatedLagrangianUPressureElement::AxisymUpdatedLagrangianUPressureElement()
      : UpdatedLagrangianUPressureElement()
   {
      //DO NOT CALL IT: only needed for Register and Serialization!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPressureElement::AxisymUpdatedLagrangianUPressureElement( IndexType NewId, GeometryType::Pointer pGeometry )
      : UpdatedLagrangianUPressureElement( NewId, pGeometry )
   {
      //DO NOT ADD DOFS HERE!!!
   }


   //******************************CONSTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPressureElement::AxisymUpdatedLagrangianUPressureElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
      : UpdatedLagrangianUPressureElement( NewId, pGeometry, pProperties )
   {
   }


   //******************************COPY CONSTRUCTOR**************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPressureElement::AxisymUpdatedLagrangianUPressureElement( AxisymUpdatedLagrangianUPressureElement const& rOther)
      :UpdatedLagrangianUPressureElement(rOther)
   {
   }


   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPressureElement&  AxisymUpdatedLagrangianUPressureElement::operator=(AxisymUpdatedLagrangianUPressureElement const& rOther)
   {
      UpdatedLagrangianUPressureElement::operator=(rOther);

      return *this;
   }


   //*********************************OPERATIONS*****************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUPressureElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
   {
      return Element::Pointer( new AxisymUpdatedLagrangianUPressureElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
   }


   //************************************CLONE*******************************************
   //************************************************************************************

   Element::Pointer AxisymUpdatedLagrangianUPressureElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
   {

      AxisymUpdatedLagrangianUPressureElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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

      return Element::Pointer( new AxisymUpdatedLagrangianUPressureElement(NewElement) );
   }


   //*******************************DESTRUCTOR*******************************************
   //************************************************************************************

   AxisymUpdatedLagrangianUPressureElement::~AxisymUpdatedLagrangianUPressureElement()
   {
   }



   //************************************************************************************
   //************************************************************************************

   int  AxisymUpdatedLagrangianUPressureElement::Check( const ProcessInfo& rCurrentProcessInfo )
   {
      KRATOS_TRY

      int correct = 0;

      correct = LargeDisplacementElement::Check(rCurrentProcessInfo);

      //verify compatibility with the constitutive law
      ConstitutiveLaw::Features LawFeatures;
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

      if(LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW))
         KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-P element type ", " UpdatedLagrangianUPressureElement" );

      //verify that the variables are correctly initialized

      if ( PRESSURE.Key() == 0 )
         KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" );

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



   //********************** COMPUTE LHS **********************************************************
   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
   {

      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      UpdatedLagrangianUPressureElement::CalculateAndAddLHS( rLocalSystem, rVariables, IntegrationWeight);
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
   {
      double IntegrationWeight = rIntegrationWeight * 2.0 * 3.141592654 * rVariables.CurrentRadius / GetProperties()[THICKNESS];

      UpdatedLagrangianUPressureElement::CalculateAndAddRHS( rLocalSystem, rVariables, rVolumeForce, IntegrationWeight);

   }

   //*************************************************************************
   //*************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddPressureForcesElemUP(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables,
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

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] += rVariables.N[i] * rVariables.N[j] *  Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber; //2D

         }

         rRightHandSideVector[indexp] -= rElementVariables.ElementalMeanStress * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;


         indexp += (dimension + 1);
      }


      KRATOS_CATCH( "" )

   }



   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddStabilizedPressureElemUP(VectorType& rRightHandSideVector,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      // VectorType Fh=rRightHandSideVector;
      // std::cout<<" Element "<<this->Id()<<" "<<std::endl;

      //use of this variable for the complete parameter: (deffault: 4)
      double AlphaStabilization  = 4.0;
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_P];
      AlphaStabilization *= StabilizationFactor;

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      AlphaStabilization=(AlphaStabilization/(18.0*LameMu));
      AlphaStabilization *= BulkModulus;


      if ( YoungModulus < 0.00001 )
      {
         AlphaStabilization = 4.0 * StabilizationFactor / 18.0 ;

         ProcessInfo SomeProcessInfo;
         std::vector<double> Values;
         LargeDisplacementElement::GetValueOnIntegrationPoints( SHEAR_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization /= Values[0];

         LargeDisplacementElement::GetValueOnIntegrationPoints( BULK_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization *= Values[0];

      }


      unsigned int indexp = dimension;
      double consistent = 1;

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            consistent=(-1)*AlphaStabilization;
            if(i==j)
               consistent=2*AlphaStabilization;

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;

         }


         indexp += (dimension + 1);
      }



      KRATOS_CATCH( "" )

   }

   // ******************************** KUUM **************************************************
   // ****************************************************************************************
   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddKuumElemUP ( MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //LMV
      //assemble into rk the material uu contribution:
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();


      Matrix ECConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      Matrix ConstitutiveMatrix = ZeroMatrix(rElementVariables.voigtsize,rElementVariables.voigtsize);

      ECConstitutiveMatrix = prod( rElementVariables.DeviatoricTensor, ECConstitutiveMatrix);


      // a. 2 ( pElem - p Nodal) I4S
      for ( unsigned int i = 0; i < 6; i++) {
         double voigtnumber = 1.0;
         if ( i > 2)
            voigtnumber = 0.500;
         ECConstitutiveMatrix(i,i) += 2.0*voigtnumber * ( rElementVariables.ElementalMeanStress - rElementVariables.NodalMeanStress);
      }

      // b. p 1 cross 1
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            ECConstitutiveMatrix(i,j) += rElementVariables.NodalMeanStress ;
         }
      }

      // c . - 2/3 I times sigma
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            double voigtnumber = 1.0;
            if ( j > 2)
               voigtnumber = 1.00;
            ECConstitutiveMatrix(i,j) -= (2.0/6.0) *voigtnumber* rVariables.StressVector(j);
         }
      }

      // LMV: try to put the sim part
      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 6; j++) {
           double voigtnumber = 1.0;
           if ( j > 2)
              voigtnumber = 1.00;
           ECConstitutiveMatrix(j,i) -= (2.0/6.0) *voigtnumber* rVariables.StressVector(j);
         }
      }

      for (unsigned int i = 0; i < 4; i++) {
         for (unsigned int j = 0; j < 4 ; j++)
         {
            ConstitutiveMatrix(i,j) = ECConstitutiveMatrix(i, j);
         }
      }

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
   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddKuugElemUP ( MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables,
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

   //********************* KuP ***********************************************
   //*************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddKupElemUP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         ThisElementData& rElementVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //LMV
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
               rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight; // * rVariables.detF; (aquest és el que posem a u, o sigui que no cal comentar)

               // axi term
               if ( k == 0)
                  rLeftHandSideMatrix(indexup+k, indexp) += rVariables.N[i] * rVariables.N[j] * ( 1.0 / rVariables.CurrentRadius) * rIntegrationWeight;
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

   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddKpuElemUP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         ThisElementData& rElementVariables,
         double& rIntegrationWeight)

   {
      KRATOS_TRY


      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      Matrix Kpu;
      //MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration


      Matrix ECConstitutiveMatrix = rVariables.ConstitutiveMatrix;

      Matrix Identity = ZeroMatrix(1,rElementVariables.voigtsize);
      for (unsigned int i = 0; i < 3 ; ++i) {
         Identity(0,i) = 1.0;
      }

      Matrix AuxConstitutiveVector = prod( Identity, ECConstitutiveMatrix);
      AuxConstitutiveVector *= ( 1.0/ 3.0 );

      // EXTRA TERMS.
      // 2 / 3 * Stress
      for (unsigned int i = 0; i < 6; i++) {
         double voigtnumber = 1.0;
         if ( i > 2)
            voigtnumber = 1.0; // ithink,...
         AuxConstitutiveVector(0, i) += (2.0/3.0) * voigtnumber * rVariables.StressVector(i);
      }
      // - ElemenPressure * 1
      for (unsigned int i = 0; i < 3; i++)
         AuxConstitutiveVector(0,i) -= rElementVariables.ElementalMeanStress;


      Matrix ConstitutiveVector = ZeroMatrix( 1, rElementVariables.voigtsize);
      if ( rElementVariables.voigtsize == 6)
      {
         ConstitutiveVector = AuxConstitutiveVector;
      }
      else
      {
         for (unsigned int i = 0; i < 4; i++) {
            ConstitutiveVector(0,i) = AuxConstitutiveVector(0,i);
         }
      }



      Kpu = prod( ConstitutiveVector, rVariables.B);

      Matrix Kpu2 = ZeroMatrix( number_of_nodes, number_of_nodes*dimension);
      for (unsigned int i = 0; i < number_of_nodes;  i++) {
         for (unsigned int j = 0; j < number_of_nodes*dimension; j++) {
            Kpu2 ( i, j ) = rVariables.N[i] * Kpu(0,j);
         }
      }



      Kpu2 *= rIntegrationWeight;

      unsigned int indexp = dimension;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            int indexup= dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexp,indexup+k) +=  Kpu2( i, j*2 + k) *mElementScalingNumber;
            }
         }
         indexp += (dimension + 1);
      }


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddKppElemUP (MatrixType& rLeftHandSideMatrix,
         ElementDataType& rVariables,
         ThisElementData& rElementVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //LMV
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      unsigned int indexpi = dimension;


      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {

            rLeftHandSideMatrix(indexpi,indexpj)  -= rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber; //2D

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
   void AxisymUpdatedLagrangianUPressureElement::CalculateAndAddKppStabElemUP (MatrixType& rLeftHandSideMatrix,
         ElementDataType & rVariables,
         ThisElementData& rElementVariables,
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //LMV
      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      // MatrixType Kh=rLeftHandSideMatrix;

      //contributions to stiffness matrix calculated on the reference configuration
      double consistent = 1.0;

      double AlphaStabilization  = 4.0;
      double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR_P];
      AlphaStabilization *= StabilizationFactor;
      if ( StabilizationFactor < 0.0001)
         return;

      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      //Experimental
      // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
      //   LameMu = rVariables.ConstitutiveMatrix(2,2);


      //use of this variable for the complete parameter:
      AlphaStabilization=(AlphaStabilization/(18.0*LameMu));
      AlphaStabilization *= BulkModulus;

      if ( YoungModulus < 0.00001 )
      {
         AlphaStabilization = 4.0 * StabilizationFactor / 18.0 ;

         ProcessInfo SomeProcessInfo;
         std::vector<double> Values;
         LargeDisplacementElement::GetValueOnIntegrationPoints( SHEAR_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization /= Values[0];

         LargeDisplacementElement::GetValueOnIntegrationPoints( BULK_MODULUS, Values, SomeProcessInfo);
         AlphaStabilization *= Values[0];


      }

      unsigned int indexpi = dimension;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         unsigned int indexpj = dimension;
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            consistent=(-1)*AlphaStabilization;
            if(indexpi==indexpj)
               consistent=2*AlphaStabilization;

            rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;     //2D

            indexpj += (dimension + 1);
         }

         indexpi += (dimension + 1);
      }

      // std::cout<<std::endl;
      // std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;


      KRATOS_CATCH( "" )
   }

   //*************************COMPUTE DEFORMATION GRADIENT*******************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateDeformationGradient(const Matrix& rDN_DX,
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

         KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

      }

      KRATOS_CATCH( "" )
   }




   //************************************************************************************
   //************************************************************************************


   void AxisymUpdatedLagrangianUPressureElement::CalculateDeformationMatrix(Matrix& rB,
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

   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateGreenLagrangeStrain(const Matrix& rF,
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

   void AxisymUpdatedLagrangianUPressureElement::CalculateAlmansiStrain(const Matrix& rF,
         Vector& rStrainVector )
   {
      KRATOS_TRY

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      //Left Cauchy-Green Calculation
      Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

      //Calculating the inverse of the jacobian
      Matrix InverseLeftCauchyGreen ( 3, 3 );
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


   //*************************COMPUTE AXYSIMMETRIC RADIUS********************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateRadius(double & rCurrentRadius,
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
         //std::cout<<" CurrentRadius "<<rCurrentRadius<<std::endl;
      }


      if ( dimension == 3 )
      {
         std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
      }

      KRATOS_CATCH( "" )
   }


   //************* STARTING - ENDING  METHODS
   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::Initialize()
   {
      KRATOS_TRY


      LargeDisplacementElement::Initialize();

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


      KRATOS_CATCH( "" )
   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::CalculateKinematics(ElementDataType& rVariables,
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

      //std::cout<<" detJ "<<rVariables.detJ<<" Area "<<2*GetGeometry().DomainSize()<<std::endl;

      //Compute cartesian derivatives [dN/dx_n]
      noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

      //Set Shape Functions Values for this integration point
      rVariables.N=row( Ncontainer, rPointNumber);

      //Calculate IntegrationPoint radius
      this->CalculateRadius (rVariables.CurrentRadius, rVariables.ReferenceRadius, rVariables.N);

      //Current Deformation Gradient [dx_n+1/dx_n]
      CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.DeltaPosition, rVariables.CurrentRadius, rVariables.ReferenceRadius);

      //Determinant of the deformation gradient F
      rVariables.detF  = MathUtils<double>::Det(rVariables.F);

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


   void AxisymUpdatedLagrangianUPressureElement::InitializeElementData (ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
   {

      // copy with modifications since I use everything

      const unsigned int number_of_nodes = GetGeometry().size();

      rVariables.detF  = 1;

      rVariables.detF0 = 1;

      rVariables.detH = 1;

      rVariables.B.resize( 4 , number_of_nodes * 2 );

      rVariables.F.resize( 3, 3 );

      rVariables.F0.resize( 3, 3 );

      rVariables.H.resize( 3, 3 );

      rVariables.ConstitutiveMatrix.resize( 6, 6 );

      rVariables.StrainVector.resize( 4 );

      rVariables.StressVector.resize( 6 );

      rVariables.DN_DX.resize( number_of_nodes, 2 );

      //set variables including all integration points values

      //reading shape functions
      rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

      //reading shape functions local gradients
      rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

      //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
      rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );


      //Calculate Delta Position
      ElementUtilities::CalculateDeltaPosition(rVariables.DeltaPosition,this->GetGeometry());

      //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
      rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );


   }

   void AxisymUpdatedLagrangianUPressureElement::CalculateThisElementData( ThisElementData& rElementVariables, const ElementDataType & rVariables)
   {

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();

      rElementVariables.voigtsize = 4;

      rElementVariables.NodalMeanStress = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         rElementVariables.NodalMeanStress += GetGeometry()[i].GetSolutionStepValue( PRESSURE ) * rVariables.N[i];

      rElementVariables.ElementalMeanStress = 0;
      for (unsigned int i = 0; i < 3 ; i++)
         rElementVariables.ElementalMeanStress += rVariables.StressVector(i);
      rElementVariables.ElementalMeanStress /= 3.0;

      Vector AuxStress = ZeroVector(6);
      AuxStress = rVariables.StressVector;
      for (unsigned int i = 0; i < 3; i++)
         AuxStress(i) += ( rElementVariables.NodalMeanStress - rElementVariables.ElementalMeanStress);

      rElementVariables.StressVector = ZeroVector(rElementVariables.voigtsize);
      rElementVariables.StressVector(0) = AuxStress(0);
      rElementVariables.StressVector(1) = AuxStress(1);
      rElementVariables.StressVector(2) = AuxStress(2);
      rElementVariables.StressVector(3) = AuxStress(3);

      rElementVariables.DeviatoricTensor = ZeroMatrix(6,6);
      for (unsigned int i = 0; i < 6; i++)
         rElementVariables.DeviatoricTensor(i,i) = 1;

      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            rElementVariables.DeviatoricTensor(i,j) -= 1.0/3.0;
         }
      }



   }


   //************************************************************************************
   //************************************************************************************

   void AxisymUpdatedLagrangianUPressureElement::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangianUPressureElement )
   }

   void AxisymUpdatedLagrangianUPressureElement::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangianUPressureElement )
   }



}  // Namespace Kratos
