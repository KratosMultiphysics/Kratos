//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#include "includes/model_part.h"

/* System includes */


/* External includes */


/* Project includes */
#include "custom_processes/non_local_plasticity_process.hpp"
#include "utilities/openmp_utils.h"

namespace Kratos
{

   //*********************************************************************************************
   // constructor
   NonLocalPlasticityProcess::NonLocalPlasticityProcess( ModelPart & rModelPart, Parameters rParameters)
      : Process(Flags()), mrModelPart( rModelPart)
   {
      KRATOS_TRY

      Parameters default_parameters( R"(
      {
         "echo_level": 1,
         "model_part_name": "Main_Domain",
         "characteristic_length": 0.0,
         "local_variables": ["PLASTIC_VOL_DEF", "PLASTIC_DEV_DEF"],
         "non_local_variables": ["NONLOCAL_PLASTIC_VOL_DEF"]
      } )" );


      rParameters.ValidateAndAssignDefaults( default_parameters);
      mCharacteristicLength = rParameters["characteristic_length"].GetDouble();
      if ( mCharacteristicLength == 0.0) {
         std::cout << " NonLocalPlasticityProcess:: " << std::endl;
         std::cout << "       using the default value for the characteristic length: DANGEROUS! " << std::endl;
      }


      //std::vector< std::string > var = rParameters["local_variables"];

      std::vector< std::string> localVariables;
      std::vector< std::string> nonlocalVariables;


      localVariables.push_back("PLASTIC_VOL_DEF");
      localVariables.push_back("PLASTIC_VOL_DEF_ABS");
      localVariables.push_back("PLASTIC_DEV_DEF");

      nonlocalVariables.push_back("NONLOCAL_PLASTIC_VOL_DEF");
      nonlocalVariables.push_back("NONLOCAL_PLASTIC_VOL_DEF_ABS");
      nonlocalVariables.push_back("NONLOCAL_PLASTIC_DEV_DEF");

      for (unsigned int var = 0; var < nonlocalVariables.size() ; var++) {
         const Variable<double> & rThisLocalVariable = KratosComponents< Variable<double> >::Get(localVariables[var]);
         const Variable<double> & rThisNonLocalVariable = KratosComponents< Variable<double> >::Get(nonlocalVariables[var]);
         mLocalVariables.push_back( rThisLocalVariable );
         mNonLocalVariables.push_back( rThisNonLocalVariable );
      }

      if ( mLocalVariables.size() != mNonLocalVariables.size() ) {
         KRATOS_ERROR << " NonLocalPlasticityProcess :: the number of local and nonlocal variables is not correct " << std::endl;
      }

      KRATOS_CATCH("")
   }

   //*********************************************************************************************
   // destructor
   NonLocalPlasticityProcess::~NonLocalPlasticityProcess()
   {

   }

   //*********************************************************************************************
   // Execute
   void NonLocalPlasticityProcess::Execute()
   {
      KRATOS_TRY

      const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
      double CharacteristicLength = mCharacteristicLength;

      std::vector< GaussPoint > GaussPointsVector;


      double start = OpenMPUtils::GetCurrentTime();

      this->PerformGaussPointSearch( GaussPointsVector, 3.0*CharacteristicLength);
      double start2 = OpenMPUtils::GetCurrentTime();



      const unsigned int numberGaussPoints = GaussPointsVector.size();
      std::vector< double > LocalVariableVector(numberGaussPoints);

      for ( unsigned int variable = 0; variable < mNonLocalVariables.size(); variable++) {

         const Variable<double> & rLocalVariable = mLocalVariables[variable];
         const Variable<double> & rNonLocalVariable = mNonLocalVariables[variable];

         for (unsigned int gp = 0; gp < numberGaussPoints; gp++) {
            LocalVariableVector[gp] = GaussPointsVector[gp].pConstitutiveLaw->GetValue( rLocalVariable, LocalVariableVector[gp]);
         }

         for (unsigned int gp = 0; gp < numberGaussPoints; gp++) {

            const GaussPoint & rGP = GaussPointsVector[gp];

            double numerator = 0;
            double denominator = 0;
            for (unsigned int nei = 0; nei < rGP.NeighbourWeight.size(); nei++) {
               numerator +=  LocalVariableVector[ rGP.NeighbourGP[nei] ] * rGP.NeighbourWeight[nei]; 
               denominator += rGP.NeighbourWeight[nei];
            }

            if ( denominator > 0) {
               numerator = numerator / denominator;
               rGP.pConstitutiveLaw->SetValue( rNonLocalVariable, numerator, rCurrentProcessInfo);
            } else {
               rGP.pConstitutiveLaw->SetValue( rNonLocalVariable, LocalVariableVector[gp], rCurrentProcessInfo);
            }


         }

      }


      double start3 = OpenMPUtils::GetCurrentTime();

      std::cout << "        SearchTime " << start2-start << " smooothing " << start3-start2 << " totalTime " << start3-start << std::endl;



      KRATOS_CATCH("")
   }

   //*********************************************************************************************
   // weighting function
   double & NonLocalPlasticityProcess::ComputeWeightFunction( const double & rDistance, const double & rCharacteristicLength, double & rAlpha)
   {
      KRATOS_TRY

      // Galavi & Schweiger
      rAlpha = rDistance * exp( - pow(rDistance/rCharacteristicLength, 2) );

      // Gaussian function
      //rAlpha =  exp( - pow(rDistance/rCharacteristicLength, 2) );
      return rAlpha;

      KRATOS_CATCH("")
   }


   //*********************************************************************************************
   // create  a list of the nodes that have an influence
   void NonLocalPlasticityProcess::PerformGaussPointSearch( std::vector< GaussPoint > & rGPVector, 
         const double CharacteristicLength)
   {
      KRATOS_TRY

      const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

      GeometryData::IntegrationMethod MyIntegrationMethod;
      array_1d<double,3> AuxLocalCoordinates;
      array_1d<double,3> AuxGlobalCoordinates;

      for ( ModelPart::ElementsContainerType::iterator ie = mrModelPart.ElementsBegin(); ie != mrModelPart.ElementsEnd(); ie++)
      {
         Element::GeometryType & rGeom = ie->GetGeometry();
         MyIntegrationMethod = ie->GetIntegrationMethod();
         const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
         unsigned int numberOfGP = IntegrationPoints.size();

         std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector(numberOfGP);
         ie->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW, ConstitutiveLawVector, rCurrentProcessInfo);

         for ( unsigned int nGP = 0 ; nGP < numberOfGP; nGP++) {

            for (unsigned int i = 0; i < 3; i++)
               AuxLocalCoordinates[i] = IntegrationPoints[nGP][i];

            rGeom.GlobalCoordinates(AuxGlobalCoordinates, AuxLocalCoordinates);

            rGPVector.push_back( GaussPoint( ConstitutiveLawVector[nGP], AuxGlobalCoordinates) );
         }

      }

      double alpha = 0;
      double zeroDistanceWeight;
      zeroDistanceWeight = ComputeWeightFunction( 0, mCharacteristicLength, zeroDistanceWeight);

      for (unsigned int ii = 0; ii < rGPVector.size(); ii++) {
         rGPVector[ii].AddNeighbour(ii, zeroDistanceWeight);

         const array_1d<double, 3> & rCoordII = rGPVector[ii].Coordinates;

         for (unsigned int jj = ii+1; jj < rGPVector.size(); jj++) {
            const array_1d<double, 3> & rCoordJJ = rGPVector[jj].Coordinates;

            double distance = 0;
            for (unsigned int i = 0; i < 3; i++)
               distance += pow( rCoordII[i] - rCoordJJ[i], 2);
            distance = pow(distance, 0.5);

            if ( distance < CharacteristicLength) {
               alpha = ComputeWeightFunction( distance, mCharacteristicLength, alpha); 
               rGPVector[ii].AddNeighbour( jj, alpha);
               rGPVector[jj].AddNeighbour( ii, alpha);
            }
         }
      }


      KRATOS_CATCH("")
   }



} // end namespace Kratos
