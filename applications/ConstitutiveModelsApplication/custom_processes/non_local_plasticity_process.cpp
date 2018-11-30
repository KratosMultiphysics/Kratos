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
         "characteristic_length": 1.0,
         "local_variables": ["PLASTIC_VOL_DEF", "PLASTIC_DEV_DEF"],
         "non_local_variables": ["NONLOCAL_PLASTIC_VOL_DEF"]
      } )" );


      rParameters.ValidateAndAssignDefaults( default_parameters);
      mCharacteristicLength = rParameters["characteristic_length"].GetDouble();
      if ( mCharacteristicLength == 1.0) {
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

      std::cout << " plotTheModelPart " << mrModelPart.NumberOfNodes() << std::endl;

      const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

      double CharacteristicLength = mCharacteristicLength;

      std::vector< GaussPoint > NeighbourGaussPoint(0);

      this->PerformGaussPointSearch( NeighbourGaussPoint, 3.0*CharacteristicLength);

      // And now we should transfer information

      for ( unsigned int i = 0; i < NeighbourGaussPoint.size() ; i++)
      {
         const ConstitutiveLaw::Pointer & rRecieveCL = NeighbourGaussPoint[i].pConstitutiveLaw;
         const std::vector<ConstitutiveLaw::Pointer > & rNeighLaws = NeighbourGaussPoint[i].NeighbourLaws;

         for ( unsigned int variable = 0; variable < mNonLocalVariables.size(); variable++) {

            const Variable<double> & rLocalVariable = mLocalVariables[variable];
            const Variable<double> & rNonLocalVariable = mNonLocalVariables[variable];

            double numerator = 0, denom = 0, value = 0, alpha = 0;

            for (unsigned int nei = 0; nei < rNeighLaws.size(); nei++)
            {
               value = NeighbourGaussPoint[i].NeighbourLaws[nei]->GetValue( rLocalVariable, value);
               alpha = ComputeWeightFunction( NeighbourGaussPoint[i].NeighbourDistances[nei], CharacteristicLength, alpha);

               numerator += value*alpha;
               denom += alpha;
            }

            value = numerator / denom;
            if ( denom > 0) {
               rRecieveCL->SetValue( rNonLocalVariable, value, rCurrentProcessInfo );
            } else {
               value = rRecieveCL->GetValue( rLocalVariable, value);
               rRecieveCL->SetValue( rNonLocalVariable, value, rCurrentProcessInfo);
            }

         }

      }




      KRATOS_CATCH("")
   }

   //*********************************************************************************************
   // weighting function
   double & NonLocalPlasticityProcess::ComputeWeightFunction( const double & rDistance, const double & rCharacteristicLength, double & rAlpha)
   {
      KRATOS_TRY

      rAlpha = rDistance * exp( - pow(rDistance/rCharacteristicLength, 2) );
      return rAlpha;

      KRATOS_CATCH("")
   }


   //*********************************************************************************************
   // create  a list of the nodes that have an influence
   void NonLocalPlasticityProcess::PerformGaussPointSearch( std::vector< GaussPoint > & rNeighbourGP, 
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

            rNeighbourGP.push_back( GaussPoint( ConstitutiveLawVector[nGP], AuxGlobalCoordinates) );
         }

      }

      for (unsigned int ii = 0; ii < rNeighbourGP.size(); ii++) {
         const array_1d<double, 3> & rCoordII = rNeighbourGP[ii].Coordinates;

         for (unsigned int jj = ii+1; jj < rNeighbourGP.size(); jj++) {
            const array_1d<double, 3> & rCoordJJ = rNeighbourGP[jj].Coordinates;

            double distance = 0;
            for (unsigned int i = 0; i < 3; i++)
               distance += pow( rCoordII[i] - rCoordJJ[i], 2);
            distance = pow(distance, 0.5);

            if ( distance < CharacteristicLength) {
               rNeighbourGP[ii].AddNeighbour( rNeighbourGP[jj].pConstitutiveLaw, distance);
               rNeighbourGP[jj].AddNeighbour( rNeighbourGP[ii].pConstitutiveLaw, distance);
            }
         }
         //std::cout << " jj " << ii << " has neightbours " << rNeighbourGP[ii].NeighbourLaws.size() << std::endl;
      }


      KRATOS_CATCH("")
   }



} // end namespace Kratos
