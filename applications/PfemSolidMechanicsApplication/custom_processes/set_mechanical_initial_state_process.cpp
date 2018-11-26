//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#include "includes/model_part.h"

/* System includes */


/* External includes */


/* Project includes */
#include "custom_processes/set_mechanical_initial_state_process.hpp"

namespace Kratos
{

   // constructor 1
/*   SetMechanicalInitialStateProcess::SetMechanicalInitialStateProcess(ModelPart& rModelPart)
      : mrModelPart(rModelPart)
   {

   }

   // constructor 2
   SetMechanicalInitialStateProcess::SetMechanicalInitialStateProcess(ModelPart& rModelPart, const bool rGravity, const double rSV, const double rSH )
      : mrModelPart(rModelPart)
   {
      mGravity = rGravity;
      mInitialStress.resize(0);
      mInitialStress.push_back(rSV);
      mInitialStress.push_back(rSH);
   }
*/
   SetMechanicalInitialStateProcess::SetMechanicalInitialStateProcess( ModelPart & rModelPart, Parameters rParameters)
      : Process(Flags()), mrModelPart( rModelPart)
   {
      KRATOS_TRY

      std::cout << " FINALLY CONSTRUCTED " << std::endl;
      Parameters default_parameters( R"(
      {
         "model_part_name":"MODEL_PART_NAME",
         "gravity_active": false, 
         "constant_vertical_stress": -10.0,
         "constant_horizontal_stress": -5.0,
         "constant_water_pressure" : -10.0,
         "top_surface_load_bool": false,
         "top_surface_load": 0.0,
         "top_water_pressure": 0.0
      } )" );

      mGravity = rParameters["gravity_active"].GetBool();

      mInitialStress.resize(0);
      double sv = rParameters["constant_vertical_stress"].GetDouble();
      mInitialStress.push_back( sv);
      double sh = rParameters["constant_horizontal_stress"].GetDouble();
      mInitialStress.push_back( sh);

      mInitialWaterPressure = rParameters["constant_water_pressure"].GetDouble();
      mSurfaceLoadBool = rParameters["top_surface_load_bool"].GetBool();
      mSurfaceLoad = rParameters["top_surface_load"].GetDouble();
      mWaterLoad = rParameters["top_water_pressure"].GetDouble();

      std::cout << " FINALLY CONSTRUCTED " << std::endl;
      std::cout << " WRITTING THE MODEL PART OR SOMETHING " << rModelPart << std::endl;
      KRATOS_CATCH("")
   }
   
   // destructor 
   SetMechanicalInitialStateProcess::~SetMechanicalInitialStateProcess()
   {

   }

   void SetMechanicalInitialStateProcess::Execute()
   {

      if ( mGravity ) {
         SetInitialMechanicalState( mrModelPart, 1);
      }
      else {
         // set the same stress state to all the elements of the domain (i.e. gravity == 0)
         SetInitialMechanicalStateConstant( mrModelPart, mInitialStress[0], mInitialStress[1], mInitialWaterPressure,  1);
      }


   }



   // THE FUNCTION
   void SetMechanicalInitialStateProcess::SetInitialMechanicalState(ModelPart& rModelPart, int EchoLevel)
   {

      if( EchoLevel > 0 )
         std::cout << "  [ InitialState, gravity " << std::endl;

      const unsigned int NumberOfMeshes = rModelPart.NumberOfMeshes();
      if( EchoLevel > 0 )
         std::cout << "    number of meshes: " << NumberOfMeshes << " meshes" << std::endl;


      //SEARCH for Ymax
      double Ymax = rModelPart.NodesBegin()->Y(); 
      for (ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in != rModelPart.NodesEnd(); ++in)
      {
         if ( ( Ymax < in->Y() ) && ( in->IsNot(RIGID) ) ) {
            Ymax = in->Y();
         }
      } 

      this->SetMechanicalState( rModelPart, EchoLevel, Ymax);

      if( EchoLevel > 0 )
         std::cout << "    End InitialState, gravity ]" << std::endl;

   }


   void SetMechanicalInitialStateProcess::SetInitialMechanicalStateConstant(ModelPart& rModelPart, double S1, double S2, double WaterPressure,  int EchoLevel)
   {
      if( EchoLevel > 0 )
         std::cout << "  [ InitialState, constant-state " << std::endl;

      const unsigned int NumberOfMeshes = rModelPart.NumberOfMeshes();
      if( EchoLevel > 0 )
         std::cout << "    number of meshes: " << NumberOfMeshes << " meshes" << std::endl;

      std::cout << " nEl " << rModelPart.NumberOfElements() << std::endl;

      if ( rModelPart.NumberOfElements() ) {
	ModelPart::ElementsContainerType::const_iterator FirstElement = rModelPart.ElementsBegin();

	ConstitutiveLaw::Features LawFeatures;
	FirstElement->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

	if (LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW) ) {
	  std::cout << "    begin of setting UP constant state " << std::endl;
	  this->SetMechanicalStateConstantUP( rModelPart, S1, S2, EchoLevel);
	}
	else {
	  std::cout << "    begin of setting U or UwP constant state " << std::endl;
	  this->SetMechanicalStateConstant( rModelPart, S1, S2, WaterPressure, EchoLevel);
	}
      }

      if( EchoLevel > 0 )
         std::cout << "    End InitialState, constant-state ]" << std::endl;
   }

   void SetMechanicalInitialStateProcess::SetMechanicalStateConstant(ModelPart& rModelPart, const double& rS1, const double& rS2, const double& rWaterPressure, int& EchoLevel)
   {
      std::vector<Vector> StressVector;
      Vector ThisVector = ZeroVector(6);
      ThisVector(0) = rS2; ThisVector(2) = rS2;
      ThisVector(1) = rS1;
      StressVector.push_back(ThisVector);

      ProcessInfo SomeProcessInfo;

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(); pElement != rModelPart.ElementsEnd() ; pElement++)
      {
         pElement->SetValueOnIntegrationPoints( ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS, StressVector, SomeProcessInfo); 
      }
      // AND NOW SET THE WATER PRESSURE

      for (ModelPart::NodesContainerType::const_iterator pNode = rModelPart.NodesBegin(); pNode != rModelPart.NodesEnd(); pNode++)
      {
         if ( pNode->SolutionStepsDataHas( WATER_PRESSURE)  ) 
         {
            double & rNodeWaterPressure = pNode->FastGetSolutionStepValue( WATER_PRESSURE );
            rNodeWaterPressure = rWaterPressure;
            double & rNodeWaterPressureOld = pNode->FastGetSolutionStepValue( WATER_PRESSURE , 1);
            rNodeWaterPressureOld = rWaterPressure;
         }

      }

   }


   void SetMechanicalInitialStateProcess::SetMechanicalState(ModelPart& rModelPart, int& EchoLevel, const double& rYmax)
   {

      if ( !rModelPart.NumberOfElements()) {
         if( EchoLevel > 0 )
            std::cout << "    end; no elements." << std::endl;
         return;
      }

      ModelPart::ElementsContainerType::const_iterator FirstElement = rModelPart.ElementsBegin();
      ConstitutiveLaw::Features LawFeatures;

      FirstElement->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);


      // Look if any element has water pressure degree of freedom
      bool WaterPressureDofs = false;
      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(); pElement != rModelPart.ElementsEnd(); pElement++) {
         for ( unsigned int i = 0; i < pElement->GetGeometry().size(); ++i) {
            if ( pElement->GetGeometry()[i].SolutionStepsDataHas( WATER_PRESSURE ) == true ) {
               WaterPressureDofs = true;
            }
         }
      }


      if ( LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW) ) {
         if( EchoLevel > 0 )
            std::cout << "    begin of setting UP gravity state " << std::endl;
         this->SetMechanicalStateUP( rModelPart, EchoLevel, rYmax);
      }
      else 
      {
         if ( WaterPressureDofs ) {
            if( EchoLevel > 0 )
               std::cout << "    begin of setting UwP gravity state " << std::endl;
            this->SetMechanicalStateUwP(rModelPart, EchoLevel, rYmax);
         }
         else 
         {
            if( EchoLevel > 0 )
               std::cout << "    begin of setting U gravity state " << std::endl;
            this->SetMechanicalStateU( rModelPart, EchoLevel, rYmax);

         }
         if( EchoLevel > 0 )
            std::cout << "   end with this mesh " << std::endl;
      }

   }


   void SetMechanicalInitialStateProcess::SetMechanicalStateConstantUP(ModelPart& rModelPart, const double& rS1, const double& rS2, int& EchoLevel)
   {
      double Pressure, VerticalStress, HorizontalStress;
      Pressure = (rS1 + 2.0*rS2) / 3.0;
      VerticalStress = rS1 - Pressure;
      HorizontalStress = rS2 - Pressure;

      unsigned int Properties = rModelPart.NumberOfProperties();
      Properties -= 1;
      const double Young = rModelPart.GetProperties(Properties)[ YOUNG_MODULUS ];
      const double Poisson = rModelPart.GetProperties(Properties)[ POISSON_RATIO ];

      if( EchoLevel > 0 )
         std::cout << "    YOUNG: " << Young << " Poisson " << Poisson << std::endl;

      double BulkModulus = Young;
      BulkModulus /= 3.0 * ( 1.0 - 2.0*Poisson);

      for (ModelPart::NodesContainerType::const_iterator pNode = rModelPart.NodesBegin(); pNode != rModelPart.NodesEnd(); pNode++)
      {
         double & rPressure = pNode->FastGetSolutionStepValue( PRESSURE );
         rPressure = Pressure;
      }

      double UndeformedDet = Pressure / BulkModulus;
      UndeformedDet = std::exp(UndeformedDet);
      std::vector<Vector> StressVector;
      Vector ThisVector = ZeroVector(6);
      ThisVector(0) = HorizontalStress; 
      ThisVector(1) = VerticalStress;
      ThisVector(2) = UndeformedDet;
      StressVector.push_back(ThisVector);

      ProcessInfo SomeProcessInfo;

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(); pElement != rModelPart.ElementsEnd(); pElement++)
      {
         //pElement->SetInitialMechanicalState( StressVector );  // to be ...
         pElement->SetValueOnIntegrationPoints( ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS, StressVector, SomeProcessInfo);
      }

   }



   void SetMechanicalInitialStateProcess::SetMechanicalStateUP(ModelPart& rModelPart, int& EchoLevel, const double& rYmax)
   {

      unsigned int Properties = rModelPart.NumberOfProperties();
      Properties -= 1;
      double Density = rModelPart.GetProperties(Properties)[ DENSITY ];
      const double Knot = rModelPart.GetProperties(Properties)[ K0 ];
      const double Young = rModelPart.GetProperties(Properties)[ YOUNG_MODULUS ];
      const double Poisson = rModelPart.GetProperties(Properties)[ POISSON_RATIO ];
      double Su = rModelPart.GetProperties(Properties)[ YIELD_STRESS ];
      if( EchoLevel > 0 )
         std::cout << " YOUNG: " << Young << " Poisson " << Poisson << " Density " << Density << " K0 " << Knot << std::endl;


      double BulkModulus = Young;
      BulkModulus /= 3.0 * ( 1.0 - 2.0*Poisson);

      double OverLoad = 0;
      if ( mSurfaceLoadBool == true)
      {
         OverLoad = mSurfaceLoad; 
      }

      double Pressure, VerticalStress, HorizontalStress;
      for (ModelPart::NodesContainerType::const_iterator pNode = rModelPart.NodesBegin(); pNode != rModelPart.NodesEnd() ; pNode++)
      {

         VerticalStress = 10.0 * Density * (pNode->Y() - rYmax) + OverLoad;
         if ( VerticalStress > 0.0)
            VerticalStress = 0.0;
         HorizontalStress = Knot * VerticalStress;

         if ( fabs(VerticalStress - HorizontalStress) > 2.0*Su ) {
            HorizontalStress = VerticalStress + 2.0*Su;
         }

         Pressure = ( VerticalStress + 2.0*HorizontalStress) / 3.0;

         double& rPressure = pNode->FastGetSolutionStepValue( PRESSURE );
         rPressure = Pressure;
      }


      double MeanStress;
      ProcessInfo SomeProcessInfo;

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(); pElement!=rModelPart.ElementsEnd() ; ++pElement)
      {
         Geometry<Node <3> >&  rGeom = (pElement)->GetGeometry();
         double Y = 0;
         for (unsigned int i = 0; i < rGeom.size(); ++i) {
            Y += rGeom[i].Y();
         }
         Y /= double( rGeom.size() );

         VerticalStress = 10.0*Density*( Y - rYmax) + OverLoad; 
         if ( VerticalStress > 0.0)
            VerticalStress = 0.0;
         HorizontalStress = Knot * VerticalStress;

         if ( fabs(VerticalStress - HorizontalStress) > 2.0*Su) {
            HorizontalStress = VerticalStress + 2.0*Su;
         }

         MeanStress = (VerticalStress + 2.0*HorizontalStress)/3.0;
         //UndeformedDet = MeanStress / BulkModulus;
         //UndeformedDet = std::exp(UndeformedDet);


         std::vector<Vector> StressVector;
         unsigned int NumberOfGaussPoints = 1;
         for (unsigned int i = 0; i < NumberOfGaussPoints; ++i) {
            Vector ThisVector = ZeroVector(6);
            ThisVector(0) = HorizontalStress - MeanStress;
            ThisVector(1) = VerticalStress - MeanStress;
            ThisVector(2) = ThisVector(0);
            StressVector.push_back(ThisVector);
         }

         pElement->SetValueOnIntegrationPoints( ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS, StressVector, SomeProcessInfo);
      }


   }

   void SetMechanicalInitialStateProcess::SetMechanicalStateU( ModelPart& rModelPart, int& EchoLevel, const double& rYmax)
   {
      unsigned int Properties = rModelPart.NumberOfProperties();
      Properties -= 1;
      double MixtureDensity = rModelPart.GetProperties(Properties)[ DENSITY ];
      const double Knot = rModelPart.GetProperties(Properties)[ K0 ];

      double VerticalStress, HorizontalStress;
      ProcessInfo SomeProcessInfo;

      double OverLoad = 0;
      if ( mSurfaceLoadBool == true)
      {
         OverLoad = mSurfaceLoad; 
      }

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(); pElement!=rModelPart.ElementsEnd() ; ++pElement)
      {
         Geometry<Node <3> >&  rGeom = (pElement)->GetGeometry();
         double Y = 0;
         for (unsigned int i = 0; i < rGeom.size(); ++i)
            Y += rGeom[i].Y();
         Y /= double( rGeom.size() );

         VerticalStress = 10.0*MixtureDensity*( Y - rYmax) + OverLoad; 
         if ( VerticalStress > 0.0)
            VerticalStress = 0.0;
         HorizontalStress = Knot * VerticalStress;

         std::vector<Vector> StressVector;
         unsigned int NumberOfGaussPoints = 1;
         for (unsigned int i = 0; i < NumberOfGaussPoints; ++i) {
            Vector ThisVector = ZeroVector(6);
            ThisVector(0) = HorizontalStress;
            ThisVector(1) = VerticalStress;
            ThisVector(2) = HorizontalStress;
            StressVector.push_back(ThisVector);
         }
         //pElement->SetInitialMechanicalState( StressVector ); to be ...
         pElement->SetValueOnIntegrationPoints( ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS, StressVector, SomeProcessInfo);
      }

      // THE PART TO PUT THE PRESSURE FOR THE NEW ELEMENTS.
      double Pressure; 
      for (ModelPart::NodesContainerType::const_iterator pNode = rModelPart.NodesBegin(); pNode != rModelPart.NodesEnd(); pNode++)
      {
         if ( pNode->SolutionStepsDataHas( PRESSURE ) ) 
         {
            VerticalStress = 10.0 * MixtureDensity * (pNode->Y() - rYmax);
            HorizontalStress = Knot * VerticalStress;


            Pressure = ( VerticalStress + 2.0*HorizontalStress) / 3.0;

            double& rPressure = pNode->FastGetSolutionStepValue( PRESSURE );
            rPressure = Pressure;

         }
      }



   }


   void SetMechanicalInitialStateProcess::SetMechanicalStateUwP(ModelPart& rModelPart, int& EchoLevel, const double& rYmax)
   {

      unsigned int Properties = rModelPart.NumberOfProperties();
      Properties -= 1;
      const double WaterDensity = rModelPart.GetProperties(Properties)[ DENSITY_WATER ];
      double MixtureDensity = rModelPart.GetProperties(Properties)[ DENSITY ];
      const double Knot = rModelPart.GetProperties(Properties)[ K0 ];
      MixtureDensity -= WaterDensity;
      double WaterPressure;

      if( EchoLevel > 0 )
         std::cout << " WaterDensity: " << WaterDensity<< " MixtureDensity " << MixtureDensity+WaterDensity << " K0 " << Knot << std::endl;

      double OverLoad = 0;
      double WaterOverLoad = 0;
      if ( mSurfaceLoadBool == true)
      {
         OverLoad = mSurfaceLoad; 
         WaterOverLoad = mWaterLoad;
      }


      // try to put zero water pressure below 0 stress.
      /*double Ymax2 = rYmax - OverLoad/10.0/(MixtureDensity + WaterDensity) ;
      std::cout << " YMax2 " << Ymax2 << std::endl;
      std::cout << " rYmax " << rYmax << std::endl;
      std::cout << " number " << OverLoad / 10.0 / (MixtureDensity +WaterDensity) << std::endl;
      std::cout << " olverLoad " << OverLoad << std::endl; */

      double sign = 1.0;
      if ( (rModelPart.NodesBegin() )->SolutionStepsDataHas( WATER_DISPLACEMENT_X) )
         sign = -1.0;
      for (ModelPart::NodesContainerType::const_iterator pNode = rModelPart.NodesBegin(); pNode != rModelPart.NodesEnd() ; pNode++) {

         WaterPressure = 10.0*WaterDensity * ( pNode->Y() -rYmax ) + WaterOverLoad;

         if ( WaterPressure > 0.0)
            WaterPressure = 0.0;

         if ( pNode->Y() > rYmax)
            WaterPressure = 0.0;

         double& rWaterPressure = pNode->FastGetSolutionStepValue( WATER_PRESSURE );
         rWaterPressure = sign * WaterPressure ;
         double& rWaterPressureOld = pNode->FastGetSolutionStepValue( WATER_PRESSURE, 1 );
         rWaterPressureOld = sign * WaterPressure ;
      }

      double VerticalStress, HorizontalStress;
      ProcessInfo SomeProcessInfo;

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(); pElement!=rModelPart.ElementsEnd() ; ++pElement)
      {
         Geometry<Node <3> >&  rGeom = (pElement)->GetGeometry();
         double Y = 0;
         for (unsigned int i = 0; i < rGeom.size(); ++i)
            Y += rGeom[i].Y();
         Y /= double( rGeom.size() );

         VerticalStress = 10.0*MixtureDensity*( Y - rYmax) + OverLoad; 
         if ( VerticalStress > 0.0)
            VerticalStress = 0.0;
         HorizontalStress = Knot * VerticalStress;

         std::vector<Vector> StressVector;
         unsigned int NumberOfGaussPoints = 1;
         for (unsigned int i = 0; i < NumberOfGaussPoints; ++i) {
            Vector ThisVector = ZeroVector(6);
            ThisVector(0) = HorizontalStress;
            ThisVector(1) = VerticalStress;
            ThisVector(2) = HorizontalStress;
            StressVector.push_back(ThisVector);
         }
         pElement->SetValueOnIntegrationPoints( ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS, StressVector, SomeProcessInfo);
      }


   }



} // end namespace Kratos
