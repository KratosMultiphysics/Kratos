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
   SetMechanicalInitialStateProcess::SetMechanicalInitialStateProcess(ModelPart& rModelPart)
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

   // destructor 
   SetMechanicalInitialStateProcess::~SetMechanicalInitialStateProcess()
   {

   }

   void SetMechanicalInitialStateProcess::ExecuteInitialize()
   {

      if ( mGravity ) {
         SetInitialMechanicalState( mrModelPart, 1);
      }
      else {
         // set the same stress state to all the elements of the domain (i.e. gravity == 0)
         SetInitialMechanicalStateConstant( mrModelPart, mInitialStress[0], mInitialStress[1], 1);
      }


   }


   // THIS IS NOT THE PLACE TO PROGRAM THIS. However..
   // This function removes previous boundary conditions at nodes that are now in contact (that is, removes Dirichlet water pressure conditions from contacting nodes).

   void SetMechanicalInitialStateProcess::ExecuteFinalizeSolutionStep()
   {
      std::cout << " [ Trying to remove boundary conditions " << std::endl;

      const unsigned int NumberOfMeshes = mrModelPart.NumberOfMeshes();

      if (NumberOfMeshes < 2) {
         std::cout << " Nothing To Be Done ] " << std::endl;
         return;
      }

      unsigned int start = 1, arranged = 0, allPossible = 0;
      array_1d<double, 3 > ContactForce, NeigContactForce;
      double CFModul, CFModulNeig;
      int ContactNeig;
      for (unsigned int MeshId = start; MeshId < NumberOfMeshes; MeshId++)
      {
         for ( ModelPart::NodesContainerType::const_iterator in = mrModelPart.NodesBegin(start); in != mrModelPart.NodesEnd(start); ++in)
         {

            ContactForce = in->GetSolutionStepValue(CONTACT_FORCE);
            CFModul = fabs( ContactForce[0] ) + fabs( ContactForce[1] );
            if ( CFModul > 1e-5)
            {
               ContactNeig = 0;
               WeakPointerVector<Node<3 > >  & rN = in->GetValue( NEIGHBOUR_NODES );
               for ( unsigned int neig = 0; neig < rN.size(); neig++)
               {
                  NeigContactForce = rN[neig].GetSolutionStepValue( CONTACT_FORCE);
                  CFModulNeig = fabs( NeigContactForce[0]) + fabs( NeigContactForce[1]);
                  if ( CFModulNeig > 1e-5)
                     ContactNeig += 1;
               }

               if (ContactNeig == 2)
               {
                  if ( in->SolutionStepsDataHas(LINE_LOAD) )
                  {
                     array_1d<double, 3 > & rLineLoad = in->GetSolutionStepValue( LINE_LOAD);
                     if ( fabs( rLineLoad[0]) + fabs( rLineLoad[1]) > 1e-5)
                     {
                        rLineLoad *= 0.0;
                        //in->SetSolutionStepValue( LINE_LOAD, LineLoad); // ja està, no facis coses rares....
                     }
                  }
                  else
                  {
                     std::cout << " ES RARO PQ no HAY LINE LOAD " << std::endl;
                  }

                  if ( in->SolutionStepsDataHas( WATER_PRESSURE ) )
                  {
                     if ( in->IsFixed( WATER_PRESSURE ) )
                     {
                        in->Free( WATER_PRESSURE);
                        arranged++;
                     }


                  }

                  allPossible += 1;
               }
            }


         }

      }

      std::cout << " We have Done " << arranged << " from a possible bicontacting "<< allPossible << " in the BCCorrection ]"<< std::endl;

   }

   // THE FUNCTION
   void SetMechanicalInitialStateProcess::SetInitialMechanicalState(ModelPart& rModelPart, int EchoLevel)
   {

      if( EchoLevel > 0 )
         std::cout << "  [ InitialState, gravity " << std::endl;

      const unsigned int NumberOfMeshes = rModelPart.NumberOfMeshes();
      if( EchoLevel > 0 )
         std::cout << "    number of meshes: " << NumberOfMeshes << " meshes" << std::endl;

      unsigned int start = 0;
      if (NumberOfMeshes>1)
         start = 1;


      //SEARCH for Ymax
      double Ymax = rModelPart.NodesBegin(start)->Y(); 
      for (ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(start); in != rModelPart.NodesEnd(start); ++in)
      {
         if ( Ymax < in->Y() ) {
            Ymax = in->Y();
         }
      } 

      for (unsigned int MeshId = start; MeshId < NumberOfMeshes; MeshId++) {

         this->SetMechanicalState( rModelPart, MeshId, EchoLevel, Ymax);

      }

      if( EchoLevel > 0 )
         std::cout << "    End InitialState, gravity ]" << std::endl;

   }


   void SetMechanicalInitialStateProcess::SetInitialMechanicalStateConstant(ModelPart& rModelPart, double S1, double S2, int EchoLevel)
   {
      if( EchoLevel > 0 )
         std::cout << "  [ InitialState, constant-state " << std::endl;

      const unsigned int NumberOfMeshes = rModelPart.NumberOfMeshes();
      if( EchoLevel > 0 )
         std::cout << "    number of meshes: " << NumberOfMeshes << " meshes" << std::endl;

      unsigned int start = 0;
      if (NumberOfMeshes > 1)
         start = 1;

      for ( unsigned int MeshId = start; MeshId < NumberOfMeshes; MeshId++)  {
         std::cout << "    working on mesh: " << MeshId << " of the total " << NumberOfMeshes << " number of elements: " <<  rModelPart.NumberOfElements(MeshId) << std::endl;

         if ( rModelPart.NumberOfElements(MeshId) ) {
            ModelPart::ElementsContainerType::const_iterator FirstElement = rModelPart.ElementsBegin(MeshId);

            ConstitutiveLaw::Features LawFeatures;
            FirstElement->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

            if (LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW) ) {
               std::cout << "    begin of setting UP constant state " << std::endl;
               this->SetMechanicalStateConstantUP( rModelPart, MeshId, S1, S2, EchoLevel);
            }
            else {
               std::cout << "    begin of setting U or UwP constant state " << std::endl;
               this->SetMechanicalStateConstant( rModelPart, MeshId, S1, S2, EchoLevel);
            }
         }
         std::cout << "   end with this mesh " << std::endl;
      }

      if( EchoLevel > 0 )
         std::cout << "    End InitialState, constant-state ]" << std::endl;
   }

   void SetMechanicalInitialStateProcess::SetMechanicalStateConstant(ModelPart& rModelPart, const unsigned int& MeshId, const double& rS1, const double& rS2, int& EchoLevel)
   {
      std::vector<Vector> StressVector;
      Vector ThisVector = ZeroVector(6);
      ThisVector(0) = rS2; ThisVector(2) = rS2;
      ThisVector(1) = rS1;
      StressVector.push_back(ThisVector);

      ProcessInfo SomeProcessInfo;

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(MeshId); pElement != rModelPart.ElementsEnd(MeshId) ; pElement++)
      {
         pElement->SetValueOnIntegrationPoints( ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS, StressVector, SomeProcessInfo); 
      }
   }


   void SetMechanicalInitialStateProcess::SetMechanicalState(ModelPart& rModelPart, const unsigned int& MeshId, int& EchoLevel, const double& rYmax)
   {
      if( EchoLevel > 0 )
         std::cout << "    working on mesh: " << MeshId << " of the total " << rModelPart.NumberOfElements(MeshId)  << " elements " << std::endl;

      if ( !rModelPart.NumberOfElements(MeshId)) {
         if( EchoLevel > 0 )
            std::cout << "    end; no elements." << std::endl;
         return;
      }

      ModelPart::ElementsContainerType::const_iterator FirstElement = rModelPart.ElementsBegin(MeshId);
      ConstitutiveLaw::Features LawFeatures;

      FirstElement->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);
      if ( LawFeatures.mOptions.Is(ConstitutiveLaw::U_P_LAW) ) {
         if( EchoLevel > 0 )
            std::cout << "    begin of setting UP gravity state " << std::endl;
         this->SetMechanicalStateUP( rModelPart, MeshId, EchoLevel, rYmax);
      }
      else 
      {
         // MIRAR SI ES HIDRODINAMICO
         bool WaterPressureDofs = false;
         for ( unsigned int i = 0; i < FirstElement->GetGeometry().size(); ++i) {
            if ( FirstElement->GetGeometry()[i].SolutionStepsDataHas( WATER_PRESSURE ) == true ) {
               WaterPressureDofs = true;
            }
         }
         if ( WaterPressureDofs ) {
            if( EchoLevel > 0 )
               std::cout << "    begin of setting UwP gravity state " << std::endl;
            this->SetMechanicalStateUwP(rModelPart, MeshId, EchoLevel, rYmax);
         }
         else 
         {
            if( EchoLevel > 0 )
               std::cout << "    begin of setting U gravity state " << std::endl;
            this->SetMechanicalStateU( rModelPart, MeshId, EchoLevel, rYmax);

         }
         if( EchoLevel > 0 )
            std::cout << "   end with this mesh " << std::endl;
      }

   }


   void SetMechanicalInitialStateProcess::SetMechanicalStateConstantUP(ModelPart& rModelPart, const unsigned int& MeshId, const double& rS1, const double& rS2, int& EchoLevel)
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

      for (ModelPart::NodesContainerType::const_iterator pNode = rModelPart.NodesBegin(MeshId); pNode != rModelPart.NodesEnd(MeshId); pNode++)
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

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(MeshId); pElement != rModelPart.ElementsEnd(MeshId); pElement++)
      {
         //pElement->SetInitialMechanicalState( StressVector );  // to be ...
         pElement->SetValueOnIntegrationPoints( ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS, StressVector, SomeProcessInfo);
      }

   }



   void SetMechanicalInitialStateProcess::SetMechanicalStateUP(ModelPart& rModelPart, const unsigned int& MeshId, int& EchoLevel, const double& rYmax)
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

      double Pressure, VerticalStress, HorizontalStress;
      for (ModelPart::NodesContainerType::const_iterator pNode = rModelPart.NodesBegin(MeshId); pNode != rModelPart.NodesEnd(MeshId) ; pNode++)
      {

         VerticalStress = 10.0 * Density * (pNode->Y() - rYmax);
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

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(MeshId); pElement!=rModelPart.ElementsEnd(MeshId) ; ++pElement)
      {
         Geometry<Node <3> >&  rGeom = (pElement)->GetGeometry();
         double Y = 0;
         for (unsigned int i = 0; i < rGeom.size(); ++i) {
            Y += rGeom[i].Y();
         }
         Y /= double( rGeom.size() );

         VerticalStress = 10.0*Density*( Y - rYmax) ; 
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

   void SetMechanicalInitialStateProcess::SetMechanicalStateU( ModelPart& rModelPart, const unsigned int& MeshId, int& EchoLevel, const double& rYmax)
   {
      unsigned int Properties = rModelPart.NumberOfProperties();
      Properties -= 1;
      double MixtureDensity = rModelPart.GetProperties(Properties)[ DENSITY ];
      const double Knot = rModelPart.GetProperties(Properties)[ K0 ];

      double VerticalStress, HorizontalStress;
      ProcessInfo SomeProcessInfo;

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(MeshId); pElement!=rModelPart.ElementsEnd(MeshId) ; ++pElement)
      {
         Geometry<Node <3> >&  rGeom = (pElement)->GetGeometry();
         double Y = 0;
         for (unsigned int i = 0; i < rGeom.size(); ++i)
            Y += rGeom[i].Y();
         Y /= double( rGeom.size() );

         VerticalStress = 10.0*MixtureDensity*( Y - rYmax) ; 
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

   }


   void SetMechanicalInitialStateProcess::SetMechanicalStateUwP(ModelPart& rModelPart, const unsigned int& MeshId, int& EchoLevel, const double& rYmax)
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

      for (ModelPart::NodesContainerType::const_iterator pNode = rModelPart.NodesBegin(MeshId); pNode != rModelPart.NodesEnd(MeshId) ; pNode++) {

         WaterPressure = 10.0*WaterDensity * ( pNode->Y() -rYmax );

         double& rWaterPressure = pNode->FastGetSolutionStepValue( WATER_PRESSURE );
         rWaterPressure = WaterPressure ;
      }

      double VerticalStress, HorizontalStress;
      ProcessInfo SomeProcessInfo;

      for (ModelPart::ElementsContainerType::const_iterator pElement = rModelPart.ElementsBegin(MeshId); pElement!=rModelPart.ElementsEnd(MeshId) ; ++pElement)
      {
         Geometry<Node <3> >&  rGeom = (pElement)->GetGeometry();
         double Y = 0;
         for (unsigned int i = 0; i < rGeom.size(); ++i)
            Y += rGeom[i].Y();
         Y /= double( rGeom.size() );

         VerticalStress = 10.0*MixtureDensity*( Y - rYmax) ; 
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
