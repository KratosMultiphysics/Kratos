//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:       JMCarbonell $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

// System includes


// External includes


// Project includes
#include "umat_application.h"
#include "umat_application_variables.h"


namespace Kratos {

  KratosUmatApplication::KratosUmatApplication():
    KratosApplication("UmatApplication")
    {}

  void KratosUmatApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "             _   _            _              "<< std::endl;
    std::cout << "     KRATOS | | | |_ __  __ _| |_            "<< std::endl;
    std::cout << "            | |_| | '  \\/ _` |  _|           "<< std::endl;
    std::cout << "             \\___/|_|_|_\\__,_|\\__| INTERFACE "<< std::endl;
    std::cout << "Initializing KratosConstitutiveModelsApplication... " << std::endl;
    std::cout << "Initializing KratosUmatApplication... " << std::endl;

    

    Serializer::Register( "FabricUmatSmallStrainModel", mFabricSmallStrainUmatModel);
    Serializer::Register( "SanisandSmallStrainUmatModel", mSanisandSmallStrainUmatModel);
    Serializer::Register( "VonMisesUmatSmallStrainModel", mVonMisesSmallStrainUmatModel);
    Serializer::Register( "HypoplasticUmatSmallStrainModel", mHypoplasticSmallStrainUmatModel);
    Serializer::Register( "VonMisesUmatLargeStrainModel", mVonMisesLargeStrainUmatModel);
 
   KRATOS_REGISTER_VARIABLE( ALPHA )
   KRATOS_REGISTER_VARIABLE( BETA )   
   KRATOS_REGISTER_VARIABLE( MF )   
   KRATOS_REGISTER_VARIABLE( CC )   
   KRATOS_REGISTER_VARIABLE( MM )   
   KRATOS_REGISTER_VARIABLE( KSIS )   
   KRATOS_REGISTER_VARIABLE( RHOM )   
   KRATOS_REGISTER_VARIABLE( PC0 )   
   KRATOS_REGISTER_VARIABLE( VOID_RATIO )   
   KRATOS_REGISTER_VARIABLE( PLASTIC_MULTIPLIER )   

   KRATOS_REGISTER_VARIABLE( P_ATM  )
   KRATOS_REGISTER_VARIABLE( E0 )
   KRATOS_REGISTER_VARIABLE( LAMBDA_SANISAND )
   KRATOS_REGISTER_VARIABLE( EPSI )
   KRATOS_REGISTER_VARIABLE( MC )
   KRATOS_REGISTER_VARIABLE( ME )
   KRATOS_REGISTER_VARIABLE( M_SANISAND )
   KRATOS_REGISTER_VARIABLE( G0_SANISAND )
   KRATOS_REGISTER_VARIABLE( NU_SANISAND )
   KRATOS_REGISTER_VARIABLE( H0 )
   KRATOS_REGISTER_VARIABLE( CH )
   KRATOS_REGISTER_VARIABLE( NB )
   KRATOS_REGISTER_VARIABLE( A_SANISAND )
   KRATOS_REGISTER_VARIABLE( ND )
   KRATOS_REGISTER_VARIABLE( PTMULT )

  }
}  // namespace Kratos.
