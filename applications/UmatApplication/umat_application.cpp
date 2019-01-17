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



    KRATOS_REGISTER_CONSTITUTIVE_LAW( "VonMisesUmatSmallStrainModel", mVonMisesSmallStrainUmatModel );
    KRATOS_REGISTER_CONSTITUTIVE_LAW( "HypoplasticUmatSmallStrainModel", mHypoplasticSmallStrainUmatModel );
    KRATOS_REGISTER_CONSTITUTIVE_LAW( "VonMisesUmatLargeStrainModel", mVonMisesLargeStrainUmatModel );


  }
}  // namespace Kratos.
