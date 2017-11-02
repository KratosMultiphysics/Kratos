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

  KratosUmatApplication::KratosUmatApplication() {}

  void KratosUmatApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "             _   _            _              "<< std::endl;
    std::cout << "     KRATOS | | | |_ __  __ _| |_            "<< std::endl;
    std::cout << "            | |_| | '  \\/ _` |  _|           "<< std::endl;
    std::cout << "             \\___/|_|_|_\\__,_|\\__| INTERFACE "<< std::endl;
    std::cout << "Initializing KratosConstitutiveModelsApplication... " << std::endl;
    std::cout << "Initializing KratosUmatApplication... " << std::endl;

    

    Serializer::Register( "VonMisesUmatSmallStrainModel", mVonMisesSmallStrainUmatModel);
    Serializer::Register( "HypoplasticUmatSmallStrainModel", mHypoplasticSmallStrainUmatModel);
    Serializer::Register( "VonMisesUmatLargeStrainModel", mVonMisesLargeStrainUmatModel);
    

  }
}  // namespace Kratos.
