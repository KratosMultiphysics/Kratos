// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_GENERIC_CONSTITUTIVE_LAWS_APPLICATION_UTILITIES_H_INCLUDED)
#define KRATOS_GENERIC_CONSTITUTIVE_LAWS_APPLICATION_UTILITIES_H_INCLUDED

// System includes
#include "spaces/ublas_space.h"
// Project includes

// Configures

// External includes

namespace Kratos
{

  ///@addtogroup ConstitutiveLawsApplication
  ///@{

  ///@name Kratos Classes
  ///@{

  /**
 * @class GenericConstitutiveLawsApplicationUtilities
 * @ingroup ConstitutiveLawsApplication
 * @details This class several methods required for thermomechanical calculations
 * @author Alejandro Cornejo
 */

  class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericConstitutiveLawsApplicationUtilities
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GenericConstitutiveLawsApplicationUtilities
    KRATOS_CLASS_POINTER_DEFINITION(GenericConstitutiveLawsApplicationUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenericConstitutiveLawsApplicationUtilities() {}

    /// Destructor.
    ~GenericConstitutiveLawsApplicationUtilities()
    {
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
      std::stringstream buffer;
      buffer << "GenericConstitutiveLawsApplicationUtilities";

      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const { rOStream << "GenericConstitutiveLawsApplicationUtilities"; }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const {}

    ///@}

  private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    GenericConstitutiveLawsApplicationUtilities &operator=(GenericConstitutiveLawsApplicationUtilities const &rOther)
    {
      return *this;
    }

    /// Copy constructor.
    GenericConstitutiveLawsApplicationUtilities(GenericConstitutiveLawsApplicationUtilities const &rOther)
    {
      *this = rOther;
    }

    ///@}

  }; // Class GenericConstitutiveLawsApplicationUtilities

  ///@}

  ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_GENERIC_CONSTITUTIVE_LAWS_APPLICATION_UTILITIES_H_INCLUDED  defined