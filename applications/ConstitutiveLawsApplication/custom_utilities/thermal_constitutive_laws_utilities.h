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

#if !defined(KRATOS_THERMAL_CONSTITUTIVE_LAWS_UTILITIES_H_INCLUDED)
#define KRATOS_THERMAL_CONSTITUTIVE_LAWS_UTILITIES_H_INCLUDED

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
 * @class ThermalConstitutiveLawsUtilities
 * @ingroup ConstitutiveLawsApplication
 * @details This class several methods required for thermomechanical calculations
 * @author Alejandro Cornejo
 */

  class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ThermalConstitutiveLawsUtilities
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ThermalConstitutiveLawsUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ThermalConstitutiveLawsUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ThermalConstitutiveLawsUtilities() {}

    /// Destructor.
    ~ThermalConstitutiveLawsUtilities()
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /**
       * @brief Computes and Set the value of reference temperature via SetValue()
       **/
    void ComputeAndSetReferenceTemperature(
        ModelPart &rModelPart);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
      std::stringstream buffer;
      buffer << "ThermalConstitutiveLawsUtilities";

      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "ThermalConstitutiveLawsUtilities"; }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const {}

    ///@}

  private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ThermalConstitutiveLawsUtilities &operator=(ThermalConstitutiveLawsUtilities const &rOther)
    {
      return *this;
    }

    /// Copy constructor.
    ThermalConstitutiveLawsUtilities(ThermalConstitutiveLawsUtilities const &rOther)
    {
      *this = rOther;
    }

    ///@}

  }; // Class ThermalConstitutiveLawsUtilities

  ///@}

  ///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_THERMAL_CONSTITUTIVE_LAWS_UTILITIES_H_INCLUDED  defined