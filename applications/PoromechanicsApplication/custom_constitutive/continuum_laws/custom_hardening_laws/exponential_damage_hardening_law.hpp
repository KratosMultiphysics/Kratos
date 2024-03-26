//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EXPONENTIAL_DAMAGE_HARDENING_LAW_H_INCLUDED)
#define  KRATOS_EXPONENTIAL_DAMAGE_HARDENING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/continuum_laws/custom_hardening_laws/hardening_law.hpp"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(POROMECHANICS_APPLICATION) ExponentialDamageHardeningLaw
	: public HardeningLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ExponentialDamageHardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION( ExponentialDamageHardeningLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ExponentialDamageHardeningLaw();


    /// Copy constructor.
    ExponentialDamageHardeningLaw(ExponentialDamageHardeningLaw const& rOther);

    /// Assignment operator.
    ExponentialDamageHardeningLaw& operator=(ExponentialDamageHardeningLaw const& rOther);

    /// Destructor.
    ~ExponentialDamageHardeningLaw() override;

    ///@}
    ///@name Operators
    ///@{

   /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this hardening law
     */
    HardeningLaw::Pointer Clone() const override;

    ///@}
    ///@name Operations
    ///@{

    double& CalculateHardening(double &rHardening, const Parameters& rValues) override;

    double& CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class ExponentialDamageHardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// // input stream function

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_EXPONENTIAL_DAMAGE_HARDENING_LAW_H_INCLUDED  defined


