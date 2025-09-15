//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BAKER_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_BAKER_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_hardening_laws/hardening_law.hpp"

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
class KRATOS_API(SOLID_MECHANICS_APPLICATION) BakerJohnsonCookThermalHardeningLaw
	: public HardeningLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BakerJohnsonCookThermalHardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION(BakerJohnsonCookThermalHardeningLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BakerJohnsonCookThermalHardeningLaw();


    /// Copy constructor.
    BakerJohnsonCookThermalHardeningLaw(BakerJohnsonCookThermalHardeningLaw const& rOther);

    /// Assignment operator.
    BakerJohnsonCookThermalHardeningLaw& operator=(BakerJohnsonCookThermalHardeningLaw const& rOther);

    /// Destructor.
    ~BakerJohnsonCookThermalHardeningLaw() override;

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

    double& CalculateIsotropicHardening(double &rIsotropicHardening, const Parameters& rValues) override;

    double& CalculateKinematicHardening(double &rKinematicHardening, const Parameters& rValues) override;


    double& CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues) override;


    double& CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const Parameters& rValues) override;

    double& CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening, const Parameters& rValues) override;


    double& CalculateDeltaThermalHardening(double &rDeltaThermalHardening, const Parameters& rValues) override;


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    // /// Turn back information as a string.
    // std::string Info() const;

    // /// Print information about this object.
    // void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // void PrintData(std::ostream& rOStream) const;


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

}; // Class BakerJohnsonCookThermalHardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// // input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   BakerJohnsonCookThermalHardeningLaw& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const BakerJohnsonCookThermalHardeningLaw& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BAKER_JOHNSON_COOK_THERMAL_HARDENING_LAW_H_INCLUDED  defined
