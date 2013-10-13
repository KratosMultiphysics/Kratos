//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_NON_LINEAR_ISOTROPIC_KINEMATIC_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_NON_LINEAR_ISOTROPIC_KINEMATIC_HARDENING_LAW_H_INCLUDED



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
class NonLinearIsotropicKinematicHardeningLaw 
	: public HardeningLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearIsotropicKinematicHardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION(NonLinearIsotropicKinematicHardeningLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NonLinearIsotropicKinematicHardeningLaw();


    /// Copy constructor.
    NonLinearIsotropicKinematicHardeningLaw(NonLinearIsotropicKinematicHardeningLaw const& rOther);

    /// Assignment operator.
    NonLinearIsotropicKinematicHardeningLaw& operator=(NonLinearIsotropicKinematicHardeningLaw const& rOther);

    /// Destructor.
    ~NonLinearIsotropicKinematicHardeningLaw();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual double& CalculateHardening(double &rHardening, const double &rAlpha, double rTemperature = 0);
  
    virtual double& CalculateIsotropicHardening(double &rIsotropicHardening, const double &rAlpha, double rTemperature = 0);

    virtual double& CalculateKinematicHardening(double &rKinematicHardening, const double &rAlpha, double rTemperature = 0);

    virtual double& CalculateDeltaHardening(double &rDeltaHardening, const double &rAlpha, double rTemperature = 0);

    virtual double& CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const double &rAlpha, double rTemperature = 0);

    virtual double& CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening, const double &rAlpha, double rTemperature = 0);

    virtual double CalculateThermalReferenceEffect(const double &rTemperature);

    virtual double CalculateThermalCurrentEffect(const double &rTemperature);

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


    /**
     * Pure isotropic hardening Theta=1;  pure kinematic hardening Theta= 0; combined isotropic-kinematic 0<Theta<1
     */

    double mTheta; 
	
     
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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class NonLinearIsotropicKinematicHardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// // input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   NonLinearIsotropicKinematicHardeningLaw& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const NonLinearIsotropicKinematicHardeningLaw& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NON_LINEAR_ISOTROPIC_KINEMATIC_HARDENING_LAW_H_INCLUDED  defined 


