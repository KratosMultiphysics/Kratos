//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_HARDENING_LAW_H_INCLUDED )
#define  KRATOS_HARDENING_LAW_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"


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
class HardeningLaw
{
public:
    ///@name Type Definitions
    ///@{

    typedef Properties::Pointer            PropertiesPointer;


    /// Pointer definition of HardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION(HardeningLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HardeningLaw();

    /// Copy constructor.
    HardeningLaw(HardeningLaw const& rOther);

    /// Assignment operator.
    HardeningLaw& operator=(HardeningLaw const& rOther);

    /// Destructor.
   virtual ~HardeningLaw() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    void InitializeMaterial (PropertiesPointer pProperties)
	{
		mpProperties = pProperties;
	}


    Properties& GetProperties()
	{
		return *mpProperties;
	}

    virtual double& CalculateHardening(double &Hardening, const double & rAlpha){};
  
    virtual double& CalculateIsotropicHardening(double &IsotropicHardening, const double & rAlpha){};

    virtual double& CalculateKinematicHardening(double &KinematicHardening, const double & rAlpha){};


    virtual double& CalculateDeltaHardening(double &DeltaHardening, const double & rAlpha){};

    virtual double& CalculateDeltaIsotropicHardening(double &DeltaIsotropicHardening, const double & rAlpha){};

    virtual double& CalculateDeltaKinematicHardening(double &DeltaKinematicHardening, const double & rAlpha){};

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


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

    const PropertiesPointer mpProperties;
     
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

    virtual void save(Serializer& rSerializer) const
    {
	    rSerializer.save("Properties",mpProperties);
    };

    virtual void load(Serializer& rSerializer)
    {
	    rSerializer.load("Properties",mpProperties);
    };

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class HardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  HardeningLaw& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const HardeningLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HARDENING_LAW_H_INCLUDED  defined 


