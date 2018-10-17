//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/custom_hardening_laws/hardening_law.hpp


#if !defined(KRATOS_MPM_HARDENING_LAW_H_INCLUDED)
#define      KRATOS_MPM_HARDENING_LAW_H_INCLUDED

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
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) MPMHardeningLaw
{
public:

	  ///@name Type Definitions
    ///@{

    typedef const Properties*   PropertiesPointer;


    /// Pointer definition of MPMHardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION( MPMHardeningLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMHardeningLaw(){ mpProperties = NULL; };

    /// Copy constructor.
    MPMHardeningLaw(MPMHardeningLaw const& rOther)
    :mpProperties(rOther.mpProperties)
    {};

    /// Assignment operator.
    MPMHardeningLaw& operator=(MPMHardeningLaw const& rOther)
    {
      mpProperties = rOther.mpProperties;
      return *this;
    };

    /// Destructor.
    virtual ~MPMHardeningLaw() {};


    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this hardening law
     */
    virtual MPMHardeningLaw::Pointer Clone() const
    {
      return Kratos::make_shared<MPMHardeningLaw>(*this);
    }

    ///@}
    ///@name Operations
    ///@{
    void InitializeMaterial (const Properties& rMaterialProperties)
	{
	  SetProperties(rMaterialProperties);
	}


    void SetProperties (const Properties& rMaterialProperties)
	{
	  mpProperties = (PropertiesPointer)(&rMaterialProperties);
	}


    const Properties& GetProperties()
	{
	  return *mpProperties;
	}


    virtual double& CalculateHardening(double &rHardening, const Parameters& rValues){ return rHardening; };

    virtual double& CalculateHardening(double &rHardening, const double& rAlpha, const double rTemperature = 0.0) {return rHardening; };

    virtual double& CalculateHardening(double &rHardening, const double& rAlpha, const Variable<double>& rThisVariable) {return rHardening;};
    
    virtual double& CalculateIsotropicHardening(double &rIsotropicHardening, const Parameters& rValues){ return rIsotropicHardening; };

    virtual double& CalculateKinematicHardening(double &rKinematicHardening, const Parameters& rValues){ return rKinematicHardening; };

    virtual double& CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues){ return rDeltaHardening; };

    virtual double& CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const Parameters& rValues){ return rDeltaIsotropicHardening; };

    virtual double& CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening, const Parameters& rValues){ return rDeltaKinematicHardening; };

    virtual double& CalculateDeltaThermalHardening(double &rDeltaThermalHardening, const Parameters& rValues){ return rDeltaThermalHardening; };

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

    PropertiesPointer mpProperties;

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
      //Properties can not be stored in serializer
      //because Properties have a ConstitutiveLaw pointer
      //when the constitutive law pointer is called to be saved
      //a recursive call is done if properties are saved.
      //rSerializer.save("Properties",mpProperties);
    };

    virtual void load(Serializer& rSerializer)
    {
      //Properties* pProperties;
      //rSerializer.load("Properties",pProperties);
      //mpProperties = pProperties;
    };

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class MPMHardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MPMHardeningLaw& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MPMHardeningLaw& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MPM_HARDENING_LAW_H_INCLUDED  defined
