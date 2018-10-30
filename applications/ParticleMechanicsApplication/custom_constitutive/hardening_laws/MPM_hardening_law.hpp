//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					    Kratos default license: kratos/license.txt
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
#include "input_output/logger.h"


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

    struct Parameters
    {
    private:

          double        mRateFactor;

          const double* mpDeltaTime;

    public:

          //Set Parameters
          void SetRateFactor  (double rRateFactor)         { mRateFactor = rRateFactor;     };
          void SetDeltaTime   (const double& rDeltaTime)   { mpDeltaTime = &rDeltaTime;     };

          //Get Parameters
          const double& GetRateFactor  () const { return  mRateFactor;   };
          const double& GetDeltaTime   () const { return *mpDeltaTime;   };

          void print() const
          {
            KRATOS_INFO("MPMHardeningLaw.Parameters") << " RateFactor " << mRateFactor  << std::endl;
            KRATOS_INFO("MPMHardeningLaw.Parameters") << " DeltaTime  " << *mpDeltaTime << std::endl;
          }

    };

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

    // Calculate Hardening Parameter
    virtual double& CalculateHardening(double &rHardening, const Parameters& rValues){ return rHardening; };

    virtual double& CalculateHardening(double &rHardening, const double& rAlpha, const double rTemperature = 0.0) {return rHardening; };

    virtual double& CalculateHardening(double &rHardening, const double& rAlpha, const Variable<double>& rThisVariable) {return rHardening;};

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
