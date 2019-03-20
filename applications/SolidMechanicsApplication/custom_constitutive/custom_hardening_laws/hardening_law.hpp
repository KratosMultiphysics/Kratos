//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_HARDENING_LAW_H_INCLUDED)
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
class KRATOS_API(SOLID_MECHANICS_APPLICATION) HardeningLaw
{
public:

	struct Parameters
	{
	private:

        double        mRateFactor;

        const double* mpDeltaGamma;
        const double* mpLameMu_bar;

        const double* mpDeltaTime;
        const double* mpTemperature;

        const double* mpEquivalentPlasticStrain;
        const double* mpEquivalentPlasticStrainOld;

        const double* mpCharacteristicSize;

        const Matrix* mpStrainMatrix;
        const Matrix* mpStressMatrix;

	public:

        //Set Parameters
        void SetRateFactor  (double rRateFactor)         { mRateFactor = rRateFactor;     };
        void SetDeltaGamma  (const double& rDeltaGamma)  { mpDeltaGamma = &rDeltaGamma;   };
        void SetLameMu_bar  (const double& rLameMu_bar)  { mpLameMu_bar = &rLameMu_bar;   };
        void SetDeltaTime   (const double& rDeltaTime)   { mpDeltaTime = &rDeltaTime;     };
        void SetTemperature (const double& rTemperature) { mpTemperature = &rTemperature; };

        void SetEquivalentPlasticStrain    (const double& rEquivalentPlasticStrain)    { mpEquivalentPlasticStrain = &rEquivalentPlasticStrain;       };
        void SetEquivalentPlasticStrainOld (const double& rEquivalentPlasticStrainOld) { mpEquivalentPlasticStrainOld = &rEquivalentPlasticStrainOld; };

        void SetCharacteristicSize (const double& rCharacteristicSize) {mpCharacteristicSize = &rCharacteristicSize;}

        void SetStrainMatrix (const Matrix& rStrainMatrix) {mpStrainMatrix = &rStrainMatrix;}
        void SetStressMatrix (const Matrix& rStressMatrix) {mpStressMatrix = &rStressMatrix;}

        //Get Parameters
        const double& GetRateFactor  () const { return  mRateFactor;   };
        const double& GetDeltaGamma  () const { return *mpDeltaGamma;  };
        const double& GetLameMu_bar  () const { return *mpLameMu_bar;  };
        const double& GetDeltaTime   () const { return *mpDeltaTime;   };
        const double& GetTemperature () const { return *mpTemperature; };

        const double& GetEquivalentPlasticStrain       () const { return *mpEquivalentPlasticStrain;       };
        const double& GetEquivalentPlasticStrainOld    () const { return *mpEquivalentPlasticStrainOld;    };

        const double& GetCharacteristicSize () const {return *mpCharacteristicSize; }

        const Matrix& GetStrainMatrix () const { return *mpStrainMatrix; }
        const Matrix& GetStressMatrix () const { return *mpStressMatrix; }


	void print() const
	  {
	    std::cout<<" RateFactor "<<mRateFactor<<std::endl;
	    std::cout<<" DeltaGamma "<<*mpDeltaGamma<<std::endl;
	    std::cout<<" LameMubar "<<*mpLameMu_bar<<std::endl;
	    std::cout<<" DeltaTime "<<*mpDeltaTime<<std::endl;
	    std::cout<<" Temperature "<<*mpTemperature<<std::endl;
	    std::cout<<" EquivalentPlastic "<<*mpEquivalentPlasticStrain<<std::endl;
	    std::cout<<" mpEquivalentPlasticStrainOld "<<*mpEquivalentPlasticStrainOld<<std::endl;
	    // std::cout<<" CharacteristicSize "<<*mpCharacteristicSize<<std::endl;
	    // std::cout<<" StrainMatrix "<<*mpStrainMatrix<<std::endl;
	    // std::cout<<" StressMatrix "<<*mpStressMatrix<<std::endl;
	  }

	// /// Constructor
	// Parameters(){};

        // /// Copy constructor.
	// Parameters(Parameters const& rOther)
	//   :mRateFactor(rOther.mRateFactor)
	//   ,mpDeltaGamma(rOther.mpDeltaGamma)
	//   ,mpLameMu_bar(rOther.mpLameMu_bar)
	//   ,mpDeltaTime(rOther.mpDeltaTime)
	//   ,mpTemperature(rOther.mpTemperature)
	//   ,mpEquivalentPlasticStrain(rOther.mpEquivalentPlasticStrain)
	//   ,mpEquivalentPlasticStrainOld(rOther.mpEquivalentPlasticStrainOld)
	//   ,mpCharacteristicSize(rOther.mpCharacteristicSize)
	//   ,mpStrainMatrix(rOther.mpStrainMatrix)
	//   ,mpStressMatrix(rOther.mpStressMatrix)
	//   {

	//   };

	// /// Assignment operator.
        // Parameters& operator=(Parameters const& rOther)
	// {
	//   mRateFactor  = rOther.mRateFactor;
	//   mpDeltaGamma = rOther.mpDeltaGamma;
	//   mpLameMu_bar = rOther.mpLameMu_bar;
	//   mpDeltaTime  = rOther.mpDeltaTime;
	//   mpTemperature = rOther.mpTemperature;
	//   mpEquivalentPlasticStrain = rOther.mpEquivalentPlasticStrain;
	//   mpEquivalentPlasticStrainOld = rOther.mpEquivalentPlasticStrainOld;
	//   mpCharacteristicSize = rOther.mpCharacteristicSize;
	//   mpStrainMatrix = rOther.mpStrainMatrix;
	//   mpStressMatrix = rOther.mpStressMatrix;
	//   return *this;
	// }

	// /// Destructor.
	// virtual ~Parameters() {};


	};


    ///@name Type Definitions
    ///@{

    typedef const Properties*   PropertiesPointer;


    /// Pointer definition of HardeningLaw
    KRATOS_CLASS_POINTER_DEFINITION( HardeningLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HardeningLaw(){ mpProperties = NULL; };

    /// Copy constructor.
    HardeningLaw(HardeningLaw const& rOther)
    :mpProperties(rOther.mpProperties)
    {};

    /// Assignment operator.
    HardeningLaw& operator=(HardeningLaw const& rOther)
    {
      //this assignment operator do not exists for const Properties::Pointer
      //for this reason mpProperties is a const Properties*
      mpProperties = rOther.mpProperties;
      return *this;
    };

    /// Destructor.
    virtual ~HardeningLaw() {};


    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this hardening law
     */
    virtual HardeningLaw::Pointer Clone() const
    {
      return Kratos::make_shared<HardeningLaw>(*this);
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

    virtual double& CalculateIsotropicHardening(double &rIsotropicHardening, const Parameters& rValues){ return rIsotropicHardening; };

    virtual double& CalculateKinematicHardening(double &rKinematicHardening, const Parameters& rValues){ return rKinematicHardening; };


    virtual double& CalculateDeltaHardening(double &rDeltaHardening, const Parameters& rValues){ return rDeltaHardening; };

    virtual double& CalculateDeltaIsotropicHardening(double &rDeltaIsotropicHardening, const Parameters& rValues){ return rDeltaIsotropicHardening; };

    virtual double& CalculateDeltaKinematicHardening(double &rDeltaKinematicHardening, const Parameters& rValues){ return rDeltaKinematicHardening; };

    virtual double& CalculateDeltaThermalHardening(double &rDeltaThermalHardening, const Parameters& rValues){ return rDeltaThermalHardening; };

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

    // /// Turn back information as a string.
    // virtual std::string Info() const;

    // /// Print information about this object.
    // virtual void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const;


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

    virtual double CalculateThermalReferenceEffect(const double &rTemperature){ return 1; };

    virtual double CalculateThermalCurrentEffect(const double &rTemperature){ return 1; };

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

}; // Class HardeningLaw

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   HardeningLaw& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const HardeningLaw& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HARDENING_LAW_H_INCLUDED  defined
