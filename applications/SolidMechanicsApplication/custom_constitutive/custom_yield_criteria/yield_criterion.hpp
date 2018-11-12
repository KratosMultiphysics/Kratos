//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_YIELD_CRITERION_H_INCLUDED )
#define  KRATOS_YIELD_CRITERION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
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
class KRATOS_API(SOLID_MECHANICS_APPLICATION) YieldCriterion
{
    public:


	struct Parameters
	{
	private:

	  const double* mpStressNorm;
	  HardeningLaw::Parameters HardeningParameters;

	public:
    
	  //Set Parameters
	  void SetStressNorm  (const double& rStressNorm)  { mpStressNorm = &rStressNorm; };

	  //Get Parameters
	  const double& GetStressNorm  () const { return *mpStressNorm;  };

	  //Set Hardening Parameters
	  void SetRateFactor  (double rRateFactor)         { HardeningParameters.SetRateFactor(rRateFactor);   };
	  void SetDeltaGamma  (const double& rDeltaGamma)  { HardeningParameters.SetDeltaGamma(rDeltaGamma);   };
	  void SetLameMu_bar  (const double& rLameMu_bar)  { HardeningParameters.SetLameMu_bar(rLameMu_bar);   };
	  void SetDeltaTime   (const double& rDeltaTime)   { HardeningParameters.SetDeltaTime(rDeltaTime);     };
	  void SetTemperature (const double& rTemperature) { HardeningParameters.SetTemperature(rTemperature); };
		
	  void SetEquivalentPlasticStrain       (const double& rEquivalentPlasticStrain)       {  HardeningParameters.SetEquivalentPlasticStrain(rEquivalentPlasticStrain);       };
	  void SetEquivalentPlasticStrainOld    (const double& rEquivalentPlasticStrainOld)    {  HardeningParameters.SetEquivalentPlasticStrainOld(rEquivalentPlasticStrainOld); };

	  void SetCharacteristicSize (const double& rCharacteristicSize) {HardeningParameters.SetCharacteristicSize(rCharacteristicSize);}
        
	  void SetStrainMatrix (const Matrix& rStrainMatrix) { HardeningParameters.SetStrainMatrix(rStrainMatrix); }
	  void SetStressMatrix (const Matrix& rStressMatrix) { HardeningParameters.SetStressMatrix(rStressMatrix); }

	  //Get Hardening Parameters
	  const double& GetRateFactor  () const { return HardeningParameters.GetRateFactor();   };
	  const double& GetDeltaGamma  () const { return HardeningParameters.GetDeltaGamma();   };
	  const double& GetLameMu_bar  () const { return HardeningParameters.GetLameMu_bar();   };
	  const double& GetDeltaTime   () const { return HardeningParameters.GetDeltaTime();    };
	  const double& GetTemperature () const { return HardeningParameters.GetTemperature();  };
		
	  const double& GetEquivalentPlasticStrain       () const { return HardeningParameters.GetEquivalentPlasticStrain();    };
	  const double& GetEquivalentPlasticStrainOld    () const { return HardeningParameters.GetEquivalentPlasticStrainOld(); };
		
	  const HardeningLaw::Parameters& GetHardeningParameters  () const { return HardeningParameters; };

	  const double& GetCharacteristicSize () const {return HardeningParameters.GetCharacteristicSize();}
        
	  const Matrix& GetStrainMatrix () const { return HardeningParameters.GetStrainMatrix(); }
	  const Matrix& GetStressMatrix () const { return HardeningParameters.GetStressMatrix(); }
	};

        ///@name Type Definitions
        ///@{
        typedef HardeningLaw::Pointer        HardeningLawPointer;
        typedef HardeningLaw::PlasticVariables        PlasticVariablesType;

        /// Pointer definition of YieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION( YieldCriterion );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        YieldCriterion()
	{
	  //KRATOS_THROW_ERROR( std::logic_error, "calling the default constructor in YieldCriterion ... illegal operation!!", "" )
	};


        /// Initialization constructor.
        YieldCriterion(HardeningLawPointer pHardeningLaw)
	:mpHardeningLaw(pHardeningLaw)
	{
	};

        /// Copy constructor.
        YieldCriterion(YieldCriterion const& rOther)
	:mpHardeningLaw(rOther.mpHardeningLaw)
	{

	};

        /// Assignment operator.
        YieldCriterion& operator=(YieldCriterion const& rOther)
	{
		mpHardeningLaw = rOther.mpHardeningLaw;
		return *this;
	}


        /// Destructor.
        virtual ~YieldCriterion() {};


        ///@}
        ///@name Operators
        ///@{

        /**
	 * Clone function (has to be implemented by any derived class)
	 * @return a pointer to a new instance of this yield criterion
	 */
        virtual YieldCriterion::Pointer Clone() const
        {
	  YieldCriterion::Pointer p_clone(new YieldCriterion(*this));
	  return p_clone;
	}


        ///@}
        ///@name Operations
        ///@{
        void InitializeMaterial (HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties)
	{
	        mpHardeningLaw = pHardeningLaw;
		mpHardeningLaw->InitializeMaterial(rMaterialProperties);
	}


        void SetHardeningLaw(HardeningLaw& rHardeningLaw)
        {      
	  mpHardeningLaw = (HardeningLawPointer) (&rHardeningLaw);
        }

        void pSetHardeningLaw(HardeningLawPointer& pHardeningLaw)
        {      
	  mpHardeningLaw = pHardeningLaw;
        }

        HardeningLaw& GetHardeningLaw()
        {      
	  return *mpHardeningLaw;
        }

        HardeningLawPointer pGetHardeningLaw()
        {      
	  return mpHardeningLaw;
        }

 
        virtual double& CalculateYieldCondition(double & rStateFunction, const Parameters& rVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )

		return rStateFunction;
	};

        virtual double& CalculateStateFunction(double & rStateFunction, const Parameters& rVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )

		return rStateFunction;
	};

        virtual double& CalculateDeltaStateFunction(double & rDeltaStateFunction, const Parameters& rVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )

		return rDeltaStateFunction;
	};


        virtual double& CalculatePlasticDissipation(double & rPlasticDissipation, const Parameters& rVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )

		return rPlasticDissipation;
	};


        virtual double& CalculateDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )

		return rDeltaPlasticDissipation;
	};


        virtual double& CalculateImplexPlasticDissipation(double & rPlasticDissipation, const Parameters& rVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )

		return rPlasticDissipation;
	};


        virtual double& CalculateImplexDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )

		return rDeltaPlasticDissipation;
	};


//Ll:
	virtual void CalculateYieldFunctionDerivative(const Vector& rPrincipalStress, Vector& rFirstDerivative)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )
  };

	virtual void CalculateYieldFunctionDerivative(const Vector& rPrincipalStress, Vector& rFirstDerivative, const double& rAlpha)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )
  };

	virtual double& CalculateYieldCondition(double& rStateFunction, const Vector& rPrincipalStress, const double& rAlpha)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )
  };

//Laurin:
	virtual double EvaluateThirdInvariantEffectMC( const double& rLodeAngle)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )
  };
  
  virtual double& CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )
  };
  
  virtual double& CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const PlasticVariablesType& rPlasticVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )
  };
  
  virtual void CalculateYieldFunctionDerivative(const Vector& rPrincipalStress, Vector& rFirstDerivative, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )
  };
  
  virtual void CalculateYieldFunctionDerivative(const Vector& rPrincipalStress, Vector& rFirstDerivative, const PlasticVariablesType& rPlasticVariables)
	{
		KRATOS_THROW_ERROR( std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!", "" )
  };
  
  

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
	
	HardeningLawPointer mpHardeningLaw;
	
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
	  rSerializer.save("mpHardeningLaw",mpHardeningLaw);
	};

	virtual void load(Serializer& rSerializer)
	{
	  rSerializer.load("mpHardeningLaw",mpHardeningLaw);
	};

        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        ///@}

    }; // Class YieldCriterion

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    // /// input stream function
    // inline std::istream& operator >> (std::istream& rIStream,
    //                                   YieldCriterion& rThis);

    // /// output stream function
    // inline std::ostream& operator << (std::ostream& rOStream,
    //                                   const YieldCriterion& rThis)
    // {
    //     rThis.PrintInfo(rOStream);
    //     rOStream << std::endl;
    //     rThis.PrintData(rOStream);

    //     return rOStream;
    // }
    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_YIELD_CRITERION_H_INCLUDED  defined 


