//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//			Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/custom_yield_criteria/yield_criterion.hpp


#if !defined(KRATOS_MPM_YIELD_CRITERION_H_INCLUDED)
#define      KRATOS_MPM_YIELD_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "custom_constitutive/hardening_laws/mpm_hardening_law.hpp"

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
class KRATOS_API(MPM_APPLICATION) MPMYieldCriterion
{
    public:
        struct Parameters
        {
        private:

        const double* mpStressNorm;
        MPMHardeningLaw::Parameters HardeningParameters;

        public:
	  //Set Parameters
	  void SetStressNorm  (const double& rStressNorm)  { mpStressNorm = &rStressNorm; };

	  //Get Parameters
	  const double& GetStressNorm  () const { return *mpStressNorm;  };
	  const MPMHardeningLaw::Parameters& GetHardeningParameters  () const { return HardeningParameters; };

	  //Set Hardening Parameters
	  void SetRateFactor  (double rRateFactor)         { HardeningParameters.SetRateFactor(rRateFactor);   };
	  void SetDeltaTime   (const double& rDeltaTime)   { HardeningParameters.SetDeltaTime(rDeltaTime);     };

	  //Get Hardening Parameters
	  const double& GetRateFactor  () const { return HardeningParameters.GetRateFactor();   };
	  const double& GetDeltaTime   () const { return HardeningParameters.GetDeltaTime();    };

	};

        ///@name Type Definitions
        ///@{
        typedef MPMHardeningLaw::Pointer        HardeningLawPointer;

        /// Pointer definition of MPMYieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION( MPMYieldCriterion );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        MPMYieldCriterion()
	{
        };

        /// Initialization constructor.
        MPMYieldCriterion(HardeningLawPointer pHardeningLaw)
	:mpHardeningLaw(pHardeningLaw)
	{
	};

        /// Copy constructor.
        MPMYieldCriterion(MPMYieldCriterion const& rOther)
	:mpHardeningLaw(rOther.mpHardeningLaw)
	{

	};

        /// Assignment operator.
        MPMYieldCriterion& operator=(MPMYieldCriterion const& rOther)
	{
        mpHardeningLaw = rOther.mpHardeningLaw;
        return *this;
	}

        /// Destructor.
        virtual ~MPMYieldCriterion() {};


        ///@}
        ///@name Operators
        ///@{

        /**
	 * Clone function (has to be implemented by any derived class)
	 * @return a pointer to a new instance of this yield criterion
	 */
        virtual MPMYieldCriterion::Pointer Clone() const
        {
          return Kratos::make_shared<MPMYieldCriterion>(*this);
	}


        ///@}
        ///@name Operations
        ///@{
        void InitializeMaterial (HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties)
	{
	        mpHardeningLaw = pHardeningLaw;
	}


        void SetHardeningLaw(MPMHardeningLaw& rHardeningLaw)
        {
	  mpHardeningLaw = (HardeningLawPointer) (&rHardeningLaw);
        }

        void pSetHardeningLaw(HardeningLawPointer& pHardeningLaw)
        {
	  mpHardeningLaw = pHardeningLaw;
        }

        MPMHardeningLaw& GetHardeningLaw()
        {
	  return *mpHardeningLaw;
        }

        HardeningLawPointer pGetHardeningLaw()
        {
	  return mpHardeningLaw;
        }


        /*
        * @brief This function return the yield criterion at the given principal stress condition
        * @param[in/out] rStateFunction Computed yield criterion
        * @param[in] rVariables Parameter of yields
        * @param[in] rPrincipalStress Current principal stress
        * @param[in] rAlpha Used parameters
        * @param[in] rBeta Used parameters
        * @return Yield criterion
        */
        virtual double& CalculateYieldCondition(double & rStateFunction, const Parameters& rVariables, const Properties& rProp)
	{
		KRATOS_ERROR << "Calling the base class function in MPMYieldCriterion ... illegal operation!!" << std::endl;
		return rStateFunction;
	};

        virtual double& CalculateYieldCondition(double& rStateFunction, const Vector& rPrincipalStress, const double& rAlpha, const Properties& rProp)
	{
		KRATOS_ERROR << "Calling the base class function in MPMYieldCriterion ... illegal operation!!" << std::endl;
        };

        virtual double& CalculateYieldCondition(double& rStateFunction, const Vector& rPrincipalStress, const double& rAlpha, const double& rBeta, const Properties& rProp)
	{
		KRATOS_ERROR << "Calling the base class function in MPMYieldCriterion ... illegal operation!!" << std::endl;
        };


        /*
        * @brief This function return the first derivative of yield criterion at the given principal stress condition
        * @param[in] rPrincipalStress Principal stresses
        * @param[in/out] rFirstDerivative First stress derivative value of yield function
        * @param[in] rAlpha Used parameters
        * @param[in] rBeta Used parameters
        */
	virtual void CalculateYieldFunctionDerivative(const Vector& rPrincipalStress, Vector& rFirstDerivative, const Properties& rProp)
	{
		KRATOS_ERROR << "Calling the base class function in MPMYieldCriterion ... illegal operation!!" << std::endl;

        };

	virtual void CalculateYieldFunctionDerivative(const Vector& rPrincipalStress, Vector& rFirstDerivative, const double& rAlpha, const Properties& rProp)
	{
		KRATOS_ERROR << "Calling the base class function in MPMYieldCriterion ... illegal operation!!" << std::endl;
        };

        virtual void CalculateYieldFunctionDerivative(const Vector& rPrincipalStress, Vector& rFirstDerivative, const double& rAlpha, const double& rBeta, const Properties& rProp)
	{
		KRATOS_ERROR << "Calling the base class function in MPMYieldCriterion ... illegal operation!!" << std::endl;
        };


        /*
        * @brief This function return the second derivative of yield criterion at the given principal stress condition
        * @param[in] rStressVector Principal stresses
        * @param[in/out] rSecondDerivative Second stress derivative value of yield function
        */
        virtual void CalculateYieldFunctionSecondDerivative(const Vector& rPrincipalStress, Vector& rSecondDerivative, const Properties& rProp)
	{
		KRATOS_ERROR << "Calling the base class function in MPMYieldCriterion ... illegal operation!!" << std::endl;
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

    }; // Class MPMYieldCriterion

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    // /// input stream function
    // inline std::istream& operator >> (std::istream& rIStream,
    //                                   MPMYieldCriterion& rThis);

    // /// output stream function
    // inline std::ostream& operator << (std::ostream& rOStream,
    //                                   const MPMYieldCriterion& rThis)
    // {
    //     rThis.PrintInfo(rOStream);
    //     rOStream << std::endl;
    //     rThis.PrintData(rOStream);

    //     return rOStream;
    // }
    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MPM_YIELD_CRITERION_H_INCLUDED  defined
