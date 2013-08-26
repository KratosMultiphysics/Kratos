//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
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
class YieldCriterion
{
    public:
        ///@name Type Definitions
        ///@{
        typedef HardeningLaw::Pointer        HardeningLawPointer;

        /// Pointer definition of YieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION(YieldCriterion);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        YieldCriterion()
	{
	  //KRATOS_ERROR(std::logic_error, "calling the default constructor in YieldCriterion ... illegal operation!!","");
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


        ///@}
        ///@name Operations
        ///@{
	void InitializeMaterial (HardeningLawPointer pHardeningLaw)
	{
		mpHardeningLaw = pHardeningLaw;
	}


        virtual double& CalculateYieldCondition(double & rStateFunction, const double& rNormStress, const double& rAlpha)
	{
		KRATOS_ERROR(std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!","");

		return rStateFunction;
	};

        virtual double& CalculateYieldCondition(double & rStateFunction, const Matrix& rStressMatrix, const double& rAlpha)
	{
		KRATOS_ERROR(std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!","");

		return rStateFunction;
	};


	virtual double& CalculateStateFunction(double & rStateFunction,const double& rNormStress, const double & rDeltaGamma, const double& rLameMu_bar, const double& rAlpha, const double& rAlphaOld)
	{
		KRATOS_ERROR(std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!","");

		return rStateFunction;
	};

	virtual double& CalculateDeltaStateFunction(double & rDeltaStateFunction, const double& rLameMu_bar, const double& rAlpha)
	{
		KRATOS_ERROR(std::logic_error, "calling the base class function in YieldCriterion ... illegal operation!!","");

		return rDeltaStateFunction;
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
	};

	virtual void load(Serializer& rSerializer)
	{
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


