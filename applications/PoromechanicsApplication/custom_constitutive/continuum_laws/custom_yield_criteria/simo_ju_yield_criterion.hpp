//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SIMO_JU_YIELD_CRITERION_H_INCLUDED)
#define  KRATOS_SIMO_JU_YIELD_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/continuum_laws/custom_yield_criteria/yield_criterion.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

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
class KRATOS_API(POROMECHANICS_APPLICATION) SimoJuYieldCriterion
	: public YieldCriterion
{
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of SimoJuYieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION( SimoJuYieldCriterion );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        SimoJuYieldCriterion();

        /// Initialization constructor.
        SimoJuYieldCriterion(HardeningLawPointer pHardeningLaw);

        /// Copy constructor.
        SimoJuYieldCriterion(SimoJuYieldCriterion const& rOther);

        /// Assignment operator.
        SimoJuYieldCriterion& operator=(SimoJuYieldCriterion const& rOther);


        /// Destructor.
        ~SimoJuYieldCriterion() override;


        ///@}
        ///@name Operators
        ///@{

        /**
	 * Clone function (has to be implemented by any derived class)
	 * @return a pointer to a new instance of this yield criterion
	 */
        YieldCriterion::Pointer Clone() const override;

        ///@}
        ///@name Operations
        ///@{

        double& CalculateYieldCondition(double & rStateFunction, const Parameters& rVariables) override;


        double& CalculateStateFunction(double & rStateFunction, const Parameters& rVariables) override;


        double& CalculateDeltaStateFunction(double & rDeltaStateFunction, const Parameters& rVariables) override;

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

    }; // Class SimoJuYieldCriterion

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SIMO_JU_YIELD_CRITERION_H_INCLUDED  defined
