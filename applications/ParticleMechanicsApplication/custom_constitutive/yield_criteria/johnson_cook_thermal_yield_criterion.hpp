//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson (adapted from Solid Mechanics App)
//

#if !defined(KRATOS_JOHNSON_COOK_THERMAL_YIELD_CRITERION_H_INCLUDED )
#define  KRATOS_JOHNSON_COOK_THERMAL_YIELD_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/yield_criteria/particle_yield_criterion.hpp"

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
class JohnsonCookThermalYieldCriterion
	: public ParticleYieldCriterion
{
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of JohnsonCookThermalYieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookThermalYieldCriterion);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        JohnsonCookThermalYieldCriterion();

        /// Initialization constructor.
        JohnsonCookThermalYieldCriterion(HardeningLawPointer pHardeningLaw);

        /// Copy constructor.
        JohnsonCookThermalYieldCriterion(JohnsonCookThermalYieldCriterion const& rOther);

        /// Assignment operator.
        JohnsonCookThermalYieldCriterion& operator=(JohnsonCookThermalYieldCriterion const& rOther);


        /// Destructor.
        ~JohnsonCookThermalYieldCriterion() override;


        ///@}
        ///@name Operators
        ///@{

        /**
	 * Clone function (has to be implemented by any derived class)
	 * @return a pointer to a new instance of this yield criterion
	 */
        ParticleYieldCriterion::Pointer Clone() const override;

        ///@}
        ///@name Operations
        ///@{

        double& CalculateYieldCondition(double& rStateFunction, const Parameters& rValues) override;

        double& CalculateStateFunction(double& rStateFunction, const Parameters& rValues) override;

        double& CalculateDeltaStateFunction(double& rDeltaStateFunction, const Parameters& rValues) override;

        double& CalculatePlasticDissipation(double & rPlasticDissipation, const Parameters& rValues) override;

        double& CalculateDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rValues) override;

        double& CalculateImplexPlasticDissipation(double & rPlasticDissipation, const Parameters& rValues) override;

        double& CalculateImplexDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rValues) override;

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

    }; // Class JohnsonCookThermalYieldCriterion

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{

}  // namespace Kratos.

#endif // KRATOS_JOHNSON_COOK_THERMAL_YIELD_CRITERION_H_INCLUDED  defined
