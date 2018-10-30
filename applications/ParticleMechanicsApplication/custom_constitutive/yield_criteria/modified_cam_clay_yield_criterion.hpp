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


#if !defined(KRATOS_MODIFIED_CAM_CLAY_YIELD_CRITERION_H_INCLUDED )
#define      KRATOS_MODIFIED_CAM_CLAY_YIELD_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/yield_criteria/MPM_yield_criterion.hpp"
#include "custom_constitutive/hardening_laws/cam_clay_hardening_law.hpp"

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
class ModifiedCamClayYieldCriterion
	: public MPMYieldCriterion
{
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of MisesHuberYieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION( ModifiedCamClayYieldCriterion );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        ModifiedCamClayYieldCriterion();

        /// Copy constructor.
        ModifiedCamClayYieldCriterion(ModifiedCamClayYieldCriterion const& rOther);

        /// Initialization constructor.
        ModifiedCamClayYieldCriterion(HardeningLawPointer pHardeningLaw);

        /// Assignment operator.
        ModifiedCamClayYieldCriterion& operator=(ModifiedCamClayYieldCriterion const& rOther);


        /// Destructor.
        virtual ~ModifiedCamClayYieldCriterion();


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{


        /*
        * @brief This function return the Modified cam clay yield criterion at the given principal stress condition
        * @param[in/out] rStateFunction Mohr coulomb yield criterion
        * @param[in] rStressVector Principal stresses
        * @param[in] rAlpha Plastic volumetric strain
        * @param[in] rOldPreconsolidationPressure The value of Preconsolidation Stress at the previous time step
        * @return Modified cam clay yield criterion
        */
        double& CalculateYieldCondition(double & rStateFunction, const Vector& rStressVector, const double& rAlpha, const double& rOldPreconsolidationPressure) override;


        /*
        * @brief This function return the first derivative of modified cam clay yield criterion at the given principal stress condition
        * @param[in] rStressVector Principal stresses
        * @param[in/out] rFirstDerivative First stress derivative value of yield function
        * @param[in] rAlpha Plastic volumetric strain
        * @param[in] rOldPreconsolidationPressure The value of Preconsolidation Stress at the previous time step
        * @return Modified cam clay yield criterion first derivative
        */
        void CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rFirstDerivative, const double& rAlpha, const double& rOldPreconsolidationPressure) override;


        /*
        * @brief This function return the second derivative of modified cam clay yield criterion at the given principal stress condition
        * @param[in] rStressVector Principal stresses
        * @param[in/out] rSecondDerivative Second stress derivative value of yield function
        * @return Modified cam clay yield criterion second derivative
        */
        void CalculateYieldFunctionSecondDerivative(const Vector& rStressVector, Vector& rSecondDerivative) override;

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
        // std::string Info() const override;

        // /// Print information about this object.
        // void PrintInfo(std::ostream& rOStream) const override;

        // /// Print object's data.
        // void PrintData(std::ostream& rOStream) const override;


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

        double GetPI();

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

    }; // Class MisesHuberYieldCriterion

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    // /// input stream function
    // inline std::istream& operator >> (std::istream& rIStream,
    //                                   MisesHuberYieldCriterion& rThis);

    // /// output stream function
    // inline std::ostream& operator << (std::ostream& rOStream,
    //                                   const MisesHuberYieldCriterion& rThis)
    // {
    //     rThis.PrintInfo(rOStream);
    //     rOStream << std::endl;
    //     rThis.PrintData(rOStream);

    //     return rOStream;
    // }
    ///@}

    ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODIFIED_CAM_CLAY_YIELD_CRITERION_H_INCLUDED  defined


