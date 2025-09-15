//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NEW_TRESCA_YIELD_CRITERION_H_INCLUDED )
#define      KRATOS_NEW_TRESCA_YIELD_CRITERION_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_yield_criteria/tresca_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.hpp"
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
   /* Implementation of the yield surface generalized treca yield surface (Krabbenhoft & Lyamin, 2015)
      The sharp edges of the yield surface are smoothed according to Gesto, Gens & Vaunat (2011)
      the gradient is botained with numerical differentiation (Perez-Foguet et al, 2000)
      */

/** Detail class definition.
*/
class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) NewTrescaYieldCriterion
	: public TrescaYieldCriterion 
{
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of MisesHuberYieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION( NewTrescaYieldCriterion );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        NewTrescaYieldCriterion();

        /// Initialization constructor.
        NewTrescaYieldCriterion(HardeningLawPointer pHardeningLaw);

        /// Copy constructor.
        NewTrescaYieldCriterion(NewTrescaYieldCriterion const& rOther);

        /// Assignment operator.
        NewTrescaYieldCriterion& operator=(NewTrescaYieldCriterion const& rOther);


        /// Destructor.
        virtual ~NewTrescaYieldCriterion();


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        double& CalculateYieldCondition(double & rStateFunction, const Vector& rStressVector, const double& rAlpha) override;

        void CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rFirstDerivative, const double& rAlpha) override;
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

#endif // KRATOS_NEW_TRESCA_YIELD_CRITERION_H_INCLUDED  defined 


