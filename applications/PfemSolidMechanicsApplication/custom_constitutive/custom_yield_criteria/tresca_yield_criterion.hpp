//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_TRESCA_YIELD_CRITERION_H_INCLUDED )
#define      KRATOS_TRESCA_YIELD_CRITERION_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_yield_criteria/yield_criterion.hpp"
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
   struct TrescaStressInvariants {

       double MeanStress;
       double J2InvSQ;
       double LodeAngle;

       double A;
       double B;
       double C;

   };



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
class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) TrescaYieldCriterion
	: public YieldCriterion 
{
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of MisesHuberYieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION( TrescaYieldCriterion );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        TrescaYieldCriterion();

        /// Initialization constructor.
        TrescaYieldCriterion(HardeningLawPointer pHardeningLaw);

        /// Copy constructor.
        TrescaYieldCriterion(TrescaYieldCriterion const& rOther);

        /// Assignment operator.
        TrescaYieldCriterion& operator=(TrescaYieldCriterion const& rOther);


        /// Destructor.
        virtual ~TrescaYieldCriterion();


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        double& CalculateYieldCondition(double & rStateFunction, const Vector& rStressVector, const double& rAlpha) override;

	//void CalculateYieldFunctionSecondDerivative(const Vector& rStressVector, Matrix& SecondDerivative, const double& rAlpha);
	
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

        void CalculateSmoothingInvariants(TrescaStressInvariants& rStressInvariants);

        void CalculateStressInvariants( const Vector& rStressVector, TrescaStressInvariants& rStressInvariants);

        double GetSmoothingLodeAngle();

        double GetPI();

        void ComputeC2andC3Vector( const Vector& rStressVector, TrescaStressInvariants& rStressInvariants, Vector& C2Vector, Vector& C3Vector);
        
        void ComputeC2andC3VectorDD( const Vector& rStressVector, TrescaStressInvariants& rStressInvariants, Vector& C2Vector, Vector& C3Vector, Matrix& C2Matrix, Matrix& C3Matrix, Vector & rLodeDerivative);

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

#endif // KRATOS_TRESCA_YIELD_CRITERION_H_INCLUDED  defined 


