//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CASM_CEM_YIELD_CRITERION_H_INCLUDED )
#define      KRATOS_CASM_CEM_YIELD_CRITERION_H_INCLUDED



// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_yield_criteria/yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/casm_cem_hardening_law.hpp"
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
class CasmCemYieldCriterion
	: public YieldCriterion 
{
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of MisesHuberYieldCriterion
        KRATOS_CLASS_POINTER_DEFINITION( CasmCemYieldCriterion );

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        CasmCemYieldCriterion();

        /// Copy constructor.
        CasmCemYieldCriterion(CasmCemYieldCriterion const& rOther);

        /// Initialization constructor.
        CasmCemYieldCriterion(HardeningLawPointer pHardeningLaw);

        /// Assignment operator.
        CasmCemYieldCriterion& operator=(CasmCemYieldCriterion const& rOther);


        /// Destructor.
        virtual ~CasmCemYieldCriterion();


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        //double& CalculateYieldCondition( double& rStateFunction, const Vector& rStressVector, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum);
				
				double& CalculateYieldCondition( double& rStateFunction, const Vector& rStressVector, const PlasticVariablesType& rPlasticVariables);
				
				//void CalculateYieldFunctionDerivative( const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha, const double& rBeta, const double& rAlphaCum, const double& rBetaCum);

				void CalculateYieldFunctionDerivative( const Vector& rStressVector, Vector& rYieldFunctionD, const PlasticVariablesType& rPlasticVariables);

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
        
        
        double EvaluateThirdInvariantEffectSheng( const double& rLodeAngle );

        void CalculateAndAddThirdInvDerivativeMC(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rPt);
        
        double EvaluateThirdInvariantEffectMC( const double& rLodeAngle );
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

        double GetSmoothingLodeAngle();
 
        void GetSmoothingConstants(double& rA, double& rB, const double& rLode);

        ///@}
        ///@name Private  Access
        ///@{

	
	///@}
	///@name Serialization
	///@{
	friend class Serializer;

	// A private default constructor necessary for serialization

	virtual void save(Serializer& rSerializer) const;

	virtual void load(Serializer& rSerializer);

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

#endif // KRATOS_CASM_CEM_YIELD_CRITERION_H_INCLUDED  defined 


