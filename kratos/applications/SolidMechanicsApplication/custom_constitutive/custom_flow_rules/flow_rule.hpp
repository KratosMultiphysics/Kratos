//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_FLOW_RULE_H_INCLUDED )
#define  KRATOS_FLOW_RULE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_constitutive/custom_yield_criteria/yield_criterion.hpp"
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
	class FlowRule
	{
	public:
		///@name Type Definitions
		///@{

		struct PlasticFactors
		{
			double Beta1;
			double Beta2;
			double Beta3;
			double Beta4;	   

			Matrix  Normal;
			Matrix  Dev_Normal;
		};


		struct RadialReturnVariables
		{
			double NormIsochoricStress;
			double TrialStateFunction;
			dobule DeltaGamma;
			double LameMu_bar;
		}


		/// Pointer definition of FlowRule
			KRATOS_CLASS_POINTER_DEFINITION(FlowRule);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		FlowRule();

	
		/// Destructor.
		virtual ~FlowRule();


		///@}
		///@name Operators
		///@{


		///@}
		///@name Operations
		///@{
    
		void InitializeMaterial (YieldCriterion & rYieldCriterion, HardeningLaw & rHardeningLaw, const Properties& rProperties)
		{
			mpYieldCriterion = &rYieldCriterion;
			mpHardeningLaw   = &rHardeningLaw;
			mpProperties     = &rProperties;

			mpYieldCriterion.InitializeMaterial(rHardeningLaw, rProperties);		
		};

		Properties & GetProperties()
		{
			return *mpProperties;
		}
	

		void CalculateReturnMapping(  Matrix& rIsoStressMatrix, const double& rTrace_b_bar );


		void CalculateScalingFactors( const Matrix & rIsoStressMatrix, const double& rTrace_b_bar, PlasticFactors& rScalingFactors );

		///@}
		///@name Access
		///@{


		///@}
		///@name Inquiry
		///@{


		///@}
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const;

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const;

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const;


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
    
		YieldCriterion*     mpYieldCriterion;
		HardeningLaw*       mpHardeningLaw;
		Properties*         mpProperties;
	
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
		///@name Private Inquiry
		///@{


		///@}
		///@name Un accessible methods
		///@{

		/// Assignment operator.
		FlowRule& operator=(FlowRule const& rOther);

		/// Copy constructor.
		FlowRule(FlowRule const& rOther);


		///@}

	}; // Class FlowRule

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
	inline std::istream& operator >> (std::istream& rIStream,
					  FlowRule& rThis);

/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream,
					  const FlowRule& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
///@}

///@} addtogroup block



///@}
///@ Template Operations
///@{

	bool CalculateReturnMapping( Matrix& rIsoStressMatrix, const double& rTrace_b_bar)
	{

		//0.- Initialize Variables
		bool Plasticity = false;
		YieldCriterion::InternalVariables   PlasticVariables = YieldCriterion->GetInternalVariables();

		const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
		const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

		mReturnMappingVariables.DeltaGamma  = 0;		
		mReturnMappingVariables.LameMu_bar  = YoungModulus/(2*(1+PoissonCoefficient));
		mReturnMappingVariables.LameMu_bar *= (rTrace_b_bar/3.0);    
		
		//1.-Isochoric stress norm
		mReturnMappingVariables.NormIsochoricStress = sqrt(IsoStressMatrix( 0 , 0 )*IsoStressMatrix( 0 , 0 )+
								   IsoStressMatrix( 1 , 1 )*IsoStressMatrix( 1 , 1 )+
								   IsoStressMatrix( 2 , 2 )*IsoStressMatrix( 2 , 2 )+
								   2.0 * IsoStressMatrix( 0 , 1 )*IsoStressMatrix( 0 , 1 ) );
		//2.- Check yield condition
		mReturnMappingVariables.TrialStateFunction = YieldCriterion->CheckYieldCondition( PlasticVariables, mReturnMappingVariables, rIsoStressMatrix )

			if( TrialStateFunction > 0 )
			{
				rPlasticVariables.DeltaPlasticStrain = 0;
				Plasticity = false;
			}
			else
			{
				bool converged = this->CalculateConsistencyCondition( PlasticVariables );

				if(!converged)
					std::cout<<" ConstitutiveLaw did not converge "<<std::endl;

				//3.- Update back stress, plastic strain and stress
				this->Update( PlasticVariables, mReturnMappingVariables, IsoStressMatrix );

	
				Plasticity = true;
			}


		return Plasticity;
	};

	bool CalculateConstistencyCondition( YieldCriterion::InternalVariables& rPlasticVariables )
	{
		//Set convergence parameters
		unsigned int iter    = 0;
		double Tolerance     = 1e-5;
		double MaxIterations = 50;

		//Initialize Parameters used in the determination of the Delta Plastic Strain
		double KinematicHardening = 0;
		double IsotropicHardening = 0;

		double DeltaKinematicHardening = 0;
		double DeltaIsotropicHardening = 0;

		//start
		double DeltaFunction = 0;
		mReturnMappingVariables.DeltaGamma    = 0;

		double StateFunction                  = mReturnMappingVariables.TrialStateFunction;
		double InitialEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain ;
		double InitialKinematicHardening      = rHardeningLaw.CalculateKinematicHardening( InitialEquivalentPlasticStrain );

   
		while ( fabs(Function)>=Tolerance && iter<=MaxIterations)
		{
			//Calculate DeltaFunction:
			DeltaKinematicHardening = mpHardeningLaw->CalculateDeltaKinematicHardening( rPlasticVariables.EquivalentPlasticStrain, rProperties );
			DeltaIsotropicHardening = mpHardeningLaw->CalculateDeltaIsotropicHardening( rPlasticVariables.EquivalentPlasticStrain, rProperties );

			DeltaFunction  = 2.0 * mReturnMappingVariables.LameMu_bar;
			DeltaFunction += (DeltaKinematicHardening + DeltaIsotropicHardening) * (2.0/3.0);

			//Calculate DeltaGamma:
			DeltaDeltaGamma  = StateFunction/DeltaFunction;
			mReturnMappingVariables.DeltaGamma += DeltaDeltaGamma;
	       
			//Update Equivalent Plastic Strain:
			rPlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * mReturnMappingVariables.DeltaGamma;
			rPlasticVariables.EquivalentPlasticStrain  = InitialEquivalentPlasticStrain + rPlasticVariables.DeltaPlasticStrain;
	       
			//Calculate Function:
			KinematicHardening = mpHardeningLaw->CalculateKinematicHardening( rPlasticVariables.EquivalentPlasticStrain );
			IsotropicHardening = mpHardeningLaw->CalculateIsotropicHardening( rPlasticVariables.EquivalentPlasticStrain );

			StateFunction  = mReturnMappingVariables.NormIsochoricStress;
			StateFunction -= 2.0 * LameMu_bar * mReturnMappingVariables.DeltaGamma;
			StateFunction -= sqrt(2.0/3.0) * ( IsotropicHardening + ( KinematicHardening - InitialKinematicHardening ));

			iter++;
		}
	   

		if(iter>MaxIterations)
			return false;

		return true;	
	}


	void Update( Matrix & rIsoStressMatrix )
	{
		//Stress Update: 
		double Auxiliar  = 2.0 * mReturnMappingVariables.LameMu_bar * mReturnMappingVariables.DeltaGamma;
		Matrix Normal    = rIsoStressMatrix * ( 1.0 / mReturnMappingVariables.NormIsochoricStress );
		IsoStressMatrix -= ( Normal * Auxiliar );
	}


	bool InternalVariablesUpdate( const Matrix & rIsoStressMatrix )
	{
		YieldCriterion::InternalVariables&  PlasticVariables = YieldCriterion->FastGetInternalVariables();
		
		PlasticVariables.DeltaPlasticStrain       = sqrt(2.0/3.0) * mReturnMappingVariables.DeltaGamma;
		PlasticVariables.EquivalentPlasticStrain += PlasticVariables.DeltaPlasticStrain;

		//To update in the constitutive law, ElasticLeftCauchyGreen is not an internal variable.

		Matrix IdentityMatrix  = identity_matrix<double> (3);

		const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
		const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

		LameMu  = YoungModulus/(2*(1+PoissonCoefficient));

		PlasticVariables.ElasticLeftCauchyGreen   = ( IsoStressMatrix * ( 1.0 / LameMu ) );
		PlasticVariables.ElasticLeftCauchyGreen  += ( (rTrace_b_bar/3.0) * IdentityMatrix );	
	}

	void CalculateScalingFactors( const Matrix & rIsoStressMatrix, PlasticFactors& rScalingFactors )
	{
	
		//1.-Identity build
		Matrix IdentityMatrix  = identity_matrix<double> (3);

		//3.-Particular Parameters
		double YieldStress           =  GetProperties()[YIELD_STRESS];
		double KinematicHardening    =  GetProperties()[KINEMATIC_HARDENING];
		double Delta                 =  GetProperties()[HARDENING_EXPONENT];
		double K_reference           =  GetProperties()[REFERENCE_HARDENING];
		double K_infinity            =  GetProperties()[INFINITY_HARDENING];

		const double& YoungModulus       = GetProperties()[YOUNG_MODULUS];
		const double& PoissonCoefficient = GetProperties()[POISSON_RATIO];
		
		//1.-Auxiliar matrices
		rScalingFactors.Normal      = rIsoStressMatrix * ( 1.0 / mReturnMappingVariables.NormIsochoricStress );

		Matrix Norm_Normal          = prod( rScalingFactors.Normal, trans(rScalingFactors.Normal) );
		double Trace_Norm_Normal    = Norm_Normal( 0, 0 ) + Norm_Normal( 1, 1 ) + Norm_Normal( 2, 2 );

		rScalingFactors.Dev_Normal  = Norm_Normal;
		rScalingFactors.Dev_Normal -= (1.0/3.0) * Trace_Norm_Normal * IdentityMatrix;


		//2.-Auxiliar constants
		const YieldCriterion::InternalVariables&  HistoricPlasticVariables = YieldCriterion->FastGetInternalVariables();
		double EquivalentPlasticStrain = HistoricPlasticVariables.EquivalentPlasticStrain + sqrt(2.0/3.0) * mReturnMappingVariables.DeltaGamma;
		double Hardening               = mpHardeningLaw->CalculateHardening( EquivalentPlasticStrain );

		rScalingFactors.Beta0 = 1.0 + Hardening/(3.0 * mReturnMappingVariables.LameMu_bar);
		
		rScalingFactors.Beta1 = 2.0 * LameMu_bar * DeltaGamma / rParameters.NormIsochoricStress;
		
		rScalingFactors.Beta2 = ( ( 1.0 - ( 1.0 / Beta0 ) ) * (2.0/3.0) * mReturnMappingVariables.NormIsochoricStress * mReturnMappingVariables.DeltaGamma )/(mReturnMappingVariables.LameMu_bar) ;
		
		rScalingFactors.Beta3 = ( ( 1.0 / Beta0 ) - Beta1 + Beta2 );
		
		rScalingFactors.Beta4 = ( ( 1.0 / Beta0 ) - Beta1 ) * mReturnMappingVariables.NormIsochoricStress / ( mReturnMappingVariables.LameMu_bar ) ;
	
	};

///@}


}  // namespace Kratos.

#endif // KRATOS_FLOW_RULE_H_INCLUDED  defined 


