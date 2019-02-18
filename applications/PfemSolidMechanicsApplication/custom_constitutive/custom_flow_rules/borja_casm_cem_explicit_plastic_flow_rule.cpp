//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "custom_constitutive/custom_flow_rules/borja_casm_cem_explicit_plastic_flow_rule.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "pfem_solid_mechanics_application_variables.h"

#include "custom_utilities/stress_invariants_utilities.hpp"

namespace Kratos
{
	// Version of Casm-like hiperelastic model. This is the Hyperelsatic model of Borja, that includes the Houlsby and the constant shear modulus as special cases.

	//************ CONSTRUCTOR ***********
	BorjaCasmCemExplicitFlowRule::BorjaCasmCemExplicitFlowRule()
		:BorjaCasmExplicitFlowRule()
	{
	}

	//*****************************INITIALIZATION CONSTRUCTOR*****************************
	//************************************************************************************

	BorjaCasmCemExplicitFlowRule::BorjaCasmCemExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
		:BorjaCasmExplicitFlowRule(pYieldCriterion)
	{
		std::cout<<"   CASM-CEM FLOW RULE constructed "<<std::endl;
	}

  //********* ASSIGMENT OPERATOR
  BorjaCasmCemExplicitFlowRule& BorjaCasmCemExplicitFlowRule::operator=(BorjaCasmCemExplicitFlowRule const& rOther)
  {
		BorjaCasmExplicitFlowRule::operator=(rOther);
		mPlasticVariables = rOther.mPlasticVariables;
    return *this;
  }

  //********** COPY CONSTRUCTOR *********
  BorjaCasmCemExplicitFlowRule::BorjaCasmCemExplicitFlowRule(BorjaCasmCemExplicitFlowRule const& rOther)
		:BorjaCasmExplicitFlowRule(rOther)
		//mPlasticVariables(rOther.mPlasticVariables)
  {
  }

  //*******   CLONE ********
  FlowRule::Pointer BorjaCasmCemExplicitFlowRule::Clone() const
  {
		FlowRule::Pointer p_clone(new BorjaCasmCemExplicitFlowRule(*this));
		return p_clone;
  }

  // ********** DESTRUCTOR **************
  BorjaCasmCemExplicitFlowRule::~BorjaCasmCemExplicitFlowRule()
  {
  }


	// *****************************************************************
	// ***** EXPLICIT RETURN MAPPING ALGORITHM (Sloan et al. 2002) *****
	// *****************************************************************
  bool BorjaCasmCemExplicitFlowRule::CalculateReturnMappingExpl(RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
  {
		//std::cout<< "  " << std::endl;
		//std::cout<< "  ..................................... " << std::endl;
		//std::cout<< "  Explicit return mapping CASM cemented " << std::endl;
		//std::cout<< "  ..................................... " << std::endl;
		//std::cout<< "    F_incr:        " << rIncrementalDeformationGradient << std::endl;
		//std::cout<< "    b_n+1:         " << rNewElasticLeftCauchyGreen << std::endl;
		//std::cout<< "    StressMatrix:  " << rStressMatrix << std::endl;
		//std::cout<< "    ReturnMappingVariables: " << rReturnMappingVariables << std::endl;
		
    // I. Initialize some variables
//std::cout<< "  ... I. Initialization ... " << std::endl;
    bool PlasticityActive = false;
    rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);
		//get internal plastic variables
		PlasticVariablesType IntPlasticVariables = mPlasticVariables;
//mPlasticVariables.print();

		//reset IncrementalPlasticStrain and Total VolumetricStrain
		rReturnMappingVariables.PlasticVariables.EquivalentPlasticStrain		= IntPlasticVariables.EquivalentPlasticStrain;
		rReturnMappingVariables.PlasticVariables.DeltaEqPlasticStrain				= 0;
		rReturnMappingVariables.PlasticVariables.PlasticShearStrain					= IntPlasticVariables.PlasticShearStrain;
		rReturnMappingVariables.PlasticVariables.DeltaPlasticShearStrain		= 0;
		rReturnMappingVariables.PlasticVariables.PreconsolidationPressure		= IntPlasticVariables.PreconsolidationPressure;
		rReturnMappingVariables.PlasticVariables.Bonding										= IntPlasticVariables.Bonding;

		Matrix PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
		Vector PreviousStressVector;
		Vector NewStressVector;
		double Tolerance = 1.0E-4; // LMV
		
//std::cout<< "  ... I. Initialization completed ... " << std::endl<<std::endl;

		// II. Check for Yield Condition (at the beginning) and if necessary do yield surface correction
//std::cout<< "  ... II. Initial yield condition ... " << std::endl;
		this->CalculateKirchhoffStressVector( PreviousElasticLeftCauchyGreen, PreviousStressVector);
//std::cout<< "    b_n:   " << PreviousElasticLeftCauchyGreen << std::endl;
//std::cout<< "    tau_n: " << PreviousStressVector << std::endl;
		rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PreviousStressVector, rReturnMappingVariables.PlasticVariables);

		if (rReturnMappingVariables.TrialStateFunction > Tolerance) {
			this->ReturnStressToYieldSurface( rReturnMappingVariables, PreviousElasticLeftCauchyGreen, PreviousStressVector, rReturnMappingVariables.TrialStateFunction, Tolerance);
			//this->ReturnStressToYieldSurfaceNormal( rReturnMappingVariables, PreviousElasticLeftCauchyGreen, PreviousStressVector, rReturnMappingVariables.TrialStateFunction, Tolerance);
			//THE INITIAL TENSIONAL STATE MAY BE OUTSIDE THE YIELD SURFACE?? (yes and no: no if you don't do nothing strange, however, if you apply LSInterpolation or change constant parameters during the simulation YES)
		}
		rReturnMappingVariables.DeltaBeta = 0;

//std::cout<< "  ... II. Initial yield condition checked ... " << std::endl<<std::endl;


		// III. Compute an ELASTIC TRIAL STEP
//std::cout<< "  ... III. Elastic trial step ... " << std::endl;
		ExplicitStressUpdateInformation  StressUpdateInformation;
		this->CalculateOneExplicitStep(rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, false, StressUpdateInformation);
		
		double ElasticTrialStateFunction;
//IntPlasticVariables.print();
		ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, NewStressVector, IntPlasticVariables);
//std::cout<< "  ... III. Elastic trial step computed: tau_trial = "<< NewStressVector << " f = " << ElasticTrialStateFunction <<" ... "<< std::endl << std::endl;


		// IVa. PURELY ELASTIC STEP
		if ( ElasticTrialStateFunction <= Tolerance)   {
			PlasticityActive = false;
			
//std::cout<< "  ... IVa. Purely elastic step computed " << std::endl;
		}
		// IVb. SOME PART OF THE STEP IS PLASTIC
		else {
			PlasticityActive = true; 

			if ( ( rReturnMappingVariables.TrialStateFunction < -Tolerance) && (ElasticTrialStateFunction > Tolerance) ) {
				//2b.1 ELASTIC & ELASTO-PLASTIC STEP 
//std::cout<< "  ... IVb.1 Elastic & elasto-plastic step " << std::endl;
				this->CalculateExplicitSolutionWithChange( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, Tolerance);
//rReturnMappingVariables.PlasticVariables.print();
//std::cout<< "  ... IVb.1 Elastic & elasto-plastic step computed " << std::endl<<std::endl;
			}
			else {         
				bool UnloadingCondition;
				UnloadingCondition =  this->EvaluateElastoPlasticUnloadingCondition( UnloadingCondition, PreviousElasticLeftCauchyGreen, rIncrementalDeformationGradient, IntPlasticVariables, Tolerance); 

				if (UnloadingCondition ) {
					//2c.1 ELASTO-PLASTIC UNLOADING
//std::cout<< "  ... IVb.2 Elasto-plastic unloading step " << std::endl;
					this->CalculateExplicitSolutionWithChange( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, Tolerance);
//rReturnMappingVariables.PlasticVariables.print();
//std::cout<< "  ... IVb.2 Elasto-plastic unloading step computed " << std::endl<<std::endl;
				}
				else {
					//2c. 2 ELASTO-PLASTIC STEP
//std::cout<< "  ... IVb.3 Elasto-plastic loading step " << std::endl;
					this->CalculateExplicitSolution( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, true, Tolerance);
//rReturnMappingVariables.PlasticVariables.print();
//std::cout<< "  ... IVb.3 Elasto-plastic loading step computed " << std::endl<<std::endl;
				}
			}
		}


		// V. YIELD SURFACE DRIFT CORRECTION AFTER STRESS INTEGRATION
		if (PlasticityActive) {

			double DriftViolation;

			DriftViolation = mpYieldCriterion->CalculateYieldCondition(DriftViolation, NewStressVector, rReturnMappingVariables.PlasticVariables);
			//rReturnMappingVariables.EigenVectors = NewElasticLeftCauchyGreen; 
			double Beta = rReturnMappingVariables.DeltaBeta;
			if ( fabs(DriftViolation) > Tolerance ) {
//std::cout<< "  ... IVb. Yield surface drift violation: "<< fabs(DriftViolation) << std::endl;
				this->ReturnStressToYieldSurface( rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, DriftViolation, Tolerance);
//rReturnMappingVariables.PlasticVariables.print();
			}
			// LMV
			rReturnMappingVariables.DeltaBeta = Beta;
			
//std::cout<< "  ... IVb. Yield surface drift corrected " << std::endl;
		}

		//set flags and compute rStressMatrix
		rStressMatrix = MathUtils<double>::StressVectorToTensor(NewStressVector);

		rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
		rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

		return PlasticityActive;
  }


  // ***********************************************
	// ***** Compute plastic hardening modulus H *****
	// ***********************************************
	void BorjaCasmCemExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const PlasticVariablesType& rPlasticVariables, double& rH)
	{
//std::cout<<"BorjaCasmExplicitFlowRule::ComputePlasticHardeningParameter"<<std::endl;

		//get material constants
		const double ShearM 					= mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		const double SpacingR 				= mpYieldCriterion->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN 					= mpYieldCriterion->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		const double SwellingSlope		= mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
		const double OtherSlope      	= mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];
		const double H1 							= mpYieldCriterion->GetHardeningLaw().GetProperties()[DEGRADATION_RATE_COMPRESSION];
		const double H2 							= mpYieldCriterion->GetHardeningLaw().GetProperties()[DEGRADATION_RATE_SHEAR];
		const double Omega 						= mpYieldCriterion->GetHardeningLaw().GetProperties()[PLASTIC_DEVIATORIC_STRAIN_HARDENING];
		const double AlphaTensile 		= mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_TENSILE];

		//calculate hardening parameters
		double Pc, Pt;
		Pc = rPlasticVariables.PreconsolidationPressure*(1+rPlasticVariables.Bonding);
		Pt = rPlasticVariables.PreconsolidationPressure*(AlphaTensile*rPlasticVariables.Bonding);

//std::cout<<"casm-cem Pc: "<<Pc<<std::endl;
//std::cout<<"casm-cem Pt: "<<Pt<<std::endl;

		//calcualte Kirchhoff stress vector and invariants
		Vector StressVector = ZeroVector(6);
		this->CalculateKirchhoffStressVector(rHenckyStrainVector, StressVector);
		double MeanStress, J2, LodeAngle;
		StressInvariantsUtilities::CalculateStressInvariants( StressVector, MeanStress, J2, LodeAngle);
		
		//calculate third invariant effect 
		double ThirdInvariantEffect = 1.0;
		ThirdInvariantEffect = mpYieldCriterion->EvaluateThirdInvariantEffectMC( LodeAngle);

//std::cout<<"casm-cem P: "<<MeanStress<<std::endl;
//std::cout<<"casm-cem J2: "<<J2<<std::endl;
//std::cout<<"casm-cem Theta: "<<LodeAngle<<std::endl;
		
		//calculate derivatives with respect to P & J2
		double FirstDerivativeP = 0;
		double FirstDerivativeJ2 = 0;
		CalculatePlasticPotentialDerivativesPJ2(StressVector, FirstDerivativeP, FirstDerivativeJ2, rPlasticVariables);

//std::cout<<"casm-cem dP: "<<FirstDerivativeP<<std::endl;
//std::cout<<"casm-cem dJ2: "<<FirstDerivativeJ2<<std::endl;

		
		//calculate d_b/d_h * d_h/d_gamma
		double dBdGamma = 0.0;
		dBdGamma = -rPlasticVariables.Bonding*( H1*fabs( FirstDerivativeP ) + H2*fabs( FirstDerivativeJ2 ) );		

//std::cout<<"casm-cem dBdGamma: "<<dBdGamma<<std::endl;

		
		//calculate d_P0/d_gamma
		double dP0dGamma = 0.0;
		dP0dGamma = ( -FirstDerivativeP + Omega * FirstDerivativeJ2 ) * rPlasticVariables.PreconsolidationPressure/(OtherSlope-SwellingSlope);

//std::cout<<"casm-cem dP0dGamma: "<<dP0dGamma<<std::endl;

		
		//calculate hardening modulus H
		//rH = 1/((Pc-Pt)*log(SpacingR)) * ( (1+rPlasticVariables.Bonding)*dP0dGamma + rPlasticVariables.PreconsolidationPressure*dBdGamma);
		double RH1 = 1/((Pc-Pt)*log(SpacingR)) * ( (1+rPlasticVariables.Bonding)*dP0dGamma + rPlasticVariables.PreconsolidationPressure*dBdGamma);
//std::cout<<"casm-cem rh: "<<rH<<std::endl;
		//rH -= ( (ShapeN*pow(J2*pow(3.0,0.5),ShapeN))/( pow(ShearM/ThirdInvariantEffect,ShapeN)*pow(-(Pc-Pt),ShapeN+1) ) + (Pc-MeanStress)/(log(SpacingR)*(MeanStress+Pt)*(Pc+Pt))) * ( (AlphaTensile*rPlasticVariables.Bonding)*dP0dGamma + rPlasticVariables.PreconsolidationPressure*AlphaTensile*dBdGamma);
		double RH2 = ( (ShapeN*pow(J2*pow(3.0,0.5),ShapeN))/( pow(ShearM/ThirdInvariantEffect,ShapeN)*pow(-(Pc-Pt),ShapeN+1) ) + (Pc-MeanStress)/(log(SpacingR)*(MeanStress+Pt)*(Pc+Pt))) * ( (AlphaTensile*rPlasticVariables.Bonding)*dP0dGamma + rPlasticVariables.PreconsolidationPressure*AlphaTensile*dBdGamma);
		rH = RH1 - RH2;
//std::cout<<"casm-cem rh: "<<rH<<std::endl;

		if ( std::isnan(rH) )
	    {
	    	std::cout<<std::endl<<"NaN in hardening modulus: "<< rH <<std::endl;
	    	std::cout<<"  P0: "<< rPlasticVariables.PreconsolidationPressure <<", b: "<< rPlasticVariables.Bonding <<std::endl;
	    	std::cout<<"  Pc: "<< Pc <<", Pt: "<< Pt <<std::endl;
	    	std::cout<<"  dP0dGamma: "<< dP0dGamma <<", dBdGamma: "<< dBdGamma <<std::endl;
	    	std::cout<<"  dGdP: "<< FirstDerivativeP <<std::endl;
	    	std::cout<<"  dGdJ2: "<< FirstDerivativeJ2 <<std::endl;
	    	std::cout<<"  RH1: "<< RH1 <<std::endl;
	    	std::cout<<"  RH2: "<< RH2 <<std::endl;
	    	std::cout<<"  RH2a: "<< ( (ShapeN*pow(J2*pow(3.0,0.5),ShapeN))/( pow(ShearM/ThirdInvariantEffect,ShapeN)*pow(-(Pc-Pt),ShapeN+1) ) + (Pc-MeanStress)/(log(SpacingR)*(MeanStress+Pt)*(Pc+Pt))) <<std::endl;
	    	std::cout<<"  RH2a1: "<< ( ShapeN*pow(J2*pow(3.0,0.5),ShapeN) )  <<std::endl;
	    	std::cout<<"  RH2a2: "<< ( pow(ShearM/ThirdInvariantEffect,ShapeN)*pow(-(Pc-Pt),ShapeN+1) ) <<std::endl;
	    	std::cout<<"  RH2a21: "<< ( pow(ShearM/ThirdInvariantEffect,ShapeN) ) <<std::endl;
	    	std::cout<<"  ShearM: "<< ( ShearM ) <<std::endl;
	    	std::cout<<"  ThirdInvariantEffect: "<< ( ThirdInvariantEffect ) <<std::endl;
	    	std::cout<<"  RH2a22: "<< ( pow(-(Pc-Pt),ShapeN+1) ) <<std::endl;
	    	std::cout<<"  -(Pc-Pt): "<< (-(Pc-Pt) ) <<std::endl;
	    	std::cout<<"  RH2a3: "<< ( (Pc-MeanStress)/(log(SpacingR)*(MeanStress+Pt)*(Pc+Pt)) ) <<std::endl;
	    	std::cout<<"  RH2b: "<< ( (AlphaTensile*rPlasticVariables.Bonding)*dP0dGamma + rPlasticVariables.PreconsolidationPressure*AlphaTensile*dBdGamma) <<std::endl;
	    }

	}
	
	
	// *************************************************************************
	// ***** Calculate plastic potential derivative with respect to P & J2 *****
	// *************************************************************************
	void BorjaCasmCemExplicitFlowRule::CalculatePlasticPotentialDerivativesPJ2( const Vector& rStressVector, double& rFirstDerivativeP, double& rFirstDerivativeJ2, const PlasticVariablesType& rPlasticVariables)
	{
		//get inclination of CS-line
		const double ShearM = mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
		const double Lambda = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];
		const double Kappa = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
		const double SpacingR = mpYieldCriterion->GetHardeningLaw().GetProperties()[SPACING_RATIO];
		const double ShapeN = mpYieldCriterion->GetHardeningLaw().GetProperties()[SHAPE_PARAMETER];
		const double AlphaTensile = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_TENSILE];
		
		//calcualte m-factor
		const double FactorM = 2.0/3.0 * (pow(ShearM*(6-ShearM),ShapeN) - pow(3.0*ShearM,ShapeN))/( (1-Kappa/Lambda)*(6.0-ShearM)*pow(3.0*ShearM,ShapeN-1) );
		
		//stress invariants & invariants derivatives
		double MeanStress, J2, LodeAngle;
		StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
		
		//calculate third invariant effect 
		double ThirdInvariantEffect = 1.0;
		ThirdInvariantEffect = mpYieldCriterion->EvaluateThirdInvariantEffectMC( LodeAngle);
		
		// calculate hardening parameters p0, b -> pt, pc
		double Pt = rPlasticVariables.PreconsolidationPressure*(AlphaTensile*rPlasticVariables.Bonding);
		
		const double flowID = 5;
		if (flowID == 1)
		{
			//Rowe
			rFirstDerivativeP = 27.0*( ShearM/ThirdInvariantEffect*MeanStress + pow(3.0,0.5)*J2 ) / ( (3.0*MeanStress-2.0*pow(3.0,0.5)*J2) * (3.0*MeanStress+pow(3.0,0.5)*J2) );
			rFirstDerivativeJ2 = 3.0*pow(3.0,0.5)*(-9.0*MeanStress - ShearM/ThirdInvariantEffect*(3.0*MeanStress + 2.0*pow(3.0,0.5)*J2)) / ( (3.0*MeanStress-2.0*pow(3.0,0.5)*J2) * (3.0*MeanStress+pow(3.0,0.5)*J2) );
std::cout<<"casm-cem ROWE: dp = "<<rFirstDerivativeP<<" dJ = "<<rFirstDerivativeJ2<<std::endl;
		}
		else if (flowID == 2)
		{
			//Yu
			rFirstDerivativeP = ( ShapeN*(FactorM-1.0)/MeanStress + ( ShapeN*FactorM*(FactorM-1.0)*pow(pow(3.0,0.5)*J2/(ShearM/ThirdInvariantEffect*(-MeanStress)),ShapeN) )/( -MeanStress*(1.0+(FactorM-1.0)* pow(pow(3.0,0.5)*J2/(ShearM/ThirdInvariantEffect*(-MeanStress)),ShapeN)) ) );
			rFirstDerivativeJ2 = ( ShapeN*FactorM*(FactorM-1.0)*pow(pow(3.0,0.5)*J2/(ShearM/ThirdInvariantEffect*(-MeanStress)),ShapeN) )/( J2*(1.0+(FactorM-1.0)* pow(pow(3.0,0.5)*J2/(ShearM/ThirdInvariantEffect*(-MeanStress)),ShapeN)) );
		}
		else if (flowID == 3)
		{
			//alternative
			const double Beta = 1.00;
			rFirstDerivativeP = 1/( (MeanStress+Pt) * log(SpacingR) ) + ( ShapeN * pow( Beta*pow(3,1/2)*J2 , ShapeN) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-(MeanStress+Pt),ShapeN+1) );
			rFirstDerivativeJ2 = ( ShapeN * pow(3*Beta,ShapeN/2) * pow(J2,ShapeN-1) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-(MeanStress+Pt),ShapeN) );
		}
		else if (flowID == 4)
		{
			//alternative
			rFirstDerivativeP = 1/( (MeanStress+Pt) * log(SpacingR) ) + ( ShapeN * pow( pow(3,1/2)*J2 , ShapeN) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-(MeanStress+Pt),ShapeN+1) );
			rFirstDerivativeJ2 = ( ShapeN * pow(3,ShapeN/2) * pow(J2,ShapeN-1) )/( pow(ShearM/ThirdInvariantEffect,ShapeN) * pow(-(MeanStress+Pt),ShapeN) );
		}
		else if (flowID == 5)
		{
			//MCC
			double chiMCC = MeanStress + 3.0/MeanStress*pow( J2*ThirdInvariantEffect/ShearM, 2 );
			rFirstDerivativeP = 2.0*MeanStress - chiMCC;
			rFirstDerivativeJ2 = 2.0 * 3.0 * pow( ThirdInvariantEffect / ShearM, 2) * J2;
//std::cout<<"casm-cem MCC: dp = "<<rFirstDerivativeP<<" dJ = "<<rFirstDerivativeJ2<<std::endl;
		}
	}

	// *************************************************************
	// ***** Set initial bonding and preconsolidation pressure *****
	//
	void BorjaCasmCemExplicitFlowRule::SetPlasticVariables( const double& rInitialPreconPressure, const double& rInitialBonding)
	{
		mPlasticVariables.Bonding = rInitialBonding;
		mPlasticVariables.PreconsolidationPressure = rInitialPreconPressure;
	}

	// *******************************
	// ***** Set initial bonding *****
	//
	void BorjaCasmCemExplicitFlowRule::SetBonding( const double& rInitialBonding)
	{
		mPlasticVariables.Bonding = rInitialBonding;
	}


	// ********** SAVE FUNCTIONS **********
	// ************************************
  	void BorjaCasmCemExplicitFlowRule::save( Serializer& rSerializer) const 
  	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BorjaCasmExplicitFlowRule )
		rSerializer.save("mPlasticVariables",mPlasticVariables);
  	}

  	void BorjaCasmCemExplicitFlowRule::load( Serializer& rSerializer)
  	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BorjaCasmExplicitFlowRule )
		rSerializer.load("mPlasticVariables",mPlasticVariables);
  	}
  	
} //end namespace kratos
